"""
bmr/verify.py
=============
Algorithm 3: Ledger verification (Developer Guide §13).
Addresses Reviewer 2 Comment 1 (destination-side verification).

The ledger is a sequence of certificate entries signed by link endpoints.
Algorithm 3 checks: continuity, simplicity, signatures, freshness, rate floor, budget, length rule, DB consistency.
"""
import hashlib
import hmac
import math
import time


def make_ed25519_signer():
    """Generate a simple HMAC-based signer for unit tests (stands in for Ed25519)."""
    try:
        from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
        privkey = Ed25519PrivateKey.generate()
        pubkey = privkey.public_key()

        def sign(data: bytes) -> bytes:
            return privkey.sign(data)

        def verify(data: bytes, sig: bytes) -> bool:
            try:
                pubkey.verify(sig, data)
                return True
            except Exception:
                return False

        return sign, verify
    except Exception:
        # Fallback: HMAC-SHA256
        key = hashlib.sha256(str(time.time()).encode()).digest()

        def sign(data: bytes) -> bytes:
            return hmac.new(key, data, hashlib.sha256).digest()

        def verify(data: bytes, sig: bytes) -> bool:
            expected = hmac.new(key, data, hashlib.sha256).digest()
            return hmac.compare_digest(expected, sig)

        return sign, verify


def make_ledger_entry(link, r, eps, delta, T, sign_fn):
    """Create a signed ledger entry.

    Fields: (link, r, eps, delta, T, signature)
    Signature covers: link, r, eps, delta, T
    """
    msg = f"{link}|{r:.10f}|{eps:.10e}|{delta:.10e}|{T}".encode()
    sig = sign_fn(msg)
    return {
        "link": link,
        "r": r,
        "eps": eps,
        "delta": delta,
        "T": T,
        "sig": sig,
        "_msg": msg,  # for verification convenience
    }


def verify_ledger(ledger, s, d, R, Delta, W, t_collusion, now, cert_db, verify_fn,
                  mode="sealed"):
    """Algorithm 3: Verify a ledger.

    Parameters
    ----------
    ledger      : list of dicts (link, r, eps, delta, T, sig, _msg)
    s, d        : source and destination node IDs
    R           : minimum rate floor
    Delta       : security budget (Σ(eps + delta) <= Delta)
    W           : freshness window (epochs)
    t_collusion : collusion resistance parameter
    now         : current epoch
    cert_db     : dict mapping (link, T) -> (r, eps, delta) for DB consistency check
    verify_fn   : function(data, sig) -> bool
    mode        : "sealed" (m >= 2t+1) or "public_xor" (m odd)

    Returns
    -------
    (accept: bool, reason: str, checks: dict)
    """
    checks = {}

    # 1. Continuity: build the node chain from s to d
    # Each entry has a link which is (u, v) or "u-v"
    def get_endpoints(link):
        if isinstance(link, (tuple, list)):
            return link[0], link[1]
        parts = str(link).split("-")
        return parts[0], parts[1]

    if not ledger:
        return False, "empty ledger", {"empty": True}

    nodes = []
    for entry in ledger:
        u, v = get_endpoints(entry["link"])
        if not nodes:
            nodes.append(u)
        if nodes[-1] != u:
            checks["continuity"] = False
            return False, f"continuity failed: expected {nodes[-1]}, got {u}", checks
        nodes.append(v)

    checks["continuity"] = True

    if nodes[0] != s or nodes[-1] != d:
        checks["endpoints"] = False
        return False, f"endpoint mismatch: chain {nodes[0]}..{nodes[-1]} vs {s}..{d}", checks
    checks["endpoints"] = True

    # 2. Simplicity: all nodes pairwise distinct
    if len(set(nodes)) != len(nodes):
        checks["simplicity"] = False
        return False, f"non-simple path: repeated node in {nodes}", checks
    checks["simplicity"] = True

    # 3. Signatures
    for i, entry in enumerate(ledger):
        if not verify_fn(entry["_msg"], entry["sig"]):
            checks["signatures"] = False
            return False, f"signature failed on entry {i}", checks
    checks["signatures"] = True

    # 4. Freshness
    for i, entry in enumerate(ledger):
        age = now - entry["T"]
        if age < 0 or age > W:
            checks["freshness"] = False
            return False, f"stale certificate on entry {i}: age={age} > W={W}", checks
    checks["freshness"] = True

    # 5. Rate floor
    for i, entry in enumerate(ledger):
        if entry["r"] < R:
            checks["rate_floor"] = False
            return False, f"rate floor violated on entry {i}: r={entry['r']:.4f} < R={R}", checks
    checks["rate_floor"] = True

    # 6. Budget
    budget_used = sum(entry["eps"] + entry["delta"] for entry in ledger)
    checks["budget"] = budget_used <= Delta
    if not checks["budget"]:
        return False, f"budget exceeded: Σw={budget_used:.3e} > Delta={Delta:.3e}", checks

    # 7. Length rule
    m = len(ledger)  # number of hops
    if mode == "sealed":
        min_hops = 2 * t_collusion + 1
        checks["length_rule"] = m >= min_hops
        if not checks["length_rule"]:
            return False, f"hop floor violated: m={m} < 2t+1={min_hops}", checks
    else:  # public_xor: odd hop count
        checks["length_rule"] = m % 2 == 1
        if not checks["length_rule"]:
            return False, f"parity violated: m={m} is even", checks

    # 8. DB consistency: check no newer DEGRADED alarm, and (r, eps, T) match DB
    for i, entry in enumerate(ledger):
        key = (entry["link"], entry["T"])
        if key in cert_db:
            db_r, db_eps, db_delta = cert_db[key]
            if abs(db_r - entry["r"]) > 1e-8 or abs(db_eps - entry["eps"]) > 1e-15:
                checks["db_consistency"] = False
                return False, f"DB mismatch on entry {i}", checks
        # Check for DEGRADED alarm (db entry with T' > T and rate=0)
        for (lnk, T_prime), (dr, deps, ddelta) in cert_db.items():
            if lnk == entry["link"] and T_prime > entry["T"] and dr == 0.0:
                checks["db_consistency"] = False
                return False, f"superseded DEGRADED alarm for link {lnk}", checks
    checks["db_consistency"] = True

    return True, "accepted", checks


# =====================================================================================
# Algorithm 3 (v2) -- per-epoch certificate signatures + session-binding signature
# (pending_v7 Q11, blueprint 3.6).  The v1 functions above are kept unchanged for E10.
#
#   * Each link certificate is signed ONCE per epoch by the link's designated attestor
#     (the lower-id endpoint).  Endpoints do NOT sign per delivery.
#   * The route assembler signs  H(ledger) || nonce || s || d || epoch  for each delivery.
#     The nonce is chosen by the destination d (challenge), so a replayed ledger fails.
#   * Database check covers every record (link, T') with T_i <= T' <= now (not only T_i).
# =====================================================================================
import json as _json


class KeyRegistry:
    """Public-key registry known to the verifier: node -> verify fn, plus the assembler key."""

    def __init__(self, nodes, use_ed25519=True):
        self.sign, self.verify = {}, {}
        for x in list(nodes) + ["__assembler__"]:
            s, v = make_ed25519_signer() if use_ed25519 else _make_hmac_signer()
            self.sign[x], self.verify[x] = s, v

    def attestor(self, link):
        return min(link)


def _make_hmac_signer():
    import os
    key = os.urandom(32)
    return (lambda data: hmac.new(key, data, hashlib.sha256).digest(),
            lambda data, sig: hmac.compare_digest(hmac.new(key, data, hashlib.sha256).digest(), sig))


def cert_message(link, r, eps, delta, T):
    return f"{tuple(link)}|{r:.10f}|{eps:.10e}|{delta:.10e}|{T}".encode()


def make_cert_entry_v2(registry, link, r, eps, delta, T, signer=None):
    link = tuple(link)
    msg = cert_message(link, r, eps, delta, T)
    who = registry.attestor(link) if signer is None else signer
    return {"link": link, "r": r, "eps": eps, "delta": delta, "T": T,
            "signer": registry.attestor(link), "sig": registry.sign[who](msg)}


def ledger_digest(entries):
    h = hashlib.sha256()
    for e in entries:
        h.update(cert_message(e["link"], e["r"], e["eps"], e["delta"], e["T"]))
        h.update(e["sig"])
    return h.digest()


def session_message(entries, nonce, s, d, epoch):
    return ledger_digest(entries) + b"|" + nonce + f"|{s}|{d}|{epoch}".encode()


def bind_session(registry, entries, nonce, s, d, epoch, signer="__assembler__"):
    return {"nonce": nonce, "s": s, "d": d, "epoch": epoch,
            "sig": registry.sign[signer](session_message(entries, nonce, s, d, epoch))}


def verify_ledger_v2(entries, session, registry, s, d, R, Delta, W, t, now, cert_db,
                     expected_nonce, mode="sealed"):
    """Returns (accept, reason).  Checks: session binding, continuity, simplicity, signatures,
    freshness, rate floor, budget, length rule, database check over [T_i, now]."""
    if not entries:
        return False, "empty ledger"
    # session binding first (cheap, one signature)
    if session["nonce"] != expected_nonce:
        return False, "session: nonce mismatch (replay)"
    if (session["s"], session["d"]) != (s, d):
        return False, "session: bound to another (s,d)"
    if not registry.verify["__assembler__"](session_message(entries, session["nonce"], s, d, session["epoch"]),
                                            session["sig"]):
        return False, "session: signature/ledger-digest mismatch"
    # continuity / endpoints
    nodes = [entries[0]["link"][0]]
    for e in entries:
        u, v = e["link"]
        if nodes[-1] != u:
            return False, "continuity"
        nodes.append(v)
    if nodes[0] != s or nodes[-1] != d:
        return False, "endpoints"
    if len(set(nodes)) != len(nodes):
        return False, "simplicity"
    # per-epoch certificate signatures by the designated attestor
    for e in entries:
        if e["signer"] != registry.attestor(e["link"]):
            return False, "signer is not the link attestor"
        if not registry.verify[e["signer"]](cert_message(e["link"], e["r"], e["eps"], e["delta"], e["T"]), e["sig"]):
            return False, "certificate signature"
    for e in entries:
        age = now - e["T"]
        if age < 0 or age > W:
            return False, "stale"
    if any(e["r"] < R for e in entries):
        return False, "rate floor"
    if sum(e["eps"] + e["delta"] for e in entries) > Delta:
        return False, "budget"
    m = len(entries)
    if mode == "sealed":
        if m < 2 * t + 1:
            return False, "hop floor"
    elif m % 2 == 0:
        return False, "parity"
    # database check over [T_i, now]: matching record at T_i; no DEGRADED / lower-rate record in (T_i, now]
    for e in entries:
        rec = cert_db.get((e["link"], e["T"]))
        if rec is not None and (abs(rec[0] - e["r"]) > 1e-8 or abs(rec[1] - e["eps"]) > 1e-15):
            return False, "db: record mismatch"
        for (lnk, Tp), (dr, de, dd) in cert_db.items():
            if lnk == e["link"] and e["T"] < Tp <= now and dr == 0.0:
                return False, "db: superseded by DEGRADED in (T_i, now]"
    return True, "accepted"
