"""
Q11: Algorithm 3 redesign -- tests T0-T12 and timing (mean +- SD).
T0 benign accept; T1-T8 as in E10 (re-implemented for the v2 signature design); 
T9  valid old ledger (valid session sig for an OLD nonce) replayed for a new delivery -> reject
T10 ledger (with valid session signature) issued for another request (s',d') -> reject
T11 path substitution: entries swapped for entries of another valid ledger (all cert sigs valid) -> reject
T12 lying endpoint signs an inflated certificate with ITS OWN valid key: ledger-only check ACCEPTS (documented limit);
    a DB that stores the independent audit record catches it (db: record mismatch) -- reported separately.
"""
import sys, os, json, time, platform, hashlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from bmr import verify as V
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print("code sha256:", hashlib.sha256(open(__file__, "rb").read()).hexdigest()[:16])
rng = np.random.default_rng(2026)
NT = 2000
NOW, R, DELTA, W, T_COLL = 10, 0.05, 1e-2, 8, 1
NODES = list(range(12))
reg = V.KeyRegistry(NODES)
forger = V.KeyRegistry(NODES)             # attacker's own keys (not in verifier's registry)


def mk_path(nodes_seq, now=NOW):
    ents = []
    for a, b in zip(nodes_seq[:-1], nodes_seq[1:]):
        link = (a, b)
        ents.append(V.make_cert_entry_v2(reg, link, 0.11 + float(rng.uniform(0, .1)), 1e-7, 1e-6,
                                         now - int(rng.integers(0, W))))
    return ents


def db_of(ents):
    return {(e["link"], e["T"]): (e["r"], e["eps"], e["delta"]) for e in ents}


def session(ents, s, d, nonce, epoch=NOW):
    return V.bind_session(reg, ents, nonce, s, d, epoch)


def run(tamper):
    acc = 0; reasons = {}
    for _ in range(NT):
        nonce = bytes(rng.integers(0, 256, 16, dtype=np.uint8))
        seq = [0, 1, 2, 3, 4]
        ents = mk_path(seq)
        db = db_of(ents)
        ents2, sess, db2, s, d, exp_nonce = tamper(ents, db, nonce, seq)
        ok, why = V.verify_ledger_v2(ents2, sess, reg, s, d, R, DELTA, W, T_COLL, NOW, db2, exp_nonce)
        acc += ok; reasons[why] = reasons.get(why, 0) + 1
    return acc / NT, reasons


def benign(ents, db, nonce, seq): return ents, session(ents, 0, 4, nonce), db, 0, 4, nonce
def T1(ents, db, nonce, seq):
    e = [dict(x) for x in ents]; i = int(rng.integers(len(e))); e[i]["r"] += .5
    return e, session(ents, 0, 4, nonce), db, 0, 4, nonce          # session digest of original -> reject
def T1b(ents, db, nonce, seq):                                      # tamper AND re-bind session with forger key
    e = [dict(x) for x in ents]; i = int(rng.integers(len(e))); e[i]["r"] += .5
    return e, V.bind_session(forger, e, nonce, 0, 4, NOW), db, 0, 4, nonce
def T2(ents, db, nonce, seq):                                       # stale: re-sign with the REAL attestor key, old T
    e = [dict(x) for x in ents]; i = int(rng.integers(len(e)))
    e[i] = V.make_cert_entry_v2(reg, e[i]["link"], e[i]["r"], e[i]["eps"], e[i]["delta"], NOW - 100)
    return e, session(e, 0, 4, nonce), db, 0, 4, nonce
def T3(ents, db, nonce, seq):
    e = [ents[0], ents[-1]]; return e, session(e, 0, 4, nonce), db, 0, 4, nonce
def T4(ents, db, nonce, seq):                                       # loop: add back-edge (1,0)
    e = list(ents); e.insert(1, V.make_cert_entry_v2(reg, (0, 1), .2, 1e-7, 1e-6, NOW - 1)); e = [e[0], e[1]] + e[2:]
    loop = [ents[0], V.make_cert_entry_v2(reg, (1, 0), .2, 1e-7, 1e-6, NOW - 1)] + ents[0:]
    return loop, session(loop, 0, 4, nonce), db, 0, 4, nonce
def T5(ents, db, nonce, seq):
    e = [dict(x) for x in ents]; i = int(rng.integers(len(e))); e[i]["eps"] /= 1000
    return e, session(ents, 0, 4, nonce), db, 0, 4, nonce
def T6(ents, db, nonce, seq):
    e = list(ents); e[0], e[1] = e[1], e[0]; return e, session(e, 0, 4, nonce), db, 0, 4, nonce
def T7(ents, db, nonce, seq):
    d2 = dict(db); d2[(ents[1]["link"], NOW)] = (0.0, 0.0, 0.0)          # DEGRADED strictly after T_i (T_i<=NOW-?)
    # ensure T_i < NOW: force the entry's T below NOW
    e = [dict(x) for x in ents]
    e[1] = V.make_cert_entry_v2(reg, e[1]["link"], e[1]["r"], 1e-7, 1e-6, NOW - 2)
    d2 = db_of(e); d2[(e[1]["link"], NOW - 1)] = (0.0, 0.0, 0.0)
    return e, session(e, 0, 4, nonce), d2, 0, 4, nonce
def T7b(ents, db, nonce, seq):                                      # DEGRADED at an epoch in (T_i, now) -- v1 check also caught; new: bounded by now
    e = [dict(x) for x in ents]
    e[0] = V.make_cert_entry_v2(reg, e[0]["link"], e[0]["r"], 1e-7, 1e-6, NOW - 4)
    d2 = db_of(e); d2[(e[0]["link"], NOW - 2)] = (0.0, 0.0, 0.0)
    return e, session(e, 0, 4, nonce), d2, 0, 4, nonce
def T8(ents, db, nonce, seq):                                       # re-signed by non-attestor with own key
    e = [dict(x) for x in ents]; i = int(rng.integers(len(e))); r = e[i]["r"] + .1
    e[i] = V.make_cert_entry_v2(forger, e[i]["link"], r, e[i]["eps"], e[i]["delta"], e[i]["T"])
    return e, session(e, 0, 4, nonce), db, 0, 4, nonce
def T9(ents, db, nonce, seq):                                       # old ledger + old valid session, new delivery nonce
    old_nonce = bytes(rng.integers(0, 256, 16, dtype=np.uint8))
    return ents, session(ents, 0, 4, old_nonce), db, 0, 4, nonce
def T10(ents, db, nonce, seq):                                      # valid ledger for (0,4) presented as request (0,5)? use (s,d)=(0,4) ledger for destination 4 but session for another d'
    return ents, session(ents, 0, 3, nonce), db, 0, 4, nonce
def T10b(ents, db, nonce, seq):                                     # ledger of another request (7->11), presented to verifier of (0,4)
    other = mk_path([7, 8, 9, 10, 11])
    return other, session(other, 7, 11, nonce), db, 0, 4, nonce
def T11(ents, db, nonce, seq):                                      # swap one entry for the entry of another valid ledger (same link, different valid cert)
    alt = mk_path(seq)
    e = list(ents); i = int(rng.integers(len(e))); e[i] = alt[i]
    return e, session(ents, 0, 4, nonce), db, 0, 4, nonce
def T11b(ents, db, nonce, seq):                                     # whole different valid path s->d, session of the original
    alt = mk_path([0, 5, 6, 7, 4])
    return alt, session(ents, 0, 4, nonce), db, 0, 4, nonce
def T12(ents, db, nonce, seq):                                      # lying attestor (lower-id endpoint) signs inflated r with ITS OWN valid key
    e = [dict(x) for x in ents]; i = int(rng.integers(len(e)))
    e[i] = V.make_cert_entry_v2(reg, e[i]["link"], e[i]["r"] + 0.3, e[i]["eps"], e[i]["delta"], e[i]["T"])
    return e, session(e, 0, 4, nonce), {}, 0, 4, nonce              # ledger-only: no audit record in DB
def T12db(ents, db, nonce, seq):                                    # same lie, DB holds the independent audit record
    e, s_, _, a, b, c = T12(ents, db, nonce, seq)
    return e, s_, db, a, b, c

tests = [("T0 benign", benign, "accept"), ("T1 inflate r (session digest)", T1, "reject"),
         ("T1b inflate r + re-bind with forger key", T1b, "reject"), ("T2 stale certificate", T2, "reject"),
         ("T3 drop a hop", T3, "reject"), ("T4 loop", T4, "reject"), ("T5 lower eps", T5, "reject"),
         ("T6 reorder", T6, "reject"), ("T7 DEGRADED in (T_i,now] (late-entry)", T7, "reject"),
         ("T7b DEGRADED in (T_i,now] (early-entry)", T7b, "reject"),
         ("T8 re-signed by non-attestor", T8, "reject"), ("T9 replay old ledger+session (old nonce)", T9, "reject"),
         ("T10 session bound to another (s,d')", T10, "reject"), ("T10b ledger of another request", T10b, "reject"),
         ("T11 entry swapped for another valid cert", T11, "reject"), ("T11b substitute other valid path", T11b, "reject"),
         ("T12 lying attestor, ledger-only", T12, "ACCEPT (documented limit)"),
         ("T12-DB lying attestor, DB has audit record", T12db, "reject (DB cross-check)")]
rows = []
print(f"{'test':50s} {'expected':28s} accept%  reasons")
for name, fn, exp in tests:
    a, why = run(fn)
    ok = (a > .99) if exp == "accept" else (a < .01) if exp.startswith("reject") else (a > .99)
    rows.append({"test": name, "expected": exp, "accept_rate": a, "reasons": why, "as_expected": bool(ok)})
    print(f"{name:50s} {exp:28s} {a*100:6.2f}  {why}")
assert all(r["as_expected"] for r in rows), "unexpected result"

# ---------------- timing ----------------
def timing(use_ed, reps=3000):
    rg = V.KeyRegistry(NODES, use_ed25519=use_ed)
    out = {}
    for m in (2, 4, 8, 11):
        seq = list(range(m + 1)); ts = []
        for _ in range(reps):
            ents = [V.make_cert_entry_v2(rg, (a, b), .2, 1e-7, 1e-6, NOW) for a, b in zip(seq[:-1], seq[1:])]
            nonce = b"n" * 16
            sess = V.bind_session(rg, ents, nonce, 0, m, NOW); db = db_of(ents)
            t0 = time.perf_counter()
            ok, _ = V.verify_ledger_v2(ents, sess, rg, 0, m, R, 1.0, W, (m - 1) // 2, NOW, db, nonce)
            ts.append((time.perf_counter() - t0) * 1e6); assert ok
        ts = np.array(ts[100:])             # drop warm-up
        per = ts / m
        out[m] = {"total_us_mean": float(ts.mean()), "total_us_sd": float(ts.std(ddof=1)), "total_us_median": float(np.median(ts)),
                  "per_entry_us_mean": float(per.mean()), "per_entry_us_sd": float(per.std(ddof=1))}
    return out
import cryptography
tm_ed = timing(True); tm_hm = timing(False)
print("\nTiming (µs): Ed25519 via cryptography", cryptography.__version__, "; HMAC-SHA256 fallback; ", platform.platform(), platform.python_version())
for m in tm_ed:
    e, h = tm_ed[m], tm_hm[m]
    print(f"  m={m:2d} hops: Ed25519 total {e['total_us_mean']:.1f}±{e['total_us_sd']:.1f} (per entry {e['per_entry_us_mean']:.1f}±{e['per_entry_us_sd']:.1f}) | "
          f"HMAC total {h['total_us_mean']:.1f}±{h['total_us_sd']:.1f} (per entry {h['per_entry_us_mean']:.1f}±{h['per_entry_us_sd']:.1f})")
json.dump({"tests": rows, "n_per_test": NT, "timing_ed25519": tm_ed, "timing_hmac": tm_hm,
           "library": f"cryptography {cryptography.__version__} (OpenSSL Ed25519)", "machine": platform.platform(),
           "cpu": platform.processor() or platform.machine(), "python": platform.python_version(),
           "note": "per-entry timing includes the one session-binding verification amortised over m entries"},
          open(os.path.join(ROOT, "data", "results", "q11_alg3_v2.json"), "w"), indent=1)
