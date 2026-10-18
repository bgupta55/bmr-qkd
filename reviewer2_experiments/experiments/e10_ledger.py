"""
E10: Ledger verification and tamper test — Developer Guide §13
==============================================================
Addresses: Reviewer 2 Comment 1 (destination-side verification, Algorithm 3).
Dataset: 10,000 random ledgers per tamper type.

Tests T0-T8 from Developer Guide:
  T0: benign ledger -> accept
  T1: inflate r_i -> reject (signature)
  T2: replay stale certificate -> reject (freshness)
  T3: drop a hop -> reject (continuity)
  T4: insert repeated node (loop) -> reject (simplicity)
  T5: lower eps_i -> reject (signature)
  T6: reorder entries -> reject (continuity)
  T7: valid but superseded by newer DEGRADED alarm -> reject (DB)
  T8: re-signed forgery by non-endpoint -> reject (signature)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import math
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bmr.verify import make_ed25519_signer, make_ledger_entry, verify_ledger

MASTER_SEED = 46
N_TRIALS = 10_000
rng = np.random.default_rng(MASTER_SEED)

# Set up a signer/verifier pair
sign_fn, verify_fn = make_ed25519_signer()

# A second key pair (forger)
sign_fn_forger, verify_fn_forger = make_ed25519_signer()

def make_valid_ledger(n_hops=4, now=10, R=0.1, Delta=1e-2, W=5, t=1):
    """Create a valid n_hops-hop ledger from node 0 to node n_hops."""
    entries = []
    for i in range(n_hops):
        link = (i, i + 1)
        r = R + 0.01 + float(rng.uniform(0, 0.1))  # rate above floor
        eps = 1e-7
        delta = 1e-6
        T = now - int(rng.integers(0, W))  # fresh
        entry = make_ledger_entry(link, r, eps, delta, T, sign_fn)
        entries.append(entry)
    return entries


def build_cert_db(ledger):
    """Build a consistent certificate database from a ledger."""
    db = {}
    for e in ledger:
        db[(e["link"], e["T"])] = (e["r"], e["eps"], e["delta"])
    return db


def run_tamper_test(test_name, tamper_fn, n_trials=N_TRIALS):
    """Run tamper test: generate ledgers, apply tamper, verify."""
    n_hops = 4
    now = 10
    R = 0.05
    Delta = 1e-2
    W = 8
    t = 1
    results = []
    for _ in range(n_trials):
        ledger = make_valid_ledger(n_hops=n_hops, now=now, R=R, Delta=Delta, W=W, t=t)
        cert_db = build_cert_db(ledger)
        ledger_tampered, cert_db_modified = tamper_fn(ledger, cert_db, now)
        accepted, reason, checks = verify_ledger(
            ledger_tampered, 0, n_hops, R, Delta, W, t, now, cert_db_modified, verify_fn
        )
        results.append(accepted)
    accept_rate = sum(results) / n_trials
    reject_rate = 1.0 - accept_rate
    return accept_rate, reject_rate


# T0: Benign
def tamper_T0(ledger, db, now):
    return list(ledger), dict(db)

# T1: Inflate r_i on one entry
def tamper_T1(ledger, db, now):
    tampered = list(ledger)
    idx = int(rng.integers(len(tampered)))
    e = dict(tampered[idx])
    e["r"] = e["r"] + 0.5   # inflate rate
    # signature stays the same (old message) => mismatch
    tampered[idx] = e
    return tampered, dict(db)

# T2: Replay stale certificate (set T too old)
def tamper_T2(ledger, db, now):
    tampered = list(ledger)
    idx = int(rng.integers(len(tampered)))
    e = dict(tampered[idx])
    old_T = e["T"]
    new_T = now - 100  # way too old
    msg_new = f"{e['link']}|{e['r']:.10f}|{e['eps']:.10e}|{e['delta']:.10e}|{new_T}".encode()
    e["T"] = new_T
    e["_msg"] = msg_new
    e["sig"] = sign_fn(msg_new)  # re-sign with new T (but old T was in DB)
    tampered[idx] = e
    return tampered, dict(db)

# T3: Drop a hop
def tamper_T3(ledger, db, now):
    if len(ledger) <= 2:
        return list(ledger), dict(db)
    tampered = [ledger[0], ledger[-1]]  # drop middle entries
    return tampered, dict(db)

# T4: Insert repeated node (loop) — add an entry that creates a repeated node
def tamper_T4(ledger, db, now):
    tampered = list(ledger)
    if len(tampered) < 2:
        return tampered, dict(db)
    # Insert a back-edge creating a loop: e.g. link (2,1) after link (1,2)
    e0 = tampered[0]
    link_loop = (1, 0)  # goes back to node 0 (already visited)
    r = 0.2; eps = 1e-7; delta = 1e-6
    T = now - 1
    entry_loop = make_ledger_entry(link_loop, r, eps, delta, T, sign_fn)
    tampered.insert(1, entry_loop)
    return tampered, dict(db)

# T5: Lower eps_i (change eps, keep old signature)
def tamper_T5(ledger, db, now):
    tampered = list(ledger)
    idx = int(rng.integers(len(tampered)))
    e = dict(tampered[idx])
    e["eps"] = e["eps"] / 1000  # lower eps
    # do NOT re-sign => old sig covers old eps
    tampered[idx] = e
    return tampered, dict(db)

# T6: Reorder entries
def tamper_T6(ledger, db, now):
    if len(ledger) <= 2:
        return list(ledger), dict(db)
    tampered = list(ledger)
    tampered[0], tampered[1] = tampered[1], tampered[0]  # swap first two
    return tampered, dict(db)

# T7: Valid but superseded by newer DEGRADED alarm
def tamper_T7(ledger, db, now):
    modified_db = dict(db)
    # Add DEGRADED alarm for the first link at a later epoch
    first_link = ledger[0]["link"]
    first_T = ledger[0]["T"]
    modified_db[(first_link, now + 1)] = (0.0, 0.0, 0.0)  # DEGRADED
    return list(ledger), modified_db

# T8: Re-signed forgery by non-endpoint (using forger key)
def tamper_T8(ledger, db, now):
    tampered = list(ledger)
    idx = int(rng.integers(len(tampered)))
    e = dict(tampered[idx])
    e["r"] = e["r"] + 0.1
    msg_new = f"{e['link']}|{e['r']:.10f}|{e['eps']:.10e}|{e['delta']:.10e}|{e['T']}".encode()
    e["_msg"] = msg_new
    e["sig"] = sign_fn_forger(msg_new)  # forger signs
    tampered[idx] = e
    return tampered, dict(db)


tampers = [
    ("T0: Benign",                     tamper_T0, "accept"),
    ("T1: Inflate r_i",                tamper_T1, "reject"),
    ("T2: Replay stale cert",          tamper_T2, "reject"),
    ("T3: Drop a hop",                 tamper_T3, "reject"),
    ("T4: Insert repeated node",       tamper_T4, "reject"),
    ("T5: Lower eps_i",                tamper_T5, "reject"),
    ("T6: Reorder entries",            tamper_T6, "reject"),
    ("T7: Superseded DEGRADED alarm",  tamper_T7, "reject"),
    ("T8: Re-signed forgery",          tamper_T8, "reject"),
]

print("=" * 70)
print("E10: Ledger Verification and Tamper Tests")
print(f"      {N_TRIALS} ledgers per test type")
print("=" * 70)
print(f"{'Test':>35} {'Expected':>10} {'Accept%':>10} {'Reject%':>10} {'PASS':>6}")
print("-" * 75)

test_results = []
t_start_all = time.time()
for test_name, tamper_fn, expected in tampers:
    t0 = time.perf_counter()
    accept_rate, reject_rate = run_tamper_test(test_name, tamper_fn)
    elapsed = time.perf_counter() - t0
    if expected == "accept":
        passed = accept_rate >= 0.99
    else:
        passed = reject_rate >= 0.99
    test_results.append({
        "test": test_name, "expected": expected,
        "accept_rate": accept_rate, "reject_rate": reject_rate,
        "passed": passed, "ms": elapsed * 1000
    })
    print(f"{test_name:>35} {expected:>10} {accept_rate*100:>9.2f}% {reject_rate*100:>9.2f}%  {'PASS' if passed else 'FAIL':>5}")

total_ms = (time.time() - t_start_all) * 1000
all_pass = all(r["passed"] for r in test_results)
print(f"\nTotal time: {total_ms:.0f}ms")
print(f"All tests passed: {all_pass}")

# --- Verification timing ---
print("\n--- Verification timing (per entry) ---")
timing_results = []
for n_hops in [2, 4, 8, 16]:
    times = []
    for _ in range(1000):
        ledger = make_valid_ledger(n_hops=n_hops, now=10)
        db = build_cert_db(ledger)
        t0 = time.perf_counter()
        verify_ledger(ledger, 0, n_hops, 0.05, 1e-2, 8, 1, 10, db, verify_fn)
        elapsed = time.perf_counter() - t0
        times.append(elapsed * 1e6)  # microseconds
    med = np.median(times)
    per_entry = med / n_hops
    timing_results.append({"n_hops": n_hops, "median_us": med, "per_entry_us": per_entry})
    print(f"  n_hops={n_hops:2d}: median={med:.1f}µs, per-entry={per_entry:.1f}µs")

# Figure
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Acceptance/rejection rates
tests = [r["test"].split(":")[0] for r in test_results]
accept_rates = [r["accept_rate"] for r in test_results]
reject_rates = [r["reject_rate"] for r in test_results]
x = np.arange(len(tests))
ax = axes[0]
bars = ax.bar(x - 0.2, accept_rates, 0.35, label="Accept rate", color="green", alpha=0.7)
bars2 = ax.bar(x + 0.2, reject_rates, 0.35, label="Reject rate", color="red", alpha=0.7)
ax.set_xticks(x)
ax.set_xticklabels(tests, rotation=30, ha='right', fontsize=8)
ax.set_ylabel("Rate")
ax.set_title("E10: Ledger Tamper Test Results\n(T0 should accept; T1-T8 should reject)")
ax.legend()
ax.axhline(1.0, color='black', linestyle=':', alpha=0.5)
ax.set_ylim([0, 1.1])
ax.grid(True, alpha=0.3)

# Timing
ax2 = axes[1]
nhops_t = [r["n_hops"] for r in timing_results]
per_entry_t = [r["per_entry_us"] for r in timing_results]
ax2.bar(range(len(nhops_t)), per_entry_t, color='steelblue', alpha=0.7, edgecolor='black')
ax2.set_xticks(range(len(nhops_t)))
ax2.set_xticklabels([str(n) for n in nhops_t])
ax2.set_xlabel("Number of hops")
ax2.set_ylabel("Verification time per entry (µs)")
ax2.set_title("E10: Verification cost per ledger entry")
ax2.grid(True, alpha=0.3)

plt.tight_layout()
outfig = os.path.join(os.path.dirname(__file__), "..", "figures", "E10_ledger_tamper.png")
plt.savefig(outfig, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved figure: {outfig}")

out = {
    "experiment": "E10",
    "description": "Ledger verification and tamper tests (Reviewer 2 Comment 1)",
    "n_trials_per_test": N_TRIALS,
    "test_results": test_results,
    "all_passed": all_pass,
    "timing_per_entry": timing_results,
}
outpath = os.path.join(os.path.dirname(__file__), "..", "data", "results", "e10_ledger.json")
with open(outpath, "w") as f:
    json.dump(out, f, indent=2)
print(f"Saved results: {outpath}")
