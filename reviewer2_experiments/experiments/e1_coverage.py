"""
E1: eta bounds coverage (Monte Carlo) — Developer Guide §3
==========================================================
Addresses: Reviewer 2 Comment 2 (derivation of Hoeffding bound, n=100 insufficient).
Dataset: Monte Carlo (T=200,000 trials), quantum correlators of honest source S=2.687.

Verifies that both eta_tight and eta_union are valid: empirical failure rate <= delta.
Expected: 0 failures for both bounds at every (n, delta) row.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import math
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bmr.bounds import eta_tight, eta_union

MASTER_SEED = 42
T_TRIALS = 200_000
S_true = 2.687
# CHSH: S = E11 + E12 + E21 - E22, correlator signs: +1, +1, +1, -1
# E_ab = <AB> for the Phi+ state measured at angles giving max CHSH
# For binomial model: P(same|s) = (1 + s * S_true/4) / 2
signs = [+1, +1, +1, -1]

rng = np.random.default_rng(MASTER_SEED)

def draw_S_hat(n, T=T_TRIALS):
    """Draw T samples of S_hat for n rounds per setting."""
    S = np.zeros(T)
    for s in signs:
        p = (1 + s * S_true / 4) / 2
        k = rng.binomial(n, p, size=T)
        S += s * (2 * k - n) / n
    return S

ns = [100, 235, 500, 1000, 2500, 10000]
deltas = [1e-1, 1e-2, 1e-3, 1e-6]

rows = []
print(f"E1: eta bounds coverage (T={T_TRIALS} trials per n)")
print(f"{'n':>8} {'delta':>10} {'fail_tight':>12} {'fail_union':>12} {'std_Shat':>10} {'eta_tight':>10} {'eta_union':>10}")
print("-" * 80)

for n in ns:
    Sh = draw_S_hat(n)
    std_Sh = float(Sh.std())
    for d in deltas:
        et = eta_tight(n, d)
        eu = eta_union(n, d)
        fail_t = float(np.mean(Sh < S_true - et))
        fail_u = float(np.mean(Sh < S_true - eu))
        row = dict(n=n, delta=d,
                   fail_tight=fail_t,
                   fail_union=fail_u,
                   std_Shat=std_Sh,
                   eta_tight=et,
                   eta_union=eu)
        rows.append(row)
        print(f"{n:>8d} {d:>10.0e} {fail_t:>12.6f} {fail_u:>12.6f} {std_Sh:>10.5f} {et:>10.5f} {eu:>10.5f}")

# Check all failure rates are zero (or at least <= delta with T trials)
# valid bound <=> empirical failure rate <= delta (the old criterion "== 0" was wrong: at delta=0.1 a rate of 0.002 is valid)
all_pass_tight = all(r["fail_tight"] <= r["delta"] for r in rows)
all_pass_union = all(r["fail_union"] <= r["delta"] for r in rows)
print(f"\nAll tight-bound empirical failure rates <= delta: {all_pass_tight}")
print(f"All union-bound empirical failure rates <= delta: {all_pass_union}")

# Verify golden values from the Developer Guide
print("\nGolden value checks:")
from bmr.bounds import h_pironio, n_min
print(f"  h_pironio(2.687) = {h_pironio(2.687):.5f}  (expected 0.47233)")
print(f"  h_pironio(2.418) = {h_pironio(2.418):.5f}  (expected 0.20614)")
print(f"  eta_union(100, 1.0) = {eta_union(100, 1.0):.4f}  (expected 0.8157)")
print(f"  n_min(2.687, 1e-6, 'tight') = {n_min(2.687, 1e-6, 'tight')}  (expected 235)")
print(f"  n_min(2.687, 1e-6, 'union') = {n_min(2.687, 1e-6, 'union')}  (expected 1078)")

# Figure F-A: empirical failure rate vs delta (log-log), with line y = delta
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, bound_key, bound_fn, title in [
    (axes[0], "fail_tight", eta_tight, "Tight bound (η_tight)"),
    (axes[1], "fail_union", eta_union, "Union bound (η_union)"),
]:
    for n in ns:
        subset = [r for r in rows if r["n"] == n]
        ds_plot = [r["delta"] for r in subset]
        fails = [r[bound_key] for r in subset]
        ax.loglog(ds_plot, fails, 'o-', label=f"n={n}", markersize=4)
    ds_ref = [1e-6, 1e-3, 1e-1]
    ax.loglog(ds_ref, ds_ref, 'k--', lw=2, label="y = δ")
    ax.set_xlabel("δ (failure probability)")
    ax.set_ylabel("Empirical failure rate")
    ax.set_title(f"Figure F-A: {title}\nEmpirical failure rate vs δ")
    ax.legend(fontsize=7, ncol=2)
    ax.set_ylim([1e-8, 1.0])
    ax.grid(True, which="both", alpha=0.3)

plt.tight_layout()
outfig = os.path.join(os.path.dirname(__file__), "..", "figures", "E1_coverage_bounds.png")
plt.savefig(outfig, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved figure: {outfig}")

# Save results
out = {
    "experiment": "E1",
    "description": "eta bounds coverage check (Reviewer 2 Comment 2)",
    "T_trials": T_TRIALS,
    "S_true": S_true,
    "master_seed": MASTER_SEED,
    "rows": rows,
    "all_pass_tight": all_pass_tight,
    "all_pass_union": all_pass_union,
    "golden_values": {
        "h_pironio_2687": h_pironio(2.687),
        "h_pironio_2418": h_pironio(2.418),
        "eta_union_100_1": eta_union(100, 1.0),
        "n_min_tight_2687_1e-6": n_min(2.687, 1e-6, "tight"),
        "n_min_union_2687_1e-6": n_min(2.687, 1e-6, "union"),
    }
}
outpath = os.path.join(os.path.dirname(__file__), "..", "data", "results", "e1_coverage.json")
with open(outpath, "w") as f:
    json.dump(out, f, indent=2)
print(f"Saved results: {outpath}")
