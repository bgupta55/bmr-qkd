"""
E2: Finite-sample table and audit grid — Developer Guide §4
============================================================
Addresses: Reviewer 2 Comment 2 (n=100 insufficient, new audit grid).
Dataset: Analytic computation (no random draws).

Generates Table T-A: for n in {100,250,500,1000,2500,10000},
delta in {1e-3,1e-6,1e-9}, S in {2.687, 2.62, 2.574, 2.418}:
  eta (tight and union), pass/fail of S - eta > 2, h(S-eta).
Also generates example certificates (Table T-B).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import math
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bmr.bounds import eta_tight, eta_union, h_pironio, n_min
from bmr.certificate import build_certificate

ns = [100, 250, 500, 1000, 2500, 10000]
deltas = [1e-3, 1e-6, 1e-9]
S_vals = [2.687, 2.62, 2.574, 2.418]
S_labels = ["S0=2.687", "S=2.62", "S1=2.574", "S1=2.418"]

print("=" * 90)
print("Table T-A: Finite-sample audit table (tight bound)")
print("=" * 90)
print(f"{'n':>8} {'delta':>8} {'eta_t':>8} {'eta_u':>8}  " +
      " | ".join([f"{'S='+str(s):>14}" for s in S_vals]))
print("-" * 100)

ta_rows = []
for n in ns:
    for d in deltas:
        et = eta_tight(n, d)
        eu = eta_union(n, d)
        row = {"n": n, "delta": d, "eta_tight": et, "eta_union": eu, "S_results": {}}
        cells = []
        for S in S_vals:
            S_lb = S - et
            if S_lb <= 2.0:
                h = 0.0
                status = "FAIL"
            else:
                h = h_pironio(S_lb)
                status = f"h={h:.4f}"
            row["S_results"][str(S)] = {"S_lb": S_lb, "h": h, "status": status}
            cells.append(f"{status:>14}")
        ta_rows.append(row)
        print(f"{n:>8d} {d:>8.0e} {et:>8.4f} {eu:>8.4f}  " + " | ".join(cells))

# Print expected values from the guide (tight bound, delta=1e-6)
print("\nExpected values from Developer Guide (tight, delta=1e-6):")
print(f"{'n':>8} {'eta':>8} {'S0=2.687':>14} {'S=2.62':>14} {'S1=2.418':>14}")
expected = [
    (100,    1.051, "FAIL",         "FAIL",         "FAIL"),
    (250,    0.665, "h=0.008",      "FAIL",         "FAIL"),
    (500,    0.470, "h=0.091",      "h=0.060",      "FAIL (alarm)"),
    (1000,   0.332, "h=0.165",      "h=0.127",      "h=0.033"),
    (2500,   0.210, "h=0.249",      "h=0.201",      "h=0.086"),
    (10000,  0.105, "h=0.342",      "h=0.280",      "h=0.141"),
]
for row in expected:
    print(f"{row[0]:>8d} {row[1]:>8.3f} {row[2]:>14} {row[3]:>14} {row[4]:>14}")

# --- Example certificates (Table T-B) ---
print("\n" + "=" * 90)
print("Table T-B: Example certificates (ell=10000, delta=1e-6, S=2.687)")
print("=" * 90)
print(f"{'n':>8} {'phi':>6} {'k':>8} {'r_e':>10} {'log2_eps':>12} {'log2_gamma':>12}")
print("-" * 60)
tb_rows = []
for n, phi in [(2500, 0.99), (2500, 0.90), (2500, 0.30), (10000, 0.99)]:
    c = build_certificate("link_A", 2.687, n, 1e-6, 10_000, phi, 0, 7)
    if c:
        print(f"{n:>8d} {phi:>6.2f} {c.k:>8d} {c.r:>10.4f} {c.log2_eps:>12.2f} {c.log2_gamma:>12.1f}")
        tb_rows.append({"n": n, "phi": phi, "k": c.k, "r": c.r,
                        "log2_eps": c.log2_eps, "log2_gamma": c.log2_gamma})
    else:
        print(f"{n:>8d} {phi:>6.2f}  DEGRADED")

# --- Design window ---
print("\nDesign window analysis (tight, delta=1e-6):")
print(f"  n_min(tight)  for S0=2.687: {n_min(2.687, 1e-6, 'tight')}")
print(f"  n_min(union)  for S0=2.687: {n_min(2.687, 1e-6, 'union')}")
print(f"  n_min(tight)  for S1=2.418: {n_min(2.418, 1e-6, 'tight')}")
print(f"  At n=500: S0=2.687 certifies (h={h_pironio(2.687-eta_tight(500,1e-6)):.4f}), S1=2.418 FAILS ({2.418-eta_tight(500,1e-6):.4f} vs 2.0)")
print(f"  At n=1000: S1=2.418 certifies (h={h_pironio(2.418-eta_tight(1000,1e-6)):.4f})")

# Figure: h(S-eta) vs n for each S, delta=1e-6
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ns_plot = list(range(100, 11000, 50))
d = 1e-6
for S, label in zip(S_vals, S_labels):
    hs = [h_pironio(S - eta_tight(n, d)) if S - eta_tight(n, d) > 2 else 0 for n in ns_plot]
    ax.plot(ns_plot, hs, label=label)
ax.set_xlabel("n (rounds per setting)")
ax.set_ylabel("h(S − η_tight) — certified min-entropy rate")
ax.set_title("E2 Table T-A: Certified rate vs n\n(tight bound, δ = 10⁻⁶)")
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_ylim([0, 0.55])
ax.axvline(235, color='gray', linestyle=':', label='n_min=235')
ax.axvline(1078, color='brown', linestyle=':', label='n_min union=1078')

ax2 = axes[1]
# Compare tight vs union at S0=2.687
S0 = 2.687
ht = [h_pironio(S0 - eta_tight(n, d)) if S0 - eta_tight(n, d) > 2 else 0 for n in ns_plot]
hu = [h_pironio(S0 - eta_union(n, d)) if S0 - eta_union(n, d) > 2 else 0 for n in ns_plot]
ax2.plot(ns_plot, ht, 'b-', label='η_tight (new bound)')
ax2.plot(ns_plot, hu, 'r--', label='η_union (old bound)')
ax2.set_xlabel("n (rounds per setting)")
ax2.set_ylabel("h(S₀ − η)")
ax2.set_title("E2: Tight vs Union bound\nS₀ = 2.687, δ = 10⁻⁶")
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
outfig = os.path.join(os.path.dirname(__file__), "..", "figures", "E2_finite_sample_table.png")
plt.savefig(outfig, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved figure: {outfig}")

# Save results
out = {
    "experiment": "E2",
    "description": "Finite-sample table + audit grid (Reviewer 2 Comment 2)",
    "table_TA": ta_rows,
    "table_TB": tb_rows,
    "n_min_tight_2687_1e6": n_min(2.687, 1e-6, "tight"),
    "n_min_union_2687_1e6": n_min(2.687, 1e-6, "union"),
    "n_min_tight_2418_1e6": n_min(2.418, 1e-6, "tight"),
}
outpath = os.path.join(os.path.dirname(__file__), "..", "data", "results", "e2_table.json")
with open(outpath, "w") as f:
    json.dump(out, f, indent=2)
print(f"Saved results: {outpath}")
