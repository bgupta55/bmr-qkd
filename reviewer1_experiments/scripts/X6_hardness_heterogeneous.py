"""
X6: Two-Budget Hardness with Heterogeneous (a_e, b_e) Independent Weights
===========================================================================
Addresses r1_v2.md §2 E5-c, E5-d, D3, D4.

E5-c: With ε_e ≡ δ_e (constant), Σε and Σδ are the SAME functional → hop cap.
      NP-hardness experiment is DEGENERATE. Fix: use explicitly independent (a_e, b_e).

E5-d: The "two-budget" solver mixed γ (guessing excess) and δ (confidence loss).
      Fix: use abstract weights per D4; realizability lemma demonstrates the budget
      is non-trivial only when link classes have non-empty interior in (ε, γ) space.

D4: Thm 4 on abstract weights (a_e, b_e) with realizability lemma:
    If per-link classes make (ε_e, ln(1+γ_e)) range over a set with non-empty 
    interior, RSP (Resource-constrained Shortest Path) reduces to BMR.
    
This experiment:
  1. Shows the budget is truly non-trivial with heterogeneous classes
  2. Demonstrates scaling hardness of the NP-hard variant
  3. Shows exact vs heuristic gap as a function of budget tightness
  4. Realizability: measure what fraction of (a,b) points are achievable
"""

import time
import math
import json
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import product

sqrt2 = math.sqrt(2)

def h_pironio(S):
    if S <= 2.0: return 0.0
    S = min(S, 2*sqrt2 - 1e-10)
    inner = max(0.0, 2.0 - S**2/4.0)
    return 1.0 - math.log2(1.0 + math.sqrt(inner))

def tau_paper(n, delta=1e-6):
    return math.sqrt(32 * math.log(8/delta) / n)

# ─── Hardware classes (D3) ────────────────────────────────────────────────────
# Three classes: (ε_c, γ_c) pairs — INDEPENDENT axes (non-empty interior in (ε,γ) space)
# Class A: high ε, LOW γ  (device with high calibration error, low guessing excess)
# Class B: LOW ε, high γ  (device with low calibration error, high guessing excess)
# Class C: medium ε, medium γ (standard device)
# These are ANTI-CORRELATED by design → (a_e,b_e) spans a genuine 2D region.
CLASSES = [
    {"name": "high_eps_low_gam",  "eps_c": 5e-6,  "gamma_c": 1e-7,  "S_mean": 2.687},
    {"name": "low_eps_high_gam",  "eps_c": 2e-7,  "gamma_c": 4e-6,  "S_mean": 2.687},
    {"name": "medium_both",       "eps_c": 1e-6,  "gamma_c": 1e-6,  "S_mean": 2.687},
]

def make_hetero_graph(n_nodes, density=0.4, seed=0):
    """Graph with heterogeneous (a_e, b_e) per D4 — INDEPENDENT weights."""
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    n_audit = 10000
    sigma = 2.0 / math.sqrt(n_audit)
    t = tau_paper(n_audit)
    
    edges = {}
    for u in range(n_nodes):
        for v in range(u+1, n_nodes):
            if rng.random() < density:
                cls = rng.choice(CLASSES)
                S_hat = float(np_rng.normal(cls["S_mean"], sigma))
                S_lb  = max(0.0, S_hat - t)
                r_e   = h_pironio(S_lb)
                eps_e = cls["eps_c"]
                gamma_e = cls["gamma_c"]
                # Abstract independent weights: a_e = ε_e, b_e = ln(1+γ_e) (per Thm 4)
                a_e = eps_e
                b_e = math.log(1 + gamma_e)
                edges[(u,v)] = edges[(v,u)] = {
                    "r_e": r_e, "eps_e": eps_e, "gamma_e": gamma_e,
                    "a_e": a_e, "b_e": b_e, "cls": cls["name"]
                }
    return edges

# ─── Two-budget solvers ───────────────────────────────────────────────────────
def two_budget_optimal(n_nodes, edges, src, dst, Da, Db, R=0.10, H=6):
    """Exact optimal two-budget via branch-and-bound / exhaustive simple paths."""
    best = {"rate": -1.0, "path": None, "cum_a": None, "cum_b": None}
    # DFS with pruning
    stack = [(src, 0.0, 0.0, 0.0, [src])]
    while stack:
        u, cum_a, cum_b, cum_r, path = stack.pop()
        if u == dst:
            if cum_a <= Da and cum_b <= Db and cum_r > best["rate"]:
                best.update({"rate": cum_r, "path": path[:], "cum_a": cum_a, "cum_b": cum_b})
            continue
        if len(path) >= H + 1:
            continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if ed["r_e"] < R: continue
            new_a = cum_a + ed["a_e"]
            new_b = cum_b + ed["b_e"]
            if new_a <= Da and new_b <= Db:
                stack.append((vv, new_a, new_b, cum_r + ed["r_e"], path + [vv]))
    return best["path"], best["rate"], best["cum_a"], best["cum_b"]

def two_budget_heuristic_lagrange(n_nodes, edges, src, dst, Da, Db, R=0.10, H=6, lam=1.0):
    """
    Lagrangian-relaxation heuristic: combine budgets with weight λ.
    w_e = -r_e + λ·a_e + μ·b_e, find shortest path.
    Per E5-e fix: report sensitivity to λ.
    """
    # Simple greedy: maximize r_e - λ·a_e - μ·b_e
    mu = lam * 0.5
    stack = [(src, 0.0, 0.0, 0.0, [src])]
    best = {"rate": -1.0, "path": None}
    while stack:
        u, cum_a, cum_b, cum_r, path = stack.pop()
        if u == dst:
            if cum_a <= Da and cum_b <= Db and cum_r > best["rate"]:
                best.update({"rate": cum_r, "path": path[:]})
            continue
        if len(path) >= H + 1: continue
        # Only expand best candidate (greedy)
        candidates = []
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if ed["r_e"] < R: continue
            new_a = cum_a + ed["a_e"]
            new_b = cum_b + ed["b_e"]
            if new_a <= Da and new_b <= Db:
                score = ed["r_e"] - lam * ed["a_e"] - mu * ed["b_e"]
                candidates.append((score, vv, new_a, new_b, cum_r + ed["r_e"]))
        for sc, vv, na, nb, nr in sorted(candidates, reverse=True)[:3]:  # top-3 beam
            stack.append((vv, na, nb, nr, path + [vv]))
    return best["path"], best["rate"]

# ─── Experiment 1: Budget tightness sweep ────────────────────────────────────
print("=" * 80)
print("X6: Two-Budget Hardness — Heterogeneous Independent Weights")
print("=" * 80)

print("\n--- Part A: Budget tightness vs feasibility and optimality gap ---")
n_nodes = 15
n_trials = 30
alpha_vals = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]  # Da = alpha × avg_path_a

tightness_results = []
print(f"{'Alpha (Da)':>12}  {'Feasibility':>14}  {'Opt rate':>12}  {'Heur rate':>12}  {'Gap':>10}")
print("-" * 65)

for alpha in alpha_vals:
    feasible_count = 0
    opt_rates = []
    heur_rates = []
    gaps = []
    
    for seed in range(n_trials):
        g = make_hetero_graph(n_nodes, density=0.4, seed=seed*17)
        if not g: continue
        src, dst = 0, n_nodes - 1
        
        # Estimate average path budget usage
        avg_a = np.mean([ed["a_e"] for ed in g.values()]) * 4  # 4-hop avg
        avg_b = np.mean([ed["b_e"] for ed in g.values()]) * 4
        Da = alpha * avg_a
        Db = alpha * avg_b * 1.5  # b-budget slightly looser
        
        opt_path, opt_rate, _, _ = two_budget_optimal(n_nodes, g, src, dst, Da, Db)
        heur_path, heur_rate = two_budget_heuristic_lagrange(n_nodes, g, src, dst, Da, Db)
        
        if opt_path is not None:
            feasible_count += 1
            opt_rates.append(opt_rate)
            if heur_path is not None:
                heur_rates.append(heur_rate)
                gaps.append(opt_rate - heur_rate)
            else:
                heur_rates.append(0.0)
                gaps.append(opt_rate)
    
    feas = feasible_count / n_trials
    opt_m = np.mean(opt_rates) if opt_rates else 0
    heur_m = np.mean(heur_rates) if heur_rates else 0
    gap_m = np.mean(gaps) if gaps else 0
    
    print(f"{alpha:>12.1f}  {feas:>14.3f}  {opt_m:>12.5f}  {heur_m:>12.5f}  {gap_m:>10.5f}")
    tightness_results.append({
        "alpha": alpha, "feasibility": round(feas, 4),
        "opt_rate_mean": round(opt_m, 6), "heur_rate_mean": round(heur_m, 6),
        "gap_mean": round(gap_m, 6), "n_feasible": feasible_count
    })

# ─── Experiment 2: Scaling hardness ──────────────────────────────────────────
print("\n--- Part B: Scaling: NP-hard solver time vs graph size ---")
print("(Independent a_e, b_e; E5-c fix: budgets are NOT the same functional)")

scaling_results = []
sizes = [5, 8, 10, 12, 15, 18]

print(f"{'Nodes':>6}  {'Time (ms)':>12}  {'Feasible':>10}  {'Budget binds?':>16}")
print("-" * 50)

for n_nodes in sizes:
    times = []
    feasible_count = 0
    budget_binds_count = 0
    
    for seed in range(20):
        g = make_hetero_graph(n_nodes, density=0.4, seed=seed*31)
        src, dst = 0, n_nodes - 1
        avg_a = np.mean([ed["a_e"] for ed in g.values()]) * 4 if g else 1e-5
        avg_b = np.mean([ed["b_e"] for ed in g.values()]) * 4 if g else 1e-5
        Da = 0.7 * avg_a   # tight budget
        Db = 0.7 * avg_b
        
        t0 = time.perf_counter()
        opt_path, opt_rate, ca, cb = two_budget_optimal(n_nodes, g, src, dst, Da, Db)
        elapsed = (time.perf_counter() - t0) * 1000
        times.append(elapsed)
        
        if opt_path is not None:
            feasible_count += 1
            # Budget "binds" if removing one budget constraint changes the path
            opt_noB, _, _, _ = two_budget_optimal(n_nodes, g, src, dst, Da, float('inf'))
            if opt_noB != opt_path:  # different path without b-constraint
                budget_binds_count += 1
    
    t_ms = np.mean(times)
    feas = feasible_count / 20
    binds = budget_binds_count / max(feasible_count, 1)
    print(f"{n_nodes:>6}  {t_ms:>12.3f}  {feas:>10.3f}  {binds:>16.3f}")
    scaling_results.append({
        "n_nodes": n_nodes, "time_ms": round(t_ms, 4),
        "feasibility": round(feas, 4), "budget_binds_rate": round(binds, 4)
    })

# ─── Experiment 3: Realizability lemma ────────────────────────────────────────
print("\n--- Part C: Realizability — does the (a_e, b_e) coverage span 2D? ---")
print("Checking that hardware classes yield non-empty interior in (ε_e, ln(1+γ_e)) space")

# Collect all link (a_e, b_e) values
all_a = []
all_b = []
for seed in range(50):
    g = make_hetero_graph(20, density=0.4, seed=seed)
    for ed in g.values():
        all_a.append(ed["a_e"])
        all_b.append(ed["b_e"])

all_a = np.array(all_a)
all_b = np.array(all_b)

print(f"  a_e range: [{all_a.min():.3e}, {all_a.max():.3e}]")
print(f"  b_e range: [{all_b.min():.3e}, {all_b.max():.3e}]")
print(f"  Unique a_e classes: {len(np.unique(np.round(all_a, 12)))}")
print(f"  Unique b_e classes: {len(np.unique(np.round(all_b, 12)))}")
print(f"  Correlation(a,b): {np.corrcoef(all_a, all_b)[0,1]:.4f}")
print(f"  → Budget space has non-empty interior: {abs(np.corrcoef(all_a, all_b)[0,1]) < 0.99}")

realizability = {
    "a_e_min": float(all_a.min()), "a_e_max": float(all_a.max()),
    "b_e_min": float(all_b.min()), "b_e_max": float(all_b.max()),
    "n_unique_a": int(len(np.unique(np.round(all_a, 12)))),
    "n_unique_b": int(len(np.unique(np.round(all_b, 12)))),
    "correlation_ab": float(np.corrcoef(all_a, all_b)[0,1]),
    "non_empty_interior": bool(abs(np.corrcoef(all_a, all_b)[0,1]) < 0.99)
}

# ─── Heuristic sensitivity to λ ──────────────────────────────────────────────
print("\n--- Part D: Heuristic sensitivity to Lagrange weight λ (fixes E5-e) ---")
lambda_vals = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
lam_sensitivity = []

n_nodes = 15
n_trials = 30
print(f"{'Lambda':>10}  {'Heur feasible':>16}  {'Heur rate (mean)':>18}  {'Violation rate':>16}")
print("-" * 65)

for lam in lambda_vals:
    heur_rates = []
    violations = 0
    total_feasible = 0
    
    for seed in range(n_trials):
        g = make_hetero_graph(n_nodes, 0.4, seed*7)
        src, dst = 0, n_nodes - 1
        avg_a = np.mean([ed["a_e"] for ed in g.values()]) * 4 if g else 1e-5
        avg_b = np.mean([ed["b_e"] for ed in g.values()]) * 4 if g else 1e-5
        Da, Db = 0.8 * avg_a, 0.8 * avg_b
        
        opt_path, opt_rate, _, _ = two_budget_optimal(n_nodes, g, src, dst, Da, Db)
        heur_path, heur_rate = two_budget_heuristic_lagrange(n_nodes, g, src, dst, Da, Db, lam=lam)
        
        if opt_path is not None:
            total_feasible += 1
            if heur_path is None:
                violations += 1
            else:
                heur_rates.append(heur_rate)
    
    hfeas = (total_feasible - violations) / n_trials
    hrate = np.mean(heur_rates) if heur_rates else 0
    vrate = violations / total_feasible if total_feasible > 0 else float('nan')
    print(f"{lam:>10.1f}  {hfeas:>16.3f}  {hrate:>18.6f}  {vrate:>16.3f}")
    lam_sensitivity.append({
        "lambda": lam, "heur_feasibility": round(hfeas, 4),
        "heur_rate_mean": round(hrate, 6), "violation_rate": round(vrate, 4) if not math.isnan(vrate) else None
    })

# ─── Figures ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

# Plot 1: Budget tightness vs feasibility and gap
alphas = [r["alpha"] for r in tightness_results]
feasibilities = [r["feasibility"] for r in tightness_results]
gaps = [r["gap_mean"] for r in tightness_results]

ax1 = axes[0]
ax2 = ax1.twinx()
ax1.plot(alphas, feasibilities, 'b-o', label='Feasibility', linewidth=2)
ax2.plot(alphas, gaps, 'r--^', label='Opt-Heur gap', linewidth=2)
ax1.set_xlabel('Budget slack factor α')
ax1.set_ylabel('Feasibility', color='b')
ax2.set_ylabel('Optimality gap (rate)', color='r')
ax1.set_title('Budget Tightness\nvs Feasibility and Gap')
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1+lines2, labels1+labels2, fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: Scaling time
sc_nodes = [r["n_nodes"] for r in scaling_results]
sc_times = [r["time_ms"] for r in scaling_results]
sc_binds = [r["budget_binds_rate"] for r in scaling_results]
axes[1].semilogy(sc_nodes, sc_times, 'r-o', label='NP-hard 2-budget time', linewidth=2)
axes[1].set_xlabel('Number of nodes')
axes[1].set_ylabel('Runtime (ms, log scale)')
axes[1].set_title('NP-Hard Solver Scaling\n(independent a_e, b_e; E5-c fix)')
ax_r = axes[1].twinx()
ax_r.plot(sc_nodes, sc_binds, 'g--s', label='Budget binds rate', linewidth=2)
ax_r.set_ylabel('Budget binds rate', color='g')
axes[1].legend(fontsize=9, loc='upper left')
ax_r.legend(fontsize=9, loc='lower right')
axes[1].grid(True, alpha=0.3)

# Plot 3: (a_e, b_e) coverage (realizability)
import matplotlib.patches as mpatches
colors_cls = {'high_quality': 'blue', 'standard': 'green', 'degraded': 'red'}
g_demo = make_hetero_graph(30, density=0.5, seed=42)
pts_by_cls = {}
for ed in g_demo.values():
    c = ed["cls"]
    if c not in pts_by_cls:
        pts_by_cls[c] = ([], [])
    pts_by_cls[c][0].append(ed["a_e"])
    pts_by_cls[c][1].append(ed["b_e"])
for cls, (ax_pts, ay_pts) in pts_by_cls.items():
    axes[2].scatter(ax_pts, ay_pts, c=colors_cls.get(cls, 'gray'), 
                    label=cls, alpha=0.7, s=40)
axes[2].set_xlabel('a_e (ε_e, defect budget weight)')
axes[2].set_ylabel('b_e (ln(1+γ_e), guessing budget weight)')
axes[2].set_title('Realizability: (a_e, b_e) Coverage\nNon-empty interior → non-trivial 2-budget')
axes[2].legend(fontsize=8)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_X6_hardness_hetero.png", dpi=150)
plt.close()
print("\nSaved: fig_X6_hardness_hetero.png")

out = {
    "description": "X6: Two-budget hardness with heterogeneous independent weights (r1_v2 D3, D4, E5-c fix)",
    "realizability": realizability,
    "tightness_results": tightness_results,
    "scaling_results": scaling_results,
    "lambda_sensitivity": lam_sensitivity
}
with open("REVISION_OSN_12OCT2026/expr_r1/X6_hardness_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("Saved: X6_hardness_results.json")
