"""
X7: Baseline Comparison — Werner-Corrected Cross-Metric Evaluation
===================================================================
Addresses r1_v2.md §2 E3-a, E4-b, D8, D9.

E3-a: E3 checked the formula but did NOT recompute any baseline.
      Comment 4 unanswered until Q-PATH/Q-CAST variants re-run on corrected quantity.

E4-b: Baselines were "unspecified surrogates" evaluated on BMR's own objective.
      Circular evaluation: BMR "wins by construction" on min r_e metric.
      Fix: cross-metric evaluation — also score on each baseline's own objective.

D8: Werner fidelity F = (1+3v)/4, NOT (v+1)/2.
D9: Rename baselines to "Q-PATH-proxy" and "Q-CAST-proxy".
    Cross-metric table mandatory.

This experiment:
  1. Corrects Werner formula error (old: F=(v+1)/2; correct: F=(1+3v)/4)
  2. Implements Q-PATH-proxy (maximize fidelity = ∏F_e)
  3. Implements Q-CAST-proxy (maximize entanglement = min E_e)
  4. Cross-metric table: evaluate each routing strategy on all three metrics:
     (a) BMR objective: Σr_e (certified rate sum)
     (b) Q-PATH objective: ∏F_e (path fidelity product)
     (c) Q-CAST objective: min E_e (bottleneck entanglement)
  5. Quantify error from old ∏F_e computation (wrong formula)
"""

import math
import json
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sqrt2 = math.sqrt(2)

def h_pironio(S):
    if S <= 2.0: return 0.0
    S = min(S, 2*sqrt2 - 1e-10)
    inner = max(0.0, 2.0 - S**2/4.0)
    return 1.0 - math.log2(1.0 + math.sqrt(inner))

def tau_paper(n, delta=1e-6):
    return math.sqrt(32 * math.log(8/delta) / n)

def sigma_physical(n):
    return 2.0 / math.sqrt(n)

# ─── Werner fidelity formulas (D8) ────────────────────────────────────────────
def fidelity_werner_correct(v):
    """Correct Werner fidelity: F = (1+3v)/4  per r1_v2.md D8."""
    return (1 + 3*v) / 4

def fidelity_werner_wrong(v):
    """OLD (incorrect) formula used in paper: F = (v+1)/2."""
    return (v + 1) / 2

def chsh_to_visibility(S):
    """v = S / (2√2)  for Werner/depolarizing state."""
    return S / (2 * sqrt2)

def fidelity_path_correct(path_S_vals):
    """
    Multi-hop path fidelity using CORRECT formula (D8).
    Swapped Werner: v_P = ∏v_e, F_P = (1+3v_P)/4.
    """
    v_prod = 1.0
    for S in path_S_vals:
        v_prod *= chsh_to_visibility(S)
    return fidelity_werner_correct(v_prod)

def fidelity_path_wrong(path_S_vals):
    """OLD incorrect path fidelity ∏F_e with wrong per-link formula."""
    f_prod = 1.0
    for S in path_S_vals:
        v = chsh_to_visibility(S)
        f_prod *= fidelity_werner_wrong(v)
    return f_prod

def fidelity_path_old_product(path_S_vals):
    """Paper's original ∏F_e with paper's Bell-diagonal formula."""
    f_prod = 1.0
    for S in path_S_vals:
        v = chsh_to_visibility(S)
        # Paper claimed F = (v+1)/2 but this is dephased Bell-diagonal
        f_prod *= (v + 1) / 2
    return f_prod

# ─── Hardware classes (D3) ────────────────────────────────────────────────────
CLASSES = [
    {"name": "high",  "eps_c": 2e-7, "S_mean": 2.720, "prob": 0.3},
    {"name": "std",   "eps_c": 1e-6, "S_mean": 2.687, "prob": 0.5},
    {"name": "deg",   "eps_c": 5e-6, "S_mean": 2.640, "prob": 0.2},
]

def make_graph(n_nodes, density=0.4, seed=0, n_audit=10000):
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    sigma = sigma_physical(n_audit)
    t = tau_paper(n_audit)
    R = 0.10
    edges = {}
    for u in range(n_nodes):
        for v in range(u+1, n_nodes):
            if rng.random() < density:
                cls = rng.choices(CLASSES, weights=[c["prob"] for c in CLASSES])[0]
                S_hat = float(np_rng.normal(cls["S_mean"], sigma))
                S_lb  = max(0.0, S_hat - t)
                r_e   = h_pironio(S_lb)
                v_val = chsh_to_visibility(S_hat)
                F_correct = fidelity_werner_correct(v_val)
                E_e = S_hat / (2*sqrt2)  # entanglement proxy ∝ visibility
                eps_e = cls["eps_c"]
                edges[(u,v)] = edges[(v,u)] = {
                    "S_hat": S_hat, "S_lb": S_lb, "r_e": r_e,
                    "v_e": v_val, "F_correct": F_correct, "E_e": E_e,
                    "eps_e": eps_e, "cls": cls["name"]
                }
    return edges

# ─── Routing algorithms ───────────────────────────────────────────────────────
def bmr_route(n_nodes, edges, src, dst, Delta, R=0.10, H=6):
    """BMR: maximize Σr_e s.t. Σε_e ≤ Δ, r_e ≥ R, simple path ≤ H hops."""
    best = {"rate": -1.0, "path": None}
    stack = [(src, 0.0, 0.0, [src])]
    while stack:
        u, cum_eps, cum_r, path = stack.pop()
        if u == dst and cum_eps <= Delta:
            if cum_r > best["rate"]:
                best.update({"rate": cum_r, "path": path[:]})
            continue
        if len(path) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if ed["r_e"] < R: continue
            new_eps = cum_eps + ed["eps_e"]
            if new_eps <= Delta:
                stack.append((vv, new_eps, cum_r + ed["r_e"], path + [vv]))
    return best["path"]

def qpath_proxy_route(n_nodes, edges, src, dst, R=0.10, H=6):
    """
    Q-PATH-proxy (D9): maximize path fidelity = F_P = (1+3·∏v_e)/4.
    Uses CORRECT Werner formula (D8).
    """
    best = {"fidelity": -1.0, "path": None}
    stack = [(src, 1.0, [src])]  # (node, cum_v_product, path)
    while stack:
        u, v_prod, path = stack.pop()
        if u == dst:
            F = fidelity_werner_correct(v_prod)
            if F > best["fidelity"]:
                best.update({"fidelity": F, "path": path[:]})
            continue
        if len(path) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if ed["r_e"] < R: continue  # minimum rate floor
            stack.append((vv, v_prod * ed["v_e"], path + [vv]))
    return best["path"]

def qcast_proxy_route(n_nodes, edges, src, dst, R=0.10, H=6):
    """
    Q-CAST-proxy (D9): maximize min entanglement E_e along path (bottleneck).
    """
    best = {"ent": -1.0, "path": None}
    stack = [(src, float('inf'), [src])]  # (node, min_E, path)
    while stack:
        u, min_E, path = stack.pop()
        if u == dst:
            if min_E > best["ent"]:
                best.update({"ent": min_E, "path": path[:]})
            continue
        if len(path) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if ed["r_e"] < R: continue
            stack.append((vv, min(min_E, ed["E_e"]), path + [vv]))
    return best["path"]

def shortest_path_fidelity(n_nodes, edges, src, dst, H=6):
    """SPF: shortest path (fewest hops), tie-break by fidelity."""
    from collections import deque
    queue = deque([(src, [src])])
    best = {"hops": float('inf'), "fidelity": -1.0, "path": None}
    while queue:
        u, path = queue.popleft()
        if u == dst:
            F = fidelity_path_correct([edges.get((path[i], path[i+1]), edges.get((path[i+1], path[i]), {})).get("S_hat", 2.0) for i in range(len(path)-1)])
            hops = len(path) - 1
            if hops < best["hops"] or (hops == best["hops"] and F > best["fidelity"]):
                best.update({"hops": hops, "fidelity": F, "path": path[:]})
            continue
        if len(path) > H: continue
        for (uu, vv), ed in edges.items():
            if uu == u and vv not in path:
                queue.append((vv, path + [vv]))
    return best["path"]

def evaluate_path(path, edges):
    """Compute all metrics for a path."""
    if path is None or len(path) < 2:
        return None
    S_vals = []
    eps_sum = 0.0
    r_sum = 0.0
    v_prod = 1.0
    min_E = float('inf')
    for i in range(len(path)-1):
        u, v = path[i], path[i+1]
        ed = edges.get((u, v), edges.get((v, u), {}))
        if not ed: return None
        S_vals.append(ed["S_hat"])
        eps_sum += ed["eps_e"]
        r_sum   += ed["r_e"]
        v_prod  *= ed["v_e"]
        min_E    = min(min_E, ed["E_e"])
    return {
        "hops":      len(path)-1,
        "r_sum":     r_sum,
        "F_correct": fidelity_werner_correct(v_prod),
        "F_wrong":   fidelity_path_wrong(S_vals),
        "min_E":     min_E,
        "eps_sum":   eps_sum,
    }

# ─── Main experiment ──────────────────────────────────────────────────────────
print("=" * 80)
print("X7: Baseline Comparison — Werner-Corrected Cross-Metric Table")
print("=" * 80)
print("\nD8: Using correct Werner F=(1+3v)/4")
print("D9: Q-PATH-proxy (maximize fidelity), Q-CAST-proxy (maximize min entanglement)")
print("E3-a fix: All baselines re-run on corrected fidelity formula")
print()

n_nodes = 20
Delta = 3e-5
R = 0.10
H = 6
N_GRAPHS = 100  # G=100 per D7

strategies = ["BMR", "Q-PATH-proxy", "Q-CAST-proxy", "SPF"]
metrics = ["r_sum", "F_correct", "min_E", "hops", "eps_sum"]

# Accumulate results across graphs
accum = {s: {m: [] for m in metrics} for s in strategies}
accum["Q-PATH-proxy"]["F_correct_vs_wrong"] = []  # extra: fidelity formula error

n_feasible = {s: 0 for s in strategies}

for seed in range(N_GRAPHS):
    g = make_graph(n_nodes, density=0.4, seed=seed*11, n_audit=10000)
    src, dst = 0, n_nodes - 1
    
    paths = {
        "BMR":         bmr_route(n_nodes, g, src, dst, Delta, R, H),
        "Q-PATH-proxy": qpath_proxy_route(n_nodes, g, src, dst, R, H),
        "Q-CAST-proxy": qcast_proxy_route(n_nodes, g, src, dst, R, H),
        "SPF":          shortest_path_fidelity(n_nodes, g, src, dst, H),
    }
    
    for s, path in paths.items():
        ev = evaluate_path(path, g)
        if ev is not None:
            n_feasible[s] += 1
            for m in metrics:
                if m in ev:
                    accum[s][m].append(ev[m])
            if s == "Q-PATH-proxy":
                accum[s]["F_correct_vs_wrong"].append(ev["F_correct"] - ev["F_wrong"])

# ─── Cross-metric table (D9 mandatory) ───────────────────────────────────────
print("Cross-Metric Table (G=100 random graphs, n=20 nodes)")
print("(Each strategy evaluated on ALL objectives, not just its own)")
print()
print(f"{'Strategy':<18}  {'Feasible':>10}  {'Mean Σr_e':>12}  {'Mean F_P':>12}  {'Mean min_E':>12}  {'Mean hops':>10}")
print("-" * 82)

cross_metric = []
for s in strategies:
    n_f = n_feasible[s]
    r_m  = np.mean(accum[s]["r_sum"]) if accum[s]["r_sum"] else 0
    f_m  = np.mean(accum[s]["F_correct"]) if accum[s]["F_correct"] else 0
    e_m  = np.mean(accum[s]["min_E"]) if accum[s]["min_E"] else 0
    h_m  = np.mean(accum[s]["hops"]) if accum[s]["hops"] else 0
    
    print(f"{s:<18}  {n_f:>10}  {r_m:>12.5f}  {f_m:>12.5f}  {e_m:>12.5f}  {h_m:>10.2f}")
    cross_metric.append({
        "strategy": s, "n_feasible": n_f,
        "r_sum_mean": round(r_m, 6), "F_correct_mean": round(f_m, 6),
        "min_E_mean": round(e_m, 6), "hops_mean": round(h_m, 4)
    })

# ─── Werner formula error quantification (E3-b) ──────────────────────────────
print("\n--- Werner formula error (E3-b) ---")
print(f"Old formula: F = (v+1)/2  (Bell-diagonal dephased, NOT Werner)")
print(f"Correct:     F = (1+3v)/4  (Werner depolarizing, per D8)")
print()
print("For specific v values:")
print(f"  {'v':>6}  {'F_wrong=(v+1)/2':>18}  {'F_correct=(1+3v)/4':>20}  {'Error (abs)':>14}")
print("-" * 68)
for v in [0.95, 0.90, 0.85, 0.80, 0.70, 0.60]:
    fw = fidelity_werner_wrong(v)
    fc = fidelity_werner_correct(v)
    print(f"  {v:>6.2f}  {fw:>18.5f}  {fc:>20.5f}  {fw-fc:>14.5f}")

print()
print("4-hop path fidelity comparison:")
print(f"  {'Scenario':>30}  {'F_wrong_product':>18}  {'F_correct_Werner':>18}  {'Ratio':>8}")
print("-" * 80)
for S_per_hop in [2.74, 2.687, 2.640, 2.55]:
    v_per_hop = chsh_to_visibility(S_per_hop)
    fw_4 = fidelity_werner_wrong(v_per_hop) ** 4
    fc_4 = fidelity_path_correct([S_per_hop]*4)
    print(f"  {'4 hops S='+str(S_per_hop):>30}  {fw_4:>18.5f}  {fc_4:>18.5f}  {fw_4/fc_4:>8.4f}")

werner_error_check = []
for v in [0.95, 0.90, 0.85, 0.80, 0.70, 0.60]:
    werner_error_check.append({
        "v": v, "F_wrong": round(fidelity_werner_wrong(v), 6),
        "F_correct": round(fidelity_werner_correct(v), 6)
    })
path_fidelity_check = []
for S_per_hop in [2.74, 2.687, 2.640, 2.55]:
    v_per_hop = chsh_to_visibility(S_per_hop)
    fw_4 = fidelity_werner_wrong(v_per_hop) ** 4
    fc_4 = fidelity_path_correct([S_per_hop]*4)
    path_fidelity_check.append({
        "S_per_hop": S_per_hop, "F_wrong_4hop": round(fw_4, 6),
        "F_correct_4hop_Werner": round(fc_4, 6)
    })

# ─── BMR vs Q-PATH: rate advantage and fidelity advantage ─────────────────────
print("\n--- BMR rate advantage over Q-PATH-proxy (cross-metric) ---")
# Compare pair by pair
n_pairs = 0
bmr_wins_rate = 0
qpath_wins_fid = 0
tied = 0
rate_diffs = []
fid_diffs  = []

for seed in range(N_GRAPHS):
    g = make_graph(n_nodes, density=0.4, seed=seed*11, n_audit=10000)
    src, dst = 0, n_nodes - 1
    
    bmr_path   = bmr_route(n_nodes, g, src, dst, Delta, R, H)
    qpath_path = qpath_proxy_route(n_nodes, g, src, dst, R, H)
    
    ev_bmr   = evaluate_path(bmr_path, g)
    ev_qpath = evaluate_path(qpath_path, g)
    
    if ev_bmr and ev_qpath:
        n_pairs += 1
        rd = ev_bmr["r_sum"] - ev_qpath["r_sum"]
        fd = ev_qpath["F_correct"] - ev_bmr["F_correct"]  # qpath wins on fidelity
        rate_diffs.append(rd)
        fid_diffs.append(fd)
        if rd > 0.001: bmr_wins_rate += 1
        if fd > 0.001: qpath_wins_fid += 1
        if abs(rd) <= 0.001 and abs(fd) <= 0.001: tied += 1

print(f"  Pairs compared: {n_pairs}")
print(f"  BMR wins on rate (Σr_e): {bmr_wins_rate}/{n_pairs} ({100*bmr_wins_rate/n_pairs:.1f}%)")
print(f"  Q-PATH wins on fidelity: {qpath_wins_fid}/{n_pairs} ({100*qpath_wins_fid/n_pairs:.1f}%)")
print(f"  Mean rate advantage of BMR: {np.mean(rate_diffs):.5f} ± {np.std(rate_diffs):.5f}")
print(f"  Mean fidelity advantage of Q-PATH: {np.mean(fid_diffs):.5f} ± {np.std(fid_diffs):.5f}")

comparison = {
    "n_pairs": n_pairs,
    "bmr_wins_rate": bmr_wins_rate, "qpath_wins_fid": qpath_wins_fid,
    "rate_diff_mean": round(float(np.mean(rate_diffs)), 6),
    "rate_diff_sd": round(float(np.std(rate_diffs)), 6),
    "fidelity_diff_mean": round(float(np.mean(fid_diffs)), 6)
}

# ─── Figures ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(17, 4.5))

# Plot 1: Cross-metric radar / grouped bar
strat_names = [r["strategy"] for r in cross_metric]
r_vals = [r["r_sum_mean"] for r in cross_metric]
F_vals = [r["F_correct_mean"] for r in cross_metric]
E_vals = [r["min_E_mean"] for r in cross_metric]

x = np.arange(len(strat_names))
w = 0.25
axes[0].bar(x - w, r_vals, w, label='Mean Σr_e (BMR obj)', color='steelblue', alpha=0.85)
axes[0].bar(x, F_vals, w, label='Mean F_P correct (Q-PATH obj)', color='orange', alpha=0.85)
axes[0].bar(x + w, E_vals, w, label='Mean min E (Q-CAST obj)', color='green', alpha=0.85)
axes[0].set_xticks(x)
axes[0].set_xticklabels(strat_names, rotation=15, fontsize=8)
axes[0].set_ylabel('Metric value')
axes[0].set_title('Cross-Metric Table\n(each strategy on all 3 objectives)')
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3, axis='y')

# Plot 2: Werner formula error
v_vals = [r["v"] for r in werner_error_check]
fw_vals = [r["F_wrong"] for r in werner_error_check]
fc_vals = [r["F_correct"] for r in werner_error_check]
axes[1].plot(v_vals, fw_vals, 'r--o', label='F_wrong = (v+1)/2', linewidth=2)
axes[1].plot(v_vals, fc_vals, 'b-s', label='F_correct = (1+3v)/4', linewidth=2)
axes[1].set_xlabel('Werner visibility v')
axes[1].set_ylabel('Fidelity F')
axes[1].set_title('Werner Fidelity Formula Error\n(D8: correct formula = (1+3v)/4)')
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

# Plot 3: Rate diff distribution BMR vs Q-PATH
axes[2].hist(rate_diffs, bins=20, color='steelblue', alpha=0.8, edgecolor='white')
axes[2].axvline(0, color='red', linestyle='--', linewidth=2, label='No difference')
axes[2].axvline(np.mean(rate_diffs), color='orange', linestyle='-', linewidth=2, 
                 label=f'Mean={np.mean(rate_diffs):.4f}')
axes[2].set_xlabel('BMR rate advantage over Q-PATH-proxy (Σr_e)')
axes[2].set_ylabel('Count (out of 100 graphs)')
axes[2].set_title('BMR vs Q-PATH-proxy\nRate Advantage Distribution')
axes[2].legend(fontsize=9)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_X7_baseline_cross_metric.png", dpi=150)
plt.close()
print("\nSaved: fig_X7_baseline_cross_metric.png")

out = {
    "description": "X7: Werner-corrected baseline comparison, cross-metric table (r1_v2 D8, D9, E3-a, E4-b)",
    "design": {
        "N_GRAPHS": N_GRAPHS, "n_nodes": n_nodes, "Delta": Delta, "R": R, "H": H,
        "D8": "Werner F=(1+3v)/4 (corrected)", "D9": "Q-PATH-proxy, Q-CAST-proxy (renamed)"
    },
    "cross_metric_table": cross_metric,
    "werner_formula_error": werner_error_check,
    "path_fidelity_comparison": path_fidelity_check,
    "bmr_vs_qpath_comparison": comparison
}
with open("REVISION_OSN_12OCT2026/expr_r1/X7_baseline_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("Saved: X7_baseline_results.json")
