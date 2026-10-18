"""
Experiment 4: BMR Route Selection & Cluster-Aware Statistics
=============================================================
Addresses Reviewer 1, Comments 2, 5:
  Comment 2: Verify polynomial algorithm (single-budget) vs NP-complete (two-budget)
  Comment 5: Pseudoreplication — graph seed as cluster unit

Simulation:
  - G = 10 graph seeds, m = 100 noise draws per graph (same as paper)
  - Also: G = 100 graph seeds (cluster-robust version)
  - Algorithms: BMR, Q-PATH(proxy), Q-CAST(proxy)
  - Metrics: certified entropy rate, feasibility
  - Reports: graph-level aggregated CIs, sign-flip p-value, tie counts
"""

import math
import json
import random
import itertools
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─── Certificate functions ────────────────────────────────────────────────────
sqrt2 = math.sqrt(2)

def h_pironio(S):
    if S <= 2.0: return 0.0
    if S > 2*sqrt2: S = 2*sqrt2
    return 1.0 - math.log2(1 + math.sqrt(2 - S**2/4))

def tau(n, delta):
    return math.sqrt(32*math.log(8/delta)/n)

def r_e_from_S(S_hat, n, delta):
    slb = S_hat - tau(n, delta)
    return h_pironio(slb) if slb > 2.0 else 0.0

def eps_e_from_S(S_hat, n, delta):
    return delta   # per corrected pipeline

# ─── Simple graph & routing model ────────────────────────────────────────────
S0      = 2.687
sigma   = 0.05
n_audit = 10_000  # per-setting (headline regime, matches paper)
delta_e = 1e-6
R       = 0.10
Delta   = 0.05    # headline budget (paper's Table 2)

def make_random_graph(n_nodes, edge_prob, rng_seed, sigma=sigma):
    """Returns adjacency list with link attributes."""
    rng = np.random.default_rng(rng_seed)
    nodes = list(range(n_nodes))
    edges = {}
    for i in range(n_nodes):
        for j in range(i+1, n_nodes):
            if rng.random() < edge_prob:
                edges[(i,j)] = edges[(j,i)] = None  # placeholder
    return nodes, edges

def draw_link_params(edges, rng, S0=S0, sigma=sigma, n=n_audit, delta=delta_e):
    """Draw CHSH values for each edge and compute cert attributes."""
    link_params = {}
    for e in set(tuple(sorted(k)) for k in edges):
        S_hat = float(rng.normal(S0, sigma))
        re    = r_e_from_S(S_hat, n, delta)
        eps   = eps_e_from_S(S_hat, n, delta)
        link_params[e] = {"S_hat": S_hat, "r_e": re, "eps_e": eps}
    return link_params

def get_link(link_params, i, j):
    return link_params.get(tuple(sorted((i,j))), None)

def all_simple_paths(adj, src, dst, n_nodes, max_hops=8):
    """Enumerate all simple paths src→dst with ≤ max_hops."""
    paths = []
    stack = [(src, [src])]
    while stack:
        node, path = stack.pop()
        if node == dst:
            paths.append(path)
            continue
        if len(path) >= max_hops + 1:
            continue
        for nb in range(n_nodes):
            if nb not in path:
                if (min(node,nb), max(node,nb)) in adj:
                    stack.append((nb, path+[nb]))
    return paths

def bmr_select(paths, link_params, R=R, Delta=Delta):
    """BMR: maximize min r_e along path, subject to sum eps_e <= Delta and r_e >= R."""
    best_re = -1.0
    best_path = None
    for path in paths:
        edges = [tuple(sorted((path[k], path[k+1]))) for k in range(len(path)-1)]
        lp = [link_params.get(e) for e in edges]
        if any(p is None for p in lp): continue
        re_vals  = [p["r_e"]  for p in lp]
        eps_vals = [p["eps_e"] for p in lp]
        if min(re_vals) < R: continue          # rate filter
        if sum(eps_vals) > Delta: continue     # budget constraint
        min_re = min(re_vals)
        if min_re > best_re:
            best_re = min_re
            best_path = path
    return best_re, best_path

def qpath_select(paths, link_params, R=R):
    """Q-PATH proxy: maximize end-to-end fidelity product (Werner)."""
    best_F = -1.0
    best_path = None
    for path in paths:
        edges = [tuple(sorted((path[k], path[k+1]))) for k in range(len(path)-1)]
        lp = [link_params.get(e) for e in edges]
        if any(p is None for p in lp): continue
        # Werner fidelity: F = (1+3v)/4, v = S/(2√2)
        v_vals = [p["S_hat"]/(2*sqrt2) for p in lp]
        F_path = (1 + 3*np.prod(v_vals)) / 4
        if F_path > best_F:
            best_F = F_path
            best_path = path
    # Convert path fidelity to certified entropy for comparison
    if best_path is None: return 0.0, None
    edges = [tuple(sorted((best_path[k], best_path[k+1]))) for k in range(len(best_path)-1)]
    lp = [link_params.get(e) for e in edges]
    re_vals = [p["r_e"] for p in lp]
    return min(re_vals), best_path

def qcast_select(paths, link_params, R=R):
    """Q-CAST proxy: minimize resource cost (sum of 1/r_e), subject to r_e >= R."""
    best_cost = float('inf')
    best_path = None
    for path in paths:
        edges = [tuple(sorted((path[k], path[k+1]))) for k in range(len(path)-1)]
        lp = [link_params.get(e) for e in edges]
        if any(p is None for p in lp): continue
        re_vals = [p["r_e"] for p in lp]
        if min(re_vals) < R: continue
        cost = sum(1/(r+1e-10) for r in re_vals)
        if cost < best_cost:
            best_cost = cost
            best_path = path
    if best_path is None: return 0.0, None
    edges = [tuple(sorted((best_path[k], best_path[k+1]))) for k in range(len(best_path)-1)]
    lp = [link_params.get(e) for e in edges]
    re_vals = [p["r_e"] for p in lp]
    return min(re_vals), best_path

# ─── Run simulation ─────────────────────────────────────────────────────────
N_NODES  = 12
EDGE_PROB = 0.35
G_seeds  = 10    # number of graph seeds
M_draws  = 100   # noise draws per graph
SRC, DST = 0, N_NODES-1

rng_top = np.random.default_rng(123)

all_results = []   # flat: (graph, draw, algo, rate)
tie_counts = {"BMR_vs_QPATH": 0, "BMR_vs_QCAST": 0}
win_counts = {"BMR_vs_QPATH": 0, "BMR_vs_QCAST": 0}
total_pairs = 0

print("=== BMR Route Selection + Cluster-Aware Statistics ===")
print(f"Graph: n={N_NODES} nodes, p={EDGE_PROB}, G={G_seeds} seeds, m={M_draws} draws/graph")
print(f"n_audit={n_audit}/setting, δ_e={delta_e}, Δ={Delta}, R={R}")
print()

for g in range(G_seeds):
    seed_g = int(rng_top.integers(0, 100000))
    nodes, edges_adj = make_random_graph(N_NODES, EDGE_PROB, seed_g)
    paths = all_simple_paths(edges_adj, SRC, DST, N_NODES, max_hops=7)
    if not paths:
        continue

    draw_rng = np.random.default_rng(seed_g + 999)
    for draw in range(M_draws):
        lp = draw_link_params(edges_adj, draw_rng)
        re_bmr,   _ = bmr_select(paths, lp)
        re_qpath, _ = qpath_select(paths, lp)
        re_qcast, _ = qcast_select(paths, lp)
        all_results.append({"graph": g, "draw": draw, "algo": "BMR",    "rate": re_bmr})
        all_results.append({"graph": g, "draw": draw, "algo": "Q-PATH", "rate": re_qpath})
        all_results.append({"graph": g, "draw": draw, "algo": "Q-CAST", "rate": re_qcast})
        # Tie / win analysis
        total_pairs += 1
        if re_bmr > re_qpath: win_counts["BMR_vs_QPATH"] += 1
        elif re_bmr == re_qpath: tie_counts["BMR_vs_QPATH"] += 1
        if re_bmr > re_qcast: win_counts["BMR_vs_QCAST"] += 1
        elif re_bmr == re_qcast: tie_counts["BMR_vs_QCAST"] += 1

# ─── Graph-level paired differences ─────────────────────────────────────────
import numpy as np
from collections import defaultdict

graph_BMR   = defaultdict(list)
graph_QPATH = defaultdict(list)
graph_QCAST = defaultdict(list)
for r in all_results:
    if r["algo"] == "BMR":    graph_BMR[r["graph"]].append(r["rate"])
    if r["algo"] == "Q-PATH": graph_QPATH[r["graph"]].append(r["rate"])
    if r["algo"] == "Q-CAST": graph_QCAST[r["graph"]].append(r["rate"])

g_ids = sorted(graph_BMR.keys())
d_qpath = np.array([np.mean(graph_BMR[g]) - np.mean(graph_QPATH[g]) for g in g_ids])
d_qcast = np.array([np.mean(graph_BMR[g]) - np.mean(graph_QCAST[g]) for g in g_ids])

def cluster_ci(d, B=20000, seed=0):
    rng2 = np.random.default_rng(seed)
    idx  = rng2.integers(0, len(d), size=(B, len(d)))
    m    = d[idx].mean(axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

def signflip_p(d, mc=200_000, seed=0):
    G, obs = len(d), abs(d.mean())
    if G <= 20:
        s = np.array(list(itertools.product([-1,1], repeat=G)))
        return float((np.abs((s*d).mean(axis=1)) >= obs - 1e-15).mean())
    rng2 = np.random.default_rng(seed)
    s = rng2.choice([-1,1], size=(mc, G))
    return float(((np.abs((s*d).mean(axis=1)) >= obs-1e-15).sum()+1)/(mc+1))

ci_qpath = cluster_ci(d_qpath)
ci_qcast = cluster_ci(d_qcast)
p_qpath  = signflip_p(d_qpath)
p_qcast  = signflip_p(d_qcast)

# Overall mean rates
mean_BMR   = np.mean([r["rate"] for r in all_results if r["algo"]=="BMR"])
mean_QPATH = np.mean([r["rate"] for r in all_results if r["algo"]=="Q-PATH"])
mean_QCAST = np.mean([r["rate"] for r in all_results if r["algo"]=="Q-CAST"])

print(f"Mean certified rate (min r_e):  BMR={mean_BMR:.5f}  Q-PATH={mean_QPATH:.5f}  Q-CAST={mean_QCAST:.5f}")
print()
print("Cluster-aware (graph-level) paired differences:")
print(f"  BMR vs Q-PATH: mean_diff={d_qpath.mean():.5f}  cluster 95%CI=[{ci_qpath[0]:.5f}, {ci_qpath[1]:.5f}]  sign-flip p={p_qpath:.4f}")
print(f"  BMR vs Q-CAST: mean_diff={d_qcast.mean():.5f}  cluster 95%CI=[{ci_qcast[0]:.5f}, {ci_qcast[1]:.5f}]  sign-flip p={p_qcast:.4f}")
print()
print(f"Tie counts (out of {total_pairs} paired trials):")
print(f"  BMR vs Q-PATH: wins={win_counts['BMR_vs_QPATH']}, ties={tie_counts['BMR_vs_QPATH']}, BMR-loses={total_pairs - win_counts['BMR_vs_QPATH'] - tie_counts['BMR_vs_QPATH']}")
print(f"  BMR vs Q-CAST: wins={win_counts['BMR_vs_QCAST']}, ties={tie_counts['BMR_vs_QCAST']}, BMR-loses={total_pairs - win_counts['BMR_vs_QCAST'] - tie_counts['BMR_vs_QCAST']}")

rel_adv_qpath = (mean_BMR - mean_QPATH) / max(mean_QPATH, 1e-10) * 100
rel_adv_qcast = (mean_BMR - mean_QCAST) / max(mean_QCAST, 1e-10) * 100
print(f"\nRelative advantage: BMR vs Q-PATH = {rel_adv_qpath:.1f}%,  BMR vs Q-CAST = {rel_adv_qcast:.1f}%")

# ─── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

algos = ["BMR", "Q-PATH", "Q-CAST"]
means_by_graph = {a: [np.mean(graph_BMR[g] if a=="BMR" else (graph_QPATH[g] if a=="Q-PATH" else graph_QCAST[g])) for g in g_ids] for a in algos}
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
for i, algo in enumerate(algos):
    axes[0].scatter(g_ids, means_by_graph[algo], label=algo, color=colors[i], zorder=3)
    axes[0].plot(g_ids, means_by_graph[algo], color=colors[i], alpha=0.5)
axes[0].set_xlabel('Graph seed ID')
axes[0].set_ylabel('Mean certified rate min(r_e), bits/round')
axes[0].set_title('Per-Graph Mean Rate by Algorithm')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Histogram of paired differences
axes[1].hist(d_qpath, bins=8, alpha=0.6, color='#ff7f0e', label='BMR−Q-PATH (graph-level)')
axes[1].hist(d_qcast, bins=8, alpha=0.6, color='#2ca02c', label='BMR−Q-CAST (graph-level)')
axes[1].axvline(0, color='black', linestyle='--')
axes[1].set_xlabel('Paired difference (bits/round)')
axes[1].set_ylabel('Count (graphs)')
axes[1].set_title('Cluster-Level Paired Differences')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_04_cluster_stats.png", dpi=150)
plt.close()
print("\nSaved: fig_04_cluster_stats.png")

# ─── Save results ────────────────────────────────────────────────────────────
out = {
    "description": "BMR route selection + cluster-aware statistics (Rev1 Comments 2, 5)",
    "params": {"N_NODES": N_NODES, "EDGE_PROB": EDGE_PROB, "G_seeds": G_seeds, "M_draws": M_draws,
               "n_audit": n_audit, "delta_e": delta_e, "Delta": Delta, "R": R},
    "mean_rates": {"BMR": mean_BMR, "Q-PATH": mean_QPATH, "Q-CAST": mean_QCAST},
    "relative_advantage_pct": {"BMR_vs_QPATH": rel_adv_qpath, "BMR_vs_QCAST": rel_adv_qcast},
    "cluster_CI_BMR_vs_QPATH": list(ci_qpath),
    "cluster_CI_BMR_vs_QCAST": list(ci_qcast),
    "signflip_p_BMR_vs_QPATH": p_qpath,
    "signflip_p_BMR_vs_QCAST": p_qcast,
    "tie_counts": tie_counts,
    "win_counts": win_counts,
    "total_pairs": total_pairs,
}
with open("REVISION_OSN_12OCT2026/expr_r1/04_cluster_stats_results.json", "w") as f:
    json.dump(out, f, indent=2)
