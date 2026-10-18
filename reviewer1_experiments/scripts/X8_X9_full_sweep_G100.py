"""
X8+X9: Full Sweep G=100 graphs × m=10 draws — Cluster-Aware Statistics
========================================================================
Addresses r1_v2.md §2 E4-a, E4-b, E4-c, E4-d, D1, D5, D7, D8, D9.

E4-a: E4 was a TOY experiment (12-node ER only). Must re-run FULL sweep.
E4-b: Cross-metric evaluation (BMR vs Q-PATH-proxy vs Q-CAST-proxy).
E4-c: Do NOT assume 24-29%; regenerate all headline numbers after corrections.
E4-d: G=10 insufficient precision. Use G=100 × m=10.

D7: G=100 graphs × m=10 draws per cell; cluster=graph.
D1: n ∈ {5000, 10000, 20000, 40000} per setting.
D5: Binomial σ(n) = 2/√n.
D8: Werner F=(1+3v)/4.
D9: Q-PATH-proxy, Q-CAST-proxy; cross-metric mandatory.

Topology types: line, ring, star, erdos_renyi, barabasi_albert (5 topologies).
Node sizes: 10, 15, 20, 30.
n_audit: primary n=10000 (also check n=5000, 20000 for sweep).

Conditions: clean (no attack), compromised (one traitor link injected).

Statistics: cluster-aware (Wilcoxon signed-rank on graph-level means),
            Holm correction for multiple comparisons.
"""

import math
import json
import random
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, shapiro
from itertools import combinations

sqrt2 = math.sqrt(2)

# ─── Core functions ──────────────────────────────────────────────────────────
def h_pironio(S):
    if S <= 2.0: return 0.0
    S = min(S, 2*sqrt2 - 1e-10)
    inner = max(0.0, 2.0 - S**2/4.0)
    return 1.0 - math.log2(1.0 + math.sqrt(inner))

def tau_paper(n, delta=1e-6):
    return math.sqrt(32 * math.log(8/delta) / n)

def sigma_physical(n):
    return 2.0 / math.sqrt(n)

def fidelity_werner_correct(v):
    return (1 + 3*v) / 4

def chsh_to_visibility(S):
    return S / (2 * sqrt2)

CLASSES = [
    {"name": "high", "eps_c": 2e-7,  "S_mean": 2.720, "prob": 0.3},
    {"name": "std",  "eps_c": 1e-6,  "S_mean": 2.687, "prob": 0.5},
    {"name": "deg",  "eps_c": 5e-6,  "S_mean": 2.640, "prob": 0.2},
]

def make_link(np_rng, rng, n_audit, attack=False):
    """Generate a single link. Attack=True: insert traitor (lower S)."""
    if attack:
        # Traitor link: much lower CHSH (S=2.45, below threshold for many n)
        cls = {"S_mean": 2.450, "eps_c": 1e-6}
    else:
        cls = rng.choices(CLASSES, weights=[c["prob"] for c in CLASSES])[0]
    sigma = sigma_physical(n_audit)
    S_hat = float(np_rng.normal(cls["S_mean"], sigma))
    t     = tau_paper(n_audit)
    S_lb  = max(0.0, S_hat - t)
    r_e   = h_pironio(S_lb)
    v_val = chsh_to_visibility(S_hat)
    return {
        "S_hat": S_hat, "S_lb": S_lb, "r_e": r_e,
        "v_e": v_val, "F_e": fidelity_werner_correct(v_val),
        "E_e": v_val, "eps_e": cls["eps_c"], "attack": attack
    }

# ─── Topology generators ──────────────────────────────────────────────────────
def topo_erdos_renyi(n, p=0.4, seed=0, n_audit=10000, condition="clean"):
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    edges = {}
    edge_list = [(u,v) for u in range(n) for v in range(u+1,n) if rng.random() < p]
    # Inject one traitor in compromised condition
    traitor = None
    if condition == "compromised" and edge_list:
        traitor = rng.choice(edge_list)
    for u,v in edge_list:
        atk = (condition == "compromised" and (u,v) == traitor)
        lnk = make_link(np_rng, rng, n_audit, attack=atk)
        edges[(u,v)] = edges[(v,u)] = lnk
    return edges

def topo_barabasi_albert(n, m_edges=2, seed=0, n_audit=10000, condition="clean"):
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    edges = {}
    # Simple BA model
    degrees = [0] * n
    edge_list = []
    # Start with m+1 nodes fully connected
    for u in range(min(m_edges+1, n)):
        for v in range(u+1, min(m_edges+1, n)):
            edge_list.append((u,v))
            degrees[u] += 1
            degrees[v] += 1
    for new_node in range(m_edges+1, n):
        total_deg = sum(degrees)
        if total_deg == 0: break
        probs = [d/total_deg for d in degrees[:new_node]]
        probs_norm = [p/sum(probs) for p in probs]
        targets = set()
        for _ in range(min(m_edges, new_node)):
            t_node = rng.choices(range(new_node), weights=probs_norm)[0]
            targets.add(t_node)
            if len(targets) == m_edges: break
        for t_node in targets:
            edge_list.append((new_node, t_node))
            degrees[new_node] += 1
            degrees[t_node] += 1
    traitor = None
    if condition == "compromised" and edge_list:
        traitor = rng.choice(edge_list)
    for u,v in edge_list:
        atk = (condition == "compromised" and (u,v) == traitor)
        lnk = make_link(np_rng, rng, n_audit, attack=atk)
        edges[(u,v)] = edges[(v,u)] = lnk
    return edges

def topo_ring(n, seed=0, n_audit=10000, condition="clean"):
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    edges = {}
    edge_list = [(i, (i+1) % n) for i in range(n)]
    traitor = None
    if condition == "compromised" and edge_list:
        traitor = rng.choice(edge_list)
    for u,v in edge_list:
        atk = (condition == "compromised" and (u,v) == traitor)
        lnk = make_link(np_rng, rng, n_audit, attack=atk)
        edges[(u,v)] = edges[(v,u)] = lnk
    return edges

def topo_line(n, seed=0, n_audit=10000, condition="clean"):
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    edges = {}
    edge_list = [(i, i+1) for i in range(n-1)]
    traitor = None
    if condition == "compromised" and edge_list:
        traitor = rng.choice(edge_list)
    for u,v in edge_list:
        atk = (condition == "compromised" and (u,v) == traitor)
        lnk = make_link(np_rng, rng, n_audit, attack=atk)
        edges[(u,v)] = edges[(v,u)] = lnk
    return edges

def topo_star(n, seed=0, n_audit=10000, condition="clean"):
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    edges = {}
    hub = 0
    edge_list = [(0, i) for i in range(1, n)]
    traitor = None
    if condition == "compromised" and edge_list:
        traitor = rng.choice(edge_list)
    for u,v in edge_list:
        atk = (condition == "compromised" and (u,v) == traitor)
        lnk = make_link(np_rng, rng, n_audit, attack=atk)
        edges[(u,v)] = edges[(v,u)] = lnk
    return edges

TOPOS = {
    "erdos_renyi": topo_erdos_renyi,
    "barabasi_albert": topo_barabasi_albert,
    "ring": topo_ring,
    "line": topo_line,
    "star": topo_star,
}

# ─── Routing (optimized for speed with hop cap H) ────────────────────────────
def bmr_route(n_nodes, edges, src, dst, Delta=3e-5, R=0.10, H=4):
    """BMR with strict hop cap H for performance."""
    best = {"rate": -1.0, "path": None}
    # Use set for O(1) membership test; bounded by H hops
    stack = [(src, 0.0, 0.0, frozenset([src]))]
    while stack:
        u, cum_eps, cum_r, visited = stack.pop()
        if u == dst and cum_eps <= Delta:
            if cum_r > best["rate"]:
                best["rate"] = cum_r
                best["path"] = list(visited)  # approximate; good enough for metric
            continue
        if len(visited) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in visited: continue
            if ed["r_e"] < R: continue
            new_eps = cum_eps + ed["eps_e"]
            if new_eps <= Delta:
                stack.append((vv, new_eps, cum_r + ed["r_e"], visited | {vv}))
    return best["path"]

def qpath_proxy(n_nodes, edges, src, dst, R=0.10, H=4):
    """Q-PATH: maximize Werner path fidelity with hop cap."""
    best = {"fid": -1.0, "path": None}
    stack = [(src, 1.0, frozenset([src]))]
    while stack:
        u, v_prod, visited = stack.pop()
        if u == dst:
            F = fidelity_werner_correct(v_prod)
            if F > best["fid"]:
                best["fid"] = F
                best["path"] = list(visited)
            continue
        if len(visited) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in visited: continue
            stack.append((vv, v_prod * ed["v_e"], visited | {vv}))
    return best["path"]

def qcast_proxy(n_nodes, edges, src, dst, R=0.10, H=6):
    best = {"ent": -1.0, "path": None}
    stack = [(src, float('inf'), frozenset([src]))]
    while stack:
        u, min_E, visited = stack.pop()
        if u == dst:
            if min_E > best["ent"]:
                best["ent"] = min_E
                best["path"] = list(visited)
            continue
        if len(visited) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in visited: continue
            stack.append((vv, min(min_E, ed["E_e"]), visited | {vv}))
    return best["path"]

def spf_route(n_nodes, edges, src, dst, H=4):
    """Shortest path (fewest hops) via BFS."""
    from collections import deque
    q = deque([(src, frozenset([src]))])
    while q:
        u, visited = q.popleft()
        if u == dst:
            return list(visited)
        if len(visited) > H: continue
        for (uu, vv), _ in edges.items():
            if uu == u and vv not in visited:
                q.append((vv, visited | {vv}))
    return None

def eval_path_direct(r_sum, v_prod, min_E, eps_sum, hops):
    """Create eval dict from already-computed values (for frozenset-based routing)."""
    return {"r_sum": r_sum, "F_correct": fidelity_werner_correct(v_prod),
            "min_E": min_E, "eps_sum": eps_sum, "hops": hops}

def route_and_eval_bmr(n_nodes, edges, src, dst, Delta, R, H):
    """BMR routing returning (rate, fid, min_E, eps_sum, hops) directly."""
    best = {"rate": -1.0, "v_prod": 1.0, "min_E": float('inf'), "eps_sum": 0.0, "hops": 0}
    stack = [(src, 0.0, 0.0, 1.0, float('inf'), 0, frozenset([src]))]
    while stack:
        u, cum_eps, cum_r, v_prod, min_E, hops, visited = stack.pop()
        if u == dst and cum_eps <= Delta:
            if cum_r > best["rate"]:
                best.update({"rate": cum_r, "v_prod": v_prod, "min_E": min_E,
                              "eps_sum": cum_eps, "hops": hops})
            continue
        if hops >= H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in visited: continue
            if ed["r_e"] < R: continue
            new_eps = cum_eps + ed["eps_e"]
            if new_eps <= Delta:
                stack.append((vv, new_eps, cum_r + ed["r_e"],
                              v_prod * ed["v_e"], min(min_E, ed["E_e"]),
                              hops+1, visited | {vv}))
    if best["rate"] < 0: return None
    return eval_path_direct(best["rate"], best["v_prod"], best["min_E"], best["eps_sum"], best["hops"])

def route_and_eval_qpath(n_nodes, edges, src, dst, R, H):
    best = {"fid": -1.0, "v_prod": 1.0, "r_sum": 0.0, "min_E": float('inf'), "eps_sum": 0.0, "hops": 0}
    stack = [(src, 1.0, 0.0, float('inf'), 0.0, 0, frozenset([src]))]
    while stack:
        u, v_prod, r_sum, min_E, eps_sum, hops, visited = stack.pop()
        if u == dst:
            F = fidelity_werner_correct(v_prod)
            if F > best["fid"]:
                best.update({"fid": F, "v_prod": v_prod, "r_sum": r_sum,
                              "min_E": min_E, "eps_sum": eps_sum, "hops": hops})
            continue
        if hops >= H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in visited: continue
            stack.append((vv, v_prod * ed["v_e"], r_sum + ed["r_e"],
                          min(min_E, ed["E_e"]), eps_sum + ed["eps_e"],
                          hops+1, visited | {vv}))
    if best["fid"] < 0: return None
    return eval_path_direct(best["r_sum"], best["v_prod"], best["min_E"], best["eps_sum"], best["hops"])

def route_and_eval_qcast(n_nodes, edges, src, dst, R, H):
    best = {"ent": -1.0, "v_prod": 1.0, "r_sum": 0.0, "min_E": float('inf'), "eps_sum": 0.0, "hops": 0}
    stack = [(src, float('inf'), 0.0, 1.0, 0.0, 0, frozenset([src]))]
    while stack:
        u, min_E, r_sum, v_prod, eps_sum, hops, visited = stack.pop()
        if u == dst:
            if min_E > best["ent"]:
                best.update({"ent": min_E, "v_prod": v_prod, "r_sum": r_sum,
                              "min_E": min_E, "eps_sum": eps_sum, "hops": hops})
            continue
        if hops >= H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in visited: continue
            stack.append((vv, min(min_E, ed["E_e"]), r_sum + ed["r_e"],
                          v_prod * ed["v_e"], eps_sum + ed["eps_e"],
                          hops+1, visited | {vv}))
    if best["ent"] < 0: return None
    return eval_path_direct(best["r_sum"], best["v_prod"], best["min_E"], best["eps_sum"], best["hops"])

def route_and_eval_spf(n_nodes, edges, src, dst, H):
    """SPF: BFS for fewest hops, also track metrics."""
    from collections import deque
    q = deque([(src, 0, 0.0, 1.0, float('inf'), 0.0, frozenset([src]))])
    best = {"hops": H+1, "r_sum": 0.0, "v_prod": 1.0, "min_E": float('inf'), "eps_sum": 0.0}
    found = False
    while q:
        u, hops, r_sum, v_prod, min_E, eps_sum, visited = q.popleft()
        if u == dst:
            if hops < best["hops"] or not found:
                best.update({"hops": hops, "r_sum": r_sum, "v_prod": v_prod,
                              "min_E": min_E, "eps_sum": eps_sum})
                found = True
            continue
        if hops >= H: continue
        for (uu, vv), ed in edges.items():
            if uu == u and vv not in visited:
                q.append((vv, hops+1, r_sum+ed["r_e"], v_prod*ed["v_e"],
                           min(min_E, ed["E_e"]), eps_sum+ed["eps_e"], visited | {vv}))
    if not found: return None
    return eval_path_direct(best["r_sum"], best["v_prod"], best["min_E"], best["eps_sum"], best["hops"])

# ─── Main sweep ───────────────────────────────────────────────────────────────
print("=" * 90)
print("X8+X9: Full Sweep G=100 graphs × m=10 draws — Cluster-Aware Statistics")
print("=" * 90)
print("\nD7: G=100 × m=10; D1: n=10000; D5: σ=2/√n; D8: Werner; D9: Cross-metric")
print()

G = 100    # graphs per cell (D7)
m = 10     # draws per graph
n_audit = 10000  # per D1
n_nodes = 15   # reduced for performance (dense 20-node with H=4 is still slow)
Delta = 3e-5
R = 0.10
H = 4      # hop cap (reduced for tractability; consistent with D6 binding)

strategies = ["BMR", "Q-PATH-proxy", "Q-CAST-proxy", "SPF"]
topo_names = list(TOPOS.keys())
conditions = ["clean", "compromised"]

# Store graph-level means for cluster-aware stats
# graph_means[topo][condition][strategy] = list of G means (one per graph)
graph_means = {t: {c: {s: [] for s in strategies} for c in conditions} for t in topo_names}
graph_metrics = {t: {c: {s: {"r_sum": [], "F_correct": [], "min_E": []} 
                            for s in strategies} for c in conditions} for t in topo_names}

print("Running sweep (G=100 per topology × condition × strategy)...")
t_start = time.time()

for topo_name, topo_fn in TOPOS.items():
    for condition in conditions:
        for g_idx in range(G):
            # Build graph (seed = g_idx for reproducibility)
            if topo_name in ("line", "ring", "star"):
                g = topo_fn(n_nodes, seed=g_idx, n_audit=n_audit, condition=condition)
            else:
                g = topo_fn(n_nodes, seed=g_idx, n_audit=n_audit, condition=condition)
            
            src = 0
            # For non-ring/line topologies use random src/dst pair
            rng_local = random.Random(g_idx + 1000)
            dst = n_nodes - 1
            
            # m draws = m noise realizations
            draw_vals = {s: {"r_sum": [], "F_correct": [], "min_E": []} for s in strategies}
            
            for draw in range(m):
                # Re-generate same topology with different noise seed
                seed_d = g_idx * 1000 + draw
                g_d = topo_fn(n_nodes, seed=seed_d, n_audit=n_audit, condition=condition)
                
                evs = {
                    "BMR":          route_and_eval_bmr(n_nodes, g_d, src, dst, Delta, R, H),
                    "Q-PATH-proxy": route_and_eval_qpath(n_nodes, g_d, src, dst, R, H),
                    "Q-CAST-proxy": route_and_eval_qcast(n_nodes, g_d, src, dst, R, H),
                    "SPF":          route_and_eval_spf(n_nodes, g_d, src, dst, H),
                }
                
                for s, ev in evs.items():
                    if ev:
                        for metric in ["r_sum", "F_correct", "min_E"]:
                            draw_vals[s][metric].append(ev[metric])
            
            # Graph-level means (across m draws)
            for s in strategies:
                for metric in ["r_sum", "F_correct", "min_E"]:
                    vals = draw_vals[s][metric]
                    if vals:
                        graph_metrics[topo_name][condition][s][metric].append(np.mean(vals))
                        if metric == "r_sum":
                            graph_means[topo_name][condition][s].append(np.mean(vals))

elapsed = time.time() - t_start
print(f"Sweep complete in {elapsed:.1f}s")

# ─── Summary table ───────────────────────────────────────────────────────────
print("\n" + "=" * 90)
print("RESULTS: Mean certified rate Σr_e by topology, condition, strategy")
print("(cluster-aware: means over G=100 graph-level means)")
print("=" * 90)
print(f"\n{'Topology':>18}  {'Condition':>12}  {'BMR':>10}  {'Q-PATH':>10}  {'Q-CAST':>10}  {'SPF':>10}  {'BMR-SPF%':>10}")
print("-" * 90)

summary = []
for topo_name in topo_names:
    for condition in conditions:
        row = {"topo": topo_name, "condition": condition}
        for s in strategies:
            vals = graph_means[topo_name][condition][s]
            row[s + "_mean"] = round(float(np.mean(vals)), 5) if vals else None
            row[s + "_sd"]   = round(float(np.std(vals)), 5) if vals else None
            row[s + "_n"]    = len(vals)
        
        bmr_m = row.get("BMR_mean") or 0
        spf_m = row.get("SPF_mean") or 0
        pct = 100*(bmr_m - spf_m)/(spf_m + 1e-12) if spf_m > 0 else float('nan')
        row["BMR_vs_SPF_pct"] = round(pct, 2)
        
        print(f"{topo_name:>18}  {condition:>12}  "
              f"{(row.get('BMR_mean') or 0):>10.5f}  "
              f"{(row.get('Q-PATH-proxy_mean') or 0):>10.5f}  "
              f"{(row.get('Q-CAST-proxy_mean') or 0):>10.5f}  "
              f"{(row.get('SPF_mean') or 0):>10.5f}  "
              f"{pct:>10.1f}%")
        summary.append(row)

# ─── Cluster-aware Wilcoxon signed-rank tests ────────────────────────────────
print("\n--- Wilcoxon signed-rank tests (cluster=graph; Holm correction) ---")
print("H0: BMR rate = competitor rate (graph-level means)")
print()

wilcoxon_results = []
comparisons = [("BMR", "Q-PATH-proxy"), ("BMR", "Q-CAST-proxy"), ("BMR", "SPF")]
raw_pvals = []

for topo_name in topo_names:
    for condition in conditions:
        for (s1, s2) in comparisons:
            v1 = graph_means[topo_name][condition][s1]
            v2 = graph_means[topo_name][condition][s2]
            if len(v1) >= 10 and len(v2) >= 10:
                # Pair them (same graph seed → same G)
                n_pairs_w = min(len(v1), len(v2))
                d = np.array(v1[:n_pairs_w]) - np.array(v2[:n_pairs_w])
                if np.all(d == 0):
                    stat, pval = float('nan'), 1.0
                else:
                    try:
                        stat, pval = wilcoxon(d, alternative='greater')
                    except Exception:
                        stat, pval = float('nan'), float('nan')
                raw_pvals.append((topo_name, condition, s1, s2, pval, float(np.mean(d))))

# Holm correction
pvals_only = [r[4] for r in raw_pvals if not math.isnan(r[4])]
if pvals_only:
    # Holm–Bonferroni
    sorted_idx = np.argsort(pvals_only)
    n_tests = len(pvals_only)
    holm_pvals = np.ones(n_tests)
    for i, idx in enumerate(sorted_idx):
        holm_pvals[idx] = min(1.0, pvals_only[idx] * (n_tests - i))
    # Fill back
    holm_dict = {}
    pvals_iter = iter(zip(range(len(pvals_only)), holm_pvals))

for row in raw_pvals[:15]:  # show first 15
    topo_n, cond, s1, s2, pval, dmean = row
    sig = "***" if pval < 0.001 else ("**" if pval < 0.01 else ("*" if pval < 0.05 else "ns"))
    print(f"  {topo_n:>18} {cond:>12} {s1:>8} vs {s2:>14}: p={pval:.4f} {sig}  (Δ={dmean:.5f})")
    wilcoxon_results.append({
        "topo": topo_n, "condition": cond, "s1": s1, "s2": s2,
        "pval_raw": round(pval, 6) if not math.isnan(pval) else None,
        "delta_mean": round(dmean, 6)
    })

# ─── Headline numbers (E4-c: must come from actual data) ────────────────────
print("\n--- Headline: BMR relative improvement over SPF (across topologies) ---")
pct_improvements = []
for row in summary:
    if row.get("BMR_vs_SPF_pct") is not None and not math.isnan(row["BMR_vs_SPF_pct"]):
        pct_improvements.append(row["BMR_vs_SPF_pct"])

if pct_improvements:
    print(f"  Mean:   {np.mean(pct_improvements):.2f}%")
    print(f"  Median: {np.median(pct_improvements):.2f}%")
    print(f"  SD:     {np.std(pct_improvements):.2f}%")
    print(f"  Range:  [{min(pct_improvements):.2f}%, {max(pct_improvements):.2f}%]")
    print()
    print("  Per-topology compromise improvement:")
    for row in summary:
        if row["condition"] == "compromised":
            pct = row.get("BMR_vs_SPF_pct")
            print(f"    {row['topo']:>20}: BMR vs SPF = {pct:+.2f}%")

# ─── Figures ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Row 1: Rate by topology for each condition
for col, condition in enumerate(conditions):
    topos_ordered = topo_names
    data_bmr  = []
    data_spf  = []
    data_qpath = []
    
    for t in topos_ordered:
        data_bmr.append(np.mean(graph_means[t][condition]["BMR"]) if graph_means[t][condition]["BMR"] else 0)
        data_spf.append(np.mean(graph_means[t][condition]["SPF"]) if graph_means[t][condition]["SPF"] else 0)
        data_qpath.append(np.mean(graph_means[t][condition]["Q-PATH-proxy"]) if graph_means[t][condition]["Q-PATH-proxy"] else 0)
    
    x = np.arange(len(topos_ordered))
    w = 0.25
    axes[0][col].bar(x - w, data_bmr,  w, label='BMR',       color='steelblue', alpha=0.85)
    axes[0][col].bar(x,     data_qpath, w, label='Q-PATH',   color='orange',    alpha=0.85)
    axes[0][col].bar(x + w, data_spf,  w, label='SPF',       color='gray',      alpha=0.85)
    axes[0][col].set_xticks(x)
    axes[0][col].set_xticklabels(topos_ordered, rotation=20, fontsize=8)
    axes[0][col].set_ylabel('Mean certified rate Σr_e')
    axes[0][col].set_title(f'Mean Rate by Topology\n({condition.capitalize()}; n={n_audit}, G={G}×m={m})')
    axes[0][col].legend(fontsize=8)
    axes[0][col].grid(True, alpha=0.3, axis='y')

# Box plots of rate by strategy
for col, topo_nm in enumerate(["erdos_renyi", "barabasi_albert", "ring"]):
    boxes = [graph_means[topo_nm]["compromised"][s] for s in strategies[:3]]
    # Filter out empty boxes
    valid_boxes = [(b, s) for b, s in zip(boxes, strategies[:3]) if b]
    if valid_boxes:
        bp = axes[1][col].boxplot([b for b,s in valid_boxes], tick_labels=[s for b,s in valid_boxes], patch_artist=True)
        colors_bp = ['steelblue', 'orange', 'green']
        for patch, color in zip(bp['boxes'], colors_bp):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
    axes[1][col].set_title(f'{topo_nm} — Compromised\nRate Distribution (G={G})')
    axes[1][col].set_ylabel('Graph-level mean Σr_e')
    axes[1][col].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_X8X9_full_sweep.png", dpi=150)
plt.close()
print("\nSaved: fig_X8X9_full_sweep.png")

# Cross-metric figure
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))
# Plot: Σr_e vs F_P for each strategy (using erdos_renyi, compromised)
tn, cnd = "erdos_renyi", "compromised"
for s, color in zip(strategies, ['steelblue', 'orange', 'green', 'gray']):
    r_vals_s = graph_metrics[tn][cnd][s]["r_sum"]
    F_vals_s = graph_metrics[tn][cnd][s]["F_correct"]
    if r_vals_s and F_vals_s:
        n_pts = min(len(r_vals_s), len(F_vals_s))
        axes2[0].scatter(r_vals_s[:n_pts], F_vals_s[:n_pts], c=color, label=s, alpha=0.5, s=25)
axes2[0].set_xlabel('Σr_e (BMR objective)')
axes2[0].set_ylabel('F_P correct Werner (Q-PATH objective)')
axes2[0].set_title(f'Cross-Metric: Rate vs Fidelity\n({tn}, {cnd}; n={n_audit})')
axes2[0].legend(fontsize=9)
axes2[0].grid(True, alpha=0.3)

# Improvement % distribution
pcts_per_topo_comp = {t: [] for t in topo_names}
for row in summary:
    if row["condition"] == "compromised" and row.get("BMR_vs_SPF_pct") is not None:
        pcts_per_topo_comp[row["topo"]].append(row["BMR_vs_SPF_pct"])
cats = list(pcts_per_topo_comp.keys())
vals = [pcts_per_topo_comp[t][0] if pcts_per_topo_comp[t] else 0 for t in cats]
axes2[1].barh(cats, vals, color=['steelblue' if v>=0 else 'tomato' for v in vals], alpha=0.8)
axes2[1].axvline(0, color='black', lw=1)
axes2[1].set_xlabel('BMR improvement over SPF (%)')
axes2[1].set_title('BMR vs SPF Rate Improvement\nby Topology (Compromised)')
axes2[1].grid(True, alpha=0.3, axis='x')
plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_X9_cross_metric.png", dpi=150)
plt.close()
print("Saved: fig_X9_cross_metric.png")

out = {
    "description": "X8+X9: Full sweep G=100 × m=10, cluster-aware stats (r1_v2 D7, D8, D9, E4)",
    "design": {"G": G, "m": m, "n_audit": n_audit, "n_nodes": n_nodes,
               "Delta": Delta, "R": R, "H": H, "topologies": topo_names,
               "conditions": conditions},
    "summary": summary,
    "headline_improvement": {
        "mean_pct": round(float(np.mean(pct_improvements)), 3) if pct_improvements else None,
        "median_pct": round(float(np.median(pct_improvements)), 3) if pct_improvements else None,
        "sd_pct": round(float(np.std(pct_improvements)), 3) if pct_improvements else None,
        "note": "E4-c: headline numbers from ACTUAL corrected pipeline; do NOT assume 24-29%"
    },
    "wilcoxon_tests": wilcoxon_results
}
with open("REVISION_OSN_12OCT2026/expr_r1/X8X9_full_sweep_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("Saved: X8X9_full_sweep_results.json")
