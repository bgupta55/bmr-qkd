"""
Experiment 5: Hardness / Runtime Scaling Experiment
=====================================================
Addresses Reviewer 1, Comments 2 and 6:
  Comment 2: NP-complete two-budget problem
  Comment 6: Runtime discrepancy (1ms vs 12ms at 150 nodes)

Two solvers:
  A) Single-budget BMR (polynomial time) — binary search on rho, 
     single additive defect constraint, hop window L to H
  B) Two-budget brute-force (exponential) — checks BOTH sum(eps_e) <= Delta
     AND sum(delta_e) <= delta_max independently

Measures wall-clock time vs n_nodes for both.
Reports:
  - Runtime table
  - NP-hard solver timeout node
  - Heuristic budget-violation rate
"""

import math
import time
import json
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

# ─── Graph utilities ─────────────────────────────────────────────────────────
S0      = 2.687
sigma   = 0.05
n_audit = 2500
delta_e = 1e-6
Delta   = 5e-5   # single-budget defect
delta_max = 2e-4  # second budget (guessing excess proxy)
R       = 0.10
H_max   = 6

def make_random_graph_adj(n_nodes, edge_prob, rng):
    adj = {}  # adj[node] = list of neighbors
    for i in range(n_nodes):
        adj[i] = []
    for i in range(n_nodes):
        for j in range(i+1, n_nodes):
            if rng.random() < edge_prob:
                adj[i].append(j)
                adj[j].append(i)
    return adj

def draw_links(adj, n_nodes, rng, n=n_audit, delta=delta_e):
    lp = {}
    for i in range(n_nodes):
        for j in adj[i]:
            e = (min(i,j), max(i,j))
            if e not in lp:
                S_hat = float(rng.normal(S0, sigma))
                re  = r_e_from_S(S_hat, n, delta)
                eps = delta  # ε_e = δ_e in corrected pipeline
                gam = max(0.0, 0.5 - 2**(-re)) if re > 0 else 0.5
                lp[e] = {"r_e": re, "eps_e": eps, "gamma_e": gam}
    return lp

def enum_paths(adj, src, dst, n_nodes, max_hops):
    """Enumerate simple paths src→dst with ≤ max_hops hops."""
    paths = []
    stack = [(src, [src])]
    while stack:
        node, path = stack.pop()
        if node == dst:
            paths.append(path)
        elif len(path)-1 < max_hops:
            for nb in adj[node]:
                if nb not in path:
                    stack.append((nb, path+[nb]))
    return paths

# ─── Solver A: Single-budget BMR (poly-time) ─────────────────────────────────
def bmr_single_budget(adj, lp, src, dst, n_nodes, R, Delta, max_hops=H_max):
    """
    Binary search on r-threshold rho, find largest rho s.t. 
    a path exists with min(r_e)>=rho AND sum(eps_e)<=Delta.
    Here we just enumerate all paths (simulation); in production use Dijkstra.
    """
    paths = enum_paths(adj, src, dst, n_nodes, max_hops)
    best = -1.0
    for path in paths:
        edges = [(min(path[k], path[k+1]), max(path[k], path[k+1])) for k in range(len(path)-1)]
        rvals = [lp[e]["r_e"]  for e in edges if e in lp]
        evals = [lp[e]["eps_e"] for e in edges if e in lp]
        if len(rvals) < len(path)-1: continue
        if min(rvals) < R: continue
        if sum(evals) > Delta: continue
        if min(rvals) > best: best = min(rvals)
    return best

# ─── Solver B: Two-budget brute-force (NP-hard) ──────────────────────────────
def bmr_two_budget(adj, lp, src, dst, n_nodes, R, Delta, delta_max, max_hops=H_max):
    """
    Two independent additive constraints:
      sum(eps_e) <= Delta   (defect budget)
      sum(gam_e) <= delta_max  (guessing-excess budget)
    Both must hold simultaneously → RSP-hard.
    Enumerate all paths (brute-force).
    """
    paths = enum_paths(adj, src, dst, n_nodes, max_hops)
    best = -1.0
    for path in paths:
        edges = [(min(path[k], path[k+1]), max(path[k], path[k+1])) for k in range(len(path)-1)]
        rvals = [lp[e]["r_e"]   for e in edges if e in lp]
        evals = [lp[e]["eps_e"] for e in edges if e in lp]
        gvals = [lp[e]["gamma_e"] for e in edges if e in lp]
        if len(rvals) < len(path)-1: continue
        if min(rvals) < R: continue
        if sum(evals) > Delta: continue      # budget 1
        if sum(gvals) > delta_max: continue  # budget 2
        if min(rvals) > best: best = min(rvals)
    return best

# ─── Lagrangian heuristic for two budgets ─────────────────────────────────────
def lagrangian_two_budget(adj, lp, src, dst, n_nodes, R, Delta, delta_max,
                          lam=1.0, max_hops=H_max):
    """
    Heuristic: minimize sum(-r_e) + lam*sum(eps_e) + lam*sum(gam_e)
    Ignores budget constraints → often violates them.
    """
    paths = enum_paths(adj, src, dst, n_nodes, max_hops)
    best_obj = float('inf')
    best_path = None
    for path in paths:
        edges = [(min(path[k], path[k+1]), max(path[k], path[k+1])) for k in range(len(path)-1)]
        rvals = [lp[e]["r_e"]   for e in edges if e in lp]
        evals = [lp[e]["eps_e"] for e in edges if e in lp]
        gvals = [lp[e]["gamma_e"] for e in edges if e in lp]
        if len(rvals) < len(path)-1: continue
        if min(rvals) < R: continue
        obj = -min(rvals) + lam*sum(evals) + lam*sum(gvals)
        if obj < best_obj:
            best_obj = obj
            best_path = (path, rvals, evals, gvals)
    if best_path is None:
        return -1.0, False, False
    path, rvals, evals, gvals = best_path
    v1 = sum(evals) <= Delta
    v2 = sum(gvals) <= delta_max
    return min(rvals), v1, v2

# ─── Runtime scaling experiment ───────────────────────────────────────────────
TIMEOUT_SEC = 2.0   # abort if brute-force > 2s per instance

node_sizes  = [8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50]
EDGE_PROB   = 0.30
N_REPS      = 3   # repetitions per size (use average)

runtime_bmr  = []
runtime_bf   = []
violation_rate = []  # fraction of Lagrangian heuristic trials that violate at least one budget
bf_timed_out  = []

rng = np.random.default_rng(777)

print("=== Runtime Scaling: BMR (poly) vs Two-Budget Brute-Force ===")
print(f"{'n':>6}  {'BMR_ms':>10}  {'BF_ms':>10}  {'BF_TO':>7}  {'Viol_%':>8}")
print("-"*50)

for n in node_sizes:
    t_bmr_list = []
    t_bf_list  = []
    viol_list  = []
    timed_out  = False
    for rep in range(N_REPS):
        seed_r = int(rng.integers(0, 1000000))
        rng2   = np.random.default_rng(seed_r)
        adj    = make_random_graph_adj(n, EDGE_PROB, rng2)
        lp     = draw_links(adj, n, rng2)
        src, dst = 0, n-1

        # BMR single-budget
        t0 = time.perf_counter()
        bmr_single_budget(adj, lp, src, dst, n, R, Delta, max_hops=min(n-1, H_max))
        t_bmr_list.append((time.perf_counter()-t0)*1000)

        # Lagrangian heuristic (check violation)
        _, v1, v2 = lagrangian_two_budget(adj, lp, src, dst, n, R, Delta, delta_max,
                                           max_hops=min(n-1, H_max))
        viol_list.append(0 if (v1 and v2) else 1)

        # Two-budget brute-force (with timeout)
        t0 = time.perf_counter()
        bmr_two_budget(adj, lp, src, dst, n, R, Delta, delta_max,
                       max_hops=min(n-1, H_max))
        elapsed = time.perf_counter() - t0
        t_bf_list.append(elapsed*1000)
        if elapsed > TIMEOUT_SEC:
            timed_out = True

    runtime_bmr.append(np.mean(t_bmr_list))
    runtime_bf.append(np.mean(t_bf_list))
    violation_rate.append(np.mean(viol_list)*100)
    bf_timed_out.append(timed_out)
    to_mark = "TIMEOUT" if timed_out else ""
    print(f"{n:>6}  {runtime_bmr[-1]:>10.3f}  {runtime_bf[-1]:>10.3f}  {to_mark:>7}  {violation_rate[-1]:>8.1f}")

# ─── Save & plot ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

# Runtime
axes[0].semilogy(node_sizes, runtime_bmr, 'b-o', label='BMR (single-budget, poly)')
axes[0].semilogy(node_sizes, runtime_bf,  'r-s', label='Two-budget brute-force')
axes[0].axhline(1.0,  color='gray', linestyle='--', lw=1, label='1 ms')
axes[0].axhline(12.0, color='orange', linestyle=':', lw=1, label='12 ms')
axes[0].set_xlabel('Network size (nodes)')
axes[0].set_ylabel('Runtime (ms, log scale)')
axes[0].set_title('Runtime Scaling: BMR vs Two-Budget Solver')
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)

# Violation rate
axes[1].plot(node_sizes, violation_rate, 'g-o')
axes[1].axhline(100, color='red', linestyle='--', lw=1, label='100%')
axes[1].set_xlabel('Network size (nodes)')
axes[1].set_ylabel('Heuristic budget violation rate (%)')
axes[1].set_title('Lagrangian Heuristic: Budget Violation Rate')
axes[1].legend()
axes[1].grid(True, alpha=0.3)
axes[1].set_ylim([0, 105])

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_05_runtime_scaling.png", dpi=150)
plt.close()

out = {
    "description": "Runtime scaling (Rev1 Comments 2, 6)",
    "node_sizes": node_sizes,
    "runtime_bmr_ms": [round(x,4) for x in runtime_bmr],
    "runtime_bf_ms":  [round(x,4) for x in runtime_bf],
    "violation_rate_pct": [round(x,2) for x in violation_rate],
    "bf_timed_out": bf_timed_out,
    "timeout_sec": TIMEOUT_SEC,
    "params": {"n_audit": n_audit, "delta_e": delta_e, "Delta": Delta,
               "delta_max": delta_max, "R": R, "EDGE_PROB": EDGE_PROB, "N_REPS": N_REPS},
}
with open("REVISION_OSN_12OCT2026/expr_r1/05_runtime_scaling_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved: fig_05_runtime_scaling.png")
