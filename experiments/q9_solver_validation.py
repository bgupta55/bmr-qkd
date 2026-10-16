"""
Q9 — Solver validation vs brute force; uncapped timing; E4b denominator fix
===========================================================================
(a) Validate bmr_route against exhaustive DFS for all n<=12 instances, t in {0,1,2,3}.
    Expect 0 mismatches.
(b) Time the full pipeline for n in {10..150}, t in {0,1,2,3,5}, 100 instances/cell, 
    1 s limit, uncapped. Report % timeouts.
(c) Report both conditional and unconditional E4b percentages.
(d) Cap analysis: how often does the cap change the result?
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import math, json, time, random
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from bmr.route import SolverTimeout
from bmr.certificate import build_certificate
from bmr.route import (bmr_route, brute_force_simple_paths, brute_force_max_rate,
                       min_weight_walk_hop_floor, is_simple)

SEED = 99

# ────────────────────────────────────────────────────────────
# Helper: make a small random graph with certificates
# ────────────────────────────────────────────────────────────
def make_instance(n, rng, p=0.4, S_mean=2.687, delta=1e-6, ell=10000, phi=0.99,
                  n_audit=500, Delta=0.5):
    """Return G_nx, s, d, Delta."""
    while True:
        G = nx.erdos_renyi_graph(n, p, seed=int(rng.integers(1e9)))
        if nx.is_connected(G): break
    for u, v in G.edges():
        S_hat = rng.normal(S_mean, 0.02)
        cert = build_certificate((u,v), S_hat, n_audit, delta, ell, phi, 0, 0)
        r_e = cert.r if cert else 0.0
        w_e = cert.w if cert else 2.0
        G[u][v]["r"] = r_e
        G[u][v]["w"] = w_e
        G[u][v]["T"] = 0
    return G, 0, n-1, Delta

# ────────────────────────────────────────────────────────────
# (a) Brute-force validation for n <= 12
# ────────────────────────────────────────────────────────────
print("="*70)
print("(a) Brute-force validation: bmr_route vs DFS enumerator (n<=12)")
print("="*70)
# Problem 1: maximise bottleneck rate over SIMPLE s-d paths with sum(w)<=Delta and
# hops>=m_min.  Ground truth = exhaustive DFS (brute_force_max_rate).  Checks per instance:
#   (i)  feasibility agreement, (ii) optimal bottleneck rate equal,
#   (iii) the path returned by bmr_route is simple, has >= m_min hops, weight <= Delta
#         and its own bottleneck equals the reported rate (independent re-verification).
# The Delta grid {0.5 (loose), tight} makes the weight constraint bind on a share of instances.
rng = np.random.default_rng(SEED)
mismatch_count = 0
total_validated = 0
truncated_excl = 0
n_feasible = 0
n_binding = 0           # instances where Delta (not just m_min) excludes some path
cells = {}
detail_rows = []
for n in range(4, 13):
    for t in [0, 1, 2, 3]:
        m_min = 2*t+1
        for k in range(1000):
            # alternate loose / tight Delta so the budget constraint is actually exercised
            G, s, d, Delta = make_instance(n, rng, p=0.5, n_audit=2500, Delta=0.5)
            if k % 2 == 1:
                ws = sorted(a["w"] for _, _, a in G.edges(data=True))
                Delta = float(sum(ws[:max(2, m_min)]) * 1.5)   # tight: ~1.5x the cheapest m_min links
            rate, bpath, bw, trunc = brute_force_max_rate(G, s, d, Delta, m_min=m_min)
            if trunc:
                truncated_excl += 1
                continue
            res = bmr_route(G, s, d, 0.0, Delta, 0, 1000, m_min=m_min)
            total_validated += 1
            c = cells.setdefault((n, t), dict(n=n, t=t, inst=0, feas=0, mism=0))
            c["inst"] += 1
            ok = True
            if (bpath is None) != (res is None):
                ok = False
            elif res is not None:
                n_feasible += 1; c["feas"] += 1
                rho, p, wt = res
                hops = len(p) - 1
                edge_r = [G[u][v]["r"] for u, v in zip(p[:-1], p[1:])]
                edge_w = sum(G[u][v]["w"] for u, v in zip(p[:-1], p[1:]))
                if abs(rho - rate) > 1e-12: ok = False
                if not is_simple(p) or hops < m_min or edge_w > Delta + 1e-12: ok = False
                if abs(min(edge_r) - rho) > 1e-12 or abs(edge_w - wt) > 1e-9: ok = False
            if not ok:
                mismatch_count += 1; c["mism"] += 1
                detail_rows.append(dict(n=n, t=t, Delta=Delta, bf=rate,
                                        res=None if res is None else res[0]))

print(f"Total validated: {total_validated}  (feasible: {n_feasible}; "
      f"excluded for DFS node-limit: {truncated_excl})")
print(f"Mismatches: {mismatch_count}")
assert mismatch_count == 0, f"MISMATCH FOUND! {detail_rows[:5]}"
print("VALIDATION PASSED — 0 mismatches\n")

# ────────────────────────────────────────────────────────────
# (b) Uncapped timing: n in {10..150}, t in {0,1,2,3,5}, 100 instances/cell
# ────────────────────────────────────────────────────────────
print("="*70)
print("(b) Uncapped timing: n={10..150}, t={0,1,2,3,5}, 100 instances/cell")
print("="*70)
sizes  = [10, 20, 30, 50, 75, 100, 150]
t_vals = [0, 1, 2, 3, 5]
N_INST = 100
TIMEOUT = 1.0   # seconds
timing_rows = []
rng2 = np.random.default_rng(SEED+1)

hdr = f"{'n':>5} {'t':>4} | {'med_ms':>8} {'p95_ms':>8} {'max_ms':>8} {'tmout%':>8} {'ok%':>6}"
print(hdr); print("-"*60)

for n in sizes:
    for t in t_vals:
        m_min = 2*t+1
        times_ms = []
        timeouts = 0
        for _ in range(N_INST):
            G, s, d, Delta = make_instance(n, rng2, p=0.3, n_audit=2500)
            t0 = time.perf_counter()
            try:
                res = bmr_route(G, s, d, 0.0, Delta, 0, 1000, m_min=m_min,
                                deadline=t0 + TIMEOUT)
            except SolverTimeout:
                timeouts += 1
            elapsed = (time.perf_counter() - t0)*1000  # ms (timed-out runs are censored at TIMEOUT)
            times_ms.append(elapsed)
        med   = float(np.median(times_ms))
        p95   = float(np.percentile(times_ms, 95))
        mx    = float(np.max(times_ms))
        tout_pct = 100.0 * timeouts / N_INST
        ok_pct  = 100.0 - tout_pct
        timing_rows.append(dict(n=n, t=t, m_min=m_min, median_ms=med, p95_ms=p95,
                                 max_ms=mx, timeout_pct=tout_pct, ok_pct=ok_pct))
        print(f"{n:>5} {t:>4} | {med:>8.3f} {p95:>8.3f} {mx:>8.3f} {tout_pct:>7.1f}% {ok_pct:>5.1f}%")
    print()

# ────────────────────────────────────────────────────────────
# (c/d) E4b denominator fix + cap analysis
# ────────────────────────────────────────────────────────────
print("="*70)
print("(c) E4b: walk vs simple path — conditional AND unconditional stats")
print("="*70)
from bmr.route import min_weight_walk_hop_floor

import itertools
def min_relay_cover(walk):
    """Smallest set of RELAY (internal) nodes whose compromise reveals every link of the
    walk (a link s-x or x-d is covered by its relay endpoint x; a link is covered if either
    endpoint is a compromised relay).  Exact: vertex cover over internal nodes, brute force
    (walk has <= ~15 distinct internal nodes).  Returns |cover|, or None if uncoverable."""
    edges = [(a, b) for a, b in zip(walk[:-1], walk[1:])]
    endpoints = {walk[0], walk[-1]}
    internal = sorted(set(walk[1:-1]))
    for size in range(0, len(internal) + 1):
        for comb in itertools.combinations(internal, size):
            cs = set(comb)
            if all((a in cs) or (b in cs) for a, b in edges):
                return size
    return None

rng3 = np.random.default_rng(SEED+2)
e4b_rows = []
N_E4B = 2200
for n in [10, 20, 30]:
    for t in [1, 2, 3]:
        K = 2*t+1
        found = ns = viol = viol_among_ns = 0
        for _ in range(N_E4B):
            G, s, d, Delta = make_instance(n, rng3, p=0.35, n_audit=2500)
            walk, _ = min_weight_walk_hop_floor(G, s, d, K)
            if walk is None:
                continue
            found += 1
            simple = is_simple(walk)
            cover = min_relay_cover(walk)
            v = (cover is not None) and (cover <= t)       # coalition of <= t relays reveals the key
            ns += (not simple); viol += v; viol_among_ns += (v and not simple)
            # sanity: a SIMPLE path with >= 2t+1 hops can never be violated (Corollary 1)
            assert simple is False or not v, "Corollary 1 violated by a simple walk!"
        row = dict(n=n, t=t, K=K, instances=N_E4B, found=found,
                   nonsimple_pct_of_all=100*ns/N_E4B,
                   nonsimple_pct_of_found=100*ns/max(found, 1),
                   violated_pct_of_all=100*viol/N_E4B,
                   violated_pct_of_found=100*viol/max(found, 1),
                   violated_pct_of_nonsimple=100*viol_among_ns/max(ns, 1))
        e4b_rows.append(row)
        print(f"n={n:3d} t={t} K={K}: found={found}/{N_E4B}  non-simple: {row['nonsimple_pct_of_all']:.1f}% of all "
              f"({row['nonsimple_pct_of_found']:.1f}% of found) | violated: {row['violated_pct_of_all']:.1f}% of all "
              f"({row['violated_pct_of_found']:.1f}% of found, {row['violated_pct_of_nonsimple']:.1f}% of non-simple)")
    print()

# (d) edge-weight / rate distribution used by (a)-(c)
rng4 = np.random.default_rng(SEED+3)
ww, rr = [], []
for _ in range(20000):
    c = build_certificate(("a", "b"), rng4.normal(2.687, 0.02), 2500, 1e-6, 10000, 0.99, 0, 0)
    ww.append(c.w); rr.append(c.r)
edge_dist = dict(S_hat="N(2.687, 0.02)", n=2500, delta=1e-6, ell=10000, phi=0.99,
                 w_mean=float(np.mean(ww)), w_sd=float(np.std(ww)), w_min=float(np.min(ww)), w_max=float(np.max(ww)),
                 r_mean=float(np.mean(rr)), r_sd=float(np.std(rr)), r_min=float(np.min(rr)), r_max=float(np.max(rr)))
print("(d) edge distribution:", edge_dist)

# ────────────────────────────────────────────────────────────
# Figure
# ────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13,5))
ax0, ax1 = axes
clrs_t = {0:"black",1:"blue",2:"green",3:"orange",5:"red"}
for t in t_vals:
    rows = [r for r in timing_rows if r["t"]==t]
    ns_t = [r["n"] for r in rows]
    meds = [r["median_ms"] for r in rows]
    p95s = [r["p95_ms"]    for r in rows]
    ax0.semilogy(ns_t, meds, "-o", color=clrs_t[t], label=f"t={t} (med)")
    ax0.semilogy(ns_t, p95s, "--", color=clrs_t[t], alpha=0.6)
ax0.set_xlabel("Network size n"); ax0.set_ylabel("Time (ms, log scale)")
ax0.set_title("Q9: Solver timing (uncapped)\nsolid=median, dashed=p95")
ax0.legend(fontsize=8); ax0.grid(True, which="both", alpha=0.3)

for t in [1,2,3]:
    rows = [r for r in e4b_rows if r["t"]==t]
    ns_e = [r["n"] for r in rows]
    vio  = [r["violated_pct_of_all"] for r in rows]
    ax1.plot(ns_e, vio, "-o", label=f"t={t} viol(all)")
ax1.set_xlabel("Network size n"); ax1.set_ylabel("Coalition violation rate (%)")
ax1.set_title("Q9(c): E4b walk-vs-path (unconditional)")
ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)
plt.tight_layout()
figpath = os.path.join(os.path.dirname(__file__),"..", "figures","Q9_solver_validation.png")
plt.savefig(figpath, dpi=150, bbox_inches='tight'); plt.close()
print(f"\nFigure saved: {figpath}")

out = {
    "experiment": "Q9",
    "description": "Solver validation vs brute force + uncapped timing + E4b denominators",
    "validation": {"total": total_validated, "feasible": n_feasible,
                   "excluded_dfs_limit": truncated_excl, "mismatches": mismatch_count,
                   "cells": [cells[k] for k in sorted(cells)],
                   "note": "ground truth = brute_force_max_rate (Problem 1). v5/earlier brute_force_simple_paths "
                           "had a bug (best_w never updated) - fixed; bmr_route itself unchanged."},
    "timing_rows": timing_rows,
    "e4b_rows": e4b_rows,
    "edge_distribution": edge_dist,
    "cap_statement": "bmr_route/oracle have NO path cap: Yen's enumeration is uncapped; a 1 s wall-clock deadline is "
                     "the only limit and is reported as 'timeout'; no 'no secure route' answer in (a) is cap-induced.",
}
outpath = os.path.join(os.path.dirname(__file__),"..","data","results","q9_solver_validation.json")
with open(outpath,"w") as f: json.dump(out, f, indent=2)
print(f"Results saved: {outpath}")
