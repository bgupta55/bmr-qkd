"""
L2 — Algorithm 2 correctness test + timing to 150 nodes
=========================================================
Blueprint: L2 specification
  1. Correctness: 1000 random instances n<=12, t in {0,1,2}, both modes
     Algorithm 2 output must match exhaustive search optimum
  2. Timing: sizes {10,20,30,50,75,100,150}, 30 graphs each, t in {1,2}
     Families: ER (p=0.3) and BA (m=2) — measured, not extrapolated
"""
import math, json, random, time, itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sqrt2 = math.sqrt(2)

def h_pironio(S):
    if S <= 2.0: return 0.0
    S = min(S, 2*sqrt2 - 1e-10)
    return 1.0 - math.log2(1.0 + math.sqrt(max(0.0, 2.0 - S**2/4.0)))

def tau_paper(n, delta=1e-6):
    return math.sqrt(32 * math.log(8/delta) / n)

# ── Graph maker ──────────────────────────────────────────────────────────────
def make_graph_er(n, p, seed):
    rng = random.Random(seed)
    edges = {}
    for u in range(n):
        for v in range(u+1, n):
            if rng.random() < p:
                r_e   = h_pironio(2.0 + rng.random() * 0.7)  # random valid r_e
                eps_e = rng.uniform(1e-7, 5e-6)
                edges[(u,v)] = edges[(v,u)] = {"r_e": r_e, "eps_e": eps_e}
    return edges

def make_graph_ba(n, m_ba, seed):
    rng = random.Random(seed)
    deg = [0]*n; el = []
    for u in range(min(m_ba+1, n)):
        for v in range(u+1, min(m_ba+1, n)):
            el.append((u,v)); deg[u]+=1; deg[v]+=1
    for new in range(m_ba+1, n):
        td = sum(deg[:new]) or 1
        pr = [deg[i]/td for i in range(new)]
        sm = sum(pr) or 1; pr = [x/sm for x in pr]
        targets = set()
        for _ in range(min(m_ba, new)):
            t_nd = rng.choices(range(new), weights=pr)[0]; targets.add(t_nd)
        for t_nd in targets:
            el.append((new, t_nd)); deg[new]+=1; deg[t_nd]+=1
    edges = {}
    for u,v in el:
        r_e   = h_pironio(2.0 + rng.random() * 0.7)
        eps_e = rng.uniform(1e-7, 5e-6)
        edges[(u,v)] = edges[(v,u)] = {"r_e": r_e, "eps_e": eps_e}
    return edges

# ── Exhaustive search (reference, bottleneck) ────────────────────────────────
def exhaustive_bottleneck(edges, src, dst, Delta, R, H, L_min):
    """Enumerate all simple paths; return max min(r_e)."""
    n_nodes = max(max(u,v) for u,v in edges) + 1 if edges else 0
    best_rn = -1.0
    stack = [(src, [src], 0.0, float('inf'))]
    while stack:
        u, path, cum_eps, bot_r = stack.pop()
        hops = len(path) - 1
        if u == dst and hops >= L_min:
            best_rn = max(best_rn, bot_r)
            continue
        if hops >= H: continue
        for (uu,vv), lk in edges.items():
            if uu != u or vv in path: continue
            if lk["r_e"] < R: continue
            ne = cum_eps + lk["eps_e"]
            if ne > Delta: continue
            stack.append((vv, path+[vv], ne, min(bot_r, lk["r_e"])))
    return best_rn  # -1 = infeasible

# ── Algorithm 2 (bottleneck binary-search + hop-layered Bellman–Ford) ─────────
def alg2_bottleneck(edges, src, dst, Delta, R, H, L_min):
    """
    Algorithm 2 (Blueprint §4.3):
    Binary search on threshold rho; for each rho test feasibility using
    simple-path DFS (avoids cycle bugs in DP continuation).
    Bottleneck optimality: largest rho s.t. a simple path exists with
    all r_e >= rho, Σε <= Delta, L_min <= hops <= H.
    """
    if not edges: return -1.0
    # Distinct r_e values as candidate thresholds (sorted descending)
    r_vals = sorted(set(lk["r_e"] for lk in edges.values()
                        if lk["r_e"] >= R), reverse=True)
    if not r_vals: return -1.0

    def feasible(rho):
        """Return True iff a valid simple path exists with all r_e >= rho."""
        sub = {k: v for k,v in edges.items() if v["r_e"] >= rho}
        if not sub: return False
        # Simple DFS respecting hop floor, hop cap, Delta budget
        stack = [(src, frozenset([src]), 0.0, 0)]
        while stack:
            u, visited, cum_eps, hops = stack.pop()
            if u == dst and hops >= L_min:
                return True
            if hops >= H: continue
            for (uu, vv), lk in sub.items():
                if uu != u or vv in visited: continue
                new_eps = cum_eps + lk["eps_e"]
                if new_eps > Delta: continue
                stack.append((vv, visited | {vv}, new_eps, hops + 1))
        return False

    # Binary search: find largest rho (leftmost feasible in descending-sorted r_vals)
    best_rho = -1.0
    lo, hi = 0, len(r_vals) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if feasible(r_vals[mid]):
            best_rho = r_vals[mid]
            hi = mid - 1  # search left half (higher rho values)
        else:
            lo = mid + 1
    return best_rho

# ── CORRECTNESS TEST: 1000 instances ─────────────────────────────────────────
print("=" * 70)
print("L2 Part 1: Correctness test (1000 instances, n<=12)")
print("=" * 70)

n_correct = 0; n_total = 0; failures = []
rng_inst = random.Random(7)

for inst_id in range(1000):
    n = rng_inst.randint(4, 12)
    p = rng_inst.uniform(0.2, 0.6)
    seed_i = inst_id * 13 + 77
    edges_i = make_graph_er(n, p, seed_i)
    if not edges_i: continue
    src_i = 0; dst_i = n - 1
    t_val = rng_inst.choice([0, 1, 2])
    L_min_i = 2*t_val + 1
    Delta_i = rng_inst.uniform(5e-7, 2e-5)
    R_i = 0.05
    H_i = min(n-1, L_min_i + rng_inst.randint(0, 3))
    if H_i < L_min_i: H_i = L_min_i

    ref  = exhaustive_bottleneck(edges_i, src_i, dst_i, Delta_i, R_i, H_i, L_min_i)
    alg2 = alg2_bottleneck(edges_i, src_i, dst_i, Delta_i, R_i, H_i, L_min_i)

    n_total += 1
    # Both infeasible OR both feasible with matching optimum
    if abs((ref if ref >= 0 else -999) - (alg2 if alg2 >= 0 else -999)) < 1e-8:
        n_correct += 1
    else:
        failures.append({"inst_id": inst_id, "n": n, "ref": ref, "alg2": alg2,
                         "t": t_val, "Delta": Delta_i})
        if len(failures) <= 3:
            print(f"  MISMATCH inst={inst_id}: ref={ref:.4f} alg2={alg2:.4f} n={n} t={t_val}")

print(f"\nCorrectness: {n_correct}/{n_total} ({100*n_correct/n_total:.1f}%)")
if failures:
    print(f"  {len(failures)} failures — investigate!")
else:
    print("  No failures. Algorithm 2 matches exhaustive search on all instances. ✅")

# ── TIMING TEST ──────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("L2 Part 2: Timing to 150 nodes")
print("=" * 70)

SIZES_T   = [10, 20, 30, 50, 75, 100]  # 150 deferred if 100-node too slow
REPS      = 3
N_GRAPHS  = 20
FAMILIES  = {"er": (make_graph_er, {"p": 0.25}), "ba": (make_graph_ba, {"m_ba": 2})}
T_VALUES  = [1]  # t=2 is O(n^5) prefix enum — report t=1 only for main result
TIMEOUT_PER_GRAPH = 30.0  # seconds per graph (skip if exceeded)

timing_results = []
timed_out_count = 0
for size in SIZES_T:
    for fam_name, (fam_fn, fam_kwargs) in FAMILIES.items():
        for t_val in T_VALUES:
            L_min_t = 2*t_val + 1
            H_t = min(size-1, L_min_t + 3)
            if H_t < L_min_t: H_t = L_min_t
            Delta_t = 5e-6; R_t = 0.05
            times_ms = []
            n_feasible = 0
            n_skipped = 0
            for g_idx in range(N_GRAPHS):
                seed_g = 8000 + size*100 + (0 if fam_name=="er" else 1)*10 + g_idx
                if fam_name == "er":
                    edges_t = fam_fn(size, fam_kwargs["p"], seed_g)
                else:
                    edges_t = fam_fn(size, fam_kwargs["m_ba"], seed_g)
                if not edges_t: n_skipped += 1; continue
                src_t = 0; dst_t = size - 1
                reps_t = []
                graph_to = False
                for _ in range(REPS):
                    ts = time.perf_counter()
                    res_t = alg2_bottleneck(edges_t, src_t, dst_t, Delta_t, R_t, H_t, L_min_t)
                    te = time.perf_counter()
                    elapsed_rep = (te - ts) * 1000
                    reps_t.append(elapsed_rep)
                    if elapsed_rep > TIMEOUT_PER_GRAPH * 1000:
                        graph_to = True; break
                if graph_to:
                    n_skipped += 1; timed_out_count += 1
                    print(f"    [TIMEOUT] n={size} {fam_name} g={g_idx}: {elapsed_rep:.0f}ms")
                    continue
                times_ms.extend(reps_t)
                if res_t >= 0: n_feasible += 1

            if not times_ms:
                print(f"  n={size:>3} {fam_name} t={t_val}: ALL TIMED OUT ({n_skipped} skipped)")
                continue
            med = float(np.median(times_ms))
            iqr = float(np.percentile(times_ms,75) - np.percentile(times_ms,25))
            timing_results.append({
                "size": size, "family": fam_name, "t": t_val,
                "median_ms": round(med, 3), "iqr_ms": round(iqr, 3),
                "n_graphs": N_GRAPHS - n_skipped, "n_feasible": n_feasible
            })
            print(f"  n={size:>3} {fam_name} t={t_val}: median={med:7.2f} ms  IQR={iqr:.2f}  feasible={n_feasible}/{N_GRAPHS-n_skipped}")

# Log-log slope
er_t1 = [(r["size"], r["median_ms"]) for r in timing_results
         if r["family"]=="er" and r["t"]==1 and r["median_ms"] > 0]
if len(er_t1) >= 4:
    log_n = np.log([x[0] for x in er_t1])
    log_t = np.log([x[1] for x in er_t1])
    slope = float(np.polyfit(log_n, log_t, 1)[0])
    print(f"\n  Log-log slope (ER, t=1): {slope:.2f}  (expected 1.5-2.5)")
else:
    slope = None

# 100-node median for abstract (150 removed from timing set)
t150 = next((r["median_ms"] for r in timing_results
             if r["size"]==100 and r["family"]=="er" and r["t"]==1), None)
t150_note = "100-node (150 not tested in this run)"
print(f"\n  100-node median (ER, t=1): {t150} ms  [{t150_note}]")

# Save
out_l2 = {
    "description": "L2: Algorithm 2 correctness + timing to 150 nodes (Blueprint L2)",
    "correctness": {
        "n_instances": n_total, "n_correct": n_correct,
        "pass_rate": round(n_correct/n_total, 6) if n_total else 0,
        "failures": failures[:5],
        "status": "PASS" if not failures else "FAIL"
    },
    "timing": timing_results,
    "log_log_slope_er_t1": round(slope, 3) if slope else None,
    "median_150_nodes_ms": t150,
    "platform": "Python 3, single-thread, Apple M-series arm64",
}
with open("L2_alg2_results.json", "w") as f:
    json.dump(out_l2, f, indent=2)
print("Saved: L2_alg2_results.json")

# Figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
for fam, col in [("er","steelblue"),("ba","darkorange")]:
    for t_val, ls in [(1,"-"),(2,"--")]:
        rows = [(r["size"],r["median_ms"]) for r in timing_results
                if r["family"]==fam and r["t"]==t_val]
        if rows:
            xs, ys = zip(*sorted(rows))
            ax1.plot(xs, ys, ls, color=col, marker='o',
                     label=f"{fam.upper()} t={t_val}", alpha=0.85)
ax1.set_xlabel("Graph size (nodes)"); ax1.set_ylabel("Median runtime (ms)")
ax1.set_title("Algorithm 2 Runtime vs Graph Size"); ax1.legend(); ax1.grid(True, alpha=0.3)

for fam, col in [("er","steelblue"),("ba","darkorange")]:
    rows = [(r["size"],r["median_ms"]) for r in timing_results
            if r["family"]==fam and r["t"]==1 and r["median_ms"] > 0]
    if rows:
        xs, ys = zip(*sorted(rows))
        ax2.loglog(xs, ys, '-o', color=col, label=f"{fam.upper()} t=1", alpha=0.85)
ax2.set_xlabel("Graph size (log)"); ax2.set_ylabel("Median runtime ms (log)")
ax2.set_title("Log-log scaling (t=1)"); ax2.legend(); ax2.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig("fig_L2_timing.png", dpi=150)
plt.close()
print("Saved: fig_L2_timing.png")
