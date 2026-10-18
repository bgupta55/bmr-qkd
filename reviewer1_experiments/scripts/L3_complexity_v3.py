"""
L3 — Complexity v3: Partition-gadget + solver suite (lightweight rewrite)
=========================================================================
Blueprint: L3 / K9
Three families:
  (a) continuous independent (a_e,b_e) ~ logU
  (c) Partition-gadget series-parallel m=5..18 stages
Solvers: exhaustive_dfs, pareto_label_dp (label-setting, correct),
         pseudo_poly_dp, larac_style, fptas_style, naive_scalarisation
All solvers use timeout guards. No scipy.milp (too slow on arm64).
"""
import math, json, random, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SEED    = 20261012
TIMEOUT = 10.0   # per-solver timeout in seconds

# ─── Graph generators ────────────────────────────────────────────────────────
def make_independent(n, seed):
    rng = np.random.default_rng(seed)
    edges = {}
    for u in range(n):
        for v in range(u+1, n):
            if rng.random() < 0.4:
                a = float(np.exp(rng.uniform(math.log(2e-7), math.log(5e-6))))
                b = float(np.exp(rng.uniform(math.log(1e-7), math.log(4e-6))))
                edges[(u,v)] = edges[(v,u)] = {"a": a, "b": b, "fr": u, "to": v}
    return edges

def make_partition_gadget(m, seed, scale=1e-7):
    """m-stage series gadget. Each stage: two parallel arcs U_i, D_i."""
    rng = random.Random(seed)
    x = [rng.randint(1, 10000) for _ in range(m)]
    X = sum(x); c = 1
    A = B = (X / 2 + m * c) * scale
    edges = {}
    for i in range(m):
        a_u = (x[i] + c) * scale; b_u = c * scale
        a_d = c * scale;           b_d = (x[i] + c) * scale
        edges[(i, i+1, "U")] = {"a": a_u, "b": b_u, "fr": i, "to": i+1}
        edges[(i, i+1, "D")] = {"a": a_d, "b": b_d, "fr": i, "to": i+1}
    return edges, A, B

def get_budgets(edges, q=0.30):
    a_vals = [lk["a"] for lk in edges.values()]
    b_vals = [lk["b"] for lk in edges.values()]
    n_hops = 3
    A = float(np.percentile(a_vals, q*100)) * n_hops
    B = float(np.percentile(b_vals, q*100)) * n_hops
    return A, B

# ─── Solvers ──────────────────────────────────────────────────────────────────
def exhaustive_dfs(edges, src, dst, A, B, H=6, L_min=1):
    """Simple DFS with budget pruning and timeout."""
    t0 = time.perf_counter()
    arcs = [(lk["fr"], lk["to"], lk["a"], lk["b"]) for lk in edges.values()
            if "fr" in lk and lk["fr"] != lk["to"]]
    stack = [(src, frozenset([src]), 0.0, 0.0, 0)]
    while stack:
        if time.perf_counter() - t0 > TIMEOUT: return None, "timeout"
        u, vis, ca, cb, h = stack.pop()
        if u == dst and h >= L_min: return True, "ok"
        if h >= H: continue
        for fr, to, ae, be in arcs:
            if fr != u or to in vis: continue
            na, nb = ca+ae, cb+be
            if na <= A and nb <= B:
                stack.append((to, vis|{to}, na, nb, h+1))
    return False, "ok"

def pareto_label_dp(edges, src, dst, A, B, H=6, L_min=1):
    """
    Correct pareto label-setting (no path tracking — only (ca,cb,hops) labels).
    Labels at each node: set of non-dominated (ca, cb) pairs at each hop depth.
    """
    t0 = time.perf_counter()
    arcs = [(lk["fr"], lk["to"], lk["a"], lk["b"]) for lk in edges.values()
            if "fr" in lk and lk["fr"] != lk["to"]]
    nodes = set(fr for fr,to,ae,be in arcs) | set(to for fr,to,ae,be in arcs)
    # dp[h][v] = list of non-dominated (ca,cb)
    INF = float('inf')
    dp = [[dict() for _ in range(H+1)] for _ in nodes]   # use dict-of-nodes
    # dp[h][v] = list of (ca,cb) Pareto-front
    Labels = {(h,v): [] for h in range(H+1) for v in nodes}
    Labels[(0, src)] = [(0.0, 0.0)]

    for h in range(H):
        if time.perf_counter() - t0 > TIMEOUT: return None, "timeout"
        for fr, to, ae, be in arcs:
            for (ca, cb) in Labels[(h, fr)]:
                na, nb = ca+ae, cb+be
                if na > A or nb > B: continue
                new_lbl = (na, nb)
                existing = Labels[(h+1, to)]
                # dominance check: skip if any existing dominates new_lbl
                if any(ex[0] <= na and ex[1] <= nb for ex in existing): continue
                # prune dominated by new_lbl
                Labels[(h+1, to)] = [ex for ex in existing
                                     if not (na <= ex[0] and nb <= ex[1])]
                Labels[(h+1, to)].append(new_lbl)

    # Check if dst is reachable in [L_min, H] hops
    for h in range(L_min, H+1):
        if Labels[(h, dst)]: return True, "ok"
    return False, "ok"

def pseudo_poly_dp(edges, src, dst, A, B, H=6, L_min=1, K=200):
    """Pseudo-polynomial DP: discretise B into K bins."""
    t0 = time.perf_counter()
    arcs = [(lk["fr"], lk["to"], lk["a"], lk["b"]) for lk in edges.values()
            if "fr" in lk and lk["fr"] != lk["to"]]
    nodes = list(set(fr for fr,to,ae,be in arcs) | set(to for fr,to,ae,be in arcs))
    if src not in nodes or dst not in nodes: return False, "ok"
    b_step = B / K if B > 0 else 1.0
    INF = float('inf')
    # dp[h][v][k] = min a-cost to reach v in h hops with b-cost <= k*b_step
    dp = {v: [[INF]*(K+1) for _ in range(H+1)] for v in nodes}
    dp[src][0][0] = 0.0
    for h in range(H):
        if time.perf_counter() - t0 > TIMEOUT: return None, "timeout"
        for fr, to, ae, be in arcs:
            b_k = min(K, int(math.ceil(be / b_step))) if b_step > 0 else 0
            for k in range(K+1 - b_k):
                old = dp[fr][h][k]
                if old < INF:
                    na = old + ae
                    if na <= A:
                        nk = k + b_k
                        if dp[to][h+1][nk] > na:
                            dp[to][h+1][nk] = na
    ok = any(dp[dst][h][k] < INF for h in range(L_min,H+1) for k in range(K+1))
    return ok, "ok"

def larac_style(edges, src, dst, A, B, H=6, L_min=1):
    """LARAC-style: Lagrangian relaxation — try several lambda values."""
    import heapq
    t0 = time.perf_counter()
    arcs = [(lk["fr"], lk["to"], lk["a"], lk["b"]) for lk in edges.values()
            if "fr" in lk and lk["fr"] != lk["to"]]
    nodes = list(set(fr for fr,to,ae,be in arcs) | set(to for fr,to,ae,be in arcs))
    if src not in nodes or dst not in nodes: return False, "ok"

    def dijkstra_weighted(lam):
        dist = {(v, h): float('inf') for v in nodes for h in range(H+1)}
        dist[(src, 0)] = 0.0
        cum_a = {}; cum_b = {}
        cum_a[(src, 0)] = 0.0; cum_b[(src, 0)] = 0.0
        pq = [(0.0, src, 0, 0.0, 0.0)]
        while pq:
            if time.perf_counter() - t0 > TIMEOUT: return None, None, None
            w, u, h, ca, cb = heapq.heappop(pq)
            if w > dist[(u, h)]: continue
            if u == dst and h >= L_min and ca <= A and cb <= B:
                return True, ca, cb
            if h >= H: continue
            for fr, to, ae, be in arcs:
                if fr != u: continue
                nw = w + ae + lam*be; na = ca+ae; nb = cb+be
                nh = h+1
                if na <= A and nb <= B and nw < dist[(to, nh)]:
                    dist[(to, nh)] = nw
                    heapq.heappush(pq, (nw, to, nh, na, nb))
        return False, None, None

    for lam in [0.01, 0.1, 1.0, 10.0, 100.0]:
        if time.perf_counter() - t0 > TIMEOUT: return None, "timeout"
        found, ca, cb = dijkstra_weighted(lam)
        if found is None: return None, "timeout"
        if found: return True, "ok"
    return False, "ok"

def fptas_style(edges, src, dst, A, B, H=6, L_min=1, eps=0.10):
    """FPTAS-style: use pseudo_poly_dp with higher K."""
    arcs = [(lk["fr"], lk["to"], lk["a"], lk["b"]) for lk in edges.values()
            if "fr" in lk and lk["fr"] != lk["to"]]
    K_fp = min(500, max(50, int(len(arcs) / eps)))
    return pseudo_poly_dp(edges, src, dst, A, B, H, L_min, K=K_fp)

def naive_scalarisation(edges, src, dst, A, B, H=6, L_min=1, lam=1.0):
    """Naive: Dijkstra on a + lam*b, no feasibility guarantee."""
    import heapq
    arcs = [(lk["fr"], lk["to"], lk["a"], lk["b"]) for lk in edges.values()
            if "fr" in lk and lk["fr"] != lk["to"]]
    nodes = list(set(fr for fr,to,ae,be in arcs) | set(to for fr,to,ae,be in arcs))
    if src not in nodes or dst not in nodes: return False, False
    dist = {(v, h): float('inf') for v in nodes for h in range(H+1)}
    dist[(src, 0)] = 0.0
    pq = [(0.0, src, 0, 0.0, 0.0)]
    while pq:
        w, u, h, ca, cb = heapq.heappop(pq)
        if w > dist[(u,h)]: continue
        if u == dst and h >= L_min:
            return True, (ca <= A and cb <= B)
        if h >= H: continue
        for fr, to, ae, be in arcs:
            if fr != u: continue
            nw = w + ae + lam*be; nh = h+1
            if nw < dist[(to, nh)]:
                dist[(to, nh)] = nw
                heapq.heappush(pq, (nw, to, nh, ca+ae, cb+be))
    return False, False

# ─── Experiments ─────────────────────────────────────────────────────────────
print("="*70); print("L3 — Complexity v3 (Blueprint K9/L3)"); print("="*70)

SIZES_L3 = [8, 12, 18]
N_INST   = 15
all_results = []

SOLVERS = [
    ("exhaustive_dfs",    lambda e,s,d,A,B,H,L: exhaustive_dfs(e,s,d,A,B,H,L)),
    ("pareto_label_dp",   lambda e,s,d,A,B,H,L: pareto_label_dp(e,s,d,A,B,H,L)),
    ("pseudo_poly_dp",    lambda e,s,d,A,B,H,L: pseudo_poly_dp(e,s,d,A,B,H,L)),
    ("larac_style",       lambda e,s,d,A,B,H,L: larac_style(e,s,d,A,B,H,L)),
    ("fptas_style",       lambda e,s,d,A,B,H,L: fptas_style(e,s,d,A,B,H,L)),
]

print("\n--- Family (a): independent log-uniform ---")
for n in SIZES_L3:
    for inst_id in range(N_INST):
        seed_i = SEED + n*1000 + inst_id
        edges_i = make_independent(n, seed_i)
        if len(edges_i) < 4: continue
        A_i, B_i = get_budgets(edges_i, 0.30)
        H_i = min(n-1, 6); L_i = 1
        row = {"family": "independent", "n": n, "inst_id": inst_id}
        exact_res = None
        for sname, sfn in SOLVERS:
            t0 = time.perf_counter()
            try:
                res, stat = sfn(edges_i, 0, n-1, A_i, B_i, H_i, L_i)
            except Exception as ex:
                res, stat = None, str(ex)
            ms = (time.perf_counter()-t0)*1000
            row[sname] = {"ms": round(ms,3), "result": res, "status": stat}
            if sname == "exhaustive_dfs": exact_res = res
        # Naive (lambda swept)
        for lam in [0.1, 1.0, 10.0]:
            t0 = time.perf_counter()
            found, feas = naive_scalarisation(edges_i, 0, n-1, A_i, B_i, H_i, L_i, lam)
            ms = (time.perf_counter()-t0)*1000
            row[f"naive_lam{lam}"] = {"ms": round(ms,3), "found": found,
                                       "feasible": feas, "violation": found and not feas}
        row["exact_feasible"] = exact_res
        all_results.append(row)
    feas_count = sum(1 for r in all_results if r["n"]==n and r.get("exact_feasible") is True)
    print(f"  n={n}: {feas_count}/{N_INST} feasible instances")

# ─── Partition-gadget ─────────────────────────────────────────────────────────
print("\n--- Family (c): Partition-gadget (hardness illustration) ---")
gadget_results = []
for m in [5, 8, 10, 12, 15, 18]:
    times_ex = []; times_ls = []; n_to_ls = 0; n_solved_ls = 0
    for inst_id in range(10):
        seed_g = SEED + m*100 + inst_id
        edges_g, A_g, B_g = make_partition_gadget(m, seed_g)
        src_g = 0; dst_g = m; H_g = m; L_g = m

        t0 = time.perf_counter()
        r_ex, s_ex = exhaustive_dfs(edges_g, src_g, dst_g, A_g, B_g, H_g, L_g)
        times_ex.append((time.perf_counter()-t0)*1000)

        t0 = time.perf_counter()
        r_ls, s_ls = pareto_label_dp(edges_g, src_g, dst_g, A_g, B_g, H_g, L_g)
        elapsed_ls = (time.perf_counter()-t0)*1000
        times_ls.append(elapsed_ls)
        if s_ls == "timeout": n_to_ls += 1
        else: n_solved_ls += 1

    med_ex = float(np.median(times_ex)); med_ls = float(np.median(times_ls))
    gadget_results.append({
        "m_stages": m, "median_ex_ms": round(med_ex,2), "median_ls_ms": round(med_ls,2),
        "n_solved_ls": n_solved_ls, "n_timeout_ls": n_to_ls
    })
    print(f"  m={m:>2}: exhaustive={med_ex:7.2f}ms  pareto_label={med_ls:7.2f}ms  "
          f"(solved {n_solved_ls}/10, timeout {n_to_ls}/10)")

# ─── Violation stats ──────────────────────────────────────────────────────────
print("\n--- Naive scalarisation violation rates ---")
naive_viol = {}
for lam in [0.1, 1.0, 10.0]:
    key = f"naive_lam{lam}"
    rows_feas = [r for r in all_results if r.get("exact_feasible") is True and key in r]
    n_viol = sum(1 for r in rows_feas if r[key]["violation"])
    pct = 100*n_viol/max(len(rows_feas),1)
    naive_viol[str(lam)] = {"n_feasible": len(rows_feas), "n_viol": n_viol, "pct": round(pct,1)}
    print(f"  lambda={lam}: {pct:.1f}% violations ({n_viol}/{len(rows_feas)})")

# ─── Solver timing summary ────────────────────────────────────────────────────
print("\n--- Solver timing (family a) ---")
solver_timing = {}
for sname, _ in SOLVERS:
    row_str = f"  {sname:>22}: "
    for n in SIZES_L3:
        rows = [r[sname]["ms"] for r in all_results if r["n"]==n and sname in r
                and r[sname]["result"] is not None]
        if rows:
            med = float(np.median(rows))
            solver_timing[f"{sname}_n{n}"] = round(med, 3)
            row_str += f"n={n}:{med:.2f}ms  "
    print(row_str)

# ─── Save ─────────────────────────────────────────────────────────────────────
out = {
    "description": "L3 complexity v3 — partition-gadget + solver suite (Blueprint K9/L3)",
    "sizes": SIZES_L3, "n_inst": N_INST,
    "n_instances_total": len(all_results),
    "solver_timing_family_a": solver_timing,
    "gadget_results": gadget_results,
    "naive_scalarisation_violations": naive_viol,
    "remark2_finite_class": (
        "With 3 discrete classes, the two-budget problem is polynomial via "
        "count-vector DP O(|E|*n^3). NP-hardness requires continuous/unbounded parameters (K9)."
    ),
}
with open("L3_complexity_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved: L3_complexity_results.json")

# ─── Figure ───────────────────────────────────────────────────────────────────
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(17, 5))

# Panel 1: solver timing family (a)
for sname, col in [("exhaustive_dfs","tomato"),("pareto_label_dp","steelblue"),
                   ("pseudo_poly_dp","green")]:
    pts = [(n, solver_timing.get(f"{sname}_n{n}")) for n in SIZES_L3]
    pts = [(x,y) for x,y in pts if y is not None]
    if pts:
        xs, ys = zip(*pts)
        ax1.semilogy(xs, ys, 'o-', label=sname.replace("_"," "), alpha=0.85)
ax1.set_xlabel("Graph size n"); ax1.set_ylabel("Median runtime (ms, log)")
ax1.set_title("Solver timing\nFamily (a) independent"); ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3, which='both')

# Panel 2: gadget blowup
ms_v = [g["m_stages"] for g in gadget_results]
ls_v = [g["median_ls_ms"] for g in gadget_results]
ex_v = [g["median_ex_ms"] for g in gadget_results]
ax2.semilogy(ms_v, ls_v, 'o-', color='purple', label="Pareto label-setting")
ax2.semilogy(ms_v, ex_v, 's--', color='gray', label="Exhaustive DFS")
ax2.set_xlabel("Gadget stages m"); ax2.set_ylabel("Median runtime (ms, log)")
ax2.set_title("Partition-gadget hardness\nlabel-setting blowup")
ax2.legend(); ax2.grid(True, alpha=0.3, which='both')

# Panel 3: naive violation rates
lam_keys = [0.1, 1.0, 10.0]
viol_pcts = [naive_viol.get(str(l),{}).get("pct",0) for l in lam_keys]
ax3.bar([str(l) for l in lam_keys], viol_pcts, color='steelblue', alpha=0.8)
ax3.set_xlabel("Lambda"); ax3.set_ylabel("Violation rate (%)")
ax3.set_title("Naive scalarisation\nbudget violations"); ax3.set_ylim(0,100)
ax3.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig("fig_L3_complexity.png", dpi=150)
plt.close()
print("Saved: fig_L3_complexity.png")
