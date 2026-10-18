"""
L2_alg2_v5b.py -- R2: Theorem 3's algorithm, correctness and timing.

Implements EXACTLY (Blueprint Appendix A.3):
  bottleneck binary search over the sorted distinct r-values >= R
    -> FEAS(rho): for each simple prefix of exactly L = 2t+1 edges from s
         (edges with r < rho removed; prefix weight w0 = sum eps; prune w0 > Delta)
         -> delete prefix vertices u_0..u_{L-1}
         -> hop-layered Bellman-Ford from u_L with weights eps, <= H-L edges
         -> accept iff w0 + D(d) <= Delta.
Both procedures are compiled with numba (same algorithmic structure; the compiled
language only changes constants, not the operation count).
Second variant ("DFS-feasibility"): bisection + simple-path DFS (the v5 procedure),
with an operation counter and a hard cap so that time-outs are reported, not hidden.
Reference: exhaustive enumeration of all simple paths (B_core.enumerate_paths).

Usage:
  python3 L2_alg2_v5b.py --mode correctness --instances 1000 --max_nodes 12 --t 0 1 2 --seed 7
  python3 L2_alg2_v5b.py --mode timing --sizes 10 20 30 50 75 100 150 --graphs 30 --families er ba --t 1 --reps 5
  python3 L2_alg2_v5b.py --mode timing --sizes 10 20 30 50 --t 2 --graphs 30
"""
import argparse, json, math, os, sys, time
import numpy as np
import numba as nb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B

INF = 1e30


# ------------------------------------------------------------------ CSR graph
def to_csr(n, edges, r, eps):
    deg = np.zeros(n, np.int64)
    for u, v in edges:
        deg[u] += 1; deg[v] += 1
    indptr = np.zeros(n + 1, np.int64)
    indptr[1:] = np.cumsum(deg)
    nbr = np.empty(indptr[-1], np.int32); eid = np.empty(indptr[-1], np.int32)
    pos = indptr[:-1].copy()
    for i, (u, v) in enumerate(edges):
        nbr[pos[u]] = v; eid[pos[u]] = i; pos[u] += 1
        nbr[pos[v]] = u; eid[pos[v]] = i; pos[v] += 1
    return n, indptr, nbr, eid, np.asarray(r, np.float64), np.asarray(eps, np.float64)


# ------------------------------------------------- Theorem 3: prefix + layered BF
@nb.njit(cache=True)
def feas_poly(n, indptr, nbr, eid, r, eps, rho, s, d, L, H, Delta, want, out_nodes, stats):
    """Decision FEAS(G[r>=rho], s, d, L, H, Delta).  If `want`, writes a witness
    (node sequence, padded with -1) to out_nodes.  stats[0] += prefixes examined,
    stats[1] += BF layers executed."""
    if H < L:
        return False
    pn = np.empty(L + 1, np.int32)           # prefix nodes
    pidx = np.zeros(L + 1, np.int64)         # next arc index per depth
    pw = np.zeros(L + 1, np.float64)         # prefix weight
    inpre = np.zeros(n, np.bool_)
    dist = np.empty(n, np.float64)
    ndist = np.empty(n, np.float64)
    K = H - L
    par = np.empty((K + 1, n), np.int32)     # parent per layer (only used if want)
    depth = 0
    pn[0] = s
    inpre[s] = True
    pidx[0] = indptr[s]
    pw[0] = 0.0
    while depth >= 0:
        u = pn[depth]
        if pidx[depth] >= indptr[u + 1]:
            inpre[u] = False
            depth -= 1
            continue
        a = pidx[depth]
        pidx[depth] += 1
        e = eid[a]
        v = nbr[a]
        if r[e] < rho or inpre[v]:
            continue
        w = pw[depth] + eps[e]
        if w > Delta:
            continue
        if v == d and depth + 1 < L:
            continue                           # d cannot be an interior prefix vertex
        if depth + 1 < L:
            depth += 1
            pn[depth] = v
            inpre[v] = True
            pidx[depth] = indptr[v]
            pw[depth] = w
            continue
        # prefix of exactly L edges ends at v = u_L, weight w
        stats[0] += 1
        if v == d:
            if want:
                for k in range(L):
                    out_nodes[k] = pn[k]
                out_nodes[L] = v
                for k in range(L + 1, out_nodes.shape[0]):
                    out_nodes[k] = -1
            return True
        if K <= 0:
            continue
        # hop-layered Bellman-Ford on G' = G minus {u_0..u_{L-1}} (prefix vertices are inpre)
        for i in range(n):
            dist[i] = INF
        dist[v] = 0.0
        if want:
            for h in range(K + 1):
                for i in range(n):
                    par[h, i] = -1
        for h in range(1, K + 1):
            stats[1] += 1
            for i in range(n):
                ndist[i] = dist[i]            # "<= h hops" (allow staying)
                if want:
                    par[h, i] = -2            # -2 : same as layer h-1
            for x in range(n):
                dx = dist[x]
                if dx >= INF or inpre[x]:
                    continue
                for b in range(indptr[x], indptr[x + 1]):
                    ee = eid[b]
                    if r[ee] < rho:
                        continue
                    y = nbr[b]
                    if inpre[y]:
                        continue
                    nd = dx + eps[ee]
                    if nd < ndist[y]:
                        ndist[y] = nd
                        if want:
                            par[h, y] = x
            for i in range(n):
                dist[i] = ndist[i]
            if w + dist[d] <= Delta:
                if want:
                    # reconstruct walk d -> v through layers
                    tmp = np.empty(H + 2, np.int32)
                    m = 0
                    cur = d
                    hh = h
                    while cur != v or hh > 0:
                        if hh == 0:
                            break
                        p_ = par[hh, cur]
                        if p_ == -2:
                            hh -= 1
                            continue
                        tmp[m] = cur
                        m += 1
                        cur = p_
                        hh -= 1
                    # cur == v here
                    for k in range(L):
                        out_nodes[k] = pn[k]
                    out_nodes[L] = v
                    pos = L + 1
                    for k in range(m - 1, -1, -1):
                        out_nodes[pos] = tmp[k]
                        pos += 1
                    for k in range(pos, out_nodes.shape[0]):
                        out_nodes[k] = -1
                return True
    return False


# ----------------------------------------------------- DFS-feasibility variant
@nb.njit(cache=True)
def feas_dfs(n, indptr, nbr, eid, r, eps, rho, s, d, L, H, Delta, op_cap, stats):
    """v5 procedure: simple-path DFS with budget pruning.  Returns 1 feasible,
    0 infeasible, -1 op-cap hit (time-out).  stats[0] += nodes expanded."""
    nodes = np.empty(H + 2, np.int32)
    idx = np.zeros(H + 2, np.int64)
    wsum = np.zeros(H + 2, np.float64)
    vis = np.zeros(n, np.bool_)
    depth = 0
    nodes[0] = s
    vis[s] = True
    idx[0] = indptr[s]
    ops = 0
    while depth >= 0:
        u = nodes[depth]
        if idx[depth] >= indptr[u + 1]:
            vis[u] = False
            depth -= 1
            continue
        a = idx[depth]
        idx[depth] += 1
        ops += 1
        if ops > op_cap:
            stats[0] += ops
            return -1
        e = eid[a]
        v = nbr[a]
        if r[e] < rho or vis[v]:
            continue
        w = wsum[depth] + eps[e]
        if w > Delta:
            continue
        hh = depth + 1
        if v == d:
            if hh >= L:
                stats[0] += ops
                return 1
            continue
        if hh < H:
            depth += 1
            nodes[depth] = v
            vis[v] = True
            idx[depth] = indptr[v]
            wsum[depth] = w
    stats[0] += ops
    return 0


def alg2(G, s, d, R, Delta, H, L, variant="poly", op_cap=2_000_000_000, want_witness=False):
    """Bottleneck binary search. Returns (rho*, witness or None, stats dict). rho*=-1 infeasible,
    rho*=None on DFS time-out."""
    n, indptr, nbr, eid, r, eps = G
    rv = np.unique(r[r >= R])[::-1]                  # descending distinct r >= R
    st = np.zeros(2, np.int64)
    if rv.size == 0 or H < L:
        return -1.0, None, dict(prefixes=0, layers=0, ops=0, timeout=False)
    out = np.full(H + 2, -1, np.int32)

    def feas(rho, want=False):
        if variant == "poly":
            return 1 if feas_poly(n, indptr, nbr, eid, r, eps, float(rho), s, d, L, H, Delta, want, out, st) else 0
        return feas_dfs(n, indptr, nbr, eid, r, eps, float(rho), s, d, L, H, Delta, op_cap, st)

    # monotone in rho: feasible(rho) => feasible(rho') for rho' < rho.  Find the first (largest) feasible.
    lo, hi, best = 0, len(rv) - 1, -1
    # quick infeasibility test at the smallest threshold
    f_last = feas(rv[hi])
    if f_last == -1:
        return None, None, dict(prefixes=int(st[0]), layers=int(st[1]), ops=int(st[0]), timeout=True)
    if f_last == 0:
        return -1.0, None, dict(prefixes=int(st[0]), layers=int(st[1]), ops=int(st[0]), timeout=False)
    best = hi
    hi -= 1
    while lo <= hi:
        mid = (lo + hi) // 2
        f = feas(rv[mid])
        if f == -1:
            return None, None, dict(prefixes=int(st[0]), layers=int(st[1]), ops=int(st[0]), timeout=True)
        if f == 1:
            best = mid; hi = mid - 1
        else:
            lo = mid + 1
    wit = None
    if want_witness and variant == "poly":
        feas(rv[best], True)
        wit = [int(x) for x in out if x >= 0]
    return float(rv[best]), wit, dict(prefixes=int(st[0]), layers=int(st[1]), ops=int(st[0]), timeout=False)


# ---------------------------------------------------------------- exhaustive ref
def exhaustive(n, edges, r, eps, s, d, R, Delta, H, L):
    if H < L:
        return -1.0
    adj = B.adj_matrix(n, edges)
    Lp, hp = B.enumerate_paths(adj, s, d, L, H, cap=4_000_000)
    if len(hp) == 0:
        return -1.0
    E = len(r)
    rr = np.concatenate([r, [np.inf]])[np.where(Lp < 0, E, Lp)]
    ee = np.concatenate([eps, [0.0]])[np.where(Lp < 0, E, Lp)]
    rmin = rr.min(1); se = ee.sum(1)
    ok = (rmin >= R) & (se <= Delta * (1 + 1e-12))
    return float(rmin[ok].max()) if ok.any() else -1.0


def check_witness(n, edges, r, eps, wit, s, d, rho, R, Delta, H, L):
    eidx = {}
    for i, (u, v) in enumerate(edges):
        eidx[(u, v)] = i; eidx[(v, u)] = i
    if wit is None or wit[0] != s or wit[-1] != d or len(set(wit)) != len(wit):
        return False
    hops = len(wit) - 1
    if not (L <= hops <= H):
        return False
    es = [eidx[(wit[i], wit[i + 1])] for i in range(hops)]
    if min(r[e] for e in es) < rho - 1e-12 or sum(eps[e] for e in es) > Delta * (1 + 1e-12):
        return False
    return abs(min(r[e] for e in es) - rho) < 1e-12 or min(r[e] for e in es) >= rho


# ---------------------------------------------------------------- generators
def gen_correctness_instance(rng):
    n = int(rng.integers(4, 13))
    p = rng.uniform(0.2, 0.6)
    iu, ju = np.triu_indices(n, 1)
    keep = rng.random(len(iu)) < p
    edges = list(zip(iu[keep].tolist(), ju[keep].tolist()))
    E = len(edges)
    if E == 0:
        return None
    r = np.array([B.h(2.0 + rng.random() * 0.7) for _ in range(E)])
    eps = rng.uniform(1e-7, 5e-6, E)
    return n, edges, r, eps


def gen_timing_graph(fam, n, seed, Delta_mode="clean"):
    rng = np.random.default_rng(np.random.SeedSequence([seed, n, 0 if fam == "er" else 1]))
    edges = B.build_graph(fam, n, rng)
    lk = B.make_links(rng, len(edges), 10000, f=0.0, noise="binomial")
    return edges, lk["r"], lk["eps"]


def run_correctness(a):
    rng = np.random.default_rng(a.seed)
    n_tot = n_ok = n_feas = n_wit = 0
    by_t = {t: dict(total=0, ok=0, feasible=0, dfs_ok=0) for t in a.t}
    fails = []
    # warm up numba
    G0 = to_csr(3, [(0, 1), (1, 2)], [0.2, 0.2], [1e-7, 1e-7])
    alg2(G0, 0, 2, 0.1, 1e-5, 2, 1)
    alg2(G0, 0, 2, 0.1, 1e-5, 2, 1, variant="dfs")
    for i in range(a.instances):
        inst = gen_correctness_instance(rng)
        if inst is None:
            continue
        n, edges, r, eps = inst
        t = int(rng.choice(a.t)); L = 2 * t + 1
        Delta = float(rng.uniform(5e-7, 2e-5)); R = 0.05
        H = min(n - 1, L + int(rng.integers(0, 4)))
        H = max(H, L)
        s, d = 0, n - 1
        ref = exhaustive(n, edges, r, eps, s, d, R, Delta, H, L)
        G = to_csr(n, edges, r, eps)
        rho, wit, _ = alg2(G, s, d, R, Delta, H, L, "poly", want_witness=True)
        rho_d, _, _ = alg2(G, s, d, R, Delta, H, L, "dfs")
        n_tot += 1; by_t[t]["total"] += 1
        same = abs(ref - rho) < 1e-12
        same_d = rho_d is not None and abs(ref - rho_d) < 1e-12
        wit_ok = True
        if ref >= 0:
            n_feas += 1; by_t[t]["feasible"] += 1
            wit_ok = check_witness(n, edges, r, eps, wit, s, d, rho, R, Delta, H, L)
            n_wit += int(wit_ok)
        if same and wit_ok:
            n_ok += 1; by_t[t]["ok"] += 1
        else:
            fails.append(dict(i=i, n=n, t=t, ref=ref, poly=rho, witness_ok=wit_ok))
        by_t[t]["dfs_ok"] += int(same_d)
    out = dict(instances=n_tot, optimum_equal=n_ok, feasible_instances=n_feas, witnesses_valid=n_wit,
               by_t=by_t, failures=fails[:5], seed=a.seed,
               note="paths may differ on ties; equality is on the optimal bottleneck value; every returned witness is "
                    "re-checked for simplicity, hop window, rate threshold and budget")
    print(json.dumps(out, indent=1))
    json.dump(out, open(a.out_corr, "w"), indent=1)


def fit_exponent(sizes, med, per_graph, B_boot=2000, seed=1):
    x = np.log(np.asarray(sizes, float)); y = np.log(np.asarray(med, float))
    slope = float(np.polyfit(x, y, 1)[0])
    rng = np.random.default_rng(seed)
    sl = []
    for _ in range(B_boot):
        yb = []
        for pg in per_graph:
            pg = np.asarray(pg); yb.append(np.median(rng.choice(pg, size=len(pg), replace=True)))
        sl.append(np.polyfit(x, np.log(np.asarray(yb)), 1)[0])
    return slope, float(np.percentile(sl, 2.5)), float(np.percentile(sl, 97.5))


def run_timing(a):
    rows = []
    # warm-up
    G0 = to_csr(4, [(0, 1), (1, 2), (2, 3)], [0.2] * 3, [1e-7] * 3)
    for var in ("poly", "dfs"):
        alg2(G0, 0, 3, 0.1, 1e-5, 4, 3, var)
    variants = a.variants
    for t in a.t:
        L = 2 * t + 1
        for fam in a.families:
            for Delta in a.Delta:
                for n in a.sizes:
                    H = a.H if a.H else L + 3 - 1      # H = L+2  (t=1 -> 5, t=2 -> 7)
                    H = min(H, n - 1)
                    for var in variants:
                        per_graph, nf, nto, prefixes, ops = [], 0, 0, [], []
                        for g in range(a.graphs):
                            edges, r, eps = gen_timing_graph(fam, n, 8000 + g)
                            if not edges:
                                continue
                            if a.stress:   # worst case: d unreachable, no budget pruning -> every prefix is tried
                                keep = [i for i, (u, v) in enumerate(edges) if u != n - 1 and v != n - 1]
                                edges = [edges[i] for i in keep]; r = r[keep]; eps = eps[keep]
                            G = to_csr(n, edges, r, eps)
                            s, d = 0, n - 1
                            ts = []
                            res = None
                            for _ in range(a.reps):
                                t0 = time.perf_counter()
                                res = alg2(G, s, d, 0.10, Delta, H, L, var, op_cap=a.op_cap)
                                ts.append((time.perf_counter() - t0) * 1000)
                                if res[0] is None:
                                    break
                            if res[0] is None:
                                nto += 1
                                continue
                            per_graph.append(float(np.median(ts)))
                            nf += int(res[0] >= 0)
                            prefixes.append(res[2]["prefixes"] if var == "poly" else res[2]["ops"])
                        if not per_graph:
                            rows.append(dict(t=t, family=fam, Delta=Delta, size=n, variant=var, graphs=0, timeouts=nto))
                            print(f"t={t} {fam} D={Delta:g} n={n:3d} {var:4s}: ALL TIMED OUT ({nto})", flush=True)
                            continue
                        pg = np.array(per_graph)
                        row = dict(t=t, family=fam, Delta=Delta, size=n, variant=var, H=H, L=L, graphs=len(pg),
                                   timeouts=nto, feasible=nf, median_ms=float(np.median(pg)),
                                   q25_ms=float(np.percentile(pg, 25)), q75_ms=float(np.percentile(pg, 75)),
                                   max_ms=float(pg.max()), median_work=float(np.median(prefixes)),
                                   per_graph_ms=[round(x, 4) for x in per_graph])
                        rows.append(row)
                        print(f"t={t} {fam} D={Delta:g} n={n:3d} {var:4s}: median={row['median_ms']:9.2f} ms "
                              f"IQR=[{row['q25_ms']:.2f},{row['q75_ms']:.2f}] feasible={nf}/{len(pg)} "
                              f"timeouts={nto} work={row['median_work']:.3g}", flush=True)
    out = dict(rows=rows, args=vars(a), platform="Apple M-series arm64, Python 3.14, numba 0.67, single thread")
    json.dump(out, open(a.out_time, "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["correctness", "timing"])
    ap.add_argument("--instances", type=int, default=1000)
    ap.add_argument("--max_nodes", type=int, default=12)
    ap.add_argument("--t", type=int, nargs="+", default=[1])
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--sizes", type=int, nargs="+", default=[10, 20, 30, 50, 75, 100, 150])
    ap.add_argument("--graphs", type=int, default=30)
    ap.add_argument("--families", nargs="+", default=["er", "ba"])
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--Delta", type=float, nargs="+", default=[5e-6])
    ap.add_argument("--H", type=int, default=0)
    ap.add_argument("--stress", action="store_true")
    ap.add_argument("--variants", nargs="+", default=["poly", "dfs"])
    ap.add_argument("--op_cap", type=int, default=3_000_000_000)
    ap.add_argument("--out_corr", default="../data/R2_correctness.json")
    ap.add_argument("--out_time", default="../data/R2_timing.json")
    a = ap.parse_args()
    run_correctness(a) if a.mode == "correctness" else run_timing(a)
