"""
L3_complexity_v4.py -- R3: corrected complexity experiment (Blueprint W9 / K9).

Every solver is instrumented with a work counter:
  enumeration    : simple s-d paths enumerated (no budget pruning, no early exit) + nodes expanded
  label-setting  : labels created (all candidate labels before dominance) and labels kept
  ILP            : HiGHS (scipy.optimize.milp) branch-and-bound nodes are not exposed by scipy -> time & status only
  pseudo-poly DP: DP cells touched
  LARAC-style    : scalar layered shortest-path calls
  FPTAS-style    : DP cells touched (rounded b, eps=0.10)
All solvers decide the same question:  exists s-d path with hops in [L,H], sum a <= A, sum b <= B ?
Weights are INTEGERS (units of 1e-9), so every solver sees an identical instance.
Walks vs paths: with L=1 a feasible walk implies a feasible simple path (shortcutting only lowers a,b), so layered
solvers (label-setting, DP, LARAC, FPTAS) are exact for L=1 on family (a); the gadget is a DAG (walks = paths).

Parts (see __main__):  fam  --n N        family (a) at one size, results -> ../data/R3_fam_n{N}.json
                       gadget            m = 5..30 all solvers            -> ../data/R3_gadget.json
                       naive             >=100 instances, lambda sweep    -> ../data/R3_naive.json
                       dp3               3-class count-vector DP timing   -> ../data/R3_dp3.json
"""
import argparse, json, math, os, sys, time
import numpy as np
import numba as nb

TIMEOUT = 60.0                 # seconds per solver call (spec)
CHUNK = 20_000_000             # numba ops between wall-clock checks
H_FAM = 6
HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------- instances
def csr(n, fr, to):
    order = np.argsort(fr, kind="stable")
    indptr = np.zeros(n + 1, np.int64)
    np.add.at(indptr, fr + 1, 1)
    indptr = np.cumsum(indptr)
    return indptr, order


def has_path(n, fr, to, s, d, H):
    seen = {s}; front = [s]
    adj = {}
    for u, v in zip(fr.tolist(), to.tolist()):
        adj.setdefault(u, []).append(v)
    for _ in range(H):
        nxt = []
        for u in front:
            for v in adj.get(u, []):
                if v not in seen:
                    seen.add(v); nxt.append(v)
        front = nxt
    return d in seen


def make_family_a(n, seed, p=0.25):
    """ER(n,p), log-uniform integer weights in units of 1e-9: a in [200,5000], b in [100,4000]."""
    k = 0
    while True:
        rng = np.random.default_rng(np.random.SeedSequence([seed, n, k])); k += 1
        iu, ju = np.triu_indices(n, 1)
        keep = rng.random(len(iu)) < p
        u, v = iu[keep], ju[keep]
        E = len(u)
        if E == 0:
            continue
        a = np.rint(np.exp(rng.uniform(math.log(200), math.log(5000), E))).astype(np.int64)
        b = np.rint(np.exp(rng.uniform(math.log(100), math.log(4000), E))).astype(np.int64)
        fr = np.concatenate([u, v]).astype(np.int32); to = np.concatenate([v, u]).astype(np.int32)
        aa = np.concatenate([a, a]); bb = np.concatenate([b, b])
        if has_path(n, fr, to, 0, n - 1, H_FAM):
            return dict(n=n, fr=fr, to=to, a=aa, b=bb, s=0, d=n - 1, H=H_FAM, L=1, a_und=a, b_und=b)


def make_gadget(m, seed, xmax=10000, c=1):
    """Partition gadget, m series stages with parallel links U_i=(x_i+c, c), D_i=(c, x_i+c); A=B=X/2+m*c.
    X is forced even, so the instance is feasible iff the multiset x has an equal-sum split (random for small m)."""
    rng = np.random.default_rng(np.random.SeedSequence([seed, m]))
    x = rng.integers(1, xmax + 1, m)
    if x.sum() % 2:
        x[-1] += 1
    X = int(x.sum()); A = B = X // 2 + m * c
    fr, to, a, b = [], [], [], []
    for i in range(m):
        fr += [i, i]; to += [i + 1, i + 1]
        a += [int(x[i]) + c, c]; b += [c, int(x[i]) + c]
    return dict(n=m + 1, fr=np.array(fr, np.int32), to=np.array(to, np.int32), a=np.array(a, np.int64),
                b=np.array(b, np.int64), s=0, d=m, H=m, L=m, A=A, B=B, x=x)


# ------------------------------------------------- solver 1: exhaustive enumeration
@nb.njit(cache=True)
def _enum_chunk(indptr, nbr, wa, wb, d, H, L, A, B, nodes, idx, ca, cb, vis, st, max_ops):
    """Resumable simple-path DFS.  st = [depth, paths, expanded, feasible_count, ops].
    NO pruning on budgets, NO early exit: every simple path with <= H hops that ends at d is counted."""
    depth = st[0]
    ops = 0
    while depth >= 0 and ops < max_ops:
        u = nodes[depth]
        if idx[depth] >= indptr[u + 1]:
            vis[u] = False
            depth -= 1
            continue
        k = idx[depth]; idx[depth] += 1
        ops += 1
        v = nbr[k]
        if vis[v]:
            continue
        na = ca[depth] + wa[k]; nbv = cb[depth] + wb[k]
        h = depth + 1
        if v == d:
            st[1] += 1
            if h >= L and na <= A and nbv <= B:
                st[3] += 1
            continue                       # a simple path cannot pass through d and come back
        if h < H:
            st[2] += 1
            depth += 1
            nodes[depth] = v; vis[v] = True; idx[depth] = indptr[v]
            ca[depth] = na; cb[depth] = nbv
    st[0] = depth
    st[4] += ops
    return depth < 0


def solve_enum(I, A, B, timeout=TIMEOUT):
    n = I["n"]; order = np.argsort(I["fr"], kind="stable")
    fr = I["fr"][order]; nbr = I["to"][order].astype(np.int32); wa = I["a"][order]; wb = I["b"][order]
    indptr = np.zeros(n + 1, np.int64); np.add.at(indptr, fr.astype(np.int64) + 1, 1); indptr = np.cumsum(indptr)
    H = I["H"]
    nodes = np.zeros(H + 2, np.int32); idx = np.zeros(H + 2, np.int64)
    ca = np.zeros(H + 2, np.int64); cb = np.zeros(H + 2, np.int64); vis = np.zeros(n, np.bool_)
    s = I["s"]; nodes[0] = s; vis[s] = True; idx[0] = indptr[s]
    st = np.zeros(5, np.int64)
    t0 = time.perf_counter(); done = False
    while not done:
        done = _enum_chunk(indptr, nbr, wa, wb, I["d"], H, I["L"], A, B, nodes, idx, ca, cb, vis, st, CHUNK)
        if not done and time.perf_counter() - t0 > timeout:
            break
    return dict(status="ok" if done else "timeout", feasible=(bool(st[3] > 0) if done else None),
                feasible_found_before_stop=bool(st[3] > 0), paths=int(st[1]), expanded=int(st[2]), ops=int(st[4]),
                ms=(time.perf_counter() - t0) * 1000)


# ------------------------------------------------- solver 2: Pareto label-setting (hop-layered)
LABEL_CAP = 20_000_000


@nb.njit(cache=True)
def _count_ext(lab_node, indptr, d):
    c = 0
    for i in range(lab_node.shape[0]):
        u = lab_node[i]
        if u != d:
            c += indptr[u + 1] - indptr[u]
    return c


@nb.njit(cache=True)
def _fill_ext(lab_node, lab_a, lab_b, indptr, nbr, wa, wb, d, A, B, out_n, out_a, out_b):
    k = 0
    for i in range(lab_node.shape[0]):
        u = lab_node[i]
        if u == d:
            continue
        for j in range(indptr[u], indptr[u + 1]):
            na = lab_a[i] + wa[j]; nb_ = lab_b[i] + wb[j]
            if na <= A and nb_ <= B:
                out_n[k] = nbr[j]; out_a[k] = na; out_b[k] = nb_; k += 1
    return k


@nb.njit(cache=True)
def _sweep(node, a, b, keep):
    m = node.shape[0]; cnt = 0
    cur = -1; minb = 1 << 62
    for i in range(m):
        if node[i] != cur:
            cur = node[i]; minb = 1 << 62
        if b[i] < minb:
            keep[i] = True; minb = b[i]; cnt += 1
        else:
            keep[i] = False
    return cnt


def _csr_of(I):
    n = I["n"]; order = np.argsort(I["fr"], kind="stable")
    fr = I["fr"][order]; nbr = I["to"][order].astype(np.int32); wa = I["a"][order]; wb = I["b"][order]
    indptr = np.zeros(n + 1, np.int64); np.add.at(indptr, fr.astype(np.int64) + 1, 1)
    return np.cumsum(indptr), nbr, wa, wb


def solve_labels(I, A, B, timeout=TIMEOUT, return_front=False):
    """Exact decision by hop-layered Pareto label-setting.  Counters: created = labels generated within both budgets,
    kept = labels surviving dominance (summed over layers), peak = largest layer."""
    indptr, nbr, wa, wb = _csr_of(I)
    n, s, d, H, L = I["n"], I["s"], I["d"], I["H"], I["L"]
    t0 = time.perf_counter()
    ln = np.array([s], np.int32); la = np.zeros(1, np.int64); lb = np.zeros(1, np.int64)
    created = kept = peak = 0; feasible = False; status = "ok"; front = None
    for h in range(1, H + 1):
        c = int(_count_ext(ln, indptr, d))
        if c > LABEL_CAP:
            status = "capped"; break
        if time.perf_counter() - t0 > timeout:
            status = "timeout"; break
        on = np.empty(c, np.int32); oa = np.empty(c, np.int64); ob = np.empty(c, np.int64)
        k = int(_fill_ext(ln, la, lb, indptr, nbr, wa, wb, d, A, B, on, oa, ob))
        on, oa, ob = on[:k], oa[:k], ob[:k]
        created += k
        if k == 0:
            ln = on; la = oa; lb = ob; break
        o = np.lexsort((ob, oa, on)); on, oa, ob = on[o], oa[o], ob[o]
        keep = np.empty(k, np.bool_); kc = int(_sweep(on, oa, ob, keep))
        ln, la, lb = on[keep], oa[keep], ob[keep]
        kept += kc; peak = max(peak, kc)
        if h >= L and np.any(ln == d):
            feasible = True
            if return_front:
                m_ = ln == d; front = (la[m_], lb[m_])
    return dict(status=status, feasible=(feasible if status == "ok" else None), created=created, kept=kept,
                peak=peak, ms=(time.perf_counter() - t0) * 1000, front=front)


# ------------------------------------------------- solver 3: pseudo-polynomial DP (exact b) and FPTAS-style
@nb.njit(cache=True)
def _dp_layers(n, fr, to, wa, wb, s, d, H, L, A, Bcap, out_minb):
    """dp[v,k] = min a-cost over walks of the current hop count whose integer b-cost equals k.  Returns cells touched.
    out_minb[0] = min k with some h in [L,H], dp[d,k] <= A (or -1)."""
    INF = 1 << 60
    K = Bcap + 1
    cur = np.full((n, K), INF, np.int64); nxt = np.empty((n, K), np.int64)
    cur[s, 0] = 0
    cells = 0; best = -1
    for h in range(1, H + 1):
        nxt[:, :] = INF
        for e in range(fr.shape[0]):
            u = fr[e]; v = to[e]
            if u == d:
                continue
            a_ = wa[e]; b_ = wb[e]
            for k in range(0, K - b_):
                x = cur[u, k]
                cells += 1
                if x < INF:
                    y = x + a_
                    if y <= A and y < nxt[v, k + b_]:
                        nxt[v, k + b_] = y
        if h >= L:
            for k in range(K):
                if nxt[d, k] <= A:
                    if best < 0 or k < best:
                        best = k
                    break
        cur, nxt = nxt, cur
    out_minb[0] = best
    return cells


def solve_dp(I, A, B, Bcap=None, timeout=TIMEOUT, cell_cap=300_000_000):
    Bc = int(B if Bcap is None else Bcap)
    if I["n"] * (Bc + 1) > cell_cap:
        return dict(status="na_size", feasible=None, cells=0, ms=0.0, table_cells=I["n"] * (Bc + 1))
    out = np.zeros(1, np.int64)
    t0 = time.perf_counter()
    cells = _dp_layers(I["n"], I["fr"], I["to"], I["a"], I["b"], I["s"], I["d"], I["H"], I["L"], int(A), Bc, out)
    return dict(status="ok", feasible=bool(out[0] >= 0 and out[0] <= B), cells=int(cells), minb=int(out[0]),
                ms=(time.perf_counter() - t0) * 1000, table_cells=I["n"] * (Bc + 1))


def solve_fptas(I, A, B, eps=0.10):
    """Round b up to multiples of theta = eps*B/H; budget B' = floor(B/theta)+H.  Accepts every truly feasible
    instance; may accept infeasible ones with true b-overshoot <= eps*B."""
    H = I["H"]; theta = eps * B / H
    b2 = np.ceil(I["b"] / theta).astype(np.int64)
    Bp = int(B // theta) + H
    J = dict(I); J["b"] = b2
    out = np.zeros(1, np.int64)
    t0 = time.perf_counter()
    cells = _dp_layers(I["n"], I["fr"], I["to"], I["a"], b2, I["s"], I["d"], I["H"], I["L"], int(A), Bp, out)
    return dict(status="ok", feasible=bool(out[0] >= 0), cells=int(cells), ms=(time.perf_counter() - t0) * 1000)


# ------------------------------------------------- solver 4: LARAC-style (Lagrangian scan) and naive scalarisation
@nb.njit(cache=True)
def _scalar_layers(n, fr, to, wa, wb, s, d, H, L, ca, cb, A, B, prune, res):
    """Hop-layered shortest path for w = ca*a + cb*b with ONE label per (hop, node).  prune: discard labels that
    exceed A or B.  res = [found, a, b, hops].  Returns number of labels relaxed."""
    INF = 1e300
    w = np.full((H + 1, n), INF); aa = np.zeros((H + 1, n), np.int64); bb = np.zeros((H + 1, n), np.int64)
    w[0, s] = 0.0
    work = 0
    for h in range(1, H + 1):
        for e in range(fr.shape[0]):
            u = fr[e]; v = to[e]
            if u == d or w[h - 1, u] >= INF:
                continue
            work += 1
            na = aa[h - 1, u] + wa[e]; nb_ = bb[h - 1, u] + wb[e]
            if prune and (na > A or nb_ > B):
                continue
            nw = ca * na + cb * nb_
            if nw < w[h, v]:
                w[h, v] = nw; aa[h, v] = na; bb[h, v] = nb_
    best = INF; res[0] = 0
    for h in range(L, H + 1):
        if w[h, d] < best:
            best = w[h, d]; res[0] = 1; res[1] = aa[h, d]; res[2] = bb[h, d]; res[3] = h
    return work


def scalar_path(I, ca, cb, A=0, B=0, prune=False):
    res = np.zeros(4, np.int64)
    work = _scalar_layers(I["n"], I["fr"], I["to"], I["a"], I["b"], I["s"], I["d"], I["H"], I["L"],
                          float(ca), float(cb), int(A), int(B), prune, res)
    return bool(res[0]), int(res[1]), int(res[2]), int(res[3]), int(work)


LAMBDAS_LARAC = [10.0 ** e for e in np.linspace(-3, 3, 31)]


def solve_larac(I, A, B):
    t0 = time.perf_counter(); calls = 0; work = 0; found = False
    for lam in LAMBDAS_LARAC:
        ok, a, b, h, w = scalar_path(I, 1.0, lam); calls += 1; work += w
        if ok and a <= A and b <= B:
            found = True; break
    return dict(status="ok", feasible=found, calls=calls, work=work, ms=(time.perf_counter() - t0) * 1000)


# ------------------------------------------------- solver 5: ILP (HiGHS via scipy.optimize.milp)
def solve_ilp(I, A, B, timeout=TIMEOUT):
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix, vstack
    n = I["n"]; fr, to = I["fr"].astype(np.int64), I["to"].astype(np.int64); m = len(fr)
    rows = np.concatenate([fr, to]); cols = np.concatenate([np.arange(m), np.arange(m)])
    vals = np.concatenate([np.ones(m), -np.ones(m)])
    Aflow = coo_matrix((vals, (rows, cols)), shape=(n, m)).tocsr()
    rhs = np.zeros(n); rhs[I["s"]] = 1; rhs[I["d"]] = -1
    mat = vstack([Aflow, coo_matrix(I["a"][None, :].astype(float)), coo_matrix(I["b"][None, :].astype(float)),
                  coo_matrix(np.ones((1, m)))]).tocsr()
    lo = np.concatenate([rhs, [-np.inf, -np.inf, I["L"]]]); hi = np.concatenate([rhs, [A, B, I["H"]]])
    t0 = time.perf_counter()
    r = milp(c=np.zeros(m), constraints=LinearConstraint(mat, lo, hi), integrality=np.ones(m),
             bounds=Bounds(0, 1), options=dict(time_limit=timeout, presolve=True))
    ms = (time.perf_counter() - t0) * 1000
    if r.status == 0:
        return dict(status="ok", feasible=True, ms=ms)
    if r.status == 2:
        return dict(status="ok", feasible=False, ms=ms)
    return dict(status="timeout", feasible=None, ms=ms, highs_status=int(r.status))


SOLVERS = ["enum", "label", "dp", "ilp", "larac", "fptas"]


def run_all(I, A, B, which=SOLVERS, timeout=TIMEOUT):
    out = {}
    if "enum" in which:  out["enum"] = solve_enum(I, A, B, timeout)
    if "label" in which: out["label"] = solve_labels(I, A, B, timeout)
    if "dp" in which:    out["dp"] = solve_dp(I, A, B)
    if "ilp" in which:   out["ilp"] = solve_ilp(I, A, B, timeout)
    if "larac" in which: out["larac"] = solve_larac(I, A, B)
    if "fptas" in which: out["fptas"] = solve_fptas(I, A, B)
    return out


# ------------------------------------------------- budgets: "30th Pareto percentile"
def pareto_budgets(I, q=0.30):
    """Unconstrained Pareto front F of s-d (a,b) pairs over hops in [L,H] (label-setting, budgets = inf).
    A = min_a + q*(max_a - min_a) over F, B likewise: the budget box covers the lowest q-quantile of the trade-off range
    in BOTH dimensions.  Returns (A, B, front a-array, b-array)."""
    big = 1 << 55
    r = solve_labels(I, big, big, timeout=TIMEOUT, return_front=True)
    # front over all layers h>=L: re-collect (solve_labels returns last layer front); gather layer-wise instead
    return r


def front_all_layers(I, timeout=TIMEOUT):
    indptr, nbr, wa, wb = _csr_of(I)
    n, s, d, H, L = I["n"], I["s"], I["d"], I["H"], I["L"]
    big = 1 << 55
    ln = np.array([s], np.int32); la = np.zeros(1, np.int64); lb = np.zeros(1, np.int64)
    fa, fb = [], []
    for h in range(1, H + 1):
        c = int(_count_ext(ln, indptr, d))
        if c == 0 or c > LABEL_CAP:
            break
        on = np.empty(c, np.int32); oa = np.empty(c, np.int64); ob = np.empty(c, np.int64)
        k = int(_fill_ext(ln, la, lb, indptr, nbr, wa, wb, d, big, big, on, oa, ob))
        on, oa, ob = on[:k], oa[:k], ob[:k]
        o = np.lexsort((ob, oa, on)); on, oa, ob = on[o], oa[o], ob[o]
        keep = np.empty(k, np.bool_); _sweep(on, oa, ob, keep)
        ln, la, lb = on[keep], oa[keep], ob[keep]
        if h >= L:
            m_ = ln == d
            fa.append(la[m_]); fb.append(lb[m_])
    if not fa:
        return None
    fa = np.concatenate(fa); fb = np.concatenate(fb)
    o = np.lexsort((fb, fa)); fa, fb = fa[o], fb[o]
    keep = np.empty(len(fa), np.bool_); minb = 1 << 62
    for i in range(len(fa)):
        keep[i] = fb[i] < minb
        minb = min(minb, fb[i])
    return fa[keep], fb[keep]


def budgets_from_front(F, q):
    """Budgets at the q-th percentile POINT of the unconstrained Pareto front (sorted by a): A=a_P, B=b_P, so P itself is
    feasible and both budgets exclude the single-objective optima whenever the front has >1 point."""
    fa, fb = F
    k = int(round(q * (len(fa) - 1)))
    return int(fa[k]), int(fb[k])


def joint_binding(F, A, B):
    """Both budgets bind: the min-a front point violates B AND the min-b front point violates A."""
    fa, fb = F
    return bool(fb[0] > B and fa[-1] > A)


# ================================================================== experiments
def planted_gadget(m, seed, xmax=10000, c=1):
    """Gadget whose item multiset has an exact equal-sum split (feasible by construction)."""
    rng = np.random.default_rng(np.random.SeedSequence([seed, m, 77]))
    while True:
        x = rng.integers(1, xmax + 1, m)
        side = rng.random(m) < 0.5
        diff = int(x[:-1][side[:-1]].sum() - x[:-1][~side[:-1]].sum())
        if diff == 0:
            continue
        x[-1] = abs(diff)
        if 1 <= x[-1] <= xmax * 3:
            break
    X = int(x.sum()); A = B = X // 2 + m * c
    G = make_gadget(m, seed)  # structure only
    fr, to, a, b = [], [], [], []
    for i in range(m):
        fr += [i, i]; to += [i + 1, i + 1]
        a += [int(x[i]) + c, c]; b += [c, int(x[i]) + c]
    return dict(n=m + 1, fr=np.array(fr, np.int32), to=np.array(to, np.int32), a=np.array(a, np.int64),
                b=np.array(b, np.int64), s=0, d=m, H=m, L=m, A=A, B=B, x=x)


def warmup():
    I = make_family_a(10, 99); F = front_all_layers(I); A, B = budgets_from_front(F, .3)
    run_all(I, A, B, timeout=5); scalar_path(I, 1, 1, A, B, True)
    G = planted_gadget(5, 1); run_all(G, G["A"], G["B"], timeout=5)


def slim(r):
    return {k: v for k, v in r.items() if k != "front"}


def part_fam(a):
    warmup()
    n = a.n; rows = []; enum_to = 0; enum_tried = 0
    for k in range(a.instances):
        I = make_family_a(n, 1000 + k)
        F = front_all_layers(I)
        if F is None:
            continue
        A, B = budgets_from_front(F, 0.30)
        which = list(SOLVERS)
        skip_enum = (enum_tried >= 5 and enum_to == enum_tried)
        if skip_enum:
            which.remove("enum")
        res = run_all(I, A, B, which, timeout=a.timeout)
        if skip_enum:
            res["enum"] = dict(status="skipped_first5_timed_out", feasible=None)
        else:
            enum_tried += 1; enum_to += int(res["enum"]["status"] == "timeout")
        row = dict(n=n, inst=k, E=len(I["a_und"]), front_size=int(len(F[0])), A=A, B=B,
                   joint_binding=joint_binding(F, A, B), exact_feasible=res["label"]["feasible"],
                   solvers={s: slim(v) for s, v in res.items()})
        rows.append(row)
        print(n, k, {s: (v.get("feasible"), round(v.get("ms", 0), 1)) for s, v in res.items()}, flush=True)
    json.dump(dict(n=n, rows=rows, timeout=a.timeout), open(f"{a.outdir}/R3_fam_n{n}.json", "w"), indent=1)


def part_gadget(a):
    warmup()
    out = []; dead = set()
    for m in a.ms:
        n_inst = a.gadget_inst if m <= 20 else a.gadget_inst_big
        rows = []
        for k in range(n_inst):
            planted = (k % 2 == 0)
            xm = 10000 if a.xmode == "fixed" else 2 ** (m + 4)
            G = planted_gadget(m, 5000 + k, xmax=xm) if planted else make_gadget(m, 5000 + k, xmax=xm)
            which = [s for s in SOLVERS if s not in dead]
            res = run_all(G, G["A"], G["B"], which, timeout=a.timeout)
            rows.append(dict(m=m, inst=k, planted=planted, paths_total=2 ** m, solvers={s: slim(v) for s, v in res.items()}))
        # kill solvers that timed out / were capped on >= half the instances at this m
        for s in SOLVERS:
            bad = sum(1 for r in rows if s in r["solvers"] and r["solvers"][s]["status"] in ("timeout", "capped", "na_size"))
            if bad * 2 >= len(rows):
                dead.add(s)
        out.append(dict(m=m, rows=rows))
        med = lambda s, f: np.median([r["solvers"][s][f] for r in rows if s in r["solvers"] and f in r["solvers"][s]]) if any(s in r["solvers"] and f in r["solvers"][s] for r in rows) else float("nan")
        print(f"m={m:2d} enum {med('enum','ms'):9.1f}ms paths={med('enum','paths'):.3g} | label {med('label','ms'):9.1f}ms "
              f"created={med('label','created'):.3g} | dp {med('dp','ms'):9.1f}ms | ilp {med('ilp','ms'):9.1f}ms | "
              f"larac {med('larac','ms'):.2f} | fptas {med('fptas','ms'):.2f} | dead={sorted(dead)}", flush=True)
        json.dump(out, open(f"{a.outdir}/R3_gadget{'' if a.xmode == 'fixed' else '_exp'}.json", "w"), indent=1)


def part_naive(a):
    warmup()
    rows = []
    sizes = a.naive_sizes
    for n in sizes:
        for k in range(a.instances):
            I = make_family_a(n, 3000 + k)
            F = front_all_layers(I)
            if F is None:
                continue
            for tag, (A, B) in (("q30", budgets_from_front(F, .30)),):
                ex = solve_labels(I, A, B)
                feas_exact = ex["feasible"]
                rec = dict(n=n, inst=k, budget=tag, exact_feasible=feas_exact, front_size=int(len(F[0])), joint_binding=joint_binding(F, A, B))
                for lam in (0.1, 1.0, 10.0, 100.0):
                    ok, aa, bb, hh, _ = scalar_path(I, 1.0, lam)                      # plain scalarisation
                    ok2, a2, b2, h2, _ = scalar_path(I, 1.0, lam, A, B, True)         # budget-pruned scalarisation
                    rec[f"lam{lam:g}"] = dict(returned=ok, violating=bool(ok and (aa > A or bb > B)),
                                              pruned_found=bool(ok2 and a2 <= A and b2 <= B))
                lar = solve_larac(I, A, B); rec["larac_feasible"] = lar["feasible"]
                fp = solve_fptas(I, A, B); rec["fptas_feasible"] = fp["feasible"]
                rows.append(rec)
    json.dump(rows, open(f"{a.outdir}/R3_naive.json", "w"), indent=1)
    print("naive done", len(rows))


@nb.njit(cache=True)
def _dp3(n, fr, to, cl, H, reach):
    """reach[v,k1,k2,k3] : exists walk s->v with k_i links of class i (k1+k2+k3<=H).  Returns states visited."""
    vis = 0
    for tot in range(H):
        for k1 in range(tot + 1):
            for k2 in range(tot + 1 - k1):
                k3 = tot - k1 - k2
                for e in range(fr.shape[0]):
                    u = fr[e]
                    if reach[u, k1, k2, k3]:
                        vis += 1
                        c = cl[e]
                        if c == 0: reach[to[e], k1 + 1, k2, k3] = True
                        elif c == 1: reach[to[e], k1, k2 + 1, k3] = True
                        else: reach[to[e], k1, k2, k3 + 1] = True
    return vis


def part_dp3(a):
    """3-class count-vector DP: n in sizes, H=6 (operating point) and H=n-1 (O(|E| n^3)); verified against label-setting."""
    warmup()
    rng = np.random.default_rng(5)
    CA = np.array([200, 1000, 5000], np.int64); CB = np.array([100, 800, 4000], np.int64)
    rows = []
    for H_mode in ("H6", "Hn-1"):
        for n in a.dp3_sizes:
            H = 6 if H_mode == "H6" else n - 1
            if n * (H + 2) ** 3 > 2.5e8:
                rows.append(dict(H_mode=H_mode, n=n, status="memory_cap")); continue
            ts = []; agree = 0; tot = 0
            for k in range(a.dp3_graphs):
                r2 = np.random.default_rng(np.random.SeedSequence([n, k, 11]))
                iu, ju = np.triu_indices(n, 1); keep = r2.random(len(iu)) < 0.25
                u, v = iu[keep], ju[keep]
                cls = r2.choice(3, size=len(u), p=[.3, .5, .2])
                fr = np.concatenate([u, v]).astype(np.int32); to = np.concatenate([v, u]).astype(np.int32)
                cl = np.concatenate([cls, cls]).astype(np.int64)
                if len(u) == 0: continue
                reach = np.zeros((n, H + 2, H + 2, H + 2), np.bool_); reach[0, 0, 0, 0] = True
                t0 = time.perf_counter()
                vis = _dp3(n, fr, to, cl, H, reach)
                # decision: exists vector at d with k.a<=A,k.b<=B (budgets = 2nd-percentile style: median of all feasible vectors' costs)
                ks = np.argwhere(reach[n - 1]); 
                ms = (time.perf_counter() - t0) * 1000
                ts.append(ms)
                if len(ks):
                    costs_a = ks @ CA; costs_b = ks @ CB; ks = ks[(ks.sum(1) >= 1)]
                    if len(ks) == 0: continue
                    ca_ = ks @ CA; cb_ = ks @ CB
                    A = int(np.percentile(ca_, 30)); B = int(np.percentile(cb_, 30))
                    dp_feas = bool(np.any((ca_ <= A) & (cb_ <= B)))
                    I = dict(n=n, fr=fr, to=to, a=CA[cl], b=CB[cl], s=0, d=n - 1, H=H, L=1)
                    if H <= 8:
                        ex = solve_labels(I, A, B)
                        tot += 1; agree += int(ex["feasible"] == dp_feas)
            rows.append(dict(H_mode=H_mode, n=n, H=H, graphs=len(ts), median_ms=float(np.median(ts)),
                             q25=float(np.percentile(ts, 25)), q75=float(np.percentile(ts, 75)),
                             state_space=int(n * math.comb(H + 3, 3)), agree=agree, agree_checked=tot))
            print(rows[-1], flush=True)
    json.dump(rows, open(f"{a.outdir}/R3_dp3.json", "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["fam", "gadget", "naive", "dp3"])
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--instances", type=int, default=30)
    ap.add_argument("--timeout", type=float, default=TIMEOUT)
    ap.add_argument("--ms", type=int, nargs="+", default=[5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
    ap.add_argument("--xmode", choices=["fixed", "exp"], default="fixed")
    ap.add_argument("--gadget_inst", type=int, default=10)
    ap.add_argument("--gadget_inst_big", type=int, default=4)
    ap.add_argument("--naive_sizes", type=int, nargs="+", default=[20, 30, 50])
    ap.add_argument("--dp3_sizes", type=int, nargs="+", default=[10, 20, 30, 50, 75, 100, 150])
    ap.add_argument("--dp3_graphs", type=int, default=10)
    ap.add_argument("--outdir", default="../data")
    a = ap.parse_args()
    {"fam": part_fam, "gadget": part_gadget, "naive": part_naive, "dp3": part_dp3}[a.part](a)
