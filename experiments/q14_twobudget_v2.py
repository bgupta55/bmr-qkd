"""
Q14: two-budget problem (Theorem 4) -- extended E11.
Problem (Theorem 4): given two additive path weights (w1 = eps+delta budget weight, w2 = second
budget weight, e.g. gamma-type), find an s-d path with W1<=B1 and W2<=B2.  Weakly NP-hard (RSP).
Because weights are non-negative, any s-d WALK satisfying both budgets can be shortcut to a simple
path that still satisfies them, so walk-based solvers decide the same feasibility question.

Solvers:  PARETO  exact label-setting;  DP  exact pseudo-polynomial (integer w1);
          FPTAS-style  scaled DP (w1 relaxed by <= eps, w2 budget kept);
          LAGR  Lagrangian-relaxation heuristic (Dijkstra on w1 + lam*w2, bisection on lam).
Families: random ER (avg deg 6) and BA(m=2), n in {10,20,30,50,100,150}, uncorrelated and
          ANTI-correlated weights;  diamond chain (hard family), k<=14, loose and tight budgets.
No number below is hand-typed; the JSON is the source.
"""
import sys, os, json, time, heapq, hashlib, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, networkx as nx
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print("code sha256:", hashlib.sha256(open(__file__, "rb").read()).hexdigest()[:16])
rng = np.random.default_rng(14)
MAXL = 300_000


def to_di(G):
    D = nx.DiGraph()
    for u, v, d in G.edges(data=True):
        D.add_edge(u, v, **d); D.add_edge(v, u, **d)
    return D


def pareto(D, s, t, B1, B2, maxl=MAXL):
    lab = {v: [] for v in D}; lab[s] = [(0, 0.0)]
    pq = [(0, 0.0, s)]; tot = 1
    while pq:
        a, b, v = heapq.heappop(pq)
        if (a, b) not in lab[v]:
            continue
        for u in D[v]:
            na, nb = a + D[v][u]["w1"], b + D[v][u]["w2"]
            if na > B1 or nb > B2 + 1e-12:
                continue
            if any(x <= na and y <= nb for x, y in lab[u]):
                continue
            lab[u] = [(x, y) for x, y in lab[u] if not (na <= x and nb <= y)] + [(na, nb)]
            tot += 1
            if tot > maxl:
                return None, tot, True
            heapq.heappush(pq, (na, nb, u))
    return len(lab[t]) > 0, tot, False


def dp_feasible(D, s, t, B1, B2, scale=None, cap=None):
    """Layered DP over integer w1-cost c (<=cap): dp[c][v]=min w2.  scale: function mapping w1->int."""
    nodes = list(D); idx = {v: i for i, v in enumerate(nodes)}; n = len(nodes)
    eu = np.array([idx[u] for u, v in D.edges()]); ev = np.array([idx[v] for u, v in D.edges()])
    ea = np.array([scale(D[u][v]["w1"]) if scale else int(D[u][v]["w1"]) for u, v in D.edges()])
    ew = np.array([D[u][v]["w2"] for u, v in D.edges()])
    cap = int(cap)
    dp = np.full((cap + 1, n), np.inf); dp[0, idx[s]] = 0.0
    z = ea == 0
    for c in range(0, cap + 1):
        ok = (ea <= c) & (ea > 0)
        if ok.any():
            cand = dp[c - ea[ok], eu[ok]] + ew[ok]
            np.minimum.at(dp[c], ev[ok], cand)
        if z.any():                       # zero-w1 edges: relax inside the layer until stable (<= n passes)
            for _ in range(n):
                before = dp[c].copy()
                np.minimum.at(dp[c], ev[z], dp[c, eu[z]] + ew[z])
                if np.array_equal(before, dp[c]):
                    break
    return bool((dp[:, idx[t]] <= B2 + 1e-12).any()), dp.size


def lagr(D, s, t, B1, B2, iters=40):
    lo, hi = 1e-6, 1e6
    seen_feasible = False
    for _ in range(iters):
        lam = math.sqrt(lo * hi)
        for u, v in D.edges():
            D[u][v]["wl"] = D[u][v]["w1"] + lam * D[u][v]["w2"]
        try:
            p = nx.dijkstra_path(D, s, t, weight="wl")
        except nx.NetworkXNoPath:
            return False
        a = sum(D[x][y]["w1"] for x, y in zip(p[:-1], p[1:])); b = sum(D[x][y]["w2"] for x, y in zip(p[:-1], p[1:]))
        if a <= B1 and b <= B2:
            return True
        if a > B1: hi = lam      # w1 too large -> penalise w1 more (smaller lam)
        else: lo = lam           # w2 too large -> penalise w2 more (larger lam)
    return False


def make_instance(kind, n, corr):
    G = nx.erdos_renyi_graph(n, min(1.0, 6.0 / (n - 1)), seed=int(rng.integers(1e9))) if kind == "ER" \
        else nx.barabasi_albert_graph(n, 2, seed=int(rng.integers(1e9)))
    if not nx.is_connected(G):
        comp = max(nx.connected_components(G), key=len); G = nx.convert_node_labels_to_integers(G.subgraph(comp))
    for u, v in G.edges():
        a = int(rng.integers(1, 21))
        b = float(rng.uniform(.02, 1.0)) if corr == "uncorr" else max(.02, (21 - a) / 20 * float(rng.uniform(.8, 1.2)))
        G[u][v]["w1"], G[u][v]["w2"] = a, b
    D = to_di(G); nodes = list(G); s, t = nodes[0], nodes[-1]
    p1 = nx.dijkstra_path(D, s, t, weight="w1"); p2 = nx.dijkstra_path(D, s, t, weight="w2")
    W = lambda p, k: sum(D[x][y][k] for x, y in zip(p[:-1], p[1:]))
    th = 0.35
    B1 = int(math.floor(W(p1, "w1") + th * (W(p2, "w1") - W(p1, "w1"))))
    B2 = W(p2, "w2") + th * (W(p1, "w2") - W(p2, "w2"))
    return D, s, t, B1, B2


rows = []
print(f"{'fam':>4} {'corr':>7} {'n':>4} | feas% | PARETO ms (to%) | DP ms | FPTAS(.1) ms viol | LAGR ms miss%")
for kind in ("ER", "BA"):
    for corr in ("uncorr", "anti"):
        for n in (10, 20, 30, 50, 100, 150):
            NI = 40 if n >= 100 else 100
            R = dict(feas=0, tp=[], td=[], tf=[], tl=[], to=0, lag_miss=0, fp_miss=0, fp_viol=[], disagree=0, labels=[])
            for _ in range(NI):
                D, s, t, B1, B2 = make_instance(kind, n, corr)
                t0 = time.perf_counter(); fp, nl, to = pareto(D, s, t, B1, B2); R["tp"].append(time.perf_counter() - t0); R["labels"].append(nl)
                t0 = time.perf_counter(); fd, _ = dp_feasible(D, s, t, B1, B2, cap=B1); R["td"].append(time.perf_counter() - t0)
                eps = 0.1; K = int(math.ceil(len(D) / eps))
                t0 = time.perf_counter()
                ff, _ = dp_feasible(D, s, t, B1, B2, scale=lambda w: int(math.ceil(w * K / max(B1, 1))), cap=K + len(D))
                R["tf"].append(time.perf_counter() - t0)
                t0 = time.perf_counter(); fl = lagr(D, s, t, B1, B2); R["tl"].append(time.perf_counter() - t0)
                if to: R["to"] += 1; fp = fd                      # take exact DP as the truth on Pareto timeouts
                if fp != fd: print('   DISAGREE', kind, corr, n, 'pareto', fp, 'dp', fd, 'B', B1, B2, flush=True)
                R["disagree"] += (fp != fd)
                R["feas"] += fd
                if fd and not ff: R["fp_miss"] += 1              # FPTAS must never miss a feasible instance
                if fd and not fl: R["lag_miss"] += 1
            rr = dict(kind=kind, corr=corr, n=n, inst=NI, feasible=R["feas"] / NI,
                      pareto_ms=1e3 * float(np.median(R["tp"])), pareto_p95_ms=1e3 * float(np.percentile(R["tp"], 95)),
                      pareto_timeout=R["to"] / NI, med_labels=float(np.median(R["labels"])),
                      dp_ms=1e3 * float(np.median(R["td"])), fptas_ms=1e3 * float(np.median(R["tf"])),
                      fptas_missed=R["fp_miss"], lagr_ms=1e3 * float(np.median(R["tl"])),
                      lagr_miss_rate=(R["lag_miss"] / R["feas"]) if R["feas"] else None,
                      pareto_vs_dp_disagree=R["disagree"])
            rows.append(rr)
            print(f"{kind:>4} {corr:>7} {n:>4} | {rr['feasible']*100:4.0f}% | {rr['pareto_ms']:7.2f} ({rr['pareto_timeout']*100:.0f}) | {rr['dp_ms']:6.2f} | {rr['fptas_ms']:6.2f} miss={rr['fptas_missed']} | {rr['lagr_ms']:6.2f} {100*(rr['lagr_miss_rate'] or 0):4.1f}%  [P/DP disagree {rr['pareto_vs_dp_disagree']}]", flush=True)

# ---------------- FPTAS violation measurement: actual relaxation actually used ----------------
# For feasible instances found ONLY thanks to relaxation (infeasible exactly, feasible at (1+eps)B1) the scheme may return YES.
rel = {"fptas_yes_exactly_infeasible": 0, "n_exactly_infeasible": 0}
for _ in range(300):
    D, s, t, B1, B2 = make_instance("ER", 30, "anti")
    fd, _ = dp_feasible(D, s, t, B1, B2, cap=B1)
    if not fd:
        rel["n_exactly_infeasible"] += 1
        eps = 0.1; K = int(math.ceil(len(D) / eps))
        ff, _ = dp_feasible(D, s, t, B1, B2, scale=lambda w: int(math.ceil(w * K / max(B1, 1))), cap=K + len(D))
        rel["fptas_yes_exactly_infeasible"] += ff
        # check that such YES is certified by a path with W1 <= (1+eps)B1
        fd2, _ = dp_feasible(D, s, t, math.floor((1 + eps) * B1), B2, cap=math.floor((1 + eps) * B1)) if ff else (True, 0)
        assert fd2, "FPTAS YES not certified by a (1+eps)B1 path"
print("\nFPTAS on exactly-infeasible ER-30 anti-correlated:", rel)

# ---------------- diamond chain (hard family) ----------------
def diamond(k, tight):
    D = nx.DiGraph(); mid = lambda i: ("m", i)
    for i in range(k):
        a, b = ("n", i), ("n", i + 1)
        D.add_edge(a, mid(i), w1=2 ** i, w2=0.0); D.add_edge(mid(i), b, w1=0, w2=0.0)
        D.add_edge(a, b, w1=0, w2=float(2 ** i))
    B1 = (2 ** (k - 1) - 1) if tight else 2 ** k; B2 = float(2 ** (k - 1)) if tight else float(2 ** k)
    return D, ("n", 0), ("n", k), B1, B2
drows = []
print("\nDiamond chain: k | PARETO labels(total) ms | DP table ms | FPTAS ms | LAGR finds feasible?  (tight budgets: B1+B2 = 2^k-1)")
for k in range(2, 15):
    D, s, t, B1, B2 = diamond(k, True)
    t0 = time.perf_counter(); fp, nl, to = pareto(D, s, t, B1, B2, maxl=2_000_000); tp = time.perf_counter() - t0
    t0 = time.perf_counter(); fd, size = dp_feasible(D, s, t, B1, B2, cap=B1); td = time.perf_counter() - t0
    eps = 0.1; K = int(math.ceil(len(D) / eps))
    t0 = time.perf_counter(); ff, fsize = dp_feasible(D, s, t, B1, B2, scale=lambda w: int(math.ceil(w * K / max(B1, 1))), cap=K + len(D)); tf = time.perf_counter() - t0
    fl = lagr(D, s, t, B1, B2)
    drows.append(dict(k=k, B1=B1, B2=B2, pareto_labels=nl, pareto_ms=tp * 1e3, pareto_timeout=to, pareto_feasible=fp,
                      dp_cells=size, dp_ms=td * 1e3, dp_feasible=fd, fptas_cells=fsize, fptas_ms=tf * 1e3, fptas_feasible=ff, lagr_feasible=fl))
    print(f"  k={k:2d}: labels={nl:8d} {tp*1e3:9.1f}ms{' TIMEOUT' if to else ''} | DP cells={size:9d} {td*1e3:8.1f}ms feas={fd} | FPTAS cells={fsize:6d} {tf*1e3:6.1f}ms feas={ff} | LAGR feas={fl}", flush=True)

fig, ax = plt.subplots(1, 3, figsize=(17, 4.5))
for (kind, corr), st in {("ER", "uncorr"): "o-", ("ER", "anti"): "s-", ("BA", "uncorr"): "o--", ("BA", "anti"): "s--"}.items():
    sub = [r for r in rows if r["kind"] == kind and r["corr"] == corr]
    ax[0].plot([r["n"] for r in sub], [r["pareto_ms"] for r in sub], st, label=f"Pareto {kind}/{corr}")
    ax[1].plot([r["n"] for r in sub], [100 * (r["lagr_miss_rate"] or 0) for r in sub], st, label=f"{kind}/{corr}")
ax[0].set_yscale("log"); ax[0].set_xlabel("n"); ax[0].set_ylabel("median ms"); ax[0].set_title("Random graphs: exact Pareto label-setting"); ax[0].legend(fontsize=7)
ax[1].set_xlabel("n"); ax[1].set_ylabel("% feasible instances missed"); ax[1].set_title("Lagrangian heuristic: misses of feasible instances"); ax[1].legend(fontsize=7)
ks = [r["k"] for r in drows]
ax[2].semilogy(ks, [r["pareto_labels"] for r in drows], "ro-", label="Pareto labels"); ax[2].semilogy(ks, [2 ** k for k in ks], "k--", label=r"$2^k$")
ax[2].semilogy(ks, [r["fptas_cells"] for r in drows], "g^-", label="FPTAS DP cells (eps=0.1)")
ax[2].set_xlabel("k diamonds"); ax[2].set_title("Diamond chain (tight budgets): worst case"); ax[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(ROOT, "figures", "Q14_twobudget.png"), dpi=140)
json.dump({"random": rows, "fptas_relaxation_check": rel, "diamond": drows, "maxl_pareto_random": MAXL,
           "eps": 0.1, "theta_budget_tightness": 0.35}, open(os.path.join(ROOT, "data", "results", "q14_twobudget_v2.json"), "w"), indent=1)
