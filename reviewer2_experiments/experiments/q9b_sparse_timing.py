"""Q9b: solver timing on SPARSE graphs (mean degree d in {4,6}); same pipeline/protocol as Q9(b)."""
import sys, os, re, json, time
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np, networkx as nx
from bmr.route import bmr_route, SolverTimeout
from bmr.certificate import build_certificate

def make_instance(n, rng, p, S_mean=2.687, delta=1e-6, ell=10000, phi=0.99, n_audit=2500, Delta=0.5):
    while True:
        G = nx.erdos_renyi_graph(n, p, seed=int(rng.integers(1e9)))
        if nx.is_connected(G): break
    for u, v in G.edges():
        c = build_certificate((u, v), rng.normal(S_mean, 0.02), n_audit, delta, ell, phi, 0, 0)
        G[u][v]["r"] = c.r if c else 0.0
        G[u][v]["w"] = c.w if c else 2.0
        G[u][v]["T"] = 0
    return G, 0, n - 1, Delta

SIZES = [10, 20, 30, 50, 75, 100, 150]; TV = [0, 1, 2, 3, 5]; N = 100; TO = 1.0
rng = np.random.default_rng(1234); rows = []
for d in (4, 6):
    for n in SIZES:
        p = min(1.0, d / (n - 1))
        for t in TV:
            ms = []; to = 0; feas = 0
            for _ in range(N):
                G, s, dd, Delta = make_instance(n, rng, p)
                t0 = time.perf_counter()
                try:
                    res = bmr_route(G, s, dd, 0.0, Delta, 0, 1000, m_min=2*t+1, deadline=t0 + TO)
                    if res is not None: feas += 1
                except SolverTimeout:
                    to += 1
                ms.append((time.perf_counter() - t0) * 1000)
            rows.append(dict(mean_deg=d, n=n, t=t, median_ms=float(np.median(ms)), p95_ms=float(np.percentile(ms, 95)),
                             max_ms=float(max(ms)), timeout_pct=100*to/N, feasible_pct=100*feas/N))
            r = rows[-1]; print(f"d={d} n={n:>3} t={t} med={r['median_ms']:.3f} p95={r['p95_ms']:.3f} max={r['max_ms']:.3f} to={r['timeout_pct']:.0f}% feas={r['feasible_pct']:.0f}%", flush=True)
json.dump(rows, open("data/results/q9b_sparse_timing.json", "w"), indent=1)
