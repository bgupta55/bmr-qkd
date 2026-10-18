"""
Q4 — Collusion experiment (pending_v7 Q4)
=========================================
For each router, over random graphs and random (s,d) requests:
 (a) worst case: smallest relay coalition covering ALL links of the delivered route
     (exact DP, endpoints s,d honest); record P(min_coalition <= t).
 (b) random: every relay dishonest independently with p in {0.05,0.10,0.20};
     record P(all links of the route covered).
 Cost: hops, secure rate, sum w.
Topologies: ER (n=30), BA (n=30), ring (n=12), hub-and-spoke (n=12, 2 hubs).
The guarantee is structural (Corollary 1); ER/BA may show a small effect only,
ring / hub-and-spoke are where the floor binds.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import math, json
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bmr.certificate import build_certificate
from bmr.route import bmr_route

SEED = 404
TRIALS = 1500
T = 1
N_AUDIT = 10_000; DELTA = 1e-3; R_MIN = 0.01; PHI = 0.99; ELL = 10_000
P_LIST = [0.05, 0.10, 0.20]


def min_cover_internal(path):
    """Exact minimum number of INTERNAL nodes covering every edge of the path
    (s,d honest => cannot be used). Returns inf if impossible (m=1)."""
    m = len(path) - 1
    if m == 1:
        return math.inf
    # DP over internal nodes v1..v_{m-1}: choose subset C, every edge (v_{i-1},v_i) needs
    # v_{i-1} in C or v_i in C, with v_0 and v_m forced NOT in C.
    INF = 10 ** 9
    # state: whether previous node chosen
    dp = {0: 0, 1: INF}               # at v_0: not chosen
    for i in range(1, m + 1):
        nd = {0: INF, 1: INF}
        for prev, c in dp.items():
            if c >= INF:
                continue
            for cur in (0, 1):
                if i == m and cur == 1:
                    continue          # v_m = d honest
                if prev == 0 and cur == 0:
                    continue          # edge (v_{i-1},v_i) uncovered
                nd[cur] = min(nd[cur], c + cur)
        dp = nd
    return min(dp.values())


def make_graph(topo, rng):
    if topo == "ER":
        while True:
            G = nx.erdos_renyi_graph(30, 4 / 29, seed=int(rng.integers(1e9)))
            if nx.is_connected(G):
                return G
    if topo == "BA":
        return nx.barabasi_albert_graph(30, 2, seed=int(rng.integers(1e9)))
    if topo == "ring":
        return nx.cycle_graph(12)
    if topo == "hub":   # two hubs, each connected to all 10 leaves, hubs connected: 2-hop paths exist
        G = nx.Graph()
        G.add_edge("H1", "H2")
        for i in range(10):
            G.add_edge("H1", i); G.add_edge("H2", i)
        return nx.convert_node_labels_to_integers(G)
    raise ValueError(topo)


def certify(G, rng):
    for u, v in G.edges():
        S_hat = rng.normal(2.687, math.sqrt((4 - 2.687 ** 2 / 4) / N_AUDIT))
        c = build_certificate((u, v), float(S_hat), N_AUDIT, 1e-6, ELL, PHI, 0, 0)
        G[u][v]["r"] = c.r if c else 0.0
        G[u][v]["w"] = c.w if c else 2.0
        G[u][v]["T"] = 0
    return G


def r_spf(G, s, d):
    return nx.shortest_path(G, s, d)


def r_rate(G, s, d):
    H = nx.Graph()
    for u, v, a in G.edges(data=True):
        if a["r"] >= R_MIN:
            H.add_edge(u, v, w=1.0 / (a["r"] + 1e-9))
    return nx.dijkstra_path(H, s, d, weight="w")


def r_bmr(t):
    def f(G, s, d):
        res = bmr_route(G, s, d, R_MIN, DELTA, 0, 100, m_min=(2 * t + 1 if t > 0 else 0))
        return None if res is None else res[1]
    return f


ROUTERS = [("SPF", r_spf), ("Rate-Greedy", r_rate), ("BMR-NoFloor", r_bmr(0)), ("BMR(t=1)", r_bmr(T))]

rng = np.random.default_rng(SEED)
rows = []
print(f"{'topo':>6} {'router':>12} | {'P(cover<=t)%':>13} {'hops':>6} {'rate':>7} {'sum_w':>9} {'%2-hop':>7} | "
      + " ".join(f"P(cov|p={p})%" for p in P_LIST) + " | N")
for topo in ["ER", "BA", "ring", "hub"]:
    acc = {r: dict(cov=[], hops=[], rate=[], w=[], two=[], rnd={p: [] for p in P_LIST}) for r, _ in ROUTERS}
    for _ in range(TRIALS):
        G = certify(make_graph(topo, rng), rng)
        nodes = list(G.nodes())
        s, d = [int(x) for x in rng.choice(len(nodes), 2, replace=False)]
        s, d = nodes[s], nodes[d]
        # draw dishonest relay indicators once per instance (paired across routers)
        dis = {p: {v: bool(rng.random() < p) for v in nodes} for p in P_LIST}
        for rn, fn in ROUTERS:
            try:
                path = fn(G, s, d)
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                path = None
            if path is None:
                continue
            m = len(path) - 1
            a = acc[rn]
            a["cov"].append(min_cover_internal(path) <= T)
            a["hops"].append(m)
            a["rate"].append(min(G[u][v]["r"] for u, v in zip(path[:-1], path[1:])))
            a["w"].append(sum(G[u][v]["w"] for u, v in zip(path[:-1], path[1:])))
            a["two"].append(m == 2)
            for p in P_LIST:
                bad = {v for v in path[1:-1] if dis[p][v]}
                covered = m >= 2 and all((path[i] in bad) or (path[i + 1] in bad) for i in range(m))
                a["rnd"][p].append(covered)
    for rn, _ in ROUTERS:
        a = acc[rn]
        N = len(a["cov"])
        row = dict(topology=topo, router=rn, N=N,
                   worst_case_pct=100 * np.mean(a["cov"]), mean_hops=float(np.mean(a["hops"])),
                   pct_2hop=100 * np.mean(a["two"]), mean_rate=float(np.mean(a["rate"])),
                   mean_sum_w=float(np.mean(a["w"])),
                   random_cover_pct={str(p): 100 * np.mean(a["rnd"][p]) for p in P_LIST})
        rows.append(row)
        print(f"{topo:>6} {rn:>12} | {row['worst_case_pct']:>13.2f} {row['mean_hops']:>6.2f} {row['mean_rate']:>7.4f} "
              f"{row['mean_sum_w']:>9.2e} {row['pct_2hop']:>7.1f} | "
              + " ".join(f"{row['random_cover_pct'][str(p)]:>12.2f}" for p in P_LIST) + f" | {N}")
    print()

# ---- structural self-checks -------------------------------------------------
for r in rows:
    if r["router"] == "BMR(t=1)":
        assert r["worst_case_pct"] == 0.0, r          # Corollary 1
print("[self-check] BMR(t=1): P(min_coalition<=t) == 0 on every topology  (Corollary 1)")
# the exact DP agrees with ceil(m/2) for m>=2
for m in range(2, 13):
    assert min_cover_internal(list(range(m + 1))) == math.ceil(m / 2) == (m - 1) // 2 + 1
print("[self-check] exact cover DP == floor((m-1)/2)+1 for m=2..12")

# ---- figure -----------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
topos = ["ER", "BA", "ring", "hub"]
names = [r for r, _ in ROUTERS]
w = 0.2
for j, rn in enumerate(names):
    ax[0].bar(np.arange(4) + (j - 1.5) * w,
              [next(r for r in rows if r["topology"] == t and r["router"] == rn)["worst_case_pct"] for t in topos],
              w, label=rn, edgecolor='k', lw=.5)
    ax[1].bar(np.arange(4) + (j - 1.5) * w,
              [next(r for r in rows if r["topology"] == t and r["router"] == rn)["random_cover_pct"]["0.1"] for t in topos],
              w, label=rn, edgecolor='k', lw=.5)
for a_, tt in zip(ax, ["Worst case: P(min coalition ≤ t=1) [%]", "Random: P(route fully covered), p=0.10 [%]"]):
    a_.set_xticks(range(4)); a_.set_xticklabels(topos); a_.set_title(tt); a_.grid(alpha=.3, axis='y'); a_.legend(fontsize=8)
plt.tight_layout()
fp = os.path.join(os.path.dirname(__file__), "..", "figures", "Q4_collusion.png")
plt.savefig(fp, dpi=150); plt.close()
out = dict(experiment="Q4", T=T, trials=TRIALS, params=dict(n_audit=N_AUDIT, Delta=DELTA, phi=PHI, ell=ELL),
           rows=rows)
json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "data", "results", "q4_collusion.json"), "w"), indent=2)
print("saved", fp)
