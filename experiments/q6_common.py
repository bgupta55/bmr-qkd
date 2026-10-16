"""
experiments/q6_common.py
========================
Shared code for the Q6 full sweep and the dynamic scenario.

Routers (all see the SAME certificate database produced by bmr.sim.LinkSim -> build_certificate):
  B0  SPF            hop-count shortest path (connectivity only)
  B1a P_lambda(0.01) Pirandola et al., PRR 8, 033221 (2026), Eq. (35), implemented AS PUBLISHED
  B1b P_lambda(1)    (hybrid rate/security Dijkstra "lambda protocol"; heuristic scalarisation;
                      no hard budget; uses (R_uv, eps_uv) = (certified rate, certified eps))
  B2  BMR-NoFloor    BMR with t = 0 (budget + bottleneck, no hop floor) = exact optimum without floor
  B3  Zhang-style    widest (max-min certified rate) path, ties -> fewest hops. RWA dropped.
  B5  Rate-greedy    Dijkstra on sum(1/r_e)  (own definition; no budget/expiry)
  B6  Q-PATH-style   Dijkstra on -log F_e   (objective-only; NOT Q-PATH)
  B7  Q-CAST-style   min hops s.t. prod F_e >= 0.75 (objective-only; NOT Q-CAST)
  B8  QBER-monitor   QBER_ub = (1 - S_hat/2sqrt2)/2 + eta/(4 sqrt2); links with QBER_ub >= 0.11 pruned;
                     widest path on 1 - 2 h2(QBER_ub) (same S_hat, same delta -> identical information)
  BMR                budget + hop floor (t = 1 -> m >= 3) + bottleneck, exact simple path

Evaluation (Definitions in pending_v7 Q1):
  feasible[i]     = exact solver (validated vs brute force, Q9) finds a route with sum_w <= Delta, hops >= 2t+1,
                    every link rate >= R_MIN     (instance property, identical for all routers)
  found[i,r]      = router returned a route
  compliant[i,r]  = found and sum_w <= Delta and hops >= 2t+1 and every link rate >= R_MIN (same as the oracle)
  secure[i,r]     = min_e r_e if compliant else 0       (mean over ALL trials)
"""
import math, heapq
import numpy as np
import networkx as nx

from bmr.route import bmr_route
from bmr.sim import LinkSim, S_HONEST, S_COMPROMISED

R_MIN = 0.01
W_FRESH = 10
T_COLL = 1
M_FLOOR = 2 * T_COLL + 1
Q_TH = 0.11
F_TH = 0.75
SQ2 = 2.0 * math.sqrt(2.0)


def h2(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


# ------------------------------------------------------------------ graphs
def make_topology(kind, n, seed):
    if kind == "ER":
        return nx.erdos_renyi_graph(n, min(1.0, 4.0 / (n - 1)), seed=seed)
    if kind == "BA":
        return nx.barabasi_albert_graph(n, 2, seed=seed)
    if kind == "WS":
        return nx.connected_watts_strogatz_graph(n, 4, 0.2, tries=50, seed=seed)
    if kind == "grid":
        rows = {10: 2, 20: 4, 30: 5, 50: 5}[n]
        G = nx.grid_2d_graph(rows, n // rows)
        return nx.convert_node_labels_to_integers(G)
    if kind == "ring":
        return nx.cycle_graph(n)
    if kind == "line":                      # added for the degenerate-topology sweep (q6_sweep_degenerate.py)
        return nx.path_graph(n)
    if kind == "star":
        return nx.star_graph(n - 1)
    raise ValueError(kind)


def load_sndlib(path):
    """Parse a SNDlib native file -> undirected nx.Graph with integer labels (names in graph attr)."""
    import re
    txt = open(path, encoding="utf-8", errors="ignore").read()
    nodes_blk = re.search(r"NODES \((.*?)\n\)", txt, re.S).group(1)
    names = [ln.split()[0] for ln in nodes_blk.strip().splitlines()
             if ln.strip() and not ln.strip().startswith("#")]
    idx = {nm: i for i, nm in enumerate(names)}
    links_blk = re.search(r"LINKS \((.*?)\n\)", txt, re.S).group(1)
    G = nx.Graph()
    G.add_nodes_from(range(len(names)))
    for ln in links_blk.strip().splitlines():
        m = re.match(r"\s*(\S+)\s*\(\s*(\S+)\s+(\S+)\s*\)", ln)
        if m:
            G.add_edge(idx[m.group(2)], idx[m.group(3)])
    G.graph["names"] = names
    return G


# ------------------------------------------------------------------ certificate database
def attach_certs(G, sim, epoch=0):
    """Publish a certificate for every link through the simulator (the only code path)."""
    for u, v in G.edges():
        c = sim.publish_certificate((u, v), epoch)
        a = G[u][v]
        s_hat = sim.last_S_hat[LinkSim._key((u, v))]
        a["cert"] = c
        a["S_hat"] = s_hat
        a["r"] = c.r if c else 0.0
        a["eps"] = c.eps if c else 1.0
        a["w"] = c.w if c else 2.0
        a["T"] = epoch
        a["F"] = (1.0 + 3.0 * min(max(s_hat, 0.0) / SQ2, 1.0)) / 4.0
        # B8 evidence: QBER measured DIRECTLY on 4n sifted rounds (trusted-device model, the most
        # favourable treatment for the baseline: a direct binomial estimate with a one-sided
        # Hoeffding radius at the same delta; it does not go through the looser CHSH radius).
        S_true = sim.true_S[LinkSim._key((u, v))]
        q_true = (1.0 - min(S_true, SQ2) / SQ2) / 2.0
        n4 = 4 * sim.n
        q_hat = sim.rng.binomial(n4, q_true) / n4
        a["q_ub"] = q_hat + math.sqrt(math.log(1.0 / sim.delta) / (2.0 * n4))
        a["r_q"] = max(0.0, 1.0 - 2.0 * h2(a["q_ub"])) if a["q_ub"] < Q_TH else 0.0


# ------------------------------------------------------------------ helpers
def _edge_attr(G, path, key):
    return [G[u][v][key] for u, v in zip(path[:-1], path[1:])]


def _widest_then_short(G, s, d, key, floor):
    """max-min `key` path, ties -> fewest hops (bottleneck threshold search + BFS)."""
    vals = sorted({a[key] for _, _, a in G.edges(data=True) if a[key] >= floor and a[key] > 0})
    lo, hi, best = 0, len(vals) - 1, None
    while lo <= hi:
        mid = (lo + hi) // 2
        H = nx.Graph()
        H.add_edges_from((u, v) for u, v, a in G.edges(data=True) if a[key] >= vals[mid])
        if s in H and d in H and nx.has_path(H, s, d):
            best = nx.shortest_path(H, s, d)
            lo = mid + 1
        else:
            hi = mid - 1
    return best


# ------------------------------------------------------------------ routers: return path or None
def r_spf(G, s, d, Delta):
    try:
        return nx.shortest_path(G, s, d)
    except nx.NetworkXNoPath:
        return None


def make_lambda(lam):
    def r_lambda(G, s, d, Delta):
        """Pirandola et al. Eq.(35): T_v <- min{T_v, max(T_u, 1/R_uv) + lam [T_u - log10(1 - eps_uv)]}."""
        INF = float("inf")
        T = {v: INF for v in G}
        T[s] = 0.0
        prev = {}
        settled = set()
        pq = [(0.0, s)]
        while pq:
            tu, u = heapq.heappop(pq)
            if u in settled:
                continue
            settled.add(u)
            if u == d:
                break
            for v, a in G[u].items():
                if v in settled or a["r"] < R_MIN:
                    continue
                eps = min(max(a["eps"], 0.0), 1.0 - 1e-15)
                cand = max(tu, 1.0 / a["r"]) + lam * (tu - math.log10(1.0 - eps))
                if cand < T[v]:
                    T[v] = cand
                    prev[v] = u
                    heapq.heappush(pq, (cand, v))
        if d not in prev:
            return None
        p = [d]
        while p[-1] != s:
            p.append(prev[p[-1]])
        return p[::-1]
    return r_lambda


def r_b3(G, s, d, Delta):
    return _widest_then_short(G, s, d, "r", R_MIN)


def r_b5(G, s, d, Delta):
    H = nx.Graph()
    for u, v, a in G.edges(data=True):
        if a["r"] >= R_MIN:
            H.add_edge(u, v, c=1.0 / a["r"])
    try:
        return nx.shortest_path(H, s, d, weight="c")
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None


def r_b6(G, s, d, Delta):
    H = nx.Graph()
    for u, v, a in G.edges(data=True):
        H.add_edge(u, v, c=-math.log(max(a["F"], 1e-12)))
    try:
        return nx.shortest_path(H, s, d, weight="c")
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None


def r_b7(G, s, d, Delta, cap=300):
    """min hops s.t. prod F >= 0.75 : enumerate simple paths by hop count (Yen), first feasible."""
    try:
        for i, p in enumerate(nx.shortest_simple_paths(G, s, d)):
            if math.prod(_edge_attr(G, p, "F")) >= F_TH:
                return p
            if i >= cap:
                return None
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None
    return None


def r_b8(G, s, d, Delta):
    return _widest_then_short(G, s, d, "r_q", 1e-9)


def r_bmr_nf(G, s, d, Delta):
    res = bmr_route(G, s, d, R_MIN, Delta, 0, W_FRESH, m_min=0)
    return res[1] if res else None


def r_bmr(G, s, d, Delta):
    res = bmr_route(G, s, d, R_MIN, Delta, 0, W_FRESH, m_min=M_FLOOR)
    return res[1] if res else None


ROUTERS = [
    ("B0-SPF", r_spf),
    ("B1a-Plambda(0.01)", make_lambda(0.01)),
    ("B1b-Plambda(1)", make_lambda(1.0)),
    ("B2-BMR-NoFloor", r_bmr_nf),
    ("B3-Zhang-style", r_b3),
    ("B5-Rate-greedy", r_b5),
    ("B6-Qpath-style", r_b6),
    ("B7-Qcast-style", r_b7),
    ("B8-QBER-monitor", r_b8),
    ("BMR", r_bmr),
]
ROUTER_NAMES = [n for n, _ in ROUTERS]


# ------------------------------------------------------------------ evaluation
def evaluate(G, path, Delta):
    """-> (found, compliant, secure_rate, hops, sum_w, compliant_budget_only, secure_rate_budget_only).
    compliant = every link rate >= R_MIN AND sum_w<=Delta AND hops>=2t+1 ; *_budget_only ignores the hop floor
    (so floor-less routers are not charged for the collusion property they do not claim)."""
    if path is None:
        return 0, 0, 0.0, 0, 0.0, 0, 0.0
    rs = _edge_attr(G, path, "r")
    ws = sum(_edge_attr(G, path, "w"))
    hops = len(path) - 1
    certified = all(r >= R_MIN for r in rs)          # same request floor as the feasibility oracle (14 Oct: N-3)
    comp = certified and ws <= Delta + 1e-15 and hops >= M_FLOOR
    comp_b = certified and ws <= Delta + 1e-15          # budget-only (floor ignored)
    return 1, int(comp), (min(rs) if comp else 0.0), hops, ws, int(comp_b), (min(rs) if comp_b else 0.0)


def min_cover(path):
    """exact min #relay nodes covering all links of a simple path = ceil(m/2) (Corollary 1)."""
    return (len(path) - 1 + 1) // 2


def feasible_oracle(G, s, d, Delta):
    res = bmr_route(G, s, d, R_MIN, Delta, 0, W_FRESH, m_min=M_FLOOR)
    return res is not None
