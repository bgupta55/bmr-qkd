"""
bmr/route.py
============
Routing solvers with hop-floor support.
Developer Guide §7.1. Addresses Reviewer 2 Comment 3.

Solvers:
  S-A: bottleneck + Dijkstra (t = 0, standard)
  S-B: layered walk (non-simple, polynomial, for comparison)
  S-C: bmr_route (exact simple path with hop floor, using Yen's k-shortest)
  S-D: brute-force DFS enumerator (ground truth, n <= 12)
"""
import networkx as nx
import math
import heapq
import itertools


def build_filtered(G, R, now, W):
    """Filter graph to edges with rate >= R and freshness within W epochs."""
    H = nx.Graph()
    H.add_edges_from(
        (u, v, a) for u, v, a in G.edges(data=True)
        if a.get("r", 0) >= R and now - a.get("T", 0) <= W
    )
    return H


class SolverTimeout(Exception):
    """Raised when a wall-clock deadline is exceeded (Q9 uncapped-timing study)."""


# Budget comparison tolerance: identical to the evaluation code (q6_common.evaluate). It absorbs floating-point
# summation-order ties such as ten hops of 1e-6 against Delta = 1e-5 (14 Oct: N-3 self-check).
_TOL = 1e-15


def oracle(Gp, s, d, Delta, m_min=0, parity=None, deadline=None):
    """Return the min-weight feasible simple path, or None.

    If m_min=0 and parity=None: standard Dijkstra (always simple for positive weights).
    Otherwise: enumerate via Yen's k-shortest simple paths until hop floor is met or budget exceeded.
    """
    if s not in Gp or d not in Gp:
        return None
    try:
        if m_min == 0 and parity is None:
            dist, path = nx.single_source_dijkstra(Gp, s, d, weight="w")
            return (path, dist) if dist <= Delta + _TOL else None
        # Yen's algorithm: paths in non-decreasing order of weight
        import time as _t
        for path in nx.shortest_simple_paths(Gp, s, d, weight="w"):
            if deadline is not None and _t.perf_counter() > deadline:
                raise SolverTimeout()
            wt = sum(Gp[u][v]["w"] for u, v in zip(path[:-1], path[1:]))
            if wt > Delta + _TOL:
                return None  # all later paths are at least as heavy
            m = len(path) - 1  # hop count
            if m >= m_min and (parity is None or m % 2 == parity):
                return (path, wt)
    except nx.NetworkXNoPath:
        return None
    return None


def bmr_route(G, s, d, R, Delta, now, W, m_min=0, parity=None, deadline=None):
    """BMR routing: binary search over bottleneck rate, exact simple path.

    Returns (bottleneck_rate, path, weight) or None if no secure route exists.
    """
    if m_min > G.number_of_nodes() - 1:      # a simple path has at most N-1 hops: floor is impossible
        return None
    H = build_filtered(G, R, now, W)
    if s not in H or d not in H:
        return None
    rhos = sorted({a["r"] for _, _, a in H.edges(data=True)})
    lo, hi, best = 0, len(rhos) - 1, None
    while lo <= hi:
        mid = (lo + hi) // 2
        Gp = nx.Graph()
        Gp.add_edges_from(
            (u, v, a) for u, v, a in H.edges(data=True) if a["r"] >= rhos[mid]
        )
        res = oracle(Gp, s, d, Delta, m_min, parity, deadline)
        if res:
            best = (rhos[mid], *res)
            lo = mid + 1
        else:
            hi = mid - 1
    return best  # (bottleneck_rate, path, weight) or None


def min_weight_walk_hop_floor(G, s, d, K):
    """Dijkstra on (node, min(hops, K)). Polynomial; may return a walk with repeats.
    Used to demonstrate that walks don't inherit the collusion-resistance bound (E4b).
    """
    INF = math.inf
    dist = {(s, 0): 0.0}
    prev = {}
    pq = [(0.0, s, 0)]
    while pq:
        dc, v, h = heapq.heappop(pq)
        if dist.get((v, h), INF) < dc:
            continue
        if v == d and h == K:
            walk, st = [], (v, h)
            while True:
                walk.append(st[0])
                if st == (s, 0):
                    break
                st = prev[st]
            return walk[::-1], dc
        for u, a in G[v].items():
            nh = min(h + 1, K)
            nd = dc + a.get("w", 1.0)
            if nd < dist.get((u, nh), INF):
                dist[(u, nh)] = nd
                prev[(u, nh)] = (v, h)
                heapq.heappush(pq, (nd, u, nh))
    return None, INF


def is_simple(walk):
    """Check if a walk is a simple path (no repeated nodes)."""
    return len(set(walk)) == len(walk)


def brute_force_simple_paths(G, s, d, Delta, m_min=0, timeout=10000):
    """Enumerate all simple paths by DFS. Ground truth for n <= 12.
    Returns the minimum-weight feasible (path, weight) or None.
    (v7 fix: the previous version compared against an immutable `best_w`
    and so returned the LAST feasible path rather than the best one.)
    """
    count = [0]
    best = [None, None]
    best_w_ref = [math.inf]

    def dfs(v, path, w):
        count[0] += 1
        if count[0] > timeout:
            return
        if v == d:
            m = len(path) - 1
            if m >= m_min and w <= Delta + _TOL and w < best_w_ref[0]:
                best[0], best[1] = list(path), w
                best_w_ref[0] = w
            return
        for u in G[v]:
            if u not in path_set:
                path_set.add(u)
                path.append(u)
                ew = G[v][u].get("w", 1.0)
                if w + ew <= Delta + _TOL:
                    dfs(u, path, w + ew)
                path.pop()
                path_set.discard(u)

    path_set = {s}
    dfs(s, [s], 0.0)
    if best[0] is None:
        return None
    return (best[0], best[1])


def brute_force_max_rate(G, s, d, Delta, m_min=0, node_limit=5_000_000):
    """Exhaustive ground truth for the BMR objective (Problem 1 of the paper):
    maximise bottleneck rate over simple s-d paths with sum(w) <= Delta and
    hops >= m_min.  Returns (rate, path, weight, truncated) or
    (None, None, None, truncated).  `truncated` is True if the node limit was
    hit (instance then excluded from the validation count)."""
    count = [0]
    best = [-1.0, None, None]
    trunc = [False]

    def dfs(v, path, w, rmin):
        count[0] += 1
        if count[0] > node_limit:
            trunc[0] = True
            return
        if v == d:
            if len(path) - 1 >= m_min and w <= Delta + _TOL and rmin > best[0]:
                best[0], best[1], best[2] = rmin, list(path), w
            return
        for u, a in G[v].items():
            if u in path_set:
                continue
            nw = w + a.get("w", 1.0)
            if nw > Delta + _TOL or min(rmin, a.get("r", 0.0)) <= best[0]:
                continue  # prune: cannot beat incumbent (or infeasible)
            path_set.add(u)
            path.append(u)
            dfs(u, path, nw, min(rmin, a.get("r", 0.0)))
            path.pop()
            path_set.discard(u)

    path_set = {s}
    dfs(s, [s], 0.0, math.inf)
    if best[1] is None:
        return None, None, None, trunc[0]
    return best[0], best[1], best[2], trunc[0]
