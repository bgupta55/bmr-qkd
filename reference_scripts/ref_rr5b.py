"""
Reference check for RR5b: three exact solvers for
   min  sum w(P)   over simple s-d paths P with 2t+1 <= |P| <= H, using only links with r >= rho
(the inner feasibility step of Algorithm 2; the outer binary search over rho is exercised too).
 (1) brute force enumeration of all simple paths
 (2) Yen-style enumeration in non-decreasing weight (networkx shortest_simple_paths)
 (3) Theorem 3d: enumerate simple prefixes of exactly L=2t+1 edges, delete prefix vertices except the
     last, hop-layered Bellman-Ford for min-weight walk with <= H-L edges (shortcut => simple path).
"""
import math, random, itertools, sys
import networkx as nx

INF = math.inf


def brute(G, s, d, L, H, rho):
    best = INF
    for p in nx.all_simple_paths(G, s, d, cutoff=H):
        m = len(p) - 1
        if m < L:
            continue
        ls = [G[p[i]][p[i + 1]] for i in range(m)]
        if min(l["r"] for l in ls) < rho:
            continue
        best = min(best, sum(l["w"] for l in ls))
    return best


def yen(G, s, d, L, H, rho):
    Gr = nx.Graph()
    Gr.add_nodes_from(G.nodes())
    for a, b, dat in G.edges(data=True):
        if dat["r"] >= rho:
            Gr.add_edge(a, b, w=dat["w"] + 0.0)
    if s not in Gr or d not in Gr or not nx.has_path(Gr, s, d):
        return INF
    for p in nx.shortest_simple_paths(Gr, s, d, weight="w"):
        m = len(p) - 1
        if m < L:
            continue
        if m > H:
            # Yen order is by weight, not hops: cannot stop; keep scanning (exact but exponential)
            continue
        return sum(Gr[p[i]][p[i + 1]]["w"] for i in range(m))
    return INF


def prefix_bf(G, s, d, L, H, rho):
    adj = {u: [(v, G[u][v]["w"]) for v in G[u] if G[u][v]["r"] >= rho] for u in G}
    best = INF

    def rec(path, wsum):
        nonlocal best
        u = path[-1]
        if len(path) - 1 == L:
            if u == d:
                best = min(best, wsum)
                return
            banned = set(path[:-1])
            # hop-layered Bellman-Ford from u to d, <= H-L edges, avoiding banned vertices
            R = H - L
            dist = {u: 0.0}
            res = INF
            for _ in range(R):
                nd = dict(dist)
                for x, dx in dist.items():
                    for y, wy in adj[x]:
                        if y in banned:
                            continue
                        if dx + wy < nd.get(y, INF):
                            nd[y] = dx + wy
                dist = nd
            res = dist.get(d, INF)
            best = min(best, wsum + res)
            return
        for v, w in adj[u]:
            if v in path or v == d and len(path) < L:  # d may only be the last vertex of the prefix
                continue
            rec(path + [v], wsum + w)

    rec([s], 0.0)
    return best


def inst(N, p, rng):
    while True:
        G = nx.gnp_random_graph(N, p, seed=rng.randrange(10**9))
        if nx.is_connected(G):
            break
    for a, b in G.edges():
        G[a][b]["r"] = rng.choice([0.05, 0.1, 0.2, 0.3, 0.34])
        G[a][b]["w"] = rng.choice([1, 2, 3, 5, 8]) * 1e-6 * rng.uniform(0.5, 1.5)
    return G


if __name__ == "__main__":
    n_inst = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    rng = random.Random(7)
    mism = {"yen": 0, "prefix": 0}; infeas = 0; tot = 0
    for _ in range(n_inst):
        N = rng.randint(5, 11); p = rng.choice([0.3, 0.5]); t = rng.choice([0, 1, 2, 3])
        G = inst(N, p, rng)
        s, d = rng.sample(range(N), 2)
        L = 2 * t + 1
        if L > N - 1:
            continue
        H = N - 1
        rho = rng.choice([0.0, 0.1, 0.2])
        b = brute(G, s, d, L, H, rho)
        y = yen(G, s, d, L, H, rho)
        q = prefix_bf(G, s, d, L, H, rho)
        tot += 1
        infeas += (b == INF)
        for name, v in (("yen", y), ("prefix", q)):
            if not ((b == INF and v == INF) or (b < INF and v < INF and abs(b - v) <= 1e-12)):
                mism[name] += 1
    print(f"instances={tot}, infeasible={infeas}, mismatches vs brute force: {mism}")
