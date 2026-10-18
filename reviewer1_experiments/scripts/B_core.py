"""
B_core.py -- shared simulation core for Section B (R1, R4, R5, R9, R10).

Design (answers to the defects D1-D5, D12 found in Section A):
  * Link noise: true binomial sampling of the four CHSH correlators
    (n shots per setting).  S_hat = 2*K/n - 4 with K ~ Binomial(4n, (1+S/4)/2).
    This has exactly the Werner variance (4 - S^2/4)/n, so the honest-link
    false-alarm rate is 1e-6 as the theory intends (fixes D1).  A legacy
    Gaussian(2/sqrt n) mode is kept only for the R1-vs-v5 comparison.
  * Exact routing: ALL simple s-d paths with L_min <= hops <= H are enumerated
    once per trial (numba DFS); every router is a mask + objective on that set,
    so every router is optimal for its own definition (no heuristic error).
  * Routers filter on validity (gate & ~alarm & fresh) where their definition
    says so (fixes D2/D3).  Clean condition f=0 is supported (fixes D4).
  * One ER generator (p=0.25) everywhere (fixes D5).  F_P, pair cost and
    throughput proxy are exported for every router (fixes D12).
"""
import math
import numpy as np
import numba as nb

SQ2 = math.sqrt(2.0)
Z_ALPHA = 4.7534            # one-sided z for alpha = 1e-6
S0_REF = 2.687              # commissioning baseline of a standard link
S_ATTACK = 2.418
EPS_ATTACK = 1e-6
# class table: S_true, eps_c, gamma_c, probability  (A1/A2 of Section A)
CLASS_S = np.array([2.720, 2.687, 2.640])
CLASS_EPS = np.array([2e-7, 1e-6, 5e-6])
CLASS_GAMMA = np.array([1e-7, 8e-7, 4e-6])
CLASS_P = np.array([0.3, 0.5, 0.2])
F_TH = 0.75
Q_SWAP = 0.9


# ------------------------------------------------------------------ formulas
def h_vec(S):
    """Pironio min-entropy bound, vectorised; 0 for S <= 2."""
    S = np.minimum(np.asarray(S, dtype=float), 2 * SQ2 - 1e-10)
    inner = np.maximum(0.0, 2.0 - S ** 2 / 4.0)
    out = 1.0 - np.log2(1.0 + np.sqrt(inner))
    return np.where(S <= 2.0, 0.0, out)


def h(S):
    return float(h_vec(S))


def tau(n, delta=1e-6):
    """Hoeffding radius (called eta in the revised manuscript)."""
    return math.sqrt(32.0 * math.log(8.0 / delta) / n)


def sigma_w(S, n):
    return np.sqrt(np.maximum(0.0, 4.0 - np.asarray(S) ** 2 / 4.0)) / math.sqrt(n)


# --------------------------------------------------------------- link model
def make_links(rng, E, n_audit, f=0.0, noise="binomial", p_exp=0.0,
               rng_comp=None, delta=1e-6):
    """Draw E links.  Returns dict of arrays (length E).
    f = fraction compromised (>=1 link if f>0).  p_exp = per-link expiry prob.
    Base draws consume `rng` identically for every f, so clean/compromised
    conditions share the same honest realisation (paired design)."""
    cls = rng.choice(3, size=E, p=CLASS_P)
    S_true = CLASS_S[cls].copy()
    eps = CLASS_EPS[cls].copy()
    gam = CLASS_GAMMA[cls].copy()
    S0_e = S_true.copy()
    u_norm = rng.standard_normal(E)
    u_bin = rng.random(E)
    u_fresh = rng.random(E)
    comp = np.zeros(E, dtype=bool)
    if f > 0 and E > 0:
        k = max(1, int(round(f * E)))
        rc = rng_comp if rng_comp is not None else rng
        comp[rc.choice(E, size=min(k, E), replace=False)] = True
        S_true[comp] = S_ATTACK
        eps[comp] = EPS_ATTACK
        gam[comp] = EPS_ATTACK
        S0_e[comp] = S0_REF
    if noise == "binomial":
        # exact Binomial(4n, (1+S/4)/2) by inverse CDF of a per-link uniform:
        # honest links are identical in the clean and compromised conditions.
        p = (1.0 + S_true / 4.0) / 2.0
        Kn = 4 * n_audit
        K = _binom_ppf(np.clip(u_bin, 1e-12, 1 - 1e-12), Kn, p)
        S_hat = 2.0 * K / n_audit - 4.0
    elif noise == "gauss02":                      # legacy v5 model (sigma = 2/sqrt n)
        S_hat = S_true + (2.0 / math.sqrt(n_audit)) * u_norm
    else:
        raise ValueError(noise)
    t = tau(n_audit, delta)
    S_lb = np.maximum(0.0, S_hat - t)
    r_e = h_vec(S_lb)
    gate = S_lb >= 2.0
    alarm = S_hat < S0_e - Z_ALPHA * sigma_w(S0_e, n_audit)
    fresh = u_fresh >= p_exp
    valid = gate & ~alarm & fresh
    v = np.clip(S_hat / (2 * SQ2), 0.0, 1.0)
    pgen = np.maximum(0.01, v)
    return dict(r=r_e, eps=eps, gamma=gam, S_hat=S_hat, S_lb=S_lb, gate=gate,
                alarm=alarm, fresh=fresh, valid=valid, comp=comp, v=v, pgen=pgen,
                S_true=S_true, cls=cls)


def _binom_ppf(u, Kn, p):
    from scipy.stats import binom
    return binom.ppf(u, Kn, p)


# -------------------------------------------------------------- graph builders
def build_graph(topo, n, rng, er_p=0.25):
    """Return list of undirected edges (u,v), u<v, and (s,d) for sealed use."""
    if topo == "er":
        iu, ju = np.triu_indices(n, 1)
        keep = rng.random(len(iu)) < er_p
        edges = list(zip(iu[keep].tolist(), ju[keep].tolist()))
    elif topo == "ba":
        import networkx as nx
        G = nx.barabasi_albert_graph(n, 2, seed=int(rng.integers(2 ** 31)))
        edges = sorted((min(a, b), max(a, b)) for a, b in G.edges())
    elif topo == "ring":
        edges = sorted((min(i, (i + 1) % n), max(i, (i + 1) % n)) for i in range(n))
    elif topo == "line":
        edges = [(i, i + 1) for i in range(n - 1)]
    elif topo == "star":
        edges = [(0, i) for i in range(1, n)]
    else:
        raise ValueError(topo)
    return edges


def choose_sd(topo, n, mode="sealed"):
    if topo == "line":
        return 0, min(3, n - 1)
    if topo == "ring":
        return 0, min(3, n // 2)
    if topo == "star":
        return (1, n - 1) if mode == "sealed" else (0, n - 1)
    return 0, n - 1


def adj_matrix(n, edges):
    adj = -np.ones((n, n), dtype=np.int32)
    for i, (u, v) in enumerate(edges):
        adj[u, v] = i
        adj[v, u] = i
    return adj


# ----------------------------------------------------------- path enumeration
@nb.njit(cache=True)
def _enum(adj, n, s, d, Lmin, H, cap, out_links, out_hops):
    cnt = 0
    nodes = np.empty(H + 2, np.int32)
    links = np.empty(H + 2, np.int32)
    nxt = np.zeros(H + 2, np.int32)
    vis = np.zeros(n, np.bool_)
    depth = 0
    nodes[0] = s
    vis[s] = True
    nxt[0] = 0
    while depth >= 0:
        u = nodes[depth]
        if nxt[depth] >= n:
            vis[u] = False
            depth -= 1
            continue
        v = nxt[depth]
        nxt[depth] += 1
        e = adj[u, v]
        if e < 0 or vis[v]:
            continue
        if v == d:
            hh = depth + 1
            if hh >= Lmin:
                if cnt >= cap:
                    return -1
                for k in range(depth):
                    out_links[cnt, k] = links[k]
                out_links[cnt, depth] = e
                out_hops[cnt] = hh
                cnt += 1
            continue
        if depth + 1 < H:
            links[depth] = e
            depth += 1
            nodes[depth] = v
            vis[v] = True
            nxt[depth] = 0
    return cnt


def enumerate_paths(adj, s, d, Lmin, H, cap=600000):
    n = adj.shape[0]
    out_links = -np.ones((cap, H), dtype=np.int32)
    out_hops = np.zeros(cap, dtype=np.int32)
    c = _enum(adj, n, s, d, Lmin, H, cap, out_links, out_hops)
    if c < 0:
        raise RuntimeError("path-enumeration capacity exceeded")
    return out_links[:c].copy(), out_hops[:c].copy()


# ------------------------------------------------------------------- routers
ROUTERS_SEALED = ["BMR", "SPF", "SPF_valid", "SPF_budget", "SPF_valid_budget",
                  "BMR_hop", "QPATHstyle", "QLEAPstyle", "QCASTstyle"]
ROUTERS_PUBLIC = ROUTERS_SEALED + ["BMR_np"]


def _pick(mask, keys, rng):
    """argmin lexicographic over `keys` (priority order) among mask, random tiebreak."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return -1
    if idx.size == 1:
        return int(idx[0])
    ks = [np.asarray(k)[idx] for k in keys] + [rng.random(idx.size)]
    order = np.lexsort(ks[::-1])
    return int(idx[order[0]])


def path_metrics(L, lk):
    """Per-path aggregate arrays from padded link-index matrix L (P,H)."""
    E = len(lk["r"])
    Lp = np.where(L < 0, E, L)

    def pad(a, fill):
        return np.concatenate([np.asarray(a), [fill]])

    r = pad(lk["r"], np.inf)[Lp]
    rmin = r.min(1)
    seps = pad(lk["eps"], 0.0)[Lp].sum(1)
    allvalid = pad(lk["valid"], True)[Lp].all(1)
    anyalarm = pad(lk["alarm"], False)[Lp].any(1)
    anyexp = (~pad(lk["fresh"], True)[Lp]).any(1)
    vprod = pad(lk["v"], 1.0)[Lp].prod(1)
    cost = pad(1.0 / lk["pgen"], 0.0)[Lp].sum(1)
    pgprod = pad(lk["pgen"], 1.0)[Lp].prod(1)
    return dict(rmin=rmin, seps=seps, allvalid=allvalid, anyalarm=anyalarm,
                anyexp=anyexp, vprod=vprod, cost=cost, pgprod=pgprod)


def select_all(L, hops, lk, Delta, R, mode, rng, routers, M=None):
    """Return {router: path index or -1}."""
    if M is None:
        M = path_metrics(L, lk)
    rmin, seps, valid = M["rmin"], M["seps"], M["allvalid"]
    odd = (hops % 2 == 1)
    bud = seps <= Delta * (1 + 1e-9)
    reqf = rmin >= R
    F = (1.0 + 3.0 * M["vprod"]) / 4.0
    tput = M["pgprod"] * np.where(hops > 1, Q_SWAP ** (hops - 1), 1.0)
    allp = np.ones(len(hops), bool)
    bmr_mask = valid & bud & reqf
    bmr_aware = bmr_mask & odd if mode == "public" else bmr_mask
    sel = {}
    sel["BMR"] = _pick(bmr_aware, [-rmin, hops], rng)
    sel["SPF"] = _pick(allp, [hops], rng)
    sel["SPF_valid"] = _pick(valid, [hops], rng)
    sel["SPF_budget"] = _pick(bud, [hops], rng)
    sel["SPF_valid_budget"] = _pick(valid & bud, [hops], rng)
    sel["BMR_hop"] = _pick(bmr_aware, [hops, -rmin], rng)
    sel["QPATHstyle"] = _pick(allp, [-M["vprod"]], rng)
    sel["QLEAPstyle"] = _pick(F >= F_TH, [M["cost"]], rng)
    sel["QCASTstyle"] = _pick(allp, [-tput], rng)
    if "BMR_np" in routers:
        sel["BMR_np"] = _pick(bmr_mask, [-rmin, hops], rng)
    return {k: sel[k] for k in routers}, M, F, tput


def trial_rows(sel, M, F, tput, L, hops, Delta, mode, Lmin):
    """Evaluate routers; returns list of dict rows (one per router)."""
    rows = []
    bmr_idx = sel.get("BMR", -1)
    for rt, i in sel.items():
        if i < 0:
            rows.append(dict(router=rt, feasible=0, path="", hops=0, r_nom=0.0,
                             r_sec=0.0, F_P=0.0, cost_pairs=0.0, throughput_proxy=0.0,
                             sum_eps=0.0, delta_violated=0, alarm_on_path=0,
                             expired_on_path=0, tie_with_BMR=int(bmr_idx < 0)))
            continue
        hp = int(hops[i])
        assert hp >= Lmin
        rn = float(M["rmin"][i])
        ok = bool(M["allvalid"][i]) and M["seps"][i] <= Delta * (1 + 1e-9)
        if mode == "public" and hp % 2 == 0:
            ok = False            # even-hop public-XOR route: zero secrecy (Prop. 3)
        rs = rn if ok else 0.0
        assert rs <= rn + 1e-12
        nodes_links = "-".join(str(int(x)) for x in L[i] if x >= 0)
        rows.append(dict(router=rt, feasible=1, path=nodes_links, hops=hp, r_nom=rn,
                         r_sec=rs, F_P=float(F[i]), cost_pairs=float(M["cost"][i]),
                         throughput_proxy=float(tput[i]), sum_eps=float(M["seps"][i]),
                         delta_violated=int(M["seps"][i] > Delta * (1 + 1e-9)),
                         alarm_on_path=int(M["anyalarm"][i]),
                         expired_on_path=int(M["anyexp"][i]),
                         tie_with_BMR=int(i == bmr_idx)))
    return rows


def check_assertions(rows, mode, R):
    """Dominance assertions requested in the PENDING file (sealed mode)."""
    if mode != "sealed":
        return
    d = {r["router"]: r for r in rows}
    b, bh = d["BMR"], d["BMR_hop"]
    assert b["feasible"] == bh["feasible"], "BMR / BMR_hop feasible sets differ"
    if bh["feasible"]:
        assert b["r_nom"] >= bh["r_nom"] - 1e-12, "BMR r_nom < BMR_hop r_nom"
    v = d["SPF_valid_budget"]
    if v["feasible"] and v["r_nom"] >= R:
        assert b["feasible"] == 1 and b["r_nom"] >= v["r_nom"] - 1e-12


def self_test():
    # 1. enumeration on a complete graph K5, s=0,d=4, Lmin=1,H=4: count simple paths
    n = 5
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    adj = adj_matrix(n, edges)
    L, hp = enumerate_paths(adj, 0, 4, 1, 4)
    # paths 0->4 in K5 with k internal nodes chosen in order: sum_k P(3,k) = 1+3+6+6
    assert len(hp) == 16, len(hp)
    # 2. binomial noise has Werner variance
    rng = np.random.default_rng(1)
    S = []
    for _ in range(40):
        lk = make_links(rng, 100, 10000)
        S.append(lk["S_hat"] - lk["S_true"])
    S = np.concatenate(S)
    sd_emp = S.std()
    sd_th = float(np.sqrt(CLASS_P @ (sigma_w(CLASS_S, 10000) ** 2)))
    assert abs(sd_emp / sd_th - 1) < 0.05, (sd_emp, sd_th)
    return dict(K5_paths=len(hp), sd_empirical=float(sd_emp), sd_theory=sd_th)


if __name__ == "__main__":
    print(self_test())
