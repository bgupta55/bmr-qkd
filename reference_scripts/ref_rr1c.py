"""
Independent reference re-implementation of the (R)-family snapshot experiment
(RR1 decomposition + ground-truth metrics RR1c).  NOT the authors' code: it exists only
to give *expected magnitudes* and pass/fail criteria for the authors' own runs.

Model (blueprint Sec. 2.1): per link Shat ~ N(S, sigma_W(S,n)), gate Shat-eta>2,
h = Pironio bound on min(Shat-eta, 2*sqrt2), H = l*h, k = floor(phi*H), r = k/l,
eps = 0.5*2^(-(H-k)/2), w = eps + delta.  Alarm: Shat < S0 - z*sigma_W(S0).
"""
import math, random, sys
import numpy as np
import networkx as nx

L = 10_000; PHI = 0.99; DELTA = 1e-6; RMIN = 0.01; T = 1; FLOOR = 2 * T + 1
S0 = 2.687; Z = 4.7534; HCAP = 8


def sigma_w(S, n):
    return math.sqrt(max(4 - S * S / 4, 1e-12) / n)


def h_of(x):
    x = min(x, 2 * math.sqrt(2))
    return 1 - math.log2(1 + math.sqrt(max(2 - x * x / 4, 0)))


def link_state(S_true, n, rng):
    eta = math.sqrt(8 * math.log(1 / DELTA) / n)
    Sh = rng.gauss(S_true, sigma_w(S_true, n))
    alarm = Sh < S0 - Z * sigma_w(S0, n)
    if Sh - eta <= 2:
        return dict(cert=False, r=0.0, w=math.inf, alarm=alarm)
    H = L * h_of(Sh - eta)
    k = math.floor(PHI * H)
    if k <= 0:
        return dict(cert=False, r=0.0, w=math.inf, alarm=alarm)
    eps = 0.5 * 2 ** (-(H - k) / 2)
    return dict(cert=True, r=k / L, w=eps + DELTA, alarm=alarm)


def trial(N, f, S_dag, delta_budget, n, rng):
    while True:
        G = nx.gnp_random_graph(N, 4 / (N - 1), seed=rng.randrange(10**9))
        if nx.is_connected(G):
            break
    E = list(G.edges())
    bad = set(rng.sample(range(len(E)), max(1, round(f * len(E)))))
    st = {}
    truly = {}
    for i, (a, b) in enumerate(E):
        S = S_dag if i in bad else S0
        key = frozenset((a, b))
        st[key] = link_state(S, n, rng)
        truly[key] = i in bad
    s, d = rng.sample(range(N), 2)
    paths = [p for p in nx.all_simple_paths(G, s, d, cutoff=HCAP)]
    if not paths:
        return None
    rng.shuffle(paths)

    def lk(p):
        return [frozenset((p[i], p[i + 1])) for i in range(len(p) - 1)]

    def compliant(p, no_alarm):
        ls = lk(p)
        if len(ls) < FLOOR:
            return False
        if any((not st[l]["cert"]) or st[l]["r"] < RMIN for l in ls):
            return False
        if no_alarm and any(st[l]["alarm"] for l in ls):
            return False
        return sum(st[l]["w"] for l in ls) <= delta_budget

    def bott(p):
        return min(st[l]["r"] for l in lk(p))

    def exposed(p):
        return any(truly[l] for l in lk(p))

    out = {}
    spf = min(paths, key=len)
    ties = [p for p in paths if len(p) == len(spf)]
    spf = ties[0]
    out["SPF"] = (spf, compliant(spf, False))
    for name, na in (("D2", False), ("D3", True)):
        cand = [p for p in sorted(paths, key=len) if compliant(p, na)]
        out[name] = (cand[0], True) if cand else (None, False)
    for name, na in (("BMR", False), ("BMR+alarm", True)):
        cand = [p for p in paths if compliant(p, na)]
        out[name] = (max(cand, key=bott), True) if cand else (None, False)
    res = {}
    feasible = any(compliant(p, False) for p in paths)
    for name, (p, comp) in out.items():
        if p is None:
            res[name] = dict(found=0, sec=0.0, tsec=0.0, exp=None)
            continue
        b = bott(p)
        sec = b if comp else 0.0
        ex = exposed(p)
        res[name] = dict(found=1, sec=sec, tsec=sec if (comp and not ex) else 0.0, exp=int(ex))
    res["_feasible"] = int(feasible)
    return res


def cell(N, f, S_dag, delta_budget, n, trials, seed):
    rng = random.Random(seed)
    names = ["SPF", "D2", "D3", "BMR", "BMR+alarm"]
    acc = {k: dict(sec=[], tsec=[], exp=[], found=[]) for k in names}
    feas = 0; m = 0
    for _ in range(trials):
        r = trial(N, f, S_dag, delta_budget, n, rng)
        if r is None:
            continue
        m += 1; feas += r["_feasible"]
        for k in names:
            acc[k]["sec"].append(r[k]["sec"]); acc[k]["tsec"].append(r[k]["tsec"])
            acc[k]["found"].append(r[k]["found"])
            if r[k]["exp"] is not None:
                acc[k]["exp"].append(r[k]["exp"])
    return acc, feas, m


if __name__ == "__main__":
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    n = 10_000
    print(f"trials/cell={trials}, N=30, f=0.20, n={n}")
    print("S_dag  Delta   feas | secure rate: SPF  D2  D3  BMR  BMR+al | true-secure: SPF D2 D3 BMR BMR+al | exposure%: SPF D2 D3 BMR BMR+al")
    for S_dag in (2.418, 2.55, 2.62):
        for D in (1e-3, 3e-5):
            acc, feas, m = cell(30, 0.20, S_dag, D, n, trials, seed=int(S_dag * 1000) + int(D * 1e6))
            names = ["SPF", "D2", "D3", "BMR", "BMR+alarm"]
            sec = [np.mean(acc[k]["sec"]) for k in names]
            tsec = [np.mean(acc[k]["tsec"]) for k in names]
            ex = [100 * np.mean(acc[k]["exp"]) if acc[k]["exp"] else float("nan") for k in names]
            vd = np.array(acc["D3"]["sec"]) - np.array(acc["D2"]["sec"])
            se = vd.std(ddof=1) / math.sqrt(len(vd))
            print(f"{S_dag:5.3f} {D:7.0e} {feas/m:5.2f} | " + " ".join(f"{x:.3f}" for x in sec) +
                  " | " + " ".join(f"{x:.3f}" for x in tsec) + " | " + " ".join(f"{x:5.1f}" for x in ex) +
                  f" | validity D3-D2 = {vd.mean():+.4f} (SE {se:.4f})")
