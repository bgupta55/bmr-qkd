"""
L7_ablation_v5b.py -- R5: ablation, Delta sweep, hop-cap sweep, F_th sweep and public-XOR parity on ONE seeded dataset.

Dataset: ER(n=30, p=0.25), G=50 graphs x m=5 link draws, n_audit = 10^4, binomial noise, s=0, d=29, Lmin=3 (t=1), H=5 unless swept.
Seeds are IDENTICAL to R1 (L1_sweep_v5b.py): the honest link realisation of a (graph,draw) is the same in the clean and the
compromised (f=0.10) condition and the same as in R1.  New here: link expiry with probability p_exp = 0.0488 per link-epoch
(R1 used p_exp = 0), so freshness is a real filter.  A dataset hash is written to ../data/R5_dataset_hash.txt.
All tables below are computed from the SAME link realisations in one pass; consistency checks are asserted at the end.

Variants (all use the same candidate set = all simple s-d paths with 3 <= hops <= H):
  full           validity(gate & no alarm & fresh) + Delta budget + rate floor R, max r_min then min hops
  no_defect      full without the budget constraint            (Delta-violation curve)
  no_freshness   full without the freshness filter
  no_alarm       full without the alarm filter
  no_alarm_R0    no alarm filter AND no rate floor (R=0)  [alarm acts alone]   ; full_R0 = full with R=0 (reference)
  hop_only       same constraints as full, objective = min hops only (random tie-break)   (true min-hop)
  hop_R0 / hop_R0_no_alarm / hop_R0_no_fresh : min-hop objective with NO rate floor (R=0) -> isolates the alarm / freshness filter
  rate_only      same constraints as full, objective = max r_min only (random tie-break)
  SPF            unconstrained shortest path (reference)
R_sec is ALWAYS evaluated against the truth: a selected path has R_sec = r_min iff every link is valid (gate, no alarm, fresh)
AND sum eps <= Delta; otherwise 0.
"""
import argparse, hashlib, json, os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B

SEED = 20261012; TOPO_ER = 1; SIZE = 30
P_EXP = 0.0488
VARIANTS = ["full", "no_defect", "no_freshness", "no_alarm", "no_alarm_R0", "full_R0", "hop_only", "rate_only", "hop_R0", "hop_R0_no_alarm", "hop_R0_no_fresh", "SPF"]
DELTA_GRID = np.logspace(-7, -3, 20)


def pathgate(L, lk):
    E = len(lk["r"]); Lp = np.where(L < 0, E, L)
    g = np.concatenate([lk["gate"], [True]])[Lp].all(1)
    c = np.concatenate([lk["comp"], [False]])[Lp].any(1)
    return g, c


def pick(mask, keys, rng):
    return B._pick(mask, keys, rng)


def select(variant, M, hops, gok, Delta, R, hcap, rng):
    seps, rmin = M["seps"], M["rmin"]
    noal = ~M["anyalarm"]; fr = ~M["anyexp"]
    hm = hops <= hcap
    bud = seps <= Delta * (1 + 1e-9)
    rf = rmin >= R
    if variant == "full":          return pick(hm & gok & noal & fr & bud & rf, [-rmin, hops], rng)
    if variant == "no_defect":     return pick(hm & gok & noal & fr & rf, [-rmin, hops], rng)
    if variant == "no_freshness":  return pick(hm & gok & noal & bud & rf, [-rmin, hops], rng)
    if variant == "no_alarm":      return pick(hm & gok & fr & bud & rf, [-rmin, hops], rng)
    if variant == "no_alarm_R0":   return pick(hm & gok & fr & bud, [-rmin, hops], rng)
    if variant == "full_R0":       return pick(hm & gok & noal & fr & bud, [-rmin, hops], rng)
    if variant == "hop_only":      return pick(hm & gok & noal & fr & bud & rf, [hops], rng)
    if variant == "rate_only":     return pick(hm & gok & noal & fr & bud & rf, [-rmin], rng)
    if variant == "hop_R0":        return pick(hm & gok & noal & fr & bud, [hops], rng)
    if variant == "hop_R0_no_alarm": return pick(hm & gok & fr & bud, [hops], rng)
    if variant == "hop_R0_no_fresh": return pick(hm & gok & noal & bud, [hops], rng)
    if variant == "SPF":           return pick(hm, [hops], rng)
    raise ValueError(variant)


def outcome(i, M, hops, comp_any, Delta, R):
    if i < 0:
        return dict(feasible=0, hops=0, r_nom=0.0, r_sec=0.0, exp_v=0, alarm_v=0, delta_v=0, sum_eps=0.0, touches_comp=0, at_cap=0)
    ok = bool(M["allvalid"][i]) and M["seps"][i] <= Delta * (1 + 1e-9)
    return dict(feasible=1, hops=int(hops[i]), r_nom=float(M["rmin"][i]), r_sec=float(M["rmin"][i]) if ok else 0.0,
                exp_v=int(M["anyexp"][i]), alarm_v=int(M["anyalarm"][i]), delta_v=int(M["seps"][i] > Delta * (1 + 1e-9)),
                sum_eps=float(M["seps"][i]), touches_comp=int(comp_any[i]), at_cap=0)


def run_graph(g):
    R = 0.10; Hmax = 7
    rg = np.random.default_rng(np.random.SeedSequence([SEED, SIZE, TOPO_ER, g, 0]))
    edges = B.build_graph("er", SIZE, rg); E = len(edges)
    adj = B.adj_matrix(SIZE, edges); s, d = 0, SIZE - 1
    L, hops = B.enumerate_paths(adj, s, d, 3, Hmax, cap=3_000_000)
    rows = []; hsh = hashlib.sha256()
    for draw in range(5):
        for f in (0.0, 0.10):
            rl = np.random.default_rng(np.random.SeedSequence([SEED, SIZE, TOPO_ER, g, draw, 1]))
            rc = np.random.default_rng(np.random.SeedSequence([SEED, SIZE, TOPO_ER, g, draw, 2, int(f * 1000)]))
            for pe_tag, pe in (("exp", P_EXP), ("noexp", 0.0)):
                lk = B.make_links(rl if pe_tag == "exp" else np.random.default_rng(np.random.SeedSequence([SEED, SIZE, TOPO_ER, g, draw, 1])),
                                  E, 10000, f=f, noise="binomial", p_exp=pe,
                                  rng_comp=rc if pe_tag == "exp" else np.random.default_rng(np.random.SeedSequence([SEED, SIZE, TOPO_ER, g, draw, 2, int(f * 1000)])))
                if pe_tag == "exp":
                    for k in ("S_hat", "r", "eps", "valid", "comp"):
                        hsh.update(np.asarray(lk[k]).tobytes())
                if len(hops) == 0:
                    continue
                M = B.path_metrics(L, lk)
                gok, cany = pathgate(L, lk)
                assert np.array_equal(M["allvalid"], gok & ~M["anyalarm"] & ~M["anyexp"])
                cond = "clean" if f == 0 else "compromised"
                rs = np.random.default_rng(np.random.SeedSequence([SEED, SIZE, TOPO_ER, g, draw, 9]))
                base = dict(graph_id=g, draw=draw, cond=cond, pexp=pe_tag)

                def emit(table, variant, Delta, hcap, i, extra=None):
                    o = outcome(i, M, hops, cany, Delta, R)
                    if i >= 0: o["at_cap"] = int(hops[i] == hcap)
                    r = dict(base); r.update(table=table, variant=variant, Delta=Delta, hcap=hcap); r.update(o)
                    if extra: r.update(extra)
                    rows.append(r)
                if pe_tag == "noexp":
                    # reproduces R1 (p_exp = 0): BMR-equivalent 'full' and SPF at Delta=5e-6, H=5
                    for v in ("full", "SPF"):
                        emit("r1check", v, 5e-6, 5, select(v, M, hops, gok, 5e-6, R, 5, rs))
                    continue
                # (b) ablation at two Deltas, H=5
                for Delta in (5e-6, 1e-6):
                    for v in VARIANTS:
                        emit("ablation", v, Delta, 5, select(v, M, hops, gok, Delta, R, 5, rs))
                # (a) Delta sweep (20 points) for full and no_defect
                for Delta in DELTA_GRID:
                    for v in ("full", "no_defect", "SPF"):
                        emit("dsweep", v, float(Delta), 5, select(v, M, hops, gok, float(Delta), R, 5, rs))
                # (c) hop-cap sweep
                for hc in (4, 5, 6, 7):
                    for Delta in (5e-6, 1e-6):
                        for v in ("full", "hop_only", "SPF"):
                            emit("hopcap", v, Delta, hc, select(v, M, hops, gok, Delta, R, hc, rs))
                        # SPF_valid_budget: shortest valid within-budget path (no rate floor) -- the strongest hop-based baseline
                        i = pick((hops <= hc) & M["allvalid"] & (M["seps"] <= Delta * (1 + 1e-9)), [hops], rs)
                        emit("hopcap", "SPF_valid_budget", Delta, hc, i)
                # (e) F_th sweep (QLEAP-style) against BMR/full, H=5, Delta=5e-6
                Fv = (1.0 + 3.0 * M["vprod"]) / 4.0
                for fth in (0.70, 0.75, 0.80):
                    i = pick((hops <= 5) & (Fv >= fth), [M["cost"]], rs)
                    emit("fth", f"QLEAP_F{fth:.2f}", 5e-6, 5, i, dict(F_th=fth))
                # (d) public-XOR parity: all routers, Delta=5e-6, H=5, Lmin=3
                m5 = hops <= 5
                L5, h5 = L[m5], hops[m5]
                M5 = B.path_metrics(L5, lk)
                sel, M2, F, tput = B.select_all(L5, h5, lk, 5e-6, R, "public", rs, B.ROUTERS_PUBLIC, M=M5)
                for r_ in B.trial_rows(sel, M2, F, tput, L5, h5, 5e-6, "public", 3):
                    rows.append(dict(base, table="parity", variant=r_["router"], Delta=5e-6, hcap=5, feasible=r_["feasible"],
                                     hops=r_["hops"], r_nom=r_["r_nom"], r_sec=r_["r_sec"], odd=int(r_["hops"] % 2 == 1) if r_["feasible"] else 0,
                                     exp_v=r_["expired_on_path"], alarm_v=r_["alarm_on_path"], delta_v=r_["delta_violated"],
                                     sum_eps=r_["sum_eps"], touches_comp=0, at_cap=0))
    return g, rows, hsh.hexdigest()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--G", type=int, default=50)
    ap.add_argument("--procs", type=int, default=10)
    ap.add_argument("--out", default="../data/R5_rows.csv.gz")
    a = ap.parse_args()
    t0 = time.time(); allrows = []; hs = {}
    with Pool(a.procs) as pool:
        for g, rows, h in pool.imap_unordered(run_graph, range(a.G)):
            allrows.extend(rows); hs[g] = h
    df = pd.DataFrame(allrows)
    df.to_csv(a.out, index=False, compression="gzip")
    dh = hashlib.sha256("".join(hs[g] for g in sorted(hs)).encode()).hexdigest()
    open(os.path.join(os.path.dirname(a.out), "R5_dataset_hash.txt"), "w").write(
        f"SHA256 of (S_hat, r, eps, valid, comp) link arrays, ER n=30 p=0.25, G={a.G}, m=5, f in {{0,0.10}}, p_exp={P_EXP}, seed={SEED}:\n{dh}\n")
    print(len(df), "rows", round(time.time() - t0, 1), "s; dataset hash", dh)
