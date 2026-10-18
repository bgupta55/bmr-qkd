"""
experiments/q6_sweep.py
=======================
Q6: Full sweep (pending_v7 Q6).  Parallel over cells (multiprocessing).

Design
  synthetic: 5 topologies {ER, BA, WS, grid, ring} x sizes {10,20,30,50} x {clean, compromised(20% links, S=2.418)}
             x n_audit {500, 10_000} x relay mode {sealed (phi=0.99), public-XOR (phi=0.30)}
             1000 trials per cell = 10 graph seeds x 100 audit-noise draws; same seeds for every router.
  real:      SNDlib Abilene(12), Polska(12), Nobel-Germany(17), GEANT(22), Germany50(50)  (1000 draws each)
  routers:   B0, B1a, B1b, B2, B3, B5, B6, B7, B8, BMR   (see q6_common.py)
  metrics:   Q1 definitions; paired comparison BMR vs each router (paired t, Wilcoxon, Holm; bootstrap CI of mean diff).
Parameters (every table): ell=1e4, phi=0.99 sealed / 0.30 public-XOR, delta=1e-6, R_MIN=0.01, W=10, t=1 (m>=3),
  Delta = 0.5 (n=500) and 1e-3 (n=10^4); Q15 shows both lie on the plateau for phi=0.99 / 0.30.
"""
import sys, os, json, math, time, hashlib, platform
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
sys.path.insert(0, os.path.dirname(__file__))
from multiprocessing import Pool
import numpy as np
import networkx as nx

from bmr.sim import LinkSim, S_COMPROMISED
from q6_common import (make_topology, load_sndlib, attach_certs, evaluate, feasible_oracle,
                       ROUTERS, ROUTER_NAMES, R_MIN, M_FLOOR)

HERE = os.path.dirname(__file__)
TRIALS = int(os.environ.get("Q6_TRIALS", 1000))
GRAPH_SEEDS = 10
DELTA_BY_N = {500: 0.5, 10_000: 1e-3}
F_COMP = 0.20
ELL, DELTA_C = 10_000, 1e-6
PHI_MODE = {"sealed": 0.99, "public-XOR": 0.30}
# Variant for Reviewer 2 Q15/D2: meaningful budget (phi = 0.9, Delta = 1e-4 at n=500, 1e-5 at n=1e4).
# Run:  Q6_VARIANT=phi09 python3 experiments/q6_sweep.py   -> data/results/q6_sweep_phi09.json
VARIANT = os.environ.get("Q6_VARIANT", "")
OUT_NAME = "q6_sweep.json"
if VARIANT == "phi09":
    DELTA_BY_N = {500: 1e-4, 10_000: 1e-5}
    PHI_MODE = {"sealed": 0.90, "public-XOR": 0.30}
    OUT_NAME = "q6_sweep_phi09.json"
SYN = ["ER", "BA", "WS", "grid", "ring"]
SIZES = [10, 20, 30, 50]
REAL = ["abilene", "polska", "nobel-germany", "geant", "germany50"]
MASTER = 61


def one_cell(spec):
    """Run one (topology, n, n_audit, condition, mode) cell -> (summary rows, per-trial arrays for stats)."""
    topo, n, n_audit, cond, mode, seed = spec
    Delta = DELTA_BY_N[n_audit]
    phi = PHI_MODE[mode]
    fcomp = F_COMP if cond == "compromised" else 0.0
    rng = np.random.default_rng(seed)
    if topo in REAL:
        gl = [load_sndlib(os.path.join(HERE, "..", "data", "topologies", topo + ".txt"))] * GRAPH_SEEDS
    else:
        gl = []
        while len(gl) < GRAPH_SEEDS:
            G = make_topology(topo, n, int(rng.integers(1e9)))
            if nx.is_connected(G):
                gl.append(G)
    per = {r: dict(rate=[], compl=[], found=[], hops=[], w=[], rate_b=[], compl_b=[]) for r in ROUTER_NAMES}
    feas = []
    npg = TRIALS // GRAPH_SEEDS
    for G0 in gl:
        edges = list(G0.edges())
        dist = dict(nx.all_pairs_shortest_path_length(G0))
        pairs = [(a, b) for a in G0 for b in G0 if a < b and dist[a][b] >= 2]
        for _ in range(npg):
            s, d = pairs[int(rng.integers(len(pairs)))]      # random request, hop distance >= 2
            sd = int(rng.integers(1e9))
            sim = LinkSim(edges, n_audit, delta=DELTA_C, ell=ELL, phi=phi, seed=sd)
            if fcomp > 0:
                k = max(1, int(round(fcomp * len(edges))))
                for i in np.random.default_rng(sd).choice(len(edges), k, replace=False):
                    sim.set_true_S(edges[i], S_COMPROMISED)
            G = G0.copy()
            attach_certs(G, sim, epoch=0)
            feas.append(int(feasible_oracle(G, s, d, Delta)))
            for rn, fn in ROUTERS:
                p = fn(G, s, d, Delta)
                f, c, rate, h, w, cb, rb_ = evaluate(G, p, Delta)
                P = per[rn]
                P["rate_b"].append(rb_); P["compl_b"].append(cb); P["rate"].append(rate); P["compl"].append(c); P["found"].append(f); P["hops"].append(h); P["w"].append(w)
    mf = float(np.mean(feas))
    rows = []
    for rn in ROUTER_NAMES:
        P = per[rn]
        hs = [h for h, f in zip(P["hops"], P["found"]) if f]
        # self-check (Q1): secure <= feasible-share * max link rate  (per-trial dominance of BMR by feasibility)
        rows.append(dict(topo=topo, n=n, n_audit=n_audit, cond=cond, mode=mode, phi=phi, Delta=Delta, router=rn,
                         mean_rate=float(np.mean(P["rate"])), mean_rate_budget_only=float(np.mean(P["rate_b"])), compl_budget_only_pct=100 * float(np.mean(P["compl_b"])), compl_pct=100 * float(np.mean(P["compl"])),
                         found_pct=100 * float(np.mean(P["found"])), refusal_pct=100 * (1 - float(np.mean(P["found"]))),
                         feasible_pct=100 * mf, mean_hops=float(np.mean(hs)) if hs else 0.0,
                         mean_sumw=float(np.mean([w for w, f in zip(P["w"], P["found"]) if f] or [0.0])),
                         trials=len(P["rate"])))
    rate_arr = {rn: per[rn]["rate"] for rn in ROUTER_NAMES}
    rate_arr_b = {rn: per[rn]["rate_b"] for rn in ROUTER_NAMES}
    # integrity assertions
    for rn in ROUTER_NAMES:
        assert all(c <= f for c, f in zip(per[rn]["compl"], per[rn]["found"]))
    assert all(c == 1 for c, f in zip(per["BMR"]["compl"], per["BMR"]["found"]) if f), "BMR delivered a non-compliant route"
    # compliant => feasible, for EVERY router (14 Oct N-3): a compliant route is a witness for the exact oracle
    for rn in ROUTER_NAMES:
        bad_i = [i for i, (c, fe) in enumerate(zip(per[rn]["compl"], feas)) if c > fe]
        assert not bad_i, f"{rn}: compliant but oracle says infeasible in {len(bad_i)} trials, cell {spec[:5]}, first trial {bad_i[0]}, w={per[rn]['w'][bad_i[0]]}, hops={per[rn]['hops'][bad_i[0]]}, Delta={Delta}"
    # BMR can only be found if feasible
    assert all(f <= fe for f, fe in zip(per["BMR"]["found"], feas))
    return rows, (rate_arr, rate_arr_b), feas


def main():
    t0 = time.time()
    rng = np.random.default_rng(MASTER)
    specs = []
    for topo in SYN:
        for n in SIZES:
            for na in (500, 10_000):
                for cond in ("clean", "compromised"):
                    for mode in ("sealed", "public-XOR"):
                        specs.append((topo, n, na, cond, mode, int(rng.integers(1e9))))
    for topo in REAL:
        for na in (500, 10_000):
            for cond in ("clean", "compromised"):
                for mode in ("sealed", "public-XOR"):
                    specs.append((topo, 0, na, cond, mode, int(rng.integers(1e9))))
    print(f"cells: {len(specs)}  trials/cell: {TRIALS}", flush=True)
    rows, rates, rates_b = [], {}, {}
    with Pool(processes=max(1, (os.cpu_count() or 4) - 1)) as pool:
        for i, (r, ra, fe) in enumerate(pool.imap(one_cell, specs, chunksize=1)):
            rows += r
            sp = specs[i]
            rates[sp[:5]] = ra[0]; rates_b[sp[:5]] = ra[1]
            if i % 10 == 0:
                print(f"  {i+1}/{len(specs)} done  ({time.time()-t0:.0f}s)", flush=True)
    # ---- paired statistics BMR vs every router over all cells (per-cell mean paired difference)
    from scipy import stats
    stat_rows = []
    pvals = []
    for metric, RT in (("rate_full(floor+budget)", rates), ("rate_budget_only", rates_b)):
        for rn in ROUTER_NAMES:
            if rn == "BMR":
                continue
            for scope in ("all", "compromised", "clean"):
                keys = [k for k in RT if scope == "all" or k[3] == scope]
                d_cell = np.array([np.mean(np.array(RT[k]["BMR"]) - np.array(RT[k][rn])) for k in keys])
                pooled = np.concatenate([np.array(RT[k]["BMR"]) - np.array(RT[k][rn]) for k in keys])
                t_stat, p_t = stats.ttest_1samp(pooled, 0.0)
                try:
                    p_w = stats.wilcoxon(pooled[pooled != 0]).pvalue if np.any(pooled != 0) else 1.0
                except ValueError:
                    p_w = 1.0
                rb = np.random.default_rng(5)
                bs = [rb.choice(d_cell, len(d_cell)).mean() for _ in range(2000)]
                stat_rows.append(dict(metric=metric, vs=rn, scope=scope, n_cells=len(d_cell), n_trials=int(len(pooled)),
                                      mean_diff=float(pooled.mean()), ci_lo=float(np.percentile(bs, 2.5)),
                                      ci_hi=float(np.percentile(bs, 97.5)), p_ttest=float(p_t), p_wilcoxon=float(p_w),
                                      frac_cells_bmr_ge=float(np.mean(d_cell >= -1e-12))))
                pvals.append(p_w)
    # Holm correction over all (router, scope) tests, Wilcoxon
    order = np.argsort(pvals); m = len(pvals)
    holm = np.empty(m)
    run = 0.0
    for rank, idx in enumerate(order):
        run = max(run, min(1.0, (m - rank) * pvals[idx]))
        holm[idx] = run
    for sr, h in zip(stat_rows, holm):
        sr["p_holm"] = float(h)
    out = dict(experiment="Q6", description="full sweep, 10 routers, synthetic + SNDlib, sealed + public-XOR",
               header=dict(trials=TRIALS, graph_seeds=GRAPH_SEEDS, delta_by_n=DELTA_BY_N, f_comp=F_COMP, ell=ELL,
                           delta_cert=DELTA_C, phi_mode=PHI_MODE, R_MIN=R_MIN, hop_floor=M_FLOOR, t=1,
                           python=platform.python_version(), seconds=round(time.time() - t0, 1),
                           code_sha256=hashlib.sha256(open(os.path.join(HERE, "q6_common.py"), "rb").read()).hexdigest()[:12]),
               rows=rows, paired_stats=stat_rows)
    with open(os.path.join(HERE, "..", "data", "results", OUT_NAME), "w") as f:
        json.dump(out, f, indent=1)
    print(f"saved {OUT_NAME}  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
