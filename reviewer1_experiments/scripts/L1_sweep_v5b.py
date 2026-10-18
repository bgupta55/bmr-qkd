"""
L1_sweep_v5b.py -- R1: core topology sweep (Section B).
Replaces v5 L1.  Exact routing on the full set of simple paths (B_core.py).

Usage (full run used for the report):
  python3 L1_sweep_v5b.py --out ../data/R1_trials.csv.gz
"""
import argparse, json, os, sys, time, gzip
from multiprocessing import Pool
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B

TOPO_ID = {"er": 1, "ba": 2, "ring": 3, "line": 4, "star": 5}


def run_cell(args):
    (topo, size, G, m, n_audit, fs, deltas, R, H, noise, seed, p_exp) = args
    rows = []
    t0 = time.time()
    for g in range(G):
        rg = np.random.default_rng(np.random.SeedSequence([seed, size, TOPO_ID[topo], g, 0]))
        edges = B.build_graph(topo, size, rg)
        E = len(edges)
        if E == 0:
            continue
        adj = B.adj_matrix(size, edges)
        # (name, Lmin, s, d, modes)
        if topo == "star":
            variants = [("sealed_t0", 1, "sealed"), ("public_t0", 1, "public")]
        else:
            variants = [("sealed", 3, "sealed"), ("public", 3, "public")]
        enum_cache = {}
        sd_of = {}
        for vname, Lmin, mode in variants:
            sd_of[vname] = B.choose_sd(topo, size, mode) if topo != "star" else \
                ((1, size - 1) if mode == "sealed" else (0, size - 1))
            enum_cache[vname] = B.enumerate_paths(adj, sd_of[vname][0], sd_of[vname][1], Lmin, H)
        for draw in range(m):
            for f in fs:
                rl = np.random.default_rng(np.random.SeedSequence([seed, size, TOPO_ID[topo], g, draw, 1]))
                rc = np.random.default_rng(np.random.SeedSequence([seed, size, TOPO_ID[topo], g, draw, 2, int(f * 1000)]))
                lk = B.make_links(rl, E, n_audit, f=f, noise=noise, p_exp=p_exp, rng_comp=rc)
                for vname, Lmin, mode in variants:
                    L, hops = enum_cache[vname]
                    if len(hops) == 0:
                        M = None
                    else:
                        M = B.path_metrics(L, lk)
                    for Delta in deltas:
                        routers = B.ROUTERS_PUBLIC if mode == "public" else B.ROUTERS_SEALED
                        rsel = np.random.default_rng(np.random.SeedSequence([seed, size, TOPO_ID[topo], g, draw, 3]))
                        if M is None:
                            rr = [dict(router=r, feasible=0, path="", hops=0, r_nom=0.0, r_sec=0.0, F_P=0.0,
                                       cost_pairs=0.0, throughput_proxy=0.0, sum_eps=0.0, delta_violated=0,
                                       alarm_on_path=0, expired_on_path=0, tie_with_BMR=1) for r in routers]
                        else:
                            sel, M2, F, tput = B.select_all(L, hops, lk, Delta, R, mode, rsel, routers, M=M)
                            rr = B.trial_rows(sel, M2, F, tput, L, hops, Delta, mode, Lmin)
                            B.check_assertions(rr, mode if vname in ("sealed", "public") and mode == "sealed" else "x", R)
                        for r in rr:
                            r.update(topology=topo, size=size, mode=vname, cond="clean" if f == 0 else "compromised",
                                     f=f, Delta=Delta, n=n_audit, graph_id=g, draw=draw, runtime_ms=0.0)
                            rows.append(r)
    el = (time.time() - t0) * 1000.0
    for r in rows:
        r["runtime_ms"] = el / max(1, len(rows))
    return (topo, size, rows, el / 1000.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="+", default=[10, 20, 30, 50])
    ap.add_argument("--topos", nargs="+", default=["er", "ba", "ring", "line", "star"])
    ap.add_argument("--G", type=int, default=50)
    ap.add_argument("--m", type=int, default=5)
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--noise", default="binomial")
    ap.add_argument("--R", type=float, default=0.10)
    ap.add_argument("--Delta", type=float, nargs="+", default=[5e-6, 1e-6])
    ap.add_argument("--H", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--out", default="../data/R1_trials.csv.gz")
    ap.add_argument("--procs", type=int, default=10)
    a = ap.parse_args()

    tasks = []
    for topo in a.topos:
        for size in a.sizes:
            fs = [0.0, 0.10]
            tasks.append((topo, size, a.G, a.m, a.n, fs, a.Delta, a.R, a.H, a.noise, a.seed, 0.0))
    # f-sweep on ER n=30 (compromised, f = 0.05, 0.2, 0.3); f=0.10 already in main run
    tasks.append(("er", 30, a.G, a.m, a.n, [0.05, 0.20, 0.30], a.Delta, a.R, a.H, a.noise, a.seed, 0.0))
    tasks.sort(key=lambda t: -(t[1] if t[0] in ("er", "ba") else 1))
    t0 = time.time()
    allrows = []
    with Pool(a.procs) as pool:
        for topo, size, rows, el in pool.imap_unordered(run_cell, tasks):
            print(f"  done {topo:5s} n={size:3d}: {len(rows):7d} rows  {el:7.1f}s", flush=True)
            allrows.extend(rows)
    df = pd.DataFrame(allrows)
    cols = ["topology", "size", "mode", "cond", "f", "Delta", "n", "graph_id", "draw", "router", "feasible",
            "path", "hops", "r_nom", "r_sec", "F_P", "cost_pairs", "throughput_proxy", "sum_eps",
            "delta_violated", "alarm_on_path", "expired_on_path", "tie_with_BMR", "runtime_ms"]
    df = df[cols]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    df.to_csv(a.out, index=False, compression="gzip")
    summ = dict(rows=len(df), runtime_s=round(time.time() - t0, 1), args=vars(a),
                assertions="r_sec<=r_nom; BMR==BMR_hop feasibility; BMR r_nom>=BMR_hop r_nom; "
                           "BMR>=SPF_valid_budget r_nom when both feasible; hops>=Lmin -- all passed")
    json.dump(summ, open(os.path.splitext(os.path.splitext(a.out)[0])[0] + "_run.json", "w"), indent=2)
    print(summ)


if __name__ == "__main__":
    main()
