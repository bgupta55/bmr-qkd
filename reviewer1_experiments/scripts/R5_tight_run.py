"""R5_tight_run.py -- repeat the full R5 pipeline (L7_ablation_v5b.run_graph, unchanged) with the tight radius
eta_tight = sqrt(8 ln(1/delta)/n) instead of eta_union = sqrt(32 ln(8/delta)/n). Same seeds, same links, same code.
Output: ../data/R5_rows_tight.csv.gz"""
import sys, os, math, time, argparse
from multiprocessing import Pool
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B, L7_ablation_v5b as L

def tau_tight(n, delta=1e-6):
    return math.sqrt(8.0 * math.log(1.0 / delta) / n)

def work(g):
    B.tau = tau_tight                      # make_links looks tau up in B at call time
    return L.run_graph(g)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--G", type=int, default=50); ap.add_argument("--procs", type=int, default=10)
    a = ap.parse_args(); t0 = time.time(); rows = []
    with Pool(a.procs) as p:
        for g, r, h in p.imap_unordered(work, range(a.G)): rows.extend(r)
    pd.DataFrame(rows).to_csv("../data/R5_rows_tight.csv.gz", index=False, compression="gzip")
    print(len(rows), "rows", round(time.time() - t0, 1), "s; eta_tight(1e4) =", tau_tight(10000))
