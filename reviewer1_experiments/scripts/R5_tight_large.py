"""R5_tight_large.py -- larger sample (G=400 graphs, same seeding rule, same code) for the alarm's false-certificate rate under eta_union and eta_tight.
Graphs 0..49 are identical to R5; 50..399 are new graphs from the same seed rule. Output: ../data/R5_large.json"""
import sys, os, math, json, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B, L7_ablation_v5b as L
def tt(n, d=1e-6): return math.sqrt(8*math.log(1/d)/n)
def work(a):
    g, rad = a
    if rad == "tight": B.tau = tt
    return L.run_graph(g)[1]
if __name__ == "__main__":
    G = 400; out = {}
    for rad in ("union", "tight"):
        with Pool(10) as p: rows = [r for rs in p.map(work, [(g, rad) for g in range(G)]) for r in rs]
        d = pd.DataFrame(rows); d = d[(d.table == "ablation") & np.isclose(d.Delta, 5e-6) & (d.cond == "compromised")]
        res = {}
        for v in ("full", "no_alarm", "SPF", "hop_R0_no_alarm"):
            s = d[d.variant == v]; f = s[s.feasible == 1]
            # per-graph false-certificate rate (cluster = graph)
            fc = f.groupby("graph_id").touches_comp.mean().reindex(range(G)).fillna(0).values
            m = fc.mean(); h = stats.t.ppf(.975, G-1)*fc.std(ddof=1)/np.sqrt(G)
            res[v] = dict(n_feasible=int(len(f)), n_false_cert=int(f.touches_comp.sum()), rate_pct=100*m, hw_pct=100*h)
        # R_sec difference full - no_alarm
        x = d[d.variant == "full"].set_index(["graph_id", "draw"]); y = d[d.variant == "no_alarm"].set_index(["graph_id", "draw"])
        dr = (x.r_sec - y.r_sec).groupby("graph_id").mean().values
        res["rsec_full_minus_noalarm"] = dr.mean()
        res["feas_full"] = float(x.feasible.mean()); res["feas_noalarm"] = float(y.feasible.mean())
        out[rad] = res; print(rad, json.dumps(res, indent=1))
    json.dump(out, open("../data/R5_large.json", "w"), indent=1)
