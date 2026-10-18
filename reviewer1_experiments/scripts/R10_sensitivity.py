"""
R10_sensitivity.py -- R10: (a) S0 sensitivity of the R1 cell {ER n=30, sealed, f in {0,0.10}} for S0 in {2.55, 2.60, 2.687}, S1 = S0 - 0.27;
(b) selected paths under r_key(S_lb) vs h(S_lb); (c) minimum honest audit length.
S0 shifts the whole honest class model: CLASS_S = (2.720, 2.687, 2.640) + (S0 - 2.687); alarm reference of each honest link = its class S;
attacked link: S_true = S0 - 0.27, alarm reference S0.  All else (graphs, links' uniforms, eps classes) identical to R1 -> paired design.
Same seeds as R1, p_exp = 0, ER p = 0.25, G = 50, m = 5, n = 10^4, Lmin = 3, H = 5.
"""
import json, math, os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B

SEED = 20261012; R = 0.10


def h2(p):
    p = np.clip(p, 1e-15, 1 - 1e-15); return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


def rkey(S):
    S = np.minimum(np.asarray(S, float), 2 * math.sqrt(2)); out = np.zeros_like(S)
    m = S > 2.0
    Q = (1 - S[m] / (2 * math.sqrt(2))) / 2
    out[m] = np.maximum(0.0, 1 - h2(Q) - h2((1 + np.sqrt((S[m] / 2) ** 2 - 1)) / 2))
    return out


def patch(S0):
    sh = S0 - 2.687
    B.CLASS_S = np.array([2.720, 2.687, 2.640]) + sh
    B.S_ATTACK = S0 - 0.27
    B.S0_REF = S0


def work(args):
    S0, g = args
    patch(S0)
    rg = np.random.default_rng(np.random.SeedSequence([SEED, 30, 1, g, 0]))
    edges = B.build_graph("er", 30, rg); E = len(edges)
    adj = B.adj_matrix(30, edges); L, hops = B.enumerate_paths(adj, 0, 29, 3, 5)
    rows = []
    if len(hops) == 0:
        return rows
    for draw in range(5):
        for f in (0.0, 0.10):
            rl = np.random.default_rng(np.random.SeedSequence([SEED, 30, 1, g, draw, 1]))
            rc = np.random.default_rng(np.random.SeedSequence([SEED, 30, 1, g, draw, 2, int(f * 1000)]))
            lk = B.make_links(rl, E, 10000, f=f, noise="binomial", p_exp=0.0, rng_comp=rc)
            M = B.path_metrics(L, lk)
            lk_key = dict(lk); lk_key["r"] = rkey(lk["S_lb"])
            Mk = B.path_metrics(L, lk_key)
            # floor: matched (same S_lb threshold) and numeric (same number R)
            S_floor = 2.0 + 0.0
            # invert h(S)=0.10 by bisection
            lo, hi = 2.0, 2.828
            for _ in range(60):
                mid = (lo + hi) / 2
                if B.h(mid) < R: lo = mid
                else: hi = mid
            S_R = hi
            for Delta in (5e-6, 1e-6):
                rs = np.random.default_rng(np.random.SeedSequence([SEED, 30, 1, g, draw, 3]))
                routers = ["BMR", "SPF", "SPF_valid_budget", "SPF_budget", "BMR_hop"]
                sel, M2, F, tput = B.select_all(L, hops, lk, Delta, R, "sealed", rs, routers, M=M)
                rr = {r["router"]: r for r in B.trial_rows(sel, M2, F, tput, L, hops, Delta, "sealed", 3)}
                # (b) BMR selection with r_key objective
                base = M["allvalid"] & (M["seps"] <= Delta * (1 + 1e-9))
                m_h = base & (M["rmin"] >= R)
                Smin = np.full(len(hops), np.inf)
                Lp = np.where(L < 0, E, L)
                Smin = np.concatenate([lk["S_lb"], [np.inf]])[Lp].min(1)
                m_match = base & (Smin >= S_R - 1e-12)                       # key objective, same S_lb floor
                m_num = base & (Mk["rmin"] >= R)                             # key objective, same numeric floor R
                rs2 = np.random.default_rng(1)
                i_h = B._pick(m_h, [-M["rmin"], hops], rs2)
                i_k = B._pick(m_match, [-Mk["rmin"], hops], np.random.default_rng(1))
                i_kn = B._pick(m_num, [-Mk["rmin"], hops], np.random.default_rng(1))
                rk_sel = float(Mk["rmin"][i_k]) if i_k >= 0 else 0.0
                rows.append(dict(S0=S0, graph_id=g, draw=draw, cond="clean" if f == 0 else "compromised", Delta=Delta,
                                 feas_BMR=rr["BMR"]["feasible"], rsec_BMR=rr["BMR"]["r_sec"], rnom_BMR=rr["BMR"]["r_nom"],
                                 feas_SVB=rr["SPF_valid_budget"]["feasible"], rsec_SVB=rr["SPF_valid_budget"]["r_sec"],
                                 rsec_SPF=rr["SPF"]["r_sec"], rsec_SPFB=rr["SPF_budget"]["r_sec"], rsec_BMRhop=rr["BMR_hop"]["r_sec"],
                                 hops_BMR=rr["BMR"]["hops"],
                                 sel_h=int(i_h), sel_key_matched=int(i_k), sel_key_numfloor=int(i_kn),
                                 same_matched=int(i_h == i_k), same_numfloor=int(i_h == i_kn),
                                 feas_h=int(i_h >= 0), feas_key_matched=int(i_k >= 0), feas_key_numfloor=int(i_kn >= 0),
                                 rkey_sel=rk_sel, rkey_sel_pos=int(rk_sel > 0),
                                 n_valid=int(base.sum())))
    return rows


if __name__ == "__main__":
    tasks = [(S0, g) for S0 in (2.687, 2.60, 2.55) for g in range(50)]
    with Pool(10) as p:
        res = p.map(work, tasks)
    df = pd.DataFrame([r for rs in res for r in rs])
    df.to_csv("../data/R10_rows.csv.gz", index=False, compression="gzip")
    md = []
    def ci(x):
        g = np.asarray(x); return g.mean(), stats.t.ppf(.975, len(g) - 1) * g.std(ddof=1) / math.sqrt(len(g))
    md.append("### T-R10-A  S₀ sensitivity of the R1 cell {ER n=30, sealed, Lmin = 3, H = 5, 50 graphs × 5 draws, n = 10⁴}\n")
    md.append("Honest class model shifted by (S₀ − 2.687); attacked link S† = S₀ − 0.27. `adv` = BMR − comparator in R_sec (cluster 95% CI over graphs). Feasibility is that of BMR.\n")
    md.append("| Δ | cond | S₀ | S† | BMR feasible | BMR R_nom\\|f | **BMR R_sec** | SPF_valid_budget R_sec | **BMR − SPF_valid_budget** [CI] | SPF R_sec | **BMR − SPF** [CI] | SPF_budget R_sec | BMR − SPF_budget [CI] |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for Delta in (5e-6, 1e-6):
        for cond in ("clean", "compromised"):
            for S0 in (2.687, 2.60, 2.55):
                s = df[(df.S0 == S0) & (df.cond == cond) & np.isclose(df.Delta, Delta)]
                gm = s.groupby("graph_id").mean(numeric_only=True)
                rnom = s[s.feas_BMR == 1].groupby("graph_id").rnom_BMR.mean().mean()
                a = ci(gm.rsec_BMR - gm.rsec_SVB); b = ci(gm.rsec_BMR - gm.rsec_SPF); c = ci(gm.rsec_BMR - gm.rsec_SPFB)
                md.append(f"| {Delta:g} | {cond} | {S0} | {S0-0.27:.3f} | {gm.feas_BMR.mean():.3f} | {rnom:.4f} | {gm.rsec_BMR.mean():.4f} | {gm.rsec_SVB.mean():.4f} | **{a[0]:+.4f} [{a[0]-a[1]:+.4f}, {a[0]+a[1]:+.4f}]** | {gm.rsec_SPF.mean():.4f} | **{b[0]:+.4f} [{b[0]-b[1]:+.4f}, {b[0]+b[1]:+.4f}]** | {gm.rsec_SPFB.mean():.4f} | {c[0]:+.4f} [{c[0]-c[1]:+.4f}, {c[0]+c[1]:+.4f}] |")
    md.append("\n### T-R10-B  Selected paths under r_key(S_lb) versus h(S_lb) (BMR objective: maximise the bottleneck; same validity, Δ-budget and tie-break)\n")
    md.append("`matched floor`: key objective with the same S_lb floor as R = 0.10 (S_lb ≥ h⁻¹(0.10)). `numeric floor`: the same number R = 0.10 applied to r_key (a different, stricter S_lb floor). Identical = same path index. `r_key>0` = share of selected key-paths with positive device-independent key rate.\n")
    md.append("| Δ | S₀ | cond | feasible (h) | feasible (key, matched) | **identical, matched floor** | feasible (key, numeric R) | identical, numeric floor | mean r_key of selected path | selected with r_key > 0 |\n|---|---|---|---|---|---|---|---|---|---|")
    for Delta in (5e-6,):
        for S0 in (2.687, 2.60, 2.55):
            for cond in ("clean", "compromised"):
                s = df[(df.S0 == S0) & (df.cond == cond) & np.isclose(df.Delta, Delta)]
                both = s[(s.feas_h == 1) | (s.feas_key_matched == 1)]
                same = both.same_matched.mean() if len(both) else float("nan")
                both2 = s[(s.feas_h == 1) | (s.feas_key_numfloor == 1)]
                same2 = both2.same_numfloor.mean() if len(both2) else float("nan")
                sel = s[s.feas_key_matched == 1]
                md.append(f"| {Delta:g} | {S0} | {cond} | {s.feas_h.mean():.3f} | {s.feas_key_matched.mean():.3f} | **{100*same:.1f}%** ({int(both.same_matched.sum())}/{len(both)}) | {s.feas_key_numfloor.mean():.3f} | {100*same2:.1f}% | {sel.rkey_sel.mean():.4f} | {100*sel.rkey_sel_pos.mean():.1f}% |")
    tau = lambda n: math.sqrt(32 * math.log(8 / 1e-6) / n)
    md.append("\n### T-R10-C  r_key(S) and the minimum honest audit length (S_lb = S₀ − τ(n) ≥ 2  ⇔  n ≥ 32 ln(8/δ)/(S₀−2)², δ = 10⁻⁶)\n")
    md.append("| S | r_key(S) | h(S) |  | S₀ | min honest n |\n|---|---|---|---|---|---|")
    for S, S0 in zip((2.687, 2.60, 2.55, 2.49, 2.418), (2.687, 2.60, 2.55, 2.49, 2.418)):
        nmin = math.ceil(32 * math.log(8 / 1e-6) / (S0 - 2) ** 2)
        md.append(f"| {S} | {float(rkey(np.array([S]))[0]):.3f} | {B.h(S):.3f} |  | {S0} | {nmin:,} |")
    open("../data/R10_tables.md", "w").write("\n".join(md) + "\n"); print("\n".join(md))
