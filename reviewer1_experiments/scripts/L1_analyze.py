"""
L1_analyze.py -- R6: graph-level (cluster) statistics on the R1 trial data.

Pipeline (K10): trial -> graph-level means (graph = independent unit)
  -> paired differences BMR - router -> cluster bootstrap 95% CI (B = 20,000)
  -> sign-flip permutation p (MC 200,000) -> d_z -> k/G wins/ties/losses
  -> Holm over the declared family -> one-way ICC / DEFF + mixed-model check.

Usage:  python3 L1_analyze.py --in ../data/R1_trials.csv.gz --outdir ../data
"""
import argparse, json, os
import numpy as np
import pandas as pd

B_BOOT = 20000
N_FLIP = 200000
TIE = 1e-9
KEY = ["topology", "size", "mode", "cond", "f", "Delta"]


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * p[i])
        adj[i] = min(1.0, run)
    return adj


def icc_oneway(x):
    """x: (G, m) array.  Returns ICC(1), DEFF = 1 + (m-1) ICC."""
    G, m = x.shape
    gm = x.mean(1)
    mu = x.mean()
    msb = m * ((gm - mu) ** 2).sum() / (G - 1)
    msw = ((x - gm[:, None]) ** 2).sum() / (G * (m - 1))
    icc = (msb - msw) / (msb + (m - 1) * msw) if (msb + (m - 1) * msw) > 0 else 0.0
    icc = max(0.0, icc)
    return float(icc), float(1 + (m - 1) * icc), float(msb), float(msw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default="../data/R1_trials.csv.gz")
    ap.add_argument("--outdir", default="../data")
    a = ap.parse_args()
    df = pd.read_csv(a.inp)
    df["r_nom_cond"] = df["r_nom"]
    # ---- graph-level means (graph = independent unit)
    gcols = KEY + ["router", "graph_id"]
    pg = (df.groupby(gcols)
            .agg(feas=("feasible", "mean"), rn=("r_nom", "mean"), rs=("r_sec", "mean"),
                 hops=("hops", "mean"), dviol=("delta_violated", "mean"),
                 alarm=("alarm_on_path", "mean"), Fp=("F_P", "mean"),
                 cost=("cost_pairs", "mean"), tput=("throughput_proxy", "mean"),
                 sum_eps=("sum_eps", "mean"), n_tr=("draw", "count"))
            .reset_index())
    os.makedirs(a.outdir, exist_ok=True)
    pg.to_csv(os.path.join(a.outdir, "R1_per_graph_means.csv.gz"), index=False, compression="gzip")

    rng_flip = np.random.default_rng(777)
    G_ref = pg.graph_id.nunique()
    SIGN = rng_flip.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(N_FLIP, G_ref))

    cell_rows, diff_rows = [], []
    cells = pg[KEY].drop_duplicates().values.tolist()
    for ci, cell in enumerate(cells):
        topo, size, mode, cond, f, Delta = cell
        sub = pg[(pg.topology == topo) & (pg["size"] == size) & (pg["mode"] == mode) &
                 (pg.cond == cond) & (pg.f == f) & (pg.Delta == Delta)]
        routers = sorted(sub.router.unique())
        G = sub.graph_id.nunique()
        gids = np.sort(sub.graph_id.unique())
        rng = np.random.default_rng(np.random.SeedSequence([20261012, ci]))
        idx = rng.integers(0, G, size=(B_BOOT, G))
        arr = {}
        for r in routers:
            x = sub[sub.router == r].set_index("graph_id").loc[gids]
            arr[r] = {c: x[c].to_numpy() for c in
                      ["feas", "rn", "rs", "hops", "dviol", "alarm", "Fp", "cost", "tput", "sum_eps"]}
        for r in routers:
            A = arr[r]
            bf = A["feas"][idx].mean(1)
            brs = A["rs"][idx].mean(1)
            brn = A["rn"][idx].sum(1) / np.maximum(A["feas"][idx].sum(1), 1e-12)
            rn_c = A["rn"].sum() / max(A["feas"].sum(), 1e-12)
            hops_c = (A["hops"].sum() / max(A["feas"].sum(), 1e-12))
            cell_rows.append(dict(
                topology=topo, size=size, mode=mode, cond=cond, f=f, Delta=Delta, router=r, G=G,
                feas=A["feas"].mean(), feas_lo=np.percentile(bf, 2.5), feas_hi=np.percentile(bf, 97.5),
                rn_given_feas=rn_c, rn_lo=np.percentile(brn, 2.5), rn_hi=np.percentile(brn, 97.5),
                rs=A["rs"].mean(), rs_lo=np.percentile(brs, 2.5), rs_hi=np.percentile(brs, 97.5),
                hops_given_feas=hops_c, delta_viol=A["dviol"].mean(),
                alarm_on_path=A["alarm"].mean(), F_P=A["Fp"].mean(), cost_pairs=A["cost"].mean(),
                throughput=A["tput"].mean(), sum_eps=A["sum_eps"].mean()))
        base = "BMR"
        for r in routers:
            if r == base:
                continue
            d = arr[base]["rs"] - arr[r]["rs"]
            bd = d[idx].mean(1)
            obs = d.mean()
            sd = d.std(ddof=1)
            if np.all(np.abs(d) < TIE):
                p = 1.0
                dz = 0.0
            else:
                stat = np.abs(SIGN[:, :G] @ d.astype(np.float32)) / G
                p = (1 + np.sum(stat >= abs(obs) - 1e-12)) / (N_FLIP + 1)
                dz = obs / sd if sd > 0 else float("inf")
            wins = int((d > TIE).sum()); losses = int((d < -TIE).sum()); ties = G - wins - losses
            # same for R_nom on common feasible?  -> report R_nom difference on all trials too
            dn = arr[base]["rn"] - arr[r]["rn"]
            diff_rows.append(dict(
                topology=topo, size=size, mode=mode, cond=cond, f=f, Delta=Delta, comparator=r, G=G,
                diff_rs=obs, lo=np.percentile(bd, 2.5), hi=np.percentile(bd, 97.5), p_signflip=p,
                d_z=dz, wins=wins, ties=ties, losses=losses,
                diff_feas=(arr[base]["feas"] - arr[r]["feas"]).mean(),
                diff_rn_uncond=dn.mean(),
                degenerate=bool(np.all(np.abs(d) < TIE))))
    cells_df = pd.DataFrame(cell_rows)
    diffs = pd.DataFrame(diff_rows)

    # ---- Holm: family = every non-degenerate sealed-mode BMR-vs-router test
    fam = (diffs["mode"] == "sealed") & (~diffs.degenerate)
    diffs["p_holm_all"] = np.nan
    diffs.loc[fam, "p_holm_all"] = holm(diffs.loc[fam, "p_signflip"].values)
    main_fam = fam & diffs.topology.isin(["er", "ba"]) & (diffs.Delta == 5e-6) & (diffs.f.isin([0.0, 0.1]))
    diffs["p_holm_main"] = np.nan
    diffs.loc[main_fam, "p_holm_main"] = holm(diffs.loc[main_fam, "p_signflip"].values)
    pub = (diffs["mode"].str.startswith("public")) & (~diffs.degenerate)
    diffs["p_holm_public"] = np.nan
    diffs.loc[pub, "p_holm_public"] = holm(diffs.loc[pub, "p_signflip"].values)

    cells_df.to_csv(os.path.join(a.outdir, "R1_cells.csv"), index=False)
    diffs.to_csv(os.path.join(a.outdir, "R1_diffs.csv"), index=False)

    # ---- decomposition (telescoping, two orders) in er/ba n=30, 50, sealed, Delta=5e-6 and 1e-6
    dec = []
    for topo in ["er", "ba"]:
        for size in [20, 30, 50]:
            for cond, f in [("clean", 0.0), ("compromised", 0.1)]:
                for Delta in [5e-6, 1e-6]:
                    sub = pg[(pg.topology == topo) & (pg["size"] == size) & (pg["mode"] == "sealed") &
                             (pg.f == f) & (pg.Delta == Delta)]
                    gids = np.sort(sub.graph_id.unique()); G = len(gids)
                    rs = {r: sub[sub.router == r].set_index("graph_id").loc[gids]["rs"].to_numpy()
                          for r in sub.router.unique()}
                    rng = np.random.default_rng(np.random.SeedSequence([5150, size, int(f * 100), int(Delta * 1e7)]))
                    idx = rng.integers(0, G, size=(B_BOOT, G))

                    def ci(v):
                        b = v[idx].mean(1)
                        return float(v.mean()), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))
                    steps = {
                        "order1_budget": rs["SPF_budget"] - rs["SPF"],
                        "order1_validity": rs["SPF_valid_budget"] - rs["SPF_budget"],
                        "order2_validity": rs["SPF_valid"] - rs["SPF"],
                        "order2_budget": rs["SPF_valid_budget"] - rs["SPF_valid"],
                        "objective_tiebreak(BMR_hop-SPFvb)": rs["BMR_hop"] - rs["SPF_valid_budget"],
                        "objective_hopwindow(BMR-BMR_hop)": rs["BMR"] - rs["BMR_hop"],
                        "objective_total(BMR-SPFvb)": rs["BMR"] - rs["SPF_valid_budget"],
                        "total(BMR-SPF)": rs["BMR"] - rs["SPF"],
                    }
                    row = dict(topology=topo, size=size, cond=cond, Delta=Delta, G=G)
                    for k, v in steps.items():
                        m_, lo, hi = ci(v)
                        row[k] = m_; row[k + "_lo"] = lo; row[k + "_hi"] = hi
                    row["validity_mean_of_orders"] = 0.5 * (row["order1_validity"] + row["order2_validity"])
                    row["budget_mean_of_orders"] = 0.5 * (row["order1_budget"] + row["order2_budget"])
                    dec.append(row)
    dec_df = pd.DataFrame(dec)
    dec_df.to_csv(os.path.join(a.outdir, "R1_decomposition.csv"), index=False)

    # ---- ICC / DEFF at trial level for two cells (BMR r_sec, SPF r_sec, BMR-SPF)
    icc = {}
    for topo in ["er", "ba"]:
        sub = df[(df.topology == topo) & (df["size"] == 30) & (df["mode"] == "sealed") &
                 (df.Delta == 5e-6) & (df.f == 0.1)]
        out = {}
        piv = {}
        for r in ["BMR", "SPF", "SPF_valid_budget"]:
            x = sub[sub.router == r].pivot(index="graph_id", columns="draw", values="r_sec").to_numpy()
            piv[r] = x
            i_, de, msb, msw = icc_oneway(x)
            out[f"{r}_r_sec"] = dict(ICC=i_, DEFF=de, MSB=msb, MSW=msw)
        x = piv["BMR"] - piv["SPF"]
        i_, de, msb, msw = icc_oneway(x)
        out["BMR-SPF_r_sec"] = dict(ICC=i_, DEFF=de, MSB=msb, MSW=msw)
        # mixed-model cross-check (random intercept per graph)
        try:
            import statsmodels.formula.api as smf
            long = pd.DataFrame({"g": np.repeat(np.arange(x.shape[0]), x.shape[1]), "y": x.ravel()})
            md = smf.mixedlm("y ~ 1", long, groups=long["g"]).fit(reml=True)
            vg = float(md.cov_re.iloc[0, 0]); ve = float(md.scale)
            out["BMR-SPF_mixedlm"] = dict(var_graph=vg, var_resid=ve, ICC=vg / (vg + ve) if vg + ve > 0 else 0.0,
                                         intercept=float(md.fe_params.iloc[0]), se=float(md.bse.iloc[0]))
        except Exception as ex:  # pragma: no cover
            out["BMR-SPF_mixedlm"] = {"error": str(ex)}
        # naive vs cluster standard error of the mean difference
        n_tr = x.size
        naive_se = float(x.std(ddof=1) / np.sqrt(n_tr))
        clus_se = float(x.mean(1).std(ddof=1) / np.sqrt(x.shape[0]))
        out["BMR-SPF_se"] = dict(naive_se_n250=naive_se, cluster_se_G50=clus_se, ratio=clus_se / naive_se)
        icc[f"{topo}30_sealed_D5e-6_f0.10"] = out
    json.dump(icc, open(os.path.join(a.outdir, "R1_icc_deff.json"), "w"), indent=2)

    # ---- false-alarm check, binomial vs legacy gaussian(2/sqrt n)
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import B_core as BC
    fa = {}
    for noise in ["binomial", "gauss02"]:
        rng = np.random.default_rng(31337)
        tot = 0; al = 0; per = {}
        for _ in range(100):
            lk = BC.make_links(rng, 2000, 10000, f=0.0, noise=noise)
            tot += len(lk["alarm"]); al += int(lk["alarm"].sum())
            for c in range(3):
                m_ = lk["cls"] == c
                per.setdefault(c, [0, 0]); per[c][0] += int(m_.sum()); per[c][1] += int(lk["alarm"][m_].sum())
        fa[noise] = dict(links=tot, alarms=al, rate=al / tot,
                         per_class={["high", "std", "deg"][c]: v[1] / v[0] for c, v in per.items()})
    # detection of the attack link (S=2.418)
    rng = np.random.default_rng(9)
    lk = BC.make_links(rng, 200000, 10000, f=1.0, noise="binomial")
    fa["attack_detection_prob_binomial"] = float(lk["alarm"].mean())
    fa["attack_gate_pass_prob_binomial"] = float(lk["gate"].mean())
    fa["attack_mean_r_e"] = float(lk["r"].mean())
    fa["attack_frac_r_e_below_R"] = float((lk["r"] < 0.10).mean())
    json.dump(fa, open(os.path.join(a.outdir, "R1_false_alarm_check.json"), "w"), indent=2)
    print("cells", len(cells_df), "diffs", len(diffs), "holm_all_family", int(fam.sum()), "main_family", int(main_fam.sum()))
    print(json.dumps(fa, indent=1))


if __name__ == "__main__":
    main()
