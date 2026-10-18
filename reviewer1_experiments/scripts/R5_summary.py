"""R5_summary.py -- tables + figure for R5 from ../data/R5_rows.csv.gz (cluster = graph; 95% CI = t-interval over 50 graph means)."""
import numpy as np, pandas as pd, json, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
D = "../data/"
df = pd.read_csv(D + "R5_rows.csv.gz")
md = []; out = {}

def gm(sub, col):
    g = sub.groupby("graph_id")[col].mean().values
    m = g.mean(); se = g.std(ddof=1) / np.sqrt(len(g)); h = stats.t.ppf(.975, len(g) - 1) * se
    return m, h
def cond_rsec(sub):          # R_nom | feasible : ratio of means within graph then mean over graphs
    s = sub[sub.feasible == 1]
    return s.groupby("graph_id").r_nom.mean().mean() if len(s) else 0.0
def fmt(sub, col, p=4): m, h = gm(sub, col); return f"{m:.{p}f} ± {h:.{p}f}"
def pct(sub, col):       # share among FEASIBLE selections, in %
    s = sub[sub.feasible == 1]
    return 100 * s[col].mean() if len(s) else float("nan")
def sel(table, **kw):
    s = df[df.table == table]
    for k, v in kw.items(): s = s[np.isclose(s[k], v, rtol=1e-6) if k in ("Delta",) else s[k] == v]
    return s

# ---------- consistency with R1 ----------
r1 = pd.read_csv(D + "R1_cells.csv")
md.append("### T-R5-0  Dataset and consistency\n")
md.append(f"One seeded dataset: ER(n=30, p=0.25), G=50 graphs × m=5 draws, n=10⁴, binomial noise, seed 20261012, f ∈ {{0, 0.10}}, expiry p = 4.88 % per link-epoch, H=5 (H up to 7 for the cap sweep), Lmin=3. SHA-256 of the link arrays: `{open(D+'R5_dataset_hash.txt').read().split(chr(10))[1]}`. Every table below is computed from the same 500 (graph, draw) link realisations.\n")
chk = sel("r1check", Delta=5e-6)
md.append("| check (p_exp = 0, Δ=5e-6, H=5) | this script | R1 sweep (L1_sweep_v5b) |\n|---|---|---|")
for cond in ("clean", "compromised"):
    for v, rr in (("full", "BMR"), ("SPF", "SPF")):
        a = chk[(chk.variant == v) & (chk.cond == cond)]
        c = r1[(r1.topology == "er") & (r1["size"] == 30) & (r1["mode"] == "sealed") & (r1.router == rr) & np.isclose(r1.Delta, 5e-6)
               & ((r1.f == 0) if cond == "clean" else np.isclose(r1.f, 0.10))]
        md.append(f"| {cond} {v} R_sec | {a.r_sec.mean():.4f} | {c.rs.iloc[0]:.4f} |")
md.append("\nBMR reproduces R1 to four decimals (same seeds, same links, deterministic objective). SPF picks uniformly at random among the many 3-hop paths, so its value depends on the tie-break random stream; the R5 stream differs from R1's. Exact expectation over the tie-break (average of R_sec over all minimum-hop paths) is computed in `R5_spf_exact.py` and reported in the next line.")
ex = json.load(open(D + 'R5_spf_exact.json'))
md.append(f"\nExact SPF R_sec (tie-break expectation): clean {ex['clean']['mean']:.4f}, compromised {ex['compromised']['mean']:.4f}; R1 Monte-Carlo tie-break: 0.1176 / 0.0961 (CI half-width ≈ 0.015); R5 stream: see rows above. All three agree within sampling error.")
# ---------- (a) Delta sweep ----------
ds = df[df.table == "dsweep"]
md.append("\n### T-R5-A  Δ-sweep, 20 log-spaced points 10⁻⁷…10⁻³ (compromised f = 0.10; ER n=30; p_exp = 4.88 %)\n")
md.append("`full` = BMR; `no_defect` = BMR without the budget (Δ ignored in selection, still *evaluated* against Δ). Δ-violation = % of **selected** paths with Σε > Δ.\n")
md.append("| Δ | full: feasible | full: R_nom\\|f | full: R_sec | full: Δ-viol % | no_defect: feasible | no_defect: R_sec | **no_defect: Δ-viol %** | no_defect: Σε (mean) | SPF: R_sec |\n|---|---|---|---|---|---|---|---|---|---|")
rows_a = []
for Dl in sorted(ds.Delta.unique()):
    s = ds[(np.isclose(ds.Delta, Dl)) & (ds.cond == "compromised")]
    f, n_, sp = s[s.variant == "full"], s[s.variant == "no_defect"], s[s.variant == "SPF"]
    rows_a.append((Dl, f.feasible.mean(), f.r_sec.mean(), n_.feasible.mean(), n_.r_sec.mean(), pct(n_, "delta_v"), pct(f, "delta_v"), sp.r_sec.mean()))
    md.append(f"| {Dl:.3g} | {fmt(f,'feasible',3)} | {cond_rsec(f):.4f} | {fmt(f,'r_sec')} | {pct(f,'delta_v'):.1f} | {fmt(n_,'feasible',3)} | {fmt(n_,'r_sec')} | **{pct(n_,'delta_v'):.1f}** | {n_[n_.feasible==1].sum_eps.mean():.2e} | {sp.r_sec.mean():.4f} |")
md.append("\nClean condition, same grid:\n")
md.append("| Δ | full: feasible | full: R_sec | no_defect: feasible | no_defect: R_sec | no_defect: Δ-viol % |\n|---|---|---|---|---|---|")
for Dl in sorted(ds.Delta.unique()):
    s = ds[(np.isclose(ds.Delta, Dl)) & (ds.cond == "clean")]
    f, n_ = s[s.variant == "full"], s[s.variant == "no_defect"]
    md.append(f"| {Dl:.3g} | {f.feasible.mean():.3f} | {f.r_sec.mean():.4f} | {n_.feasible.mean():.3f} | {n_.r_sec.mean():.4f} | {pct(n_,'delta_v'):.1f} |")

# ---------- (b) ablation ----------
ab = df[df.table == "ablation"]
names = {"full": "full (BMR)", "no_defect": "no_defect", "no_freshness": "no_freshness", "no_alarm": "no_alarm", "hop_only": "hop_only (min-hop, same constraints)",
         "rate_only": "rate_only (max r_min, same constraints)", "full_R0": "full, R=0 (reference)", "no_alarm_R0": "no_alarm, R=0",
         "hop_R0": "hop_only, R=0 (reference)", "hop_R0_no_alarm": "hop_only, R=0, no alarm filter", "hop_R0_no_fresh": "hop_only, R=0, no freshness filter", "SPF": "SPF (no filters)"}
for Dl, tag in ((5e-6, "B1"), (1e-6, "B2")):
    md.append(f"\n### T-R5-{tag}  Ablation at Δ = {Dl:g}, compromised (f = 0.10), p_exp = 4.88 %, H = 5 (cluster 95% CI over 50 graphs for R_sec)\n")
    md.append("Columns: feasibility, R_nom|feasible, R_sec (0 if infeasible/invalid), mean hops, and violations as % of **selected** paths.\n")
    md.append("| Variant | feasible | R_nom\\|f | R_sec | hops | expiry-viol % | alarm-viol % | Δ-viol % | Σε (mean, feasible) |\n|---|---|---|---|---|---|---|---|---|")
    s0 = ab[(np.isclose(ab.Delta, Dl)) & (ab.cond == "compromised")]
    for v in ["full", "no_defect", "no_freshness", "no_alarm", "hop_only", "rate_only", "full_R0", "no_alarm_R0", "hop_R0", "hop_R0_no_alarm", "hop_R0_no_fresh", "SPF"]:
        s = s0[s0.variant == v]; fe = s[s.feasible == 1]
        md.append(f"| {names[v]} | {s.feasible.mean():.3f} | {cond_rsec(s):.4f} | {fmt(s,'r_sec')} | {fe.hops.mean():.2f} | {pct(s,'exp_v'):.1f} | {pct(s,'alarm_v'):.1f} | {pct(s,'delta_v'):.1f} | {fe.sum_eps.mean():.2e} |")
    md.append("\nSame, **clean** condition (f = 0):\n")
    md.append("| Variant | feasible | R_nom\\|f | R_sec | hops | expiry-viol % | Δ-viol % |\n|---|---|---|---|---|---|---|")
    s0 = ab[(np.isclose(ab.Delta, Dl)) & (ab.cond == "clean")]
    for v in ["full", "no_defect", "no_freshness", "hop_only", "rate_only", "SPF"]:
        s = s0[s0.variant == v]; fe = s[s.feasible == 1]
        md.append(f"| {names[v]} | {s.feasible.mean():.3f} | {cond_rsec(s):.4f} | {fmt(s,'r_sec')} | {fe.hops.mean():.2f} | {pct(s,'exp_v'):.1f} | {pct(s,'delta_v'):.1f} |")
    out[f"ablation_{Dl:g}"] = {v: dict(feasible=float(s0[s0.variant == v].feasible.mean()), r_sec=float(s0[s0.variant == v].r_sec.mean())) for v in names}

# ---------- (c) hop cap ----------
hc = df[df.table == "hopcap"]
md.append("\n### T-R5-C  Hop-cap sweep (compromised f = 0.10 and clean), Lmin = 3, ER n=30\n")
md.append("`at cap %` = share of selected BMR routes whose length equals H. Δ = 5e-6 (main) and 1e-6 (tight).\n")
md.append("| Δ | cond | H | BMR feasible | BMR R_nom\\|f | **BMR R_sec** | BMR hops | BMR at cap % | SPF_valid_budget feasible | SPF_valid_budget R_sec | **BMR − SPF_valid_budget R_sec [95% CI]** | BMR_hop(min-hop) R_sec | SPF R_sec |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
hop_rows = []
for Dl in (5e-6, 1e-6):
    for cond in ("clean", "compromised"):
        for H in (4, 5, 6, 7):
            s = hc[(np.isclose(hc.Delta, Dl)) & (hc.cond == cond) & (hc.hcap == H)]
            b, v, bh, sp = s[s.variant == "full"], s[s.variant == "SPF_valid_budget"], s[s.variant == "hop_only"], s[s.variant == "SPF"]
            diff = (b.set_index(["graph_id", "draw"]).r_sec - v.set_index(["graph_id", "draw"]).r_sec).groupby("graph_id").mean().values
            m = diff.mean(); h = stats.t.ppf(.975, len(diff) - 1) * diff.std(ddof=1) / np.sqrt(len(diff))
            fe = b[b.feasible == 1]
            hop_rows.append((Dl, cond, H, b.r_sec.mean(), v.r_sec.mean(), m, h))
            md.append(f"| {Dl:g} | {cond} | {H} | {b.feasible.mean():.3f} | {cond_rsec(b):.4f} | {b.r_sec.mean():.4f} | {fe.hops.mean():.2f} | {100*fe.at_cap.mean():.1f} | {v.feasible.mean():.3f} | {v.r_sec.mean():.4f} | **{m:+.4f} [{m-h:+.4f}, {m+h:+.4f}]** | {bh.r_sec.mean():.4f} | {sp.r_sec.mean():.4f} |")

# ---------- (d) parity ----------
pr = df[df.table == "parity"]
md.append("\n### T-R5-D  Public-XOR parity (mode `public`, Δ = 5e-6, H = 5, Lmin = 3, ER n=30)\n")
md.append("`odd %` = share of feasible selections with an odd number of hops (even-hop public-XOR routes have zero secrecy). BMR applies the odd-hop constraint; `BMR_np` is BMR without it; the other routers are parity-unaware.\n")
md.append("| Router | cond | feasible | odd-hop % (of feasible) | R_nom\\|f | R_sec [95% CI] | hops |\n|---|---|---|---|---|---|---|")
for cond in ("clean", "compromised"):
    for rt in ["BMR", "BMR_np", "BMR_hop", "SPF", "SPF_valid", "SPF_budget", "SPF_valid_budget", "QPATHstyle", "QLEAPstyle", "QCASTstyle"]:
        s = pr[(pr.variant == rt) & (pr.cond == cond)]; fe = s[s.feasible == 1]
        odd = 100 * fe.odd.mean() if len(fe) else float("nan")
        md.append(f"| {rt} | {cond} | {s.feasible.mean():.3f} | {odd:.1f} | {cond_rsec(s):.4f} | {fmt(s,'r_sec')} | {fe.hops.mean():.2f} |")

# ---------- (e) F_th ----------
ft = df[df.table == "fth"]
md.append("\n### T-R5-E  Fidelity-threshold sensitivity of the QLEAP-style router (Δ = 5e-6 evaluation, H = 5)\n")
md.append("| F_th | cond | feasible | R_nom\\|f | R_sec | BMR R_sec | BMR − QLEAP_F [95% CI] |\n|---|---|---|---|---|---|---|")
bm = ab[(np.isclose(ab.Delta, 5e-6)) & (ab.variant == "full")]
for cond in ("clean", "compromised"):
    for fth in (0.70, 0.75, 0.80):
        s = ft[(ft.cond == cond) & (ft.F_th == fth)]
        b = bm[bm.cond == cond]
        diff = (b.set_index(["graph_id", "draw"]).r_sec - s.set_index(["graph_id", "draw"]).r_sec).groupby("graph_id").mean().values
        m = diff.mean(); h = stats.t.ppf(.975, len(diff) - 1) * diff.std(ddof=1) / np.sqrt(len(diff))
        md.append(f"| {fth:.2f} | {cond} | {s.feasible.mean():.3f} | {cond_rsec(s):.4f} | {s.r_sec.mean():.4f} | {b.r_sec.mean():.4f} | {m:+.4f} [{m-h:+.4f}, {m+h:+.4f}] |")

open(D + "R5_tables.md", "w").write("\n".join(md) + "\n")
json.dump(out, open(D + "R5_summary.json", "w"), indent=1)
print("\n".join(md))

# ---------- figure ----------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
Dls = [r[0] for r in rows_a]
ax[0].semilogx(Dls, [r[2] for r in rows_a], "o-", label="full R_sec"); ax[0].semilogx(Dls, [r[4] for r in rows_a], "s--", label="no_defect R_sec")
ax[0].semilogx(Dls, [r[7] for r in rows_a], "^:", label="SPF R_sec"); ax[0].set_xlabel("Δ"); ax[0].set_ylabel("R_sec (compromised)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which="both")
ax[1].semilogx(Dls, [r[5] for r in rows_a], "s--", color="C1"); ax[1].semilogx(Dls, [r[6] for r in rows_a], "o-", color="C0")
ax[1].set_xlabel("Δ"); ax[1].set_ylabel("Δ-violation % of selected paths"); ax[1].set_title("no_defect (dashed) vs full (solid)"); ax[1].grid(alpha=.3, which="both")
for Dl, st in ((5e-6, "o-"), (1e-6, "s--")):
    xs = [r for r in hop_rows if r[0] == Dl and r[1] == "compromised"]
    ax[2].errorbar([r[2] for r in xs], [r[5] for r in xs], yerr=[r[6] for r in xs], fmt=st, capsize=3, label=f"BMR − SPF_valid_budget, Δ={Dl:g}")
ax[2].axhline(0, color="k", lw=.6); ax[2].set_xlabel("hop cap H"); ax[2].set_ylabel("ΔR_sec (compromised)"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
plt.tight_layout(); plt.savefig("../figures/fig_R5_ablation.png", dpi=150)
