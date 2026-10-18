"""R2_summary.py -- exponent fits (bootstrap CI over graphs) + figure for R2."""
import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from L2_alg2_v5b import fit_exponent
D = "../data/"
def rows(fn): return json.load(open(D + fn))["rows"]
t1, t2, st = rows("R2_timing_t1.json"), rows("R2_timing_t2.json"), rows("R2_timing_stress_t1.json")
out, md = {}, []
def sel(R, fam, var): return sorted([r for r in R if r["family"] == fam and r["variant"] == var and r.get("graphs")], key=lambda r: r["size"])
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for name, R, var, fams in [("typical t=1 (Δ=5e-6)", t1, "poly", ["er", "ba"]), ("typical t=2", t2, "poly", ["er", "ba"]),
                           ("worst-case t=1 (d unreachable, Δ=1)", st, "poly", ["er", "ba"]), ("typical t=1 DFS variant", t1, "dfs", ["er", "ba"])]:
    for fam in fams:
        S = sel(R, fam, var)
        sizes = [r["size"] for r in S]; med = [r["median_ms"] for r in S]
        # fit on sizes >= 30 (below that, timer/constant overhead dominates)
        k = [i for i, s in enumerate(sizes) if s >= 30]
        sl, lo, hi = fit_exponent([sizes[i] for i in k], [med[i] for i in k], [S[i]["per_graph_ms"] for i in k])
        wk = [S[i]["median_work"] for i in k]
        slw, lw, hw = fit_exponent([sizes[i] for i in k], wk, [[w] * 2 for w in wk], B_boot=10)
        out[f"{name}|{fam}"] = dict(sizes=sizes, median_ms=med, slope=sl, ci=[lo, hi], work_slope=slw)
        md.append(f"| {name} | {fam.upper()} | {sl:.2f} [{lo:.2f}, {hi:.2f}] | {slw:.2f} | " +
                  " / ".join(f"{m:.3g}" for m in med) + " |")
        a = ax[0] if "worst" not in name else ax[1]
        if "DFS" in name or ("t=2" in name): continue
        a.loglog(sizes, med, "o-", label=f"{fam.upper()}")
ax[0].set_title("Typical operating point (t=1, Δ=5e-6)"); ax[1].set_title("Worst case: d unreachable, no pruning (one FEAS call)")
for a in ax: a.set_xlabel("nodes n"); a.set_ylabel("median ms (30/20 graphs)"); a.legend(); a.grid(True, which="both", alpha=.3)
n_ = np.array([30, 150]); ax[1].loglog(n_, 0.07 * (n_ / 30) ** 2, "k--", lw=.8, label="n^2 ref")
plt.tight_layout(); plt.savefig("../figures/fig_R2_timing.png", dpi=150)
hdr = "| Setting | Family | log-log slope (n≥30) [95% bootstrap CI] | work-count slope | median ms at n = 10/20/30/50/(75/100/150) |\n|---|---|---|---|---|\n"
open(D + "R2_tables.md", "w").write(hdr + "\n".join(md) + "\n")
json.dump(out, open(D + "R2_fits.json", "w"), indent=1); print(hdr + "\n".join(md))
