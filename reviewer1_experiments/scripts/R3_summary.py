"""R3_summary.py -- tables + figure for R3 from ../data/R3_*.json"""
import json, glob, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = "../data/"; md = []
SOL = ["enum", "label", "dp", "ilp", "larac", "fptas"]
def med(x): return float(np.median(x)) if len(x) else float("nan")
# ---------------- family (a)
md.append("### T-R3-A  Family (a): ER(p=0.25), integer log-uniform (a,b), H=6, L=1, budgets at the 30th Pareto percentile\n")
md.append("| n | inst | exact-feasible | joint-binding | front size (med) | enum ms (med) [paths enumerated, med] | label ms [labels created] | DP ms | ILP ms | LARAC-style ms / false-neg | FPTAS-style ms / false-pos | enum time-outs |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
fam = {}
for n in (10, 20, 30, 50, 75, 100, 150):
    R = json.load(open(f"{D}R3_fam_n{n}.json"))["rows"]
    feas = [r for r in R if r["exact_feasible"]]
    jb = sum(r["joint_binding"] for r in R)
    def col(s, f="ms", rows=R): return [r["solvers"][s][f] for r in rows if s in r["solvers"] and f in r["solvers"][s] and r["solvers"][s].get("status") == "ok"]
    en = col("enum"); eto = sum(1 for r in R if r["solvers"]["enum"]["status"] != "ok")
    paths = col("enum", "paths"); created = col("label", "created")
    lar_fn = sum(1 for r in feas if not r["solvers"]["larac"]["feasible"])
    fp_fp = sum(1 for r in R if (not r["exact_feasible"]) and r["solvers"]["fptas"]["feasible"])
    fam[n] = dict(n=n, inst=len(R), feas=len(feas), jb=jb, enum=med(en), label=med(col("label")), dp=med(col("dp")), ilp=med(col("ilp")),
                  larac=med(col("larac")), fptas=med(col("fptas")), paths=med(paths), created=med(created), eto=eto, lar_fn=lar_fn, fp_fp=fp_fp)
    f = fam[n]
    md.append(f"| {n} | {len(R)} | {len(feas)}/{len(R)} | {jb}/{len(R)} | {med([r['front_size'] for r in R]):.0f} | {f['enum']:.3g} [{f['paths']:.3g}] | "
              f"{f['label']:.3g} [{f['created']:.3g}] | {f['dp']:.3g} | {f['ilp']:.3g} | {f['larac']:.3g} / {lar_fn}/{len(feas)} | {f['fptas']:.3g} / {fp_fp}/{len(R)-len(feas)} | {eto}/{len(R)} |")
# ---------------- gadget
for tag, fn in (("fixed", "R3_gadget.json"), ("exp", "R3_gadget_exp.json")):
    G = json.load(open(D + fn))
    md.append(f"\n### T-R3-B{'1' if tag=='fixed' else '2'}  Partition gadget, m=5..30 stages ({'x_i <= 10^4 (fixed range)' if tag=='fixed' else 'x_i <= 2^(m+4) (weights grow with m)'}); median over instances (half planted-feasible, half random); 2^m s-d paths\n")
    md.append("| m | 2^m | enum ms [paths enumerated] | label ms [labels created / kept] | DP ms [cells] | ILP ms | LARAC-style ms (false-neg) | FPTAS-style ms (false-pos) | exact feasible |\n|---|---|---|---|---|---|---|---|---|")
    for g in G:
        R = g["rows"]; m = g["m"]
        def c(s, f):
            v = [r["solvers"][s][f] for r in R if s in r["solvers"] and r["solvers"][s].get("status") == "ok" and f in r["solvers"][s]]
            return med(v)
        def nsk(s): return sum(1 for r in R if s not in r["solvers"] or r["solvers"][s].get("status") != "ok")
        def exact(r):
            for k in ("label", "enum", "ilp", "dp"):
                v = r["solvers"].get(k)
                if v and v.get("status") == "ok" and v.get("feasible") is not None: return bool(v["feasible"])
            return False
        ex = [exact(r) for r in R]
        lfn = sum(1 for r, e in zip(R, ex) if e and not r["solvers"]["larac"]["feasible"])
        ffp = sum(1 for r, e in zip(R, ex) if (not e) and r["solvers"]["fptas"]["feasible"])
        nl = nsk("label")
        dp = f"{c('dp','ms'):.3g} [{c('dp','cells'):.2g}]" + (f" ({nsk('dp')} n/a)" if nsk("dp") else "")
        md.append(f"| {m} | {2**m:.3g} | {c('enum','ms'):.3g} [{c('enum','paths'):.3g}] | {c('label','ms'):.3g} [{c('label','created'):.3g} / {c('label','kept'):.3g}]{(' ('+str(nl)+' capped/t.o.)') if nl else ''} | {dp} | {c('ilp','ms'):.3g} | "
                  f"{c('larac','ms'):.2g} ({lfn}/{sum(ex)}) | {c('fptas','ms'):.2g} ({ffp}/{len(ex)-sum(ex)}) | {sum(ex)}/{len(ex)} |")
# ---------------- naive
N = json.load(open(D + "R3_naive.json"))
md.append("\n### T-R3-C  Naive scalarisation (min a + λ·b, no budget awareness) vs budget-pruned scalarisation vs LARAC-style vs FPTAS-style; family (a), budgets at 30th Pareto percentile, **only exact-feasible instances**\n")
md.append("| n | instances (feasible) | joint-binding | λ | returned a path | **returned-violating** | (failed-to-find-feasible of the 'pruned' variant) | LARAC-style failed | FPTAS-style failed |\n|---|---|---|---|---|---|---|---|---|")
naive = {}
for n in sorted(set(r["n"] for r in N)):
    R = [r for r in N if r["n"] == n and r["exact_feasible"]]
    for lam in ("0.1", "1", "10", "100"):
        ret = sum(r["lam" + lam]["returned"] for r in R); vio = sum(r["lam" + lam]["violating"] for r in R)
        pf = sum(1 for r in R if not r["lam" + lam]["pruned_found"])
        lf = sum(1 for r in R if not r["larac_feasible"]); ff = sum(1 for r in R if not r["fptas_feasible"])
        naive[(n, lam)] = (vio, len(R))
        md.append(f"| {n} | {len(R)} ({sum(1 for r in N if r['n']==n)} gen.) | {sum(r['joint_binding'] for r in R)}/{len(R)} | {lam} | {ret}/{len(R)} | **{vio}/{len(R)} = {100*vio/len(R):.1f}%** | {pf}/{len(R)} | {lf}/{len(R)} | {ff}/{len(R)} |")
# ---------------- dp3
P = json.load(open(D + "R3_dp3.json"))
md.append("\n### T-R3-D  3-class count-vector DP (Remark 2): reachable (vertex, k1,k2,k3) states; ER p=0.25, 10 graphs per size\n")
md.append("| H | n | state space n·C(H+3,3) | median ms [IQR] | agrees with label-setting (decision at 30th-pct budgets) |\n|---|---|---|---|---|")
for r in P:
    if "median_ms" not in r: md.append(f"| {r['H_mode']} | {r['n']} | – | skipped: {r['status']} | – |"); continue
    md.append(f"| {r['H']} | {r['n']} | {r['state_space']:.3g} | {r['median_ms']:.3g} [{r['q25']:.3g}, {r['q75']:.3g}] | {r['agree']}/{r['agree_checked']} |")
open(D + "R3_tables.md", "w").write("\n".join(md) + "\n"); print("\n".join(md))
# ---------------- figure
fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
ns = sorted(fam)
for s, lab in [("enum", "enumeration"), ("label", "label-setting"), ("dp", "pseudo-poly DP"), ("ilp", "ILP (HiGHS)"), ("larac", "LARAC-style"), ("fptas", "FPTAS-style")]:
    ax[0].loglog(ns, [max(fam[n][s], 1e-3) for n in ns], "o-", label=lab)
ax[0].set_title("Family (a): median runtime"); ax[0].set_xlabel("nodes n"); ax[0].set_ylabel("ms"); ax[0].legend(fontsize=7); ax[0].grid(True, which="both", alpha=.3)
for fn, ls, tg in (("R3_gadget.json", "-", "x≤1e4"), ("R3_gadget_exp.json", "--", "x≤2^(m+4)")):
    G = json.load(open(D + fn)); ms = [g["m"] for g in G]
    for s, c, lab in [("enum", "C0", "enumeration"), ("label", "C1", "label-setting"), ("dp", "C2", "pseudo-poly DP"), ("ilp", "C3", "ILP")]:
        y = []
        for g in G:
            v = [r["solvers"][s]["ms"] for r in g["rows"] if s in r["solvers"] and r["solvers"][s].get("status") == "ok"]
            y.append(np.median(v) if v else np.nan)
        ax[1].semilogy(ms, y, ls, color=c, marker="o", ms=3, label=f"{lab} ({tg})" if True else None)
ax[1].set_title("Partition gadget: runtime vs m"); ax[1].set_xlabel("stages m"); ax[1].set_ylabel("ms"); ax[1].legend(fontsize=6); ax[1].grid(True, which="both", alpha=.3)
G = json.load(open(D + "R3_gadget.json")); ms = [g["m"] for g in G]
ax[2].semilogy(ms, [np.median([r["solvers"]["enum"]["paths"] for r in g["rows"] if "enum" in r["solvers"] and r["solvers"]["enum"]["status"] == "ok"] or [np.nan]) for g in G], "o-", label="paths enumerated (=2^m)")
def lc(G_):
    o = []
    for g in G_:
        v = [r["solvers"]["label"]["created"] for r in g["rows"] if "label" in r["solvers"] and r["solvers"]["label"]["status"] == "ok"]
        o.append(np.median(v) if v else np.nan)
    return o
ax[2].semilogy(ms, lc(G), "s-", label="labels created (x≤1e4)")
G2 = json.load(open(D + "R3_gadget_exp.json"))
ax[2].semilogy([g["m"] for g in G2], lc(G2), "^--", label="labels created (x≤2^(m+4))")
ax[2].set_title("Work counters on the gadget"); ax[2].set_xlabel("stages m"); ax[2].legend(fontsize=7); ax[2].grid(True, which="both", alpha=.3)
plt.tight_layout(); plt.savefig("../figures/fig_R3_complexity.png", dpi=150)
