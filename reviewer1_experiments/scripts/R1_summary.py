"""
R1_summary.py -- builds R1_summary.json, markdown tables (data/R1_tables.md) and
fig_R1_sweep.png from R1_cells.csv / R1_diffs.csv / R1_decomposition.csv /
R1_trials.csv.gz.  No new randomness.
"""
import json, os, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)
c = pd.read_csv(os.path.join(D, "R1_cells.csv"))
d = pd.read_csv(os.path.join(D, "R1_diffs.csv"))
dec = pd.read_csv(os.path.join(D, "R1_decomposition.csv"))
tr = pd.read_csv(os.path.join(D, "R1_trials.csv.gz"))
md = []


def fp(x):
    return "<5e-6" if x < 5.1e-6 else f"{x:.1e}" if x < 1e-3 else f"{x:.3f}"


def cell(topo, size, mode, f, Delta, router):
    q = c[(c.topology == topo) & (c["size"] == size) & (c["mode"] == mode) & (c.f == f) &
          (c.Delta == Delta) & (c.router == router)]
    return q.iloc[0]


def diff(topo, size, mode, f, Delta, comp):
    q = d[(d.topology == topo) & (d["size"] == size) & (d["mode"] == mode) & (d.f == f) &
          (d.Delta == Delta) & (d.comparator == comp)]
    return q.iloc[0]


ROUT = ["BMR", "SPF", "SPF_valid", "SPF_budget", "SPF_valid_budget", "BMR_hop", "QPATHstyle", "QLEAPstyle", "QCASTstyle"]

# ---------- T1: main table (ER/BA, n=10..50, clean & compromised, Delta=5e-6)
md.append("### T-R1-A  Main cells, sealed t=1, Delta=5e-6 (feasibility, R_nom|feasible, R_sec; 95% cluster CI)\n")
md.append("| Topology | n | Cond | Router | Feasibility | R_nom given feasible | R_sec (zeros incl.) | hops given feas |")
md.append("|---|---|---|---|---|---|---|---|")
for topo in ["er", "ba"]:
    for size in [10, 20, 30, 50]:
        for f, cond in [(0.0, "clean"), (0.1, "f=0.10")]:
            for r in ["SPF", "SPF_valid_budget", "BMR_hop", "BMR", "QPATHstyle", "QLEAPstyle", "QCASTstyle"]:
                x = cell(topo, size, "sealed", f, 5e-6, r)
                md.append(f"| {topo.upper()} | {size} | {cond} | {r} | {x.feas:.3f} [{x.feas_lo:.3f},{x.feas_hi:.3f}] | "
                          f"{x.rn_given_feas:.4f} [{x.rn_lo:.4f},{x.rn_hi:.4f}] | {x.rs:.4f} [{x.rs_lo:.4f},{x.rs_hi:.4f}] | {x.hops_given_feas:.2f} |")

# ---------- T2: differences
md.append("\n### T-R1-B  Paired differences BMR - router on R_sec (bits/round), graph-level, sealed, Delta=5e-6\n")
md.append("| Topology | n | Cond | Comparator | Diff | 95% cluster CI | p (sign-flip) | Holm p (main family) | d_z | W/T/L (of G) |")
md.append("|---|---|---|---|---|---|---|---|---|---|")
for topo in ["er", "ba"]:
    for size in [10, 20, 30, 50]:
        for f, cond in [(0.0, "clean"), (0.1, "f=0.10")]:
            for comp in ["SPF", "SPF_valid_budget", "BMR_hop", "QPATHstyle", "QLEAPstyle", "QCASTstyle"]:
                x = diff(topo, size, "sealed", f, 5e-6, comp)
                ph = x.p_holm_main
                md.append(f"| {topo.upper()} | {size} | {cond} | {comp} | {x.diff_rs:+.4f} | [{x.lo:+.4f},{x.hi:+.4f}] | "
                          f"{fp(x.p_signflip)} | {fp(ph) if not np.isnan(ph) else 'n/a'} | {x.d_z:.2f} | {int(x.wins)}/{int(x.ties)}/{int(x.losses)} |")

# ---------- T3: f-sweep ER 30
md.append("\n### T-R1-C  f-sweep, ER n=30, sealed, Delta=5e-6 (R_sec; difference to SPF and to SPF_valid_budget)\n")
md.append("| f | SPF | SPF_budget | SPF_valid_budget | BMR | BMR - SPF [CI] | BMR - SPF_valid_budget [CI] | idealised ceiling (+%) | SPF_valid_budget vs SPF_budget (+%) | P(selected path touches compromised link): SPF / SPF_budget |")
md.append("|---|---|---|---|---|---|---|---|---|---|")
fsw = {}
E_er30 = tr[(tr.topology == "er") & (tr["size"] == 30)].shape[0]  # not used
for f in [0.0, 0.05, 0.10, 0.20, 0.30]:
    r = {k: cell("er", 30, "sealed", f, 5e-6, k) for k in ["SPF", "SPF_budget", "SPF_valid_budget", "BMR"]}
    d1 = diff("er", 30, "sealed", f, 5e-6, "SPF")
    d2 = diff("er", 30, "sealed", f, 5e-6, "SPF_valid_budget")
    # ideal ceiling: P(3-hop path touches a compromised link) = 1-(1-f)^3 ; gain over a blind router
    pt = 1 - (1 - f) ** 3
    ceil = 100 * (1 / (1 - pt) - 1) if f > 0 else 0.0
    obs = 100 * (r["SPF_valid_budget"].rs / r["SPF_budget"].rs - 1)
    tch = {k: tr[(tr.topology == "er") & (tr["size"] == 30) & (tr["mode"] == "sealed") & (tr.Delta == 5e-6)
                 & (tr.f == f) & (tr.router == k)].alarm_on_path.mean() for k in ["SPF", "SPF_budget"]}
    fsw[str(f)] = dict(SPF=r["SPF"].rs, SPF_budget=r["SPF_budget"].rs, SPF_valid_budget=r["SPF_valid_budget"].rs,
                       BMR=r["BMR"].rs, diff_SPF=d1.diff_rs, ci_SPF=[d1.lo, d1.hi],
                       diff_SVB=d2.diff_rs, ci_SVB=[d2.lo, d2.hi], ceiling_pct=ceil, p_touch=pt)
    md.append(f"| {f:.2f} | {r['SPF'].rs:.4f} | {r['SPF_budget'].rs:.4f} | {r['SPF_valid_budget'].rs:.4f} | {r['BMR'].rs:.4f} | "
              f"{d1.diff_rs:+.4f} [{d1.lo:+.4f},{d1.hi:+.4f}] | {d2.diff_rs:+.4f} [{d2.lo:+.4f},{d2.hi:+.4f}] | "
              f"{ceil:.0f} | {obs:+.0f} | {tch['SPF']:.3f} / {tch['SPF_budget']:.3f} |")
md.append("\nCeiling = 1/(1-P_touch) - 1 with P_touch = 1-(1-f)^3 (random 3-hop path, clean alternative always available). "
          "`alarm_on_path` counts a selected path as touched when it contains an alarmed link (an attacked link raises the alarm with probability 1.0). "
          "The measured exposure of SPF is close to the ceiling's P_touch; the exposure of SPF_budget is higher, because the budget filter "
          "removes paths with class-C links (eps=5e-6) while attacked links carry eps=1e-6 and pass it. "
          "The observed validity gain therefore exceeds the idealised ceiling for the budget-aware comparator; the ceiling is not a bound there.\n")

# ---------- T4: decomposition
md.append("\n### T-R1-D  Decomposition of BMR - SPF (R_sec, bits/round, graph-level CI), sealed, Delta=5e-6\n")
md.append("| Topology | n | Cond | Budget filter (SPF->SPF_budget) | Validity filter (SPF_budget->SPF_valid_budget) | Objective: min-hop tie-break (->BMR_hop) | Objective: hop window (->BMR) | Total BMR - SPF |")
md.append("|---|---|---|---|---|---|---|---|")
for topo in ["er", "ba"]:
    for size in [20, 30, 50]:
        for cond in ["clean", "compromised"]:
            x = dec[(dec.topology == topo) & (dec["size"] == size) & (dec.cond == cond) & (dec.Delta == 5e-6)].iloc[0]
            md.append(f"| {topo.upper()} | {size} | {cond} | {x.order1_budget:+.4f} [{x.order1_budget_lo:+.4f},{x.order1_budget_hi:+.4f}] | "
                      f"{x.order1_validity:+.4f} [{x.order1_validity_lo:+.4f},{x.order1_validity_hi:+.4f}] | "
                      f"{x['objective_tiebreak(BMR_hop-SPFvb)']:+.4f} | {x['objective_hopwindow(BMR-BMR_hop)']:+.4f} | "
                      f"{x['total(BMR-SPF)']:+.4f} |")
md.append("\nOrder 1 = budget first, then validity. Order 2 (validity first, then budget) is in `R1_decomposition.csv`.\n")

# ---------- T5: tight Delta
md.append("\n### T-R1-E  Tight budget Delta=1e-6 (sealed): feasibility, R_sec, BMR - SPF\n")
md.append("| Topology | n | Cond | BMR feas | SPF_budget feas | BMR R_nom|feas | BMR R_sec | SPF R_sec | SPF_valid_budget R_sec | BMR - SPF [CI] |")
md.append("|---|---|---|---|---|---|---|---|---|---|")
for topo in ["er", "ba"]:
    for size in [10, 20, 30, 50]:
        for f, cond in [(0.0, "clean"), (0.1, "f=0.10")]:
            b = cell(topo, size, "sealed", f, 1e-6, "BMR"); sb = cell(topo, size, "sealed", f, 1e-6, "SPF_budget")
            s = cell(topo, size, "sealed", f, 1e-6, "SPF"); sv = cell(topo, size, "sealed", f, 1e-6, "SPF_valid_budget")
            dd = diff(topo, size, "sealed", f, 1e-6, "SPF")
            md.append(f"| {topo.upper()} | {size} | {cond} | {b.feas:.3f} | {sb.feas:.3f} | {b.rn_given_feas:.4f} | {b.rs:.4f} | {s.rs:.4f} | {sv.rs:.4f} | "
                      f"{dd.diff_rs:+.4f} [{dd.lo:+.4f},{dd.hi:+.4f}] |")

# ---------- T6: ring/line/star
md.append("\n### T-R1-F  Ring, line, star (sealed t=1, Delta=5e-6, f=0.10): ties are ties\n")
md.append("| Topology | n | BMR feas | BMR R_sec | SPF R_sec | W/T/L BMR vs SPF | note |")
md.append("|---|---|---|---|---|---|---|")
for topo in ["ring", "line", "star"]:
    for size in [10, 20, 30, 50]:
        mode = "sealed" if topo != "star" else "sealed_t0"
        b = cell(topo, size, mode, 0.1, 5e-6, "BMR"); s = cell(topo, size, mode, 0.1, 5e-6, "SPF")
        dd = diff(topo, size, mode, 0.1, 5e-6, "SPF")
        note = {"ring": "only one simple route within H=5 (the other arc has n-3 > 5 hops)",
                "line": "exactly one simple path", "star": "t=0 reported (|P|>=1)"}[topo] if size == 10 else ""
        md.append(f"| {topo} | {size} | {b.feas:.3f} | {b.rs:.4f} | {s.rs:.4f} | {int(dd.wins)}/{int(dd.ties)}/{int(dd.losses)} | {note} |")

# ---------- T7: public mode parity
md.append("\n### T-R1-G  Public-XOR mode, odd-hop fraction of selected route, ER/BA, Delta=5e-6, f=0.10\n")
md.append("| Topology | n | BMR (parity-aware) odd | BMR_np (parity-blind) odd | SPF odd | SPF_valid_budget odd | BMR R_sec | BMR_np R_sec | SPF R_sec |")
md.append("|---|---|---|---|---|---|---|---|---|")
pub = tr[(tr["mode"] == "public") & (tr.f == 0.1) & (tr.Delta == 5e-6) & (tr.feasible == 1)].copy()
pub["odd"] = pub.hops % 2 == 1
par = {}
for topo in ["er", "ba"]:
    for size in [10, 20, 30, 50]:
        row = {}
        for r in ["BMR", "BMR_np", "SPF", "SPF_valid_budget"]:
            q = pub[(pub.topology == topo) & (pub["size"] == size) & (pub.router == r)]
            row[r] = q.odd.mean() if len(q) else float("nan")
        rs = {r: cell(topo, size, "public", 0.1, 5e-6, r).rs for r in ["BMR", "BMR_np", "SPF"]}
        par[f"{topo}{size}"] = row
        md.append(f"| {topo.upper()} | {size} | {row['BMR']:.3f} | {row['BMR_np']:.3f} | {row['SPF']:.3f} | {row['SPF_valid_budget']:.3f} | "
                  f"{rs['BMR']:.4f} | {rs['BMR_np']:.4f} | {rs['SPF']:.4f} |")

# ---------- T8: cross-metric (Q-*-style columns)
md.append("\n### T-R1-H  Cross-metric view, ER n=30 and BA n=30 (sealed, Delta=5e-6, f=0.10): each router on every metric\n")
md.append("| Topology | Router | Feas | R_nom|feas | R_sec | F_P (Werner) | pair cost (sum 1/p_gen) | throughput proxy | sum eps | budget-violating selections |")
md.append("|---|---|---|---|---|---|---|---|---|---|")
for topo in ["er", "ba"]:
    for r in ROUT:
        x = cell(topo, 30, "sealed", 0.1, 5e-6, r)
        md.append(f"| {topo.upper()} | {r} | {x.feas:.3f} | {x.rn_given_feas:.4f} | {x.rs:.4f} | {x.F_P:.4f} | {x.cost_pairs:.2f} | {x.throughput:.4f} | {x.sum_eps:.2e} | {x.delta_viol:.3f} |")
md.append("\nCost and throughput columns are unconditional means over trials (infeasible = 0). F_P, cost and throughput are the Q-*-style native objectives.\n")

open(os.path.join(D, "R1_tables.md"), "w").write("\n".join(md))

# ---------- Summary json
summ = dict(
    trials=int(len(tr)), cells=int(len(c)), graphs_per_cell=50, draws_per_graph=5,
    er30_clean_D5e6=dict(BMR_feas=cell("er", 30, "sealed", 0, 5e-6, "BMR").feas,
                         BMR_rs=cell("er", 30, "sealed", 0, 5e-6, "BMR").rs,
                         SPFvb_rs=cell("er", 30, "sealed", 0, 5e-6, "SPF_valid_budget").rs),
    f_sweep_er30=fsw, parity=par,
    false_alarm=json.load(open(os.path.join(D, "R1_false_alarm_check.json"))),
    icc=json.load(open(os.path.join(D, "R1_icc_deff.json"))),
)
json.dump(summ, open(os.path.join(D, "R1_summary.json"), "w"), indent=2, default=float)

# ---------- Figure
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
cols = {"SPF": "gray", "SPF_budget": "tab:orange", "SPF_valid_budget": "tab:green", "BMR_hop": "tab:purple", "BMR": "tab:blue"}
fs = [0.0, 0.05, 0.10, 0.20, 0.30]
for r in cols:
    y = [cell("er", 30, "sealed", f, 5e-6, r).rs for f in fs]
    lo = [cell("er", 30, "sealed", f, 5e-6, r).rs_lo for f in fs]
    hi = [cell("er", 30, "sealed", f, 5e-6, r).rs_hi for f in fs]
    ax[0].plot(fs, y, "o-", color=cols[r], label=r)
    ax[0].fill_between(fs, lo, hi, color=cols[r], alpha=0.12)
ax[0].set_xlabel("compromised fraction f"); ax[0].set_ylabel("mean R_sec (bits/round)")
ax[0].set_title("(a) f-sweep, ER n=30, sealed t=1, $\\Delta$=5e-6\n(95% graph-cluster CI)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
sizes = [10, 20, 30, 50]
for r in ["SPF", "SPF_valid_budget", "BMR"]:
    ax[1].plot(sizes, [cell("er", s, "sealed", 0.1, 5e-6, r).feas for s in sizes], "o-", color=cols[r], label=f"{r} $\\Delta$=5e-6")
ax[1].plot(sizes, [cell("er", s, "sealed", 0.1, 1e-6, "BMR").feas for s in sizes], "s--", color="tab:red", label="BMR $\\Delta$=1e-6")
ax[1].set_xlabel("nodes n"); ax[1].set_ylabel("feasibility rate")
ax[1].set_title("(b) Feasibility, ER, f=0.10"); ax[1].set_ylim(0, 1.05); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
x = np.arange(4); w = 0.2
xx = dec[(dec.topology == "er") & (dec["size"] == 30) & (dec.Delta == 5e-6)]
lab = ["budget\nfilter", "validity\nfilter", "min-hop\ntie-break", "hop\nwindow"]
for k, (cond, colr) in enumerate([("clean", "tab:cyan"), ("compromised", "tab:red")]):
    q = xx[xx.cond == cond].iloc[0]
    v = [q.order1_budget, q.order1_validity, q["objective_tiebreak(BMR_hop-SPFvb)"], q["objective_hopwindow(BMR-BMR_hop)"]]
    ax[2].bar(x + (k - 0.5) * 0.35, v, 0.35, color=colr, label=cond)
ax[2].set_xticks(x); ax[2].set_xticklabels(lab, fontsize=8)
ax[2].set_ylabel("contribution to R_sec (bits/round)")
ax[2].set_title("(c) Decomposition of BMR - SPF\nER n=30, $\\Delta$=5e-6"); ax[2].legend(); ax[2].grid(alpha=.3, axis="y")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_R1_sweep.png"), dpi=150); plt.close()
print("written R1_tables.md, R1_summary.json, fig_R1_sweep.png")
