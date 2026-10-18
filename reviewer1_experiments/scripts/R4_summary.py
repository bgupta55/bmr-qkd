"""R4_summary.py -- tables + trace figure for R4 from ../data/R4_summary_*.json and R4_trace_*.csv"""
import json, pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = "../data/"
runs = [("D5e-06_det", "Δ=5e-6, deterministic tie-break (arc A first)"), ("D5e-06_random", "Δ=5e-6, random tie-break"), ("D1e-06_det", "Δ=1e-6 (tight), deterministic tie-break")]
S = {k: json.load(open(f"{D}R4_summary_{k}.json")) for k, _ in runs}
md = []
m = S["D5e-06_det"]
md.append("### T-R4-A  Detection and false alarms (500 trials × 60 epochs, n = 10⁴, binomial noise, alarm at α = 10⁻⁶)\n")
md.append("| quantity | value |\n|---|---|")
md.append(f"| clean link-epoch decisions (8 links × 24 pre-attack epochs + 7 links × 36 post-attack epochs, × 500) | {m['clean_decisions']:,} |")
md.append(f"| expected false alarms at α = 10⁻⁶ | {m['expected_fa_at_alpha_1e-6']:.3f} |")
md.append(f"| **observed false alarms** (pre-attack / non-attacked links after the attack) | **{m['fa_pre_attack']} / {m['fa_post_attack_other_links']}** |")
md.append(f"| attacked link S† = 2.418: mean Ŝ ± SD over all attacked link-epochs | {m['S_att_hat_mean']:.4f} ± {m['S_att_hat_sd']:.4f} |")
md.append(f"| S_lb mean [min, max] | {m['S_lb_att_mean']:.4f} [{m['S_lb_att_min']:.4f}, {m['S_lb_att_max']:.4f}] |")
md.append(f"| **gate (S_lb ≥ 2) pass rate** | **{m['gate_pass_rate']:.3f}** (the gate does not catch this attack) |")
md.append(f"| r_e of the attacked link (mean / max) | {m['r_att_mean']:.4f} / {m['r_att_max']:.4f}  (< R = 0.10) |")
md.append(f"| alarm fires in the attack epoch (epoch 25) | {m['detect_epoch_hist'].get('25',0)}/500; undetected {m['undetected']} |")
md.append(f"| detection-epoch histogram | {m['detect_epoch_hist']} |")
for key, title in runs:
    s = S[key]
    md.append(f"\n### T-R4-{'B' if key=='D5e-06_det' else 'C' if key=='D5e-06_random' else 'D'}  Router performance, {title}\n")
    md.append("Feasibility, R_nom|feasible and R_sec (0 if infeasible or if the selected path contains an attacked link / violates Δ) — epochs 1–24 vs 26–60; attack epoch 25 shown separately. “on A” = share of feasible selections on the attacked arc 0-1-2-3-4.\n")
    md.append("| Router | pre feas | pre R_nom\\|f | pre R_sec | ep25 R_sec | post feas | post R_nom\\|f | post R_sec [95% CI over trials] | post on A | alarm-on-path (post) |\n|---|---|---|---|---|---|---|---|---|---|")
    for rt, r in s["table"].items():
        p, a_, q = r["pre_1_24"], r["attack_25"], r["post_26_60"]; ci = r["post_r_sec_ci95"]
        md.append(f"| {rt} | {p['feasible']:.3f} | {p['r_nom_given_feasible']:.4f} | {p['r_sec']:.4f} | {a_['r_sec']:.4f} | {q['feasible']:.3f} | {q['r_nom_given_feasible']:.4f} | {q['r_sec']:.4f} [{ci[0]:.4f}, {ci[1]:.4f}] | {q['on_attacked_arc']:.2f} | {q['alarm_on_path']:.3f} |")
md.append("\n### T-R4-E  Recovery after the attack (first epoch ≥ 25 with R_sec > 0), Δ = 5×10⁻⁶, deterministic tie-break vs random tie-break\n")
md.append("| Router | det: recovered | det: mean first secure epoch | random: recovered | random: mean first secure epoch |\n|---|---|---|---|---|")
for rt in S["D5e-06_det"]["recovery"]:
    a, b = S["D5e-06_det"]["recovery"][rt], S["D5e-06_random"]["recovery"][rt]
    f = lambda x: "–" if x is None else f"{x:.2f}"
    md.append(f"| {rt} | {a['recovered_trials']}/{a['of']} | {f(a['mean_first_secure_epoch'])} | {b['recovered_trials']}/{b['of']} | {f(b['mean_first_secure_epoch'])} |")
open(D + "R4_tables.md", "w").write("\n".join(md) + "\n"); print("\n".join(md))
# figure
fig, ax = plt.subplots(1, 3, figsize=(17, 4.4), sharey=True)
for k, (key, title) in enumerate(runs):
    df = pd.read_csv(f"{D}R4_trace_{key}.csv")
    g = df.pivot(index="epoch", columns="router", values="r_sec")
    for rt, st in [("BMR", "-"), ("SPF_valid_budget", "--"), ("SPF_valid", "--"), ("SPF", ":"), ("SPF_budget", ":"), ("QPATHstyle", "-."), ("QLEAPstyle", "-."), ("QCASTstyle", "-.")]:
        ax[k].plot(g.index, g[rt] + (0.002 if rt == "SPF_budget" else 0), st, label=rt, lw=2.2 if rt == "BMR" else 1.2)
    ax[k].axvline(25, color="k", lw=.8); ax[k].set_title(title, fontsize=9); ax[k].set_xlabel("epoch"); ax[k].grid(alpha=.3)
ax[0].set_ylabel("mean R_sec over 500 trials"); ax[0].legend(fontsize=7, ncol=2, loc="center left")
plt.tight_layout(); plt.savefig("../figures/fig_R4_trace.png", dpi=150)
