"""Q16: generate Q6 sweep and dynamic figures (was missing from q6_sweep.py)."""
import sys, os, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, numpy as np
print("code sha256:", hashlib.sha256(open(__file__, "rb").read()).hexdigest()[:16])
sw = json.load(open(os.path.join(ROOT, "data", "results", "q6_sweep.json")))
dy = json.load(open(os.path.join(ROOT, "data", "results", "q6_dynamic.json")))
ROUTERS = ["B0-SPF","B1a-Plambda(0.01)","B1b-Plambda(1)","B2-BMR-NoFloor","B3-Zhang-style",
           "B5-Rate-greedy","B6-Qpath-style","B7-Qcast-style","B8-QBER-monitor","BMR"]
COLORS = plt.cm.tab10(np.linspace(0,1,len(ROUTERS)))
R2C = {r: c for r, c in zip(ROUTERS, COLORS)}
DISP = {"B0-SPF":"SPF","B1a-Plambda(0.01)":"P_λ(0.01)","B1b-Plambda(1)":"P_λ(1)",
        "B2-BMR-NoFloor":"Pareto (BMR-NoFloor)","B3-Zhang-style":"Zhang-widest",
        "B5-Rate-greedy":"Rate-greedy","B6-Qpath-style":"Q-PATH-style",
        "B7-Qcast-style":"Q-CAST-style","B8-QBER-monitor":"QBER","BMR":"BMR"}
from collections import defaultdict
ragg = defaultdict(lambda: {"rate_f": [], "compl_f": [], "rate_bo": [], "compl_bo": [], "c_rate_f": [], "c_compl_f": []})
for row in sw["rows"]:
    k = row["router"]; r = ragg[k]
    r["rate_f"].append(row["mean_rate"]); r["compl_f"].append(row["compl_pct"])
    r["rate_bo"].append(row["mean_rate_budget_only"]); r["compl_bo"].append(row["compl_budget_only_pct"])
    if row["cond"] == "compromised":
        r["c_rate_f"].append(row["mean_rate"]); r["c_compl_f"].append(row["compl_pct"])
rns = [r for r in ROUTERS if r in ragg]
x = np.arange(len(rns)); W = 0.35
fig, ax = plt.subplots(1, 2, figsize=(14, 5))
ax[0].bar(x - W/2, [sum(ragg[r]["rate_f"])/len(ragg[r]["rate_f"]) for r in rns], W, label="Full (floor+budget)", color=[R2C[r] for r in rns], alpha=.85)
ax[0].bar(x + W/2, [sum(ragg[r]["rate_bo"])/len(ragg[r]["rate_bo"]) for r in rns], W, label="Budget-only", color=[R2C[r] for r in rns], alpha=.45, edgecolor="k", linewidth=.5)
ax[0].set_xticks(x); ax[0].set_xticklabels([DISP.get(r,r) for r in rns], rotation=30, ha="right", fontsize=8)
ax[0].set_ylabel("Mean secure rate"); ax[0].set_title("Q6 Sweep: mean rate (200 cells × 1000 trials)\nfull metric vs budget-only metric"); ax[0].legend(fontsize=9)
ax[1].bar(x - W/2, [sum(ragg[r]["compl_f"])/len(ragg[r]["compl_f"]) for r in rns], W, label="Full compliance %", color=[R2C[r] for r in rns], alpha=.85)
ax[1].bar(x + W/2, [sum(ragg[r]["c_compl_f"])/max(1,len(ragg[r]["c_compl_f"])) for r in rns], W, label="Full compliance % (compromised links)", color=[R2C[r] for r in rns], alpha=.45, edgecolor="r", linewidth=.5)
ax[1].set_xticks(x); ax[1].set_xticklabels([DISP.get(r,r) for r in rns], rotation=30, ha="right", fontsize=8)
ax[1].set_ylabel("Compliance %"); ax[1].set_title("Q6 Sweep: budget/floor compliance\nfull metric, all cells (solid) vs 20% compromised (edge)"); ax[1].legend(fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(ROOT,"figures","Q6_sweep.png"), dpi=140)
print("Q6_sweep.png saved")
# dynamic figure
from collections import defaultdict as D2
dagg = D2(list)
for row in dy["rows"]:
    key = (row["series"], row["n_audit"], row["S1"])
    dagg[key].append(row)
fig2, axes = plt.subplots(2, 3, figsize=(16, 9))
topos = ["ring","grid","ER"]
axes = axes.flatten()
for i, (topo, na, S1) in enumerate([(t, 10_000, 2.418) for t in topos] + [("grid", 500, 2.418), ("grid",10_000,2.62), ("ER",10_000,2.55)]):
    ax = axes[i]
    rows = [r for r in dy["rows"] if r.get("topo",r.get("topo_name","")) == topo and r["n_audit"] == na and abs(r["S1"]-S1)<0.001]
    if not rows: ax.text(0.5,0.5,"no data",ha="center"); ax.set_title(f"{topo} n_audit={na} S1={S1}"); continue
    series_d = {r["series"]: r for r in rows}
    for sname, row in series_d.items():
        ls = "-" if "BMR" in sname else ("--" if "B3" in sname else ":")
        ax.plot([1,25,60,80], [row.get("rate_pre"), row.get("rate_pre"), row.get("rate_post1"), row.get("rate_post2")], ls, label=sname, lw=1.5)
    ax.axvline(25,color="k",ls=":",lw=.8); ax.axvline(60,color="k",ls=":",lw=.8)
    ax.set_title(f"{topo} n_audit={na} S1={S1}"); ax.set_xlabel("epoch"); ax.set_ylabel("rate"); ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(ROOT,"figures","Q6_dynamic.png"), dpi=140)
print("Q6_dynamic.png saved")
