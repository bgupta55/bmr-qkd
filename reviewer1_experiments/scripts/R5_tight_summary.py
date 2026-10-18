"""R5_tight_summary.py -- tables for R5 under eta_union vs eta_tight (cluster = graph, t-CI over 50 graph means)."""
import pandas as pd, numpy as np, json
from scipy import stats
D="../data/"; U=pd.read_csv(D+"R5_rows.csv.gz"); T=pd.read_csv(D+"R5_rows_tight.csv.gz")
def ci(x):
    x=np.asarray(x); m=x.mean(); h=stats.t.ppf(.975,len(x)-1)*x.std(ddof=1)/np.sqrt(len(x)); return m,h
def gmean(s,col): return s.groupby("graph_id")[col].mean().values
def pct_false(s):   # selections whose path contains an alarmed link, as % of feasible selections
    f=s[s.feasible==1]; return 100*f.alarm_v.mean() if len(f) else float("nan")
md=[]; out={}
md.append("| radius | variant | feasible | R_sec [95% CI] | false-certificate % (alarmed link on certified path) | attacked link on certified path % |\n|---|---|---|---|---|---|")
for nm,d in (("union",U),("tight",T)):
    for v in ("full","no_alarm","no_freshness","hop_R0","hop_R0_no_alarm","SPF"):
        s=d[(d.table=="ablation")&np.isclose(d.Delta,5e-6)&(d.cond=="compromised")&(d.variant==v)]
        m,h=ci(gmean(s,"r_sec")); f=s[s.feasible==1]
        tc=100*f.touches_comp.mean() if len(f) else float("nan")
        md.append(f"| {nm} | {v} | {s.feasible.mean():.3f} | {m:.4f} ± {h:.4f} | {pct_false(s):.1f} | {tc:.1f} |")
        out[f"{nm}_{v}"]=dict(feasible=float(s.feasible.mean()),r_sec=m,hw=h,false_cert_pct=pct_false(s),attacked_pct=tc)
# paired difference full - no_alarm (R_sec) and rate-floor-only effect
md.append("\n| radius | contrast (compromised, Δ=5e-6) | diff in R_sec [95% CI] | diff in feasibility |\n|---|---|---|---|")
for nm,d in (("union",U),("tight",T)):
    s=d[(d.table=="ablation")&np.isclose(d.Delta,5e-6)&(d.cond=="compromised")]
    for a,b in (("full","no_alarm"),("hop_R0","hop_R0_no_alarm"),("full","SPF")):
        x=s[s.variant==a].set_index(["graph_id","draw"]); y=s[s.variant==b].set_index(["graph_id","draw"])
        dr=(x.r_sec-y.r_sec).groupby("graph_id").mean().values; m,h=ci(dr)
        md.append(f"| {nm} | {a} − {b} | {m:+.4f} [{m-h:+.4f}, {m+h:+.4f}] | {x.feasible.mean()-y.feasible.mean():+.3f} |")
        out[f"{nm}_{a}-{b}"]=dict(diff=m,lo=m-h,hi=m+h)
# BMR vs SPF_valid_budget, hop-cap H=5
md.append("\n| radius | cond | Δ | BMR − SPF_valid_budget R_sec [95% CI] (H=5) |\n|---|---|---|---|")
for nm,d in (("union",U),("tight",T)):
    for Dl in (5e-6,1e-6):
        for cond in ("clean","compromised"):
            s=d[(d.table=="hopcap")&np.isclose(d.Delta,Dl)&(d.cond==cond)&(d.hcap==5)]
            x=s[s.variant=="full"].set_index(["graph_id","draw"]).r_sec; y=s[s.variant=="SPF_valid_budget"].set_index(["graph_id","draw"]).r_sec
            m,h=ci((x-y).groupby("graph_id").mean().values)
            md.append(f"| {nm} | {cond} | {Dl:g} | {m:+.4f} [{m-h:+.4f}, {m+h:+.4f}] |")
            out[f"{nm}_hop5_{cond}_{Dl:g}"]=dict(diff=m,lo=m-h,hi=m+h)
# Delta-violation of no_defect, tight
md.append("\n| radius | Δ-violation % of no_defect selections (compromised), Δ=1e-6 / 5e-6 |\n|---|---|")
for nm,d in (("union",U),("tight",T)):
    r=[]
    for Dl in (1e-6,5e-6):
        s=d[(d.table=="ablation")&np.isclose(d.Delta,Dl)&(d.cond=="compromised")&(d.variant=="no_defect")]; f=s[s.feasible==1]
        r.append(f"{100*f.delta_v.mean():.1f}")
    md.append(f"| {nm} | {r[0]} / {r[1]} |")
open(D+"R5_tight_tables.md","w").write("\n".join(md)); json.dump(out,open(D+"R5_tight_summary.json","w"),indent=1)
print("\n".join(md))
