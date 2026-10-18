"""
Q15 — Delta sweep, regime map, robustness (P1)
=============================================
Run Delta in {1e-6 .. 5e-2} at n=500 and n=10^4.
Record feasibility and secure delivered rate vs (n, Delta, phi).
Produces: regime map showing the knee where Δ starts binding.
Also: compromised-fraction sweep (f=0..0.5).
Expected knee: phi=0.99: w ≈ 0.017 at n=500; w ≈ 4e-6 at n=1e4.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bmr.certificate import build_certificate
from bmr.route import bmr_route

SEED = 99; TRIALS = 600; R_MIN = 0.0; W = 100
N_NODES = 30; P_ER = 4/29

def run_cell(n_audit, Delta, phi, fcomp, rng):
    feasible, rates = [], []
    for _ in range(TRIALS):
        sd = int(rng.integers(1e9))
        G = nx.erdos_renyi_graph(N_NODES, P_ER, seed=sd)
        if not nx.is_connected(G): continue
        trng = np.random.default_rng(sd)
        for u,v in G.edges():
            S0 = 2.418 if (trng.random()<fcomp) else 2.687
            sh = float(trng.normal(S0, math.sqrt((4-S0**2/4)/n_audit)))
            c = build_certificate((u,v), sh, n_audit, 1e-6, 10000, phi, 0, 0)
            G[u][v]["r"] = c.r if c else 0.0
            G[u][v]["w"] = c.w if c else 2.0
            G[u][v]["T"] = 0
        res = bmr_route(G, 0, N_NODES-1, R_MIN, Delta, 0, W, m_min=3)
        feas = res is not None and res[2] <= Delta
        feasible.append(feas)
        rates.append(res[0] if feas else 0.0)
    return dict(feasible_pct=100*np.mean(feasible), mean_rate=float(np.mean(rates)))

print("="*70); print("Q15: Δ sweep regime map"); print("="*70)
rng_m = np.random.default_rng(SEED)
rows = []
for n_audit in (500, 10_000):
    for phi in (0.9, 0.99):
        # compute honest link w at this n
        c0 = build_certificate("e", 2.687, n_audit, 1e-6, 10000, phi, 0, 0)
        w0 = c0.w if c0 else 1.0
        print(f"n={n_audit} phi={phi}: honest w={w0:.4g}  (Δ knee ≈ 3*w={3*w0:.4g})")
        for Delta in [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.05]:
            r = run_cell(n_audit, Delta, phi, 0.0, np.random.default_rng(int(rng_m.integers(1e9))))
            row = dict(n_audit=n_audit, phi=phi, Delta=Delta, fcomp=0.0, **r)
            rows.append(row)
            print(f"  Delta={Delta:.1e}: feasible={r['feasible_pct']:.1f}%  rate={r['mean_rate']:.4f}")
    print()

print("Compromised fraction sweep (n=10000, Delta=1e-3 | n=500, Delta=0.5), phi=0.99")
for n_audit, Delta in ((10_000, 1e-3), (500, 0.5)):
    for fcomp in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        r = run_cell(n_audit, Delta, 0.99, fcomp, np.random.default_rng(int(rng_m.integers(1e9))))
        row = dict(n_audit=n_audit, phi=0.99, Delta=Delta, fcomp=fcomp, **r)
        rows.append(row)
        print(f"  n={n_audit} fcomp={fcomp:.2f}: feasible={r['feasible_pct']:.1f}%  rate={r['mean_rate']:.4f}")

fig, ax = plt.subplots(1, 2, figsize=(12,4))
for n_audit, phi, col, ls in [(500,0.99,'blue','-'), (10_000,0.99,'red','--'),
                              (500,0.9,'navy',':'), (10_000,0.9,'darkred',':')]:
    rr = [r for r in rows if r["n_audit"]==n_audit and r["phi"]==phi and r["fcomp"]==0.0
          and r["Delta"] in (1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,0.01,0.03,0.05)]
    ds = [r["Delta"] for r in rr]; feas = [r["feasible_pct"] for r in rr]
    rates = [r["mean_rate"] for r in rr]
    ax[0].semilogx(ds, feas, color=col, ls=ls, marker='o', label=f"n={n_audit}, phi={phi}")
    ax[1].semilogx(ds, rates, color=col, ls=ls, marker='o', label=f"n={n_audit}, phi={phi}")
for a, t, y in zip(ax, ["Feasibility %","Mean secure rate (bits/round)"],["Feasibility","Rate"]):
    a.set_xlabel("Δ (defect budget)"); a.set_ylabel(y); a.set_title(t); a.legend(); a.grid(alpha=0.3)
plt.suptitle("Q15: Δ regime map (ER n=30, phi=0.99)")
plt.tight_layout()
fp = os.path.join(os.path.dirname(__file__),"..","figures","Q15_delta_regime.png")
plt.savefig(fp, dpi=150); plt.close()
json.dump(dict(experiment="Q15", rows=rows), open(
    os.path.join(os.path.dirname(__file__),"..","data","results","q15_delta_regime.json"),"w"), indent=2)
print("saved", fp)
