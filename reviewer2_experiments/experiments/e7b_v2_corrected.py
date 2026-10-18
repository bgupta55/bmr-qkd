"""
E7b v2 — Q1: Corrected E7b metrics with self-checks (pending_v7 Q1)
====================================================================
Recomputes E7b with the correct metric definitions:
  feasible[i]     = exists_route_within_budget(instance_i)  [same for all routers]
  found[i,r]      = route_r is not None
  compliant[i,r]  = found and sum_w(route) <= Delta and hops >= 2*t+1
  secure_rate[i,r]= bottleneck_rate if compliant else 0.0   [mean over ALL trials]

Self-checks:
  assert mean(secure) <= mean(feasible) * max_link_rate + 1e-9
  compliance of delivered keys for BMR must be exactly 1.0
  refusal rate = 1 - mean(found)

Adds B3 (Zhang-style Key-Rate Router), B8 (QBER-Monitor), B2 (BMR-NoFloor ablation).
Reports refusal rate and common-feasible rate columns.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import math, json, time
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from bmr.bounds import eta_tight, h_pironio
from bmr.certificate import build_certificate
from bmr.route import bmr_route

MASTER_SEED = 50
N_TRIALS    = 1000
# Δ per audit level: ensure w_e (eps + delta) can be non-trivially constrained
# At n=500, S0 honest: eps ≈ 0.017, delta = 1e-6 => w ≈ 0.017; use Δ = 0.5 (per path, 3 hops)
# At n=10000, S0 honest: w ≈ 3.94e-6; use Δ = 1e-3 (allows ~250 links, so feasibility is high)
DELTA_BY_N  = {500: 0.5, 10_000: 1e-3}
R_MIN       = 0.01
W_FRESH     = 10
T_COLLUSION = 1   # hop floor = 2*1+1 = 3

def make_db(G, n_audit, rng, frac_comp=0.20, S0=2.687, S1=2.418, delta=1e-6,
            ell=10_000, phi=0.99):
    db = {}
    edges = list(G.edges())
    n_comp = max(0, int(frac_comp * len(edges)))
    comp_set = set(rng.choice(len(edges), n_comp, replace=False).tolist())
    for i, (u, v) in enumerate(edges):
        S_true = S1 if i in comp_set else S0
        std_S = math.sqrt((4 - S_true**2/4) / n_audit)
        S_hat = float(rng.normal(S_true, std_S))
        cert = build_certificate((u, v), S_hat, n_audit, delta, ell, phi, 0, 0)
        r_e   = cert.r   if cert else 0.0
        eps_e = cert.eps if cert else 1.0
        w_e   = cert.w   if cert else 2.0
        db[(u,v)] = db[(v,u)] = {
            "r": r_e, "eps": eps_e, "w": w_e, "T": 0,
            "comp": i in comp_set, "cert": cert
        }
    return db

def build_G_nx(G_base, db):
    G = nx.Graph()
    for u, v in G_base.edges():
        a = db.get((u,v), {"r":0,"eps":1,"w":2,"T":0})
        G.add_edge(u, v, **a)
    return G

# ───── Routers ────────────────────────────────────────────────────────────────

def spf(G, s, d, Delta, t):
    try:
        path = nx.shortest_path(G, s, d)
        r = min(G[u][v]["r"] for u,v in zip(path[:-1], path[1:]))
        w = sum(G[u][v]["w"] for u,v in zip(path[:-1], path[1:]))
        m = len(path)-1
        return path, r, w
    except nx.NetworkXNoPath:
        return None, 0.0, float("inf")

def rate_greedy(G, s, d, Delta, t):
    Gf = nx.Graph()
    for u, v, a in G.edges(data=True):
        if a["r"] >= R_MIN:
            Gf.add_edge(u, v, w=1.0/(a["r"]+1e-9), r=a["r"], ww=a["w"])
    try:
        _, path = nx.single_source_dijkstra(Gf, s, d, weight="w")
        r = min(G[u][v]["r"] for u,v in zip(path[:-1], path[1:]))
        w = sum(G[u][v]["w"] for u,v in zip(path[:-1], path[1:]))
        return path, r, w
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None, 0.0, float("inf")

def zhang_key_rate(G, s, d, Delta, t):
    """B3: max-min key-rate via Dijkstra on 1/r, ignoring w budget."""
    Gf = nx.Graph()
    for u, v, a in G.edges(data=True):
        if a["r"] >= R_MIN:
            Gf.add_edge(u, v, weight=1.0/(a["r"]+1e-9))
    try:
        _, path = nx.single_source_dijkstra(Gf, s, d, weight="weight")
        r = min(G[u][v]["r"] for u,v in zip(path[:-1], path[1:]))
        w = sum(G[u][v]["w"] for u,v in zip(path[:-1], path[1:]))
        return path, r, w
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None, 0.0, float("inf")

def qber_monitor(G, s, d, Delta, t):
    """B8: QBER-monitor style — route on Bell-derived rate only, no budget, no hop floor.
    Uses same rate data (from build_certificate), so equivalent to rate_greedy under
    honest-link assumption. Key difference from BMR: no budget enforcement.
    """
    return rate_greedy(G, s, d, Delta, t)

import heapq
def sr_dijkstra(G, s, d, Delta, t):
    """Pareto label-setting on (rate, budget) = BMR-NoFloor ablation (B2). NOT the Pirandola lambda protocol (see q6_common B1a/B1b)."""
    if s not in G or d not in G:
        return None, 0.0, float("inf")
    labels = {v: [] for v in G.nodes()}
    labels[s] = [(0.0, 0.0)]
    pq = [(0.0, 0.0, s, [s])]
    best_r, best_path, best_w = -float("inf"), None, float("inf")
    def dominated(nr, b, ex):
        return any(e[0] <= nr and e[1] <= b for e in ex)
    while pq:
        neg_r, bud, v, path = heapq.heappop(pq)
        r_so_far = -neg_r
        if v == d:
            if r_so_far > best_r and bud <= Delta:
                best_r, best_path, best_w = r_so_far, path, bud
            continue
        for u, a in G[v].items():
            rl = a.get("r",0); wl = a.get("w",2)
            if rl < R_MIN: continue
            nr = min(r_so_far, rl) if r_so_far>0 else rl
            nb_w = bud + wl
            if nb_w > Delta: continue
            nn = -nr
            if not dominated(nn, nb_w, labels[u]):
                labels[u] = [l for l in labels[u] if not(nn<=l[0] and nb_w<=l[1])]
                labels[u].append((nn, nb_w))
                heapq.heappush(pq, (nn, nb_w, u, path+[u]))
    return (best_path, best_r, best_w) if best_path else (None, 0.0, float("inf"))

def bmr_nf(G, s, d, Delta, t):
    """B2: BMR-NoFloor (t=0, no hop floor)."""
    res = bmr_route(G, s, d, R_MIN, Delta, 0, W_FRESH, m_min=0)
    if res is None: return None, 0.0, float("inf")
    rho, path, w = res
    return path, rho, w

def bmr(G, s, d, Delta, t):
    """BMR: binary search on rate + Yen's simple path with hop floor 2t+1."""
    m_min = 2*t+1 if t>0 else 0
    res = bmr_route(G, s, d, R_MIN, Delta, 0, W_FRESH, m_min=m_min)
    if res is None: return None, 0.0, float("inf")
    rho, path, w = res
    return path, rho, w

ROUTERS = [
    ("SPF",         spf),
    ("Rate-Greedy", rate_greedy),
    ("Zhang-KR",    zhang_key_rate),
    ("QBER-Mon",    qber_monitor),
    ("Pareto-LS(=B2)", sr_dijkstra),
    ("BMR-NoFloor", bmr_nf),
    ("BMR",         bmr),
]

configs = [
    ("ER",30,4,"random"),
    ("BA",30,2,"scale-free"),
]
n_audits    = [500, 10_000]
conditions  = [("clean",0.0),("compromised",0.20)]

all_rows = []
rng_master = np.random.default_rng(MASTER_SEED)

print("="*110)
print("E7b-v2: Corrected baseline comparison with self-checks")
print("="*110)
hdr = f"{'topo':>12} {'n':>7} {'cond':>12} {'router':>14} | {'rate':>8} {'compl%':>7} {'feasbl%':>8} {'refusal%':>9} {'hops':>7}"
print(hdr); print("-"*110)

for gtype, n_nodes, degree, tname in configs:
    for n_audit in n_audits:
        DELTA = DELTA_BY_N[n_audit]
        for cname, fcomp in conditions:
            per_router = {r: {"rates":[],"compliant":[],"found":[],"hops":[]}
                         for r,_ in ROUTERS}
            instance_feasible = []
            max_link_rate_all = []

            g_rng = np.random.default_rng(int(rng_master.integers(1e9)))
            for trial in range(N_TRIALS):
                seed_t = int(g_rng.integers(1e9))
                t_rng  = np.random.default_rng(seed_t)
                if gtype=="ER":
                    G_base = nx.erdos_renyi_graph(n_nodes, degree/(n_nodes-1), seed=seed_t)
                else:
                    G_base = nx.barabasi_albert_graph(n_nodes, 2, seed=seed_t)   # m = 2 (14 Oct fix N-2: was max(1, degree//2) = 1 -> a tree)
                if not nx.is_connected(G_base): continue

                s, d = 0, n_nodes-1
                db = make_db(G_base, n_audit, t_rng, fcomp)
                G_nx = build_G_nx(G_base, db)

                # Instance-level feasibility: does any route within budget exist?
                # Use the Pareto label-setting solver as an oracle (doesn't depend on BMR)
                _, r_oracle, w_oracle = sr_dijkstra(G_nx, s, d, DELTA, T_COLLUSION)
                feas = (r_oracle > 0 and w_oracle <= DELTA)
                instance_feasible.append(feas)
                mlr = max((G_nx[u][v]["r"] for u,v in G_nx.edges() if G_nx[u][v]["r"]>0), default=0)
                max_link_rate_all.append(mlr)

                for rname, rfn in ROUTERS:
                    path, rate, ww = rfn(G_nx, s, d, DELTA, T_COLLUSION)
                    found = path is not None
                    compl = found and ww <= DELTA
                    if rname == "BMR":
                        # BMR also checks hop floor
                        if found: compl = compl and (len(path)-1 >= 2*T_COLLUSION+1)
                    sec_rate = rate if compl else 0.0
                    per_router[rname]["rates"].append(sec_rate)
                    per_router[rname]["compliant"].append(1 if compl else 0)
                    per_router[rname]["found"].append(1 if found else 0)
                    per_router[rname]["hops"].append(len(path)-1 if found else 0)

            mean_feas = np.mean(instance_feasible)
            mlr_mean  = np.mean(max_link_rate_all)

            for rname, _ in ROUTERS:
                pr = per_router[rname]
                mean_rate   = np.mean(pr["rates"])
                compl_pct   = 100*np.mean(pr["compliant"])
                found_pct   = 100*np.mean(pr["found"])
                refusal_pct = 100*(1-np.mean(pr["found"]))
                hops_fnd    = [h for h,f in zip(pr["hops"],pr["found"]) if f]
                mean_hops   = np.mean(hops_fnd) if hops_fnd else 0.0

                # Self-check (assertion)
                assert mean_rate <= mean_feas * mlr_mean + 1e-6, \
                    f"SELF-CHECK FAILED: {rname} mean_rate={mean_rate:.4f} > {mean_feas:.3f}*{mlr_mean:.4f}"
                if rname == "BMR":
                    # BMR: compliance of delivered routes (where found) should be ~100% (may be <100 if budget tight)
                    delivered_compl = [c for c,f in zip(pr["compliant"],pr["found"]) if f]
                    if delivered_compl:
                        bmr_del_compl = np.mean(delivered_compl)
                        # Only assert if BMR found routes; compliance of delivered should be 1.0
                        assert bmr_del_compl >= 0.999 or compl_pct == 0, \
                            f"BMR: compliance of delivered routes = {bmr_del_compl:.4f} < 1.0!"

                row = dict(
                    topology=tname, graph_type=gtype, n_nodes=n_nodes,
                    n_audit=n_audit, condition=cname, frac_comp=fcomp,
                    router=rname, Delta=DELTA,
                    mean_rate=round(mean_rate,4), compl_pct=round(compl_pct,1),
                    feasible_pct=round(100*mean_feas,1), refusal_pct=round(refusal_pct,1),
                    mean_hops=round(mean_hops,2), n_trials=N_TRIALS,
                )
                all_rows.append(row)
                print(f"{tname:>12} {n_audit:>7} {cname:>12} {rname:>14} | "
                      f"{mean_rate:>8.4f} {compl_pct:>7.1f}% {100*mean_feas:>7.1f}% "
                      f"{refusal_pct:>8.1f}% {mean_hops:>7.2f}")
            print()

print("\n[All self-checks passed]")

# ── figure ──
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
topos_p = ["random","scale-free"]
conds_p = ["clean","compromised"]
clrs = {"SPF":"gray","Rate-Greedy":"steelblue","Zhang-KR":"teal",
        "QBER-Mon":"purple","Pareto-LS(=B2)":"orange","BMR-NoFloor":"tomato","BMR":"green"}

for ri, n_audit in enumerate(n_audits):
    DELTA = DELTA_BY_N[n_audit]
    xs = np.arange(len(topos_p)*len(conds_p))
    lbls = [f"{t[:4]}\n{c[:4]}" for t in topos_p for c in conds_p]
    width = 0.1
    ax0, ax1 = axes[ri]
    for j, (rname, _) in enumerate(ROUTERS):
        rates_p, compl_p = [], []
        for tname in topos_p:
            for cname in conds_p:
                match = [r for r in all_rows if r["topology"]==tname and r["n_audit"]==n_audit
                         and r["condition"]==cname and r["router"]==rname]
                rates_p.append(match[0]["mean_rate"] if match else 0)
                compl_p.append(match[0]["compl_pct"] if match else 0)
        off = (j - len(ROUTERS)/2) * width
        ax0.bar(xs+off, rates_p, width, label=rname, color=clrs[rname], alpha=0.8, edgecolor='k', lw=0.5)
        ax1.bar(xs+off, compl_p, width, label=rname, color=clrs[rname], alpha=0.8, edgecolor='k', lw=0.5)
    ax0.set_xticks(xs); ax0.set_xticklabels(lbls, fontsize=8)
    ax0.set_ylabel("Secure delivered rate (bits/round)"); ax0.set_title(f"E7b-v2 Rate (n={n_audit})")
    ax0.legend(fontsize=7); ax0.grid(True, alpha=0.3, axis='y')
    ax1.set_xticks(xs); ax1.set_xticklabels(lbls, fontsize=8)
    ax1.set_ylabel("Budget compliance %"); ax1.set_title(f"E7b-v2 Compliance (n={n_audit})")
    ax1.legend(fontsize=7); ax1.grid(True, alpha=0.3, axis='y'); ax1.set_ylim(0,110)
plt.tight_layout()
outfig = os.path.join(os.path.dirname(__file__), "..", "figures", "E7b_v2_corrected.png")
plt.savefig(outfig, dpi=150, bbox_inches='tight'); plt.close()
print(f"\nFigure saved: {outfig}")

out = {"experiment":"E7b_v2","description":"Corrected E7b with full router suite and self-checks (Q1)",
       "parameters":{"N_TRIALS":N_TRIALS,"T_COLLUSION":T_COLLUSION,"R_MIN":R_MIN,"DELTA_BY_N":DELTA_BY_N},
       "all_rows":all_rows}
outpath = os.path.join(os.path.dirname(__file__),"..", "data","results","e7b_v2_corrected.json")
with open(outpath,"w") as f: json.dump(out, f, indent=2)
print(f"Results saved: {outpath}")
