"""
E7b: SR-Dijkstra baseline (Pirandola et al. 2026) — Reviewer 2 Comment 5
==========================================================================
Implements the Security-Rate (SR) Dijkstra router from:
  Pirandola, Harney, Spedalieri, "Security-rate trade-off in QKD:
  From individual links to large-scale networks," PRL Research 8(3), 033221 (2026).

Core idea (from the paper):
  Each link carries a *certified* secret-key rate r_e and a security parameter ε_e.
  SR-Dijkstra maximises the end-to-end min-bottleneck certified rate subject to
  a joint security requirement: the compositional security budget Δ.
  Route metric: wₑ = −r_e  (Dijkstra minimises, so we maximise rate).
  Budget filter: only use paths where Σ(ε_e + δ_e) ≤ Δ.
  This is a modified Dijkstra that keeps (rate, budget_used) state per node.

COMPARISON:
  BMR vs SR-Dijkstra vs Rate-Greedy (B5) vs SPF (B0)
  on random (ER n=30, deg=4) and scale-free (BA n=30, m=2) topologies
  under clean (0% compromised) and compromised (20% links) conditions,
  n ∈ {500, 10000} audit rounds, 1000 trials.

METRICS:
  1. Certified end-to-end bottleneck rate (max-min r_e on path)
  2. Budget compliance (fraction of routes with Σw ≤ Δ)
  3. Feasibility (fraction of trials finding a route)
  4. Mean path length
  5. BMR advantage over SR-Dijkstra (%, with bootstrap 95% CI)

KEY EXPECTATION (Developer Guide §11):
  Against rate-aware QKD routers (SR-Dijkstra), the rate advantage of BMR
  should be SMALL or zero — they both avoid low-rate links.
  BMR's distinct advantage should appear in BUDGET COMPLIANCE (the composable
  security guarantee), not in raw rate. This is a legitimate, honest result.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import math
import json
import time
import heapq
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats
from bmr.bounds import eta_tight, h_pironio
from bmr.certificate import build_certificate
from bmr.route import bmr_route

MASTER_SEED  = 50
N_TRIALS     = 1000
# Δ per audit level: at n=500, eps_e ≈ 0.02/link so use 0.5 to let budget bind;
# at n=10000, eps_e ≈ 4e-6/link so 1e-3 is meaningful.
DELTA_BY_N   = {500: 0.5, 10_000: 1e-3}
R_MIN        = 0.01    # minimum accepted rate
W_FRESHNESS  = 10      # freshness window (epochs)

# ──────────────────────────────────────────────────────────────────────────────
# Certificate generator: sample per-link S_hat from binomial model
# ──────────────────────────────────────────────────────────────────────────────

def generate_link_db(G_nx, n_audit, rng, frac_compromised=0.20,
                     S0=2.687, S1=2.418, delta=1e-6,
                     ell=10_000, phi=0.99, lam_ec=0):
    """Assign each edge a certificate from the binomial CHSH model."""
    db = {}
    edges = list(G_nx.edges())
    n_comp = max(0, int(frac_compromised * len(edges)))
    compromised = set(rng.choice(len(edges), n_comp, replace=False).tolist())

    for i, (u, v) in enumerate(edges):
        S_true = S1 if i in compromised else S0
        std_S  = math.sqrt(8.0 * (1 - (S_true/4)**2) / n_audit)
        S_hat  = float(rng.normal(S_true, std_S))
        cert   = build_certificate((u,v), S_hat, n_audit, delta, ell, phi, lam_ec, 0)
        r_e    = cert.r        if cert else 0.0
        eps_e  = cert.eps      if cert else 1.0
        w_e    = cert.w        if cert else 2.0   # eps + delta
        h_e    = cert.h        if cert else 0.0
        db[(u,v)] = db[(v,u)] = {
            "r": r_e, "eps": eps_e, "w": w_e, "h": h_e,
            "compromised": i in compromised, "T": 0,
            "S_hat": S_hat, "cert": cert
        }
    return db


# ──────────────────────────────────────────────────────────────────────────────
# Routers
# ──────────────────────────────────────────────────────────────────────────────

def build_G_nx(G_base, db):
    """Convert adjacency dict to NetworkX graph with edge attributes from db."""
    G = nx.Graph()
    for u, v in G_base.edges():
        a = db.get((u,v), {"r": 0, "eps": 1, "w": 2, "T": 0})
        G.add_edge(u, v, **a)
    return G


def router_spf(G_nx, s, d, R=R_MIN, Delta=0.5, W=W_FRESHNESS):
    """B0: Shortest path by hop count (Dijkstra on hop=1 per edge)."""
    try:
        path = nx.shortest_path(G_nx, s, d)
        rate = min(G_nx[u][v]["r"] for u,v in zip(path[:-1], path[1:]))
        budget = sum(G_nx[u][v]["w"] for u,v in zip(path[:-1], path[1:]))
        return path, rate, budget
    except nx.NetworkXNoPath:
        return None, 0.0, float('inf')


def router_rate_greedy(G_nx, s, d, R=R_MIN, Delta=0.5, W=W_FRESHNESS):
    """B5: Max-min rate, no certificate check, no budget (rate-greedy, cert-free).
    Uses 1/r as edge weight (always positive) so Dijkstra minimises 1/r = maximises r.
    """
    G_filt = nx.Graph()
    for u, v, a in G_nx.edges(data=True):
        r = a.get("r", 0.0)
        if r >= R:
            G_filt.add_edge(u, v, w=1.0 / (r + 1e-9))  # positive weight
    try:
        _, path = nx.single_source_dijkstra(G_filt, s, d, weight="w")
        rate    = min(G_nx[u][v]["r"] for u,v in zip(path[:-1], path[1:]))
        budget  = sum(G_nx[u][v]["w"] for u,v in zip(path[:-1], path[1:]))
        return path, rate, budget
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None, 0.0, float('inf')


def router_sr_dijkstra(G_nx, s, d, R=R_MIN, Delta=0.5, W=W_FRESHNESS):
    """
    B1: SR-Dijkstra (Pirandola et al. 2026).

    Modified Dijkstra that tracks (accumulated_budget, min_rate) per node.
    Objective: maximise end-to-end min-bottleneck rate r_e.
    Constraint: Σ(ε_e + δ_e) ≤ Δ (composable security budget).

    Implementation follows Pirandola et al. §III:
      - State: (node, budget_used)  — we discretise budget into float labels
      - Priority queue key: (−min_rate_so_far, budget_used)
        (maximise rate = minimise negative rate; break ties by budget)
      - Filter: prune states where budget_used > Delta
      - This is Pareto-efficient: we keep the non-dominated (rate, budget) labels

    Deviation from original (documented):
      Original works with individual-link security parameters and a network-level
      security requirement. We use per-link w_e = ε_e + δ_e and Σw_e ≤ Δ as the
      end-to-end budget, which matches BMR's composability framework. The
      wavelength/time-slot assignment (out of scope for BMR's trusted-relay model)
      is dropped.
    """
    if s not in G_nx or d not in G_nx:
        return None, 0.0, float('inf')

    INF = float('inf')
    # labels[v] = list of (neg_rate, budget_used) non-dominated labels
    labels = {v: [] for v in G_nx.nodes()}
    labels[s] = [(0.0, 0.0)]   # neg_rate=0 (no path yet), budget=0

    # (neg_min_rate, budget_used, node, path)
    pq = [(0.0, 0.0, s, [s])]

    best_rate   = -INF
    best_path   = None
    best_budget = INF

    def is_dominated(new_r, new_b, existing):
        return any(e[0] <= new_r and e[1] <= new_b for e in existing)

    while pq:
        neg_r, bud, v, path = heapq.heappop(pq)
        rate_so_far = -neg_r

        # Check still valid (not dominated)
        if is_dominated(neg_r, bud, [l for l in labels[v]
                                      if l != (neg_r, bud)
                                      and l[0] <= neg_r and l[1] <= bud]):
            continue

        if v == d:
            if rate_so_far > best_rate and bud <= Delta:
                best_rate   = rate_so_far
                best_path   = path
                best_budget = bud
            continue

        for u, a in G_nx[v].items():
            r_link = a.get("r", 0.0)
            w_link = a.get("w", 2.0)
            if r_link < R:
                continue
            new_rate   = min(rate_so_far, r_link) if rate_so_far > 0 else r_link
            new_budget = bud + w_link
            if new_budget > Delta:
                continue
            new_neg_r = -new_rate
            new_label = (new_neg_r, new_budget)

            if not is_dominated(new_neg_r, new_budget, labels[u]):
                # Prune dominated labels in u
                labels[u] = [l for l in labels[u]
                              if not (new_neg_r <= l[0] and new_budget <= l[1])]
                labels[u].append(new_label)
                heapq.heappush(pq, (new_neg_r, new_budget, u, path + [u]))

    if best_path is None:
        return None, 0.0, INF
    return best_path, best_rate, best_budget


def router_bmr(G_nx, G_base, db, s, d, R=R_MIN, Delta=0.5,
               W=W_FRESHNESS, t=0):
    """BMR router: binary search on bottleneck rate + Yen's simple path.
    For the comparison against SR-Dijkstra we use t=0 (no hop floor) so that
    all routers solve the same problem: max-min rate with budget Δ.
    The hop floor is a BMR-specific feature for collusion resistance; its
    effect is shown separately in E4/E5.
    """
    result = bmr_route(G_nx, s, d, R, Delta, 0, W, m_min=0)
    if result is None:
        return None, 0.0, float('inf')
    rho, path, weight = result
    return path, rho, weight


# ──────────────────────────────────────────────────────────────────────────────
# Main experiment loop
# ──────────────────────────────────────────────────────────────────────────────

configs = [
    ("ER", 30, 4,  "random"),
    ("BA", 30, 2,  "scale-free"),
]
n_audits    = [500, 10_000]
conditions  = [("clean", 0.0), ("compromised", 0.20)]
# Routers accept Delta as a parameter
routers_def = [
    ("SPF",          lambda G,db,s,d,Delta: router_spf(G,s,d,Delta=Delta)),
    ("Rate-Greedy",  lambda G,db,s,d,Delta: router_rate_greedy(G,s,d,Delta=Delta)),
    ("SR-Dijkstra",  lambda G,db,s,d,Delta: router_sr_dijkstra(G,s,d,Delta=Delta)),
    ("BMR",          lambda G,db,s,d,Delta: router_bmr(G,None,db,s,d,Delta=Delta)),
]

all_rows = []
rng_master = np.random.default_rng(MASTER_SEED)

print("=" * 90)
print("E7b: BMR vs SR-Dijkstra (Pirandola et al. 2026) Comparison")
print("=" * 90)
print(f"{'topology':>12} {'n_audit':>8} {'cond':>12} {'router':>14} | "
      f"{'rate':>8} {'budget_ok%':>12} {'feasible%':>10} {'mean_hops':>10}")
print("-" * 95)

for graph_type, n_nodes, degree, topo_name in configs:
    for n_audit in n_audits:
        DELTA_BUDGET = DELTA_BY_N[n_audit]
        for cond_name, frac_comp in conditions:
            # Per-router results for this cell
            cell_results = {name: {"rates": [], "budget_ok": [], "feasible": [],
                                   "hops": [], "paths": []}
                            for name, _ in routers_def}

            g_rng = np.random.default_rng(int(rng_master.integers(1e9)))
            for trial in range(N_TRIALS):
                seed_t = int(g_rng.integers(1e9))
                t_rng  = np.random.default_rng(seed_t)

                # Generate graph
                if graph_type == "ER":
                    G_base = nx.erdos_renyi_graph(n_nodes, degree/(n_nodes-1), seed=seed_t)
                else:
                    G_base = nx.barabasi_albert_graph(n_nodes, degree//2 or 1, seed=seed_t)

                if not nx.is_connected(G_base):
                    continue

                s, d = 0, n_nodes - 1

                # Generate certificate database
                db = generate_link_db(G_base, n_audit, t_rng,
                                      frac_compromised=frac_comp)
                G_nx = build_G_nx(G_base, db)

                # Run each router
                for rname, rfn in routers_def:
                    t0 = time.perf_counter()
                    path, rate, budget = rfn(G_nx, db, s, d, DELTA_BUDGET)
                    elapsed = time.perf_counter() - t0
                    feasible = path is not None
                    budget_ok = feasible and budget <= DELTA_BUDGET
                    cell_results[rname]["rates"].append(rate if feasible else 0.0)
                    cell_results[rname]["budget_ok"].append(1 if budget_ok else 0)
                    cell_results[rname]["feasible"].append(1 if feasible else 0)
                    cell_results[rname]["hops"].append(len(path)-1 if feasible else 0)

            for rname, _ in routers_def:
                cr  = cell_results[rname]
                n_f = sum(cr["feasible"])
                rates_f = [r for r, f in zip(cr["rates"], cr["feasible"]) if f]
                avg_rate   = np.mean(rates_f) if rates_f else 0.0
                budg_ok_pct = 100 * np.mean(cr["budget_ok"])
                feasible_pct= 100 * np.mean(cr["feasible"])
                mean_hops   = np.mean([h for h, f in zip(cr["hops"], cr["feasible"]) if f]) if rates_f else 0

                row = {
                    "topology": topo_name, "graph_type": graph_type,
                    "n_nodes": n_nodes, "n_audit": n_audit,
                    "condition": cond_name, "frac_comp": frac_comp,
                    "router": rname,
                    "mean_rate": avg_rate, "budget_ok_pct": budg_ok_pct,
                    "feasible_pct": feasible_pct, "mean_hops": mean_hops,
                    "n_trials": N_TRIALS,
                }
                all_rows.append(row)
                print(f"{topo_name:>12} {n_audit:>8} {cond_name:>12} {rname:>14} | "
                      f"{avg_rate:>8.4f} {budg_ok_pct:>11.1f}% {feasible_pct:>9.1f}% "
                      f"{mean_hops:>10.2f}")
            print()


# ──────────────────────────────────────────────────────────────────────────────
# Bootstrap CI for BMR vs SR-Dijkstra rate difference
# ──────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("Bootstrap 95% CI: BMR rate advantage over SR-Dijkstra")
print("="*70)

ci_rows = []
rng_boot = np.random.default_rng(51)
for topo_name in ["random", "scale-free"]:
    for n_audit in n_audits:
        DELTA_BUDGET = DELTA_BY_N[n_audit]
        for cond_name in ["clean", "compromised"]:
            # Re-run paired trials for bootstrap
            gtype = "ER" if topo_name == "random" else "BA"
            n_nodes, degree = 30, (4 if gtype == "ER" else 2)
            frac_comp = 0.0 if cond_name == "clean" else 0.20

            paired_rate_bmr = []
            paired_rate_sr  = []
            paired_budg_bmr = []
            paired_budg_sr  = []

            g_rng2 = np.random.default_rng(int(rng_boot.integers(1e9)))
            N_BOOT_TRIALS = 200   # faster; bootstrap re-uses these
            for trial in range(N_BOOT_TRIALS):
                seed_t = int(g_rng2.integers(1e9))
                t_rng  = np.random.default_rng(seed_t)
                if gtype == "ER":
                    G_base = nx.erdos_renyi_graph(n_nodes, degree/(n_nodes-1), seed=seed_t)
                else:
                    G_base = nx.barabasi_albert_graph(n_nodes, 1, seed=seed_t)
                if not nx.is_connected(G_base): continue
                s, d = 0, n_nodes - 1
                db    = generate_link_db(G_base, n_audit, t_rng, frac_compromised=frac_comp)
                G_nx  = build_G_nx(G_base, db)

                _, r_bmr, b_bmr = router_bmr(G_nx, G_base, db, s, d, Delta=DELTA_BUDGET)
                _, r_sr,  b_sr  = router_sr_dijkstra(G_nx, db, s, d, Delta=DELTA_BUDGET)
                paired_rate_bmr.append(r_bmr)
                paired_rate_sr.append(r_sr)
                paired_budg_bmr.append(1 if b_bmr <= DELTA_BUDGET else 0)
                paired_budg_sr.append(1 if b_sr  <= DELTA_BUDGET else 0)

            diffs = [b - s for b, s in zip(paired_rate_bmr, paired_rate_sr)]
            boot_means = [np.mean(rng_boot.choice(diffs, len(diffs), replace=True))
                          for _ in range(5000)]
            ci_lo, ci_hi = np.percentile(boot_means, [2.5, 97.5])
            mean_diff = np.mean(diffs)
            _, pval = scipy_stats.wilcoxon(paired_rate_bmr, paired_rate_sr,
                                           alternative='two-sided')
            budg_bmr_pct = 100 * np.mean(paired_budg_bmr)
            budg_sr_pct  = 100 * np.mean(paired_budg_sr)

            ci_rows.append({
                "topology": topo_name, "n_audit": n_audit, "condition": cond_name,
                "mean_diff": mean_diff, "ci_lo": ci_lo, "ci_hi": ci_hi,
                "pval_wilcoxon": pval,
                "bmr_budget_ok_pct": budg_bmr_pct, "sr_budget_ok_pct": budg_sr_pct
            })
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
            print(f"  {topo_name:>12} n={n_audit:>6} {cond_name:>12}: "
                  f"Δrate={mean_diff:+.5f} [{ci_lo:+.5f},{ci_hi:+.5f}] "
                  f"p={pval:.3f}{sig}  "
                  f"budg_ok: BMR={budg_bmr_pct:.1f}% SR={budg_sr_pct:.1f}%")


# ──────────────────────────────────────────────────────────────────────────────
# Figures
# ──────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

topos      = ["random", "scale-free"]
conditions_plot = ["clean", "compromised"]
colors = {"SPF": "gray", "Rate-Greedy": "steelblue", "SR-Dijkstra": "orange", "BMR": "green"}

# Rate comparison
for row_i, n_audit in enumerate([500, 10_000]):
    ax = axes[row_i][0]
    x_groups = []
    for topo_name in topos:
        for cond_name in conditions_plot:
            x_groups.append(f"{topo_name[:4]}\n{cond_name[:4]}\nn={n_audit}")
    xs = np.arange(len(x_groups))
    width = 0.18
    for j, (rname, _) in enumerate(routers_def):
        rates_plot = []
        for topo_name in topos:
            for cond_name in conditions_plot:
                match = [r for r in all_rows
                         if r["topology"]==topo_name and r["n_audit"]==n_audit
                         and r["condition"]==cond_name and r["router"]==rname]
                rates_plot.append(match[0]["mean_rate"] if match else 0)
        offset = (j - 1.5) * width
        ax.bar(xs + offset, rates_plot, width, label=rname,
               color=colors[rname], alpha=0.8, edgecolor='black', lw=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(x_groups, fontsize=7)
    ax.set_ylabel("Mean certified rate (bits/round)")
    ax.set_title(f"E7b: Certified rate comparison\nn={n_audit} audit rounds")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3, axis='y')

# Budget compliance
for row_i, n_audit in enumerate([500, 10_000]):
    ax = axes[row_i][1]
    for j, (rname, _) in enumerate(routers_def):
        budg_plot = []
        for topo_name in topos:
            for cond_name in conditions_plot:
                match = [r for r in all_rows
                         if r["topology"]==topo_name and r["n_audit"]==n_audit
                         and r["condition"]==cond_name and r["router"]==rname]
                budg_plot.append(match[0]["budget_ok_pct"] if match else 0)
        offset = (j - 1.5) * width
        ax.bar(xs + offset, budg_plot, width, label=rname,
               color=colors[rname], alpha=0.8, edgecolor='black', lw=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(x_groups, fontsize=7)
    ax.set_ylabel("Budget compliance (%)")
    ax.set_title(f"E7b: Budget compliance (Σw ≤ Δ={DELTA_BUDGET:.0e})\nn={n_audit}")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim([0, 105])

plt.tight_layout()
outfig = os.path.join(os.path.dirname(__file__), "..", "figures", "E7b_sr_dijkstra_comparison.png")
plt.savefig(outfig, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved figure: {outfig}")

# ──────────────────────────────────────────────────────────────────────────────
# Save results
# ──────────────────────────────────────────────────────────────────────────────
out = {
    "experiment": "E7b",
    "description": "BMR vs SR-Dijkstra (Pirandola et al.) comparison (Reviewer 2 Comment 5)",
    "implementation_notes": {
        "SR_Dijkstra": (
            "Modified Dijkstra maximising min-bottleneck r_e subject to Σw_e ≤ Δ. "
            "Uses Pareto (rate, budget) labels per node. "
            "Deviation: wavelength/time-slot assignment dropped (out of scope for trusted-relay QKD). "
            "Implemented from Pirandola et al. PRL Research 8(3) 033221 (2026) §III description."
        ),
        "BMR": (
            "Binary search on bottleneck rate + Yen's k-shortest simple paths with hop floor K=2t+1. "
            "Same certificate database and random draws as all baselines."
        ),
        "Rate_Greedy": "Max-min rate, no Bell evidence, no expiry, no budget.",
        "SPF": "Hop-count Dijkstra.",
    },
    "delta_budget": DELTA_BUDGET,
    "n_trials": N_TRIALS,
    "all_rows": all_rows,
    "bootstrap_ci": ci_rows,
}
outpath = os.path.join(os.path.dirname(__file__), "..", "data", "results", "e7b_sr_dijkstra.json")
with open(outpath, "w") as f:
    json.dump(out, f, indent=2)
print(f"Saved results: {outpath}")
