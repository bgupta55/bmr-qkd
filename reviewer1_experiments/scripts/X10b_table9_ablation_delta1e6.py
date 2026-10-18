"""
X10b_table9_ablation_delta1e6.py
─────────────────────────────────
GAP D fix: produce Table 9 (ablation study) at the recommended manuscript
operating point Δ = 1e-6.

r1_v2.md GAP D: "Use Δ=1e-6 (78% feasible, hops=3.85) for main ablation
table; keep Δ sweep figure for binding demonstration."

Uses the IDENTICAL routing functions and graph generator as X10 v2 
(copy-imported) to guarantee consistency.  Only Delta changes: 5e-7 → 1e-6.
"""

import math
import json
import random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Same constants as X10 ─────────────────────────────────────────────────────
SEED          = 42
N_GRAPHS      = 100
n_nodes       = 20
n_audit       = 10_000
R             = 0.10
delta_max     = 3e-5
delta0        = 1e-6
H_req         = 6
tau_max_val   = 100
current_epoch = 5
DENSITY       = 0.4

# GAP D: Δ=1e-6 (78% feasibility operating point)
DELTA_GAP_D   = 1e-6

random.seed(SEED)
np.random.seed(SEED)

sqrt2 = math.sqrt(2)

# ── Copy-exact helpers from X10 ───────────────────────────────────────────────
CLASSES = [
    {"name": "high",  "eps_c": 2e-7,  "gamma_c": 1e-7,  "S_mean": 2.720, "prob": 0.3},
    {"name": "std",   "eps_c": 1e-6,  "gamma_c": 8e-7,  "S_mean": 2.687, "prob": 0.5},
    {"name": "deg",   "eps_c": 5e-6,  "gamma_c": 4e-6,  "S_mean": 2.640, "prob": 0.2},
]

def h_pironio(S):
    if S <= 2.0: return 0.0
    S = min(S, 2*sqrt2 - 1e-10)
    inner = max(0.0, 2.0 - S**2/4.0)
    return 1.0 - math.log2(1.0 + math.sqrt(inner))

def tau_paper(n, delta=1e-6):
    return math.sqrt(32 * math.log(8/delta) / n)

def sigma_physical(n):
    return 2.0 / math.sqrt(n)

def link_expiry_prob(epoch, tau_max):
    return 1.0 - math.exp(-epoch / tau_max)

def hop_cap(delta_max_v, delta0_v, H_req_v):
    return min(H_req_v, int(math.floor(delta_max_v / delta0_v)))

def make_graph(n_nodes_v, density=0.4, seed=0, n_audit_v=10000,
               current_epoch_v=0, tau_max_v=100):
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)
    sigma  = sigma_physical(n_audit_v)
    t_cert = tau_paper(n_audit_v)
    edges  = {}
    for u in range(n_nodes_v):
        for v in range(u+1, n_nodes_v):
            if rng.random() < density:
                cls = rng.choices(CLASSES, weights=[c["prob"] for c in CLASSES])[0]
                S_hat = float(np_rng.normal(cls["S_mean"], sigma))
                S_lb  = max(0.0, S_hat - t_cert)
                r_e   = h_pironio(S_lb)
                eps_e = cls["eps_c"]
                gamma_e = cls["gamma_c"]
                p_exp   = link_expiry_prob(current_epoch_v, tau_max_v)
                expired = rng.random() < p_exp
                gate_passed = (S_lb >= 2.0)
                edges[(u,v)] = edges[(v,u)] = {
                    "S_hat": S_hat, "S_lb": S_lb, "r_e": r_e,
                    "eps_e": eps_e, "gamma_e": gamma_e,
                    "expired": expired, "gate_passed": gate_passed,
                    "cls": cls["name"],
                }
    return edges

def route_full(n_nodes_v, edges, src, dst, Delta, R_v, H):
    best = {"rate": -1.0, "path": None, "violation": False}
    stack = [(src, 0.0, 0.0, [src], False)]
    while stack:
        u, cum_eps, cum_r, path, has_violation = stack.pop()
        if u == dst and cum_eps <= Delta:
            if cum_r > best["rate"]:
                best.update({"rate": cum_r, "path": path[:], "violation": has_violation})
            continue
        if len(path) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if not ed["gate_passed"]: continue
            if ed["expired"]: continue
            if ed["r_e"] < R_v: continue
            new_eps = cum_eps + ed["eps_e"]
            if new_eps <= Delta:
                stack.append((vv, new_eps, cum_r + ed["r_e"], path + [vv], has_violation))
    return best["path"], best["rate"], best["violation"]

def route_no_defect(n_nodes_v, edges, src, dst, R_v, H):
    best = {"rate": -1.0, "path": None}
    stack = [(src, 0.0, [src])]
    while stack:
        u, cum_r, path = stack.pop()
        if u == dst:
            if cum_r > best["rate"]:
                best.update({"rate": cum_r, "path": path[:]})
            continue
        if len(path) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if not ed["gate_passed"]: continue
            if ed["expired"]: continue
            if ed["r_e"] < R_v: continue
            stack.append((vv, cum_r + ed["r_e"], path + [vv]))
    return best["path"], best["rate"]

def route_no_freshness(n_nodes_v, edges, src, dst, Delta, R_v, H):
    best = {"rate": -1.0, "path": None, "violation": False}
    stack = [(src, 0.0, 0.0, [src], False)]
    while stack:
        u, cum_eps, cum_r, path, has_violation = stack.pop()
        if u == dst and cum_eps <= Delta:
            if cum_r > best["rate"]:
                best.update({"rate": cum_r, "path": path[:], "violation": has_violation})
            continue
        if len(path) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if not ed["gate_passed"]: continue
            if ed["r_e"] < R_v: continue
            new_eps = cum_eps + ed["eps_e"]
            if new_eps <= Delta:
                new_viol = has_violation or ed["expired"]
                stack.append((vv, new_eps, cum_r + ed["r_e"], path + [vv], new_viol))
    return best["path"], best["rate"], best["violation"]

def route_hop_only(n_nodes_v, edges, src, dst, R_v, H):
    best = {"rate": -1.0, "path": None, "violation": False}
    stack = [(src, 0.0, [src], False)]
    while stack:
        u, cum_r, path, has_violation = stack.pop()
        if u == dst:
            if cum_r > best["rate"]:
                best.update({"rate": cum_r, "path": path[:], "violation": has_violation})
            continue
        if len(path) > H: continue
        for (uu, vv), ed in edges.items():
            if uu != u or vv in path: continue
            if not ed["gate_passed"]: continue
            if ed["r_e"] < R_v: continue
            new_viol = has_violation or ed["expired"]
            stack.append((vv, cum_r + ed["r_e"], path + [vv], new_viol))
    return best["path"], best["rate"], best["violation"]

def eval_path(path, edges):
    if path is None or len(path) < 2: return None
    r_sum = 0.0; eps_sum = 0.0; hops = len(path) - 1; n_expired = 0
    for i in range(len(path) - 1):
        u, v = path[i], path[i+1]
        ed = edges.get((u,v), edges.get((v,u)))
        if not ed: return None
        r_sum   += ed["r_e"]
        eps_sum += ed["eps_e"]
        if ed["expired"]: n_expired += 1
    return {"r_sum": r_sum, "eps_sum": eps_sum, "hops": hops,
            "security_violation": n_expired > 0}


def main():
    print("X10b: Table 9 Ablation at GAP D Operating Point — Δ=1e-6")
    print("=" * 65)
    H = hop_cap(delta_max, delta0, H_req)
    print(f"  N_GRAPHS={N_GRAPHS}, n_nodes={n_nodes}, n_audit={n_audit}")
    print(f"  Δ = {DELTA_GAP_D:.1e}  ← GAP D recommended operating point (78% feasible)")
    print(f"  H = min({H_req}, floor({delta_max:.0e}/{delta0:.0e})) = {H}")
    print(f"  R={R}, epoch={current_epoch}, τ_max={tau_max_val}")
    P_exp = link_expiry_prob(current_epoch, tau_max_val)
    print(f"  P(link expired by epoch {current_epoch}) = {P_exp:.4f} ({P_exp*100:.2f}%)")
    print()

    variants  = ["full", "no_defect", "no_freshness", "hop_only"]
    accum     = {v: {"rate": [], "hops": [], "sec_violation": [],
                     "eps_sum": [], "feasible": 0} for v in variants}

    for seed in range(N_GRAPHS):
        g   = make_graph(n_nodes, DENSITY, seed, n_audit, current_epoch, tau_max_val)
        src = 0
        dst = n_nodes - 1

        p, _, _ = route_full(n_nodes, g, src, dst, DELTA_GAP_D, R, H)
        ev = eval_path(p, g)
        if ev:
            accum["full"]["feasible"] += 1
            accum["full"]["rate"].append(ev["r_sum"])
            accum["full"]["hops"].append(ev["hops"])
            accum["full"]["sec_violation"].append(float(ev["security_violation"]))
            accum["full"]["eps_sum"].append(ev["eps_sum"])

        p, _ = route_no_defect(n_nodes, g, src, dst, R, H)
        ev = eval_path(p, g)
        if ev:
            accum["no_defect"]["feasible"] += 1
            accum["no_defect"]["rate"].append(ev["r_sum"])
            accum["no_defect"]["hops"].append(ev["hops"])
            accum["no_defect"]["sec_violation"].append(float(ev["security_violation"]))
            accum["no_defect"]["eps_sum"].append(ev["eps_sum"])

        p, _, _ = route_no_freshness(n_nodes, g, src, dst, DELTA_GAP_D, R, H)
        ev = eval_path(p, g)
        if ev:
            accum["no_freshness"]["feasible"] += 1
            accum["no_freshness"]["rate"].append(ev["r_sum"])
            accum["no_freshness"]["hops"].append(ev["hops"])
            accum["no_freshness"]["sec_violation"].append(float(ev["security_violation"]))
            accum["no_freshness"]["eps_sum"].append(ev["eps_sum"])

        p, _, _ = route_hop_only(n_nodes, g, src, dst, R, H)
        ev = eval_path(p, g)
        if ev:
            accum["hop_only"]["feasible"] += 1
            accum["hop_only"]["rate"].append(ev["r_sum"])
            accum["hop_only"]["hops"].append(ev["hops"])
            accum["hop_only"]["sec_violation"].append(float(ev["security_violation"]))
            accum["hop_only"]["eps_sum"].append(ev["eps_sum"])

    # ── Print & collect results ───────────────────────────────────────────────
    print(f"{'Variant':>16}  {'Feas/100':>10}  {'Rate':>10}  "
          f"{'SD':>8}  {'Hops':>8}  {'Sec.Viol%':>12}  {'MeanΣε':>12}")
    print("-" * 85)

    ablation_results = []
    for v in variants:
        a   = accum[v]
        n_f = a["feasible"]
        r_m = float(np.mean(a["rate"]))     if a["rate"] else 0.0
        r_s = float(np.std(a["rate"]))      if a["rate"] else 0.0
        h_m = float(np.mean(a["hops"]))     if a["hops"] else 0.0
        sv  = 100*float(np.mean(a["sec_violation"])) if a["sec_violation"] else 0.0
        em  = float(np.mean(a["eps_sum"]))  if a["eps_sum"] else 0.0
        print(f"{v:>16}  {n_f:>10}  {r_m:>10.4f}  "
              f"{r_s:>8.4f}  {h_m:>8.2f}  {sv:>12.2f}  {em:>12.3e}")
        ablation_results.append({
            "variant": v, "n_feasible": n_f,
            "rate_mean": round(r_m, 6), "rate_sd": round(r_s, 6),
            "hops_mean": round(h_m, 4),
            "sec_violation_pct": round(sv, 3),
            "eps_sum_mean": round(em, 10),
        })

    full_row    = next(r for r in ablation_results if r["variant"] == "full")
    nodef_row   = next(r for r in ablation_results if r["variant"] == "no_defect")
    nofresh_row = next(r for r in ablation_results if r["variant"] == "no_freshness")
    hoponly_row = next(r for r in ablation_results if r["variant"] == "hop_only")

    feas_gap = nodef_row["n_feasible"] - full_row["n_feasible"]
    rate_gap = nodef_row["rate_mean"]  - full_row["rate_mean"]
    print()
    print(f"  Feasibility gap (no_defect − full)  : {feas_gap}pp")
    print(f"  Rate gap (no_defect − full)          : {rate_gap:+.4f} bits/round")
    print(f"  hop_only security violation rate     : {hoponly_row['sec_violation_pct']:.2f}%")
    print(f"  no_freshness security violation rate : {nofresh_row['sec_violation_pct']:.2f}%")
    print(f"  full security violation rate         : {full_row['sec_violation_pct']:.2f}% (should be 0)")

    # Verify consistency with delta_sweep row at 1e-6 from X10 results
    print()
    print("  Expected (from X10 delta_sweep): feasibility=0.78, rate≈1.02, hops≈3.85")
    print(f"  Got (full variant)             : feasibility={full_row['n_feasible']/100:.2f}, "
          f"rate={full_row['rate_mean']:.3f}, hops={full_row['hops_mean']:.2f}")

    # ── JSON ──────────────────────────────────────────────────────────────────
    output = {
        "description": (
            "X10b: Table 9 ablation at GAP D operating point Δ=1e-6. "
            "Identical routing functions/graphs as X10 v2, only Delta changed 5e-7→1e-6."
        ),
        "design": {
            "N_GRAPHS": N_GRAPHS, "n_nodes": n_nodes, "n_audit": n_audit,
            "Delta": DELTA_GAP_D,
            "R": R, "H": H, "delta_0": delta0, "delta_max": delta_max,
            "tau_max": tau_max_val, "current_epoch": current_epoch,
            "P_link_expired": round(P_exp, 6),
            "GAP_D_note": "Δ=1e-6 is the recommended Table 9 operating point per r1_v2.md",
        },
        "ablation_table9_delta1e6": ablation_results,
        "summary": {
            "feasibility_gap_pp": feas_gap,
            "rate_gap_bits_per_round": round(rate_gap, 4),
            "hop_only_sec_violation_pct": hoponly_row["sec_violation_pct"],
            "no_freshness_sec_violation_pct": nofresh_row["sec_violation_pct"],
        },
    }
    with open("X10b_table9_results.json", "w") as f:
        json.dump(output, f, indent=2)
    print("\nSaved: X10b_table9_results.json")

    # ── Figure: 3-panel ablation bar chart ───────────────────────────────────
    labels = [r["variant"].replace("_", "\n") for r in ablation_results]
    colors = ["steelblue", "coral", "goldenrod", "mediumpurple"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        f"Table 9 Ablation Study — Operating Point Δ = 1×10⁻⁶  (GAP D recommended)\n"
        f"N={N_GRAPHS} graphs, n_nodes={n_nodes}, n_audit={n_audit:,}",
        fontsize=12, fontweight="bold",
    )

    # Panel 1: Feasibility
    ax = axes[0]
    feasibles = [r["n_feasible"] for r in ablation_results]
    bars = ax.bar(labels, feasibles, color=colors, edgecolor="k", alpha=0.85)
    for bar, val in zip(bars, feasibles):
        ax.text(bar.get_x() + bar.get_width()/2, val + 1.5, str(val),
                ha="center", va="bottom", fontweight="bold", fontsize=11)
    ax.set_ylim(0, 115)
    ax.set_title("Feasibility (out of 100 graphs)", fontweight="bold")
    ax.set_ylabel("Feasible graphs")
    ax.axhline(100, color="gray", ls="--", alpha=0.4, lw=1)

    # Panel 2: Mean rate with SD
    ax = axes[1]
    rates_v = [r["rate_mean"] for r in ablation_results]
    sds_v   = [r["rate_sd"]   for r in ablation_results]
    bars = ax.bar(labels, rates_v, yerr=sds_v, color=colors,
                  edgecolor="k", alpha=0.85, capsize=5)
    for bar, val in zip(bars, rates_v):
        ax.text(bar.get_x() + bar.get_width()/2, val + max(sds_v)*0.05 + 0.01,
                f"{val:.3f}", ha="center", va="bottom", fontsize=9)
    ax.set_title("Mean Σr_e (certified rate, bits/round)", fontweight="bold")
    ax.set_ylabel("Mean Σr_e")

    # Panel 3: Security violation rate
    ax = axes[2]
    sv_v = [r["sec_violation_pct"] for r in ablation_results]
    bars = ax.bar(labels, sv_v, color=colors, edgecolor="k", alpha=0.85)
    for bar, val in zip(bars, sv_v):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width()/2, val + 0.3,
                    f"{val:.1f}%", ha="center", va="bottom",
                    fontweight="bold", fontsize=10, color="red")
        else:
            ax.text(bar.get_x() + bar.get_width()/2, 0.3,
                    "0%", ha="center", va="bottom", fontsize=9, color="green")
    ax.set_ylim(0, max(sv_v)*1.3 + 2 if max(sv_v) > 0 else 5)
    ax.set_title("Security Violation Rate (%)\n(expired cert used on path)", fontweight="bold")
    ax.set_ylabel("Violation rate (%)")

    plt.tight_layout()
    plt.savefig("fig_X10b_table9_ablation.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved: fig_X10b_table9_ablation.png")


if __name__ == "__main__":
    main()
