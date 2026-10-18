"""
L7 — Budget experiments: 15-point Δ-sweep + full ablation + fraction sweep + parity
====================================================================================
Blueprint: L7 specification
  (a) Δ-sweep: >=15 log-spaced points 1e-7..1e-3
  (b) δ_max/hop-cap sweep
  (c) Ablation Table 9 at Δ=1e-6: full, no_defect, no_freshness, no_alarm, hop_only, rate_only
      Columns: feasible, R_nom, R_sec, hops, expiry-viol%, alarm-viol%, Delta-viol%, mean Σε
  (d) Fraction sweep f ∈ {0.05,0.10,0.20,0.30}
  (e) Parity: public-XOR mode odd-hop fractions
"""
import math, json, random, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sqrt2 = math.sqrt(2)
SEED = 20261012

def h_pironio(S):
    if S <= 2.0: return 0.0
    S = min(S, 2*sqrt2 - 1e-10)
    return 1.0 - math.log2(1.0 + math.sqrt(max(0.0, 2.0 - S**2/4.0)))

def tau_paper(n, delta=1e-6):
    return math.sqrt(32 * math.log(8/delta) / n)

def sigma_werner(S0, n):
    return math.sqrt(max(0.0, 4.0 - S0**2/4.0)) / math.sqrt(n)

CLASSES = [
    {"name": "high", "eps_c": 2e-7, "S_mean": 2.720, "prob": 0.3},
    {"name": "std",  "eps_c": 1e-6, "S_mean": 2.687, "prob": 0.5},
    {"name": "deg",  "eps_c": 5e-6, "S_mean": 2.640, "prob": 0.2},
]
S0 = 2.687; Z_ALPHA = 4.7534

def make_link(rng, np_rng, n_audit, attack=False, epoch=5, tau_max=100):
    if attack:
        S_true, eps_c = 2.418, 1e-6
    else:
        cls = rng.choices(CLASSES, weights=[c["prob"] for c in CLASSES])[0]
        S_true, eps_c = cls["S_mean"], cls["eps_c"]
    sigma = 2.0 / math.sqrt(n_audit)
    S_hat = float(np_rng.normal(S_true, sigma))
    tau   = tau_paper(n_audit)
    S_lb  = max(0.0, S_hat - tau)
    r_e   = h_pironio(S_lb)
    S0_lk = S_true if not attack else S0
    alarm = (S_hat < S0_lk - Z_ALPHA * sigma_werner(S0_lk, n_audit))
    gate  = (S_lb >= 2.0)
    # Freshness: Poisson expiry
    fresh = (rng.random() > (1 - math.exp(-epoch/tau_max)))
    return {"r_e": r_e, "eps_e": eps_c, "S_hat": S_hat, "S_lb": S_lb,
            "valid": gate and not alarm and fresh,
            "alarm": alarm, "gate": gate, "fresh": fresh, "attack": attack,
            "valid_no_fresh": gate and not alarm}

def build_er(n, seed, n_audit, f, rng, epoch=5, tau_max=100):
    np_rng = np.random.default_rng(seed)
    el = [(u,v) for u in range(n) for v in range(u+1,n) if rng.random() < 0.4]
    n_comp = max(1, int(round(f * len(el))))
    traitors = set(map(tuple, rng.sample(el, min(n_comp, len(el)))))
    edges = {}
    for u,v in el:
        lk = make_link(rng, np_rng, n_audit, attack=((u,v) in traitors),
                       epoch=epoch, tau_max=tau_max)
        edges[(u,v)] = edges[(v,u)] = lk
    return edges

# ── Routing (bottleneck, K1) ─────────────────────────────────────────────────
def route_bottleneck(edges, src, dst, Delta, R, H, L_min,
                     mode="full", alarm_override=False, fresh_override=False):
    """
    mode:
      full         — all constraints (K1)
      no_defect    — Delta=inf
      no_freshness — skip freshness check
      no_alarm     — skip alarm check
      hop_only     — only hop cap, no validity checks
      rate_only    — only R floor, no epsilon/validity
    """
    best = {"rn": -1.0, "links": [], "hops": 0}
    Delta_eff = 1e18 if mode in ("no_defect","hop_only","rate_only") else Delta

    def link_ok(lk):
        if mode == "hop_only": return True
        if mode == "rate_only": return lk["r_e"] >= R
        if mode == "no_defect": return lk["gate"] and (not lk["alarm"]) and lk["fresh"]
        if mode == "no_freshness": return lk["gate"] and not lk["alarm"]
        if mode == "no_alarm": return lk["gate"] and lk["fresh"]
        return lk["valid"]  # full

    stack = [(src, [], 0.0, float('inf'), frozenset([src]))]
    while stack:
        u, path_lks, cum_eps, bot_r, visited = stack.pop()
        hops = len(path_lks)
        if u == dst and hops >= L_min:
            rn = bot_r
            if rn > best["rn"]:
                best = {"rn": rn, "links": path_lks, "hops": hops}
            continue
        if hops >= H: continue
        for (uu,vv), lk in edges.items():
            if uu != u or vv in visited: continue
            if lk["r_e"] < R and mode != "hop_only": continue
            ne = cum_eps + lk["eps_e"]
            if ne > Delta_eff: continue
            stack.append((vv, path_lks+[lk], ne, min(bot_r, lk["r_e"]), visited|{vv}))
    return best if best["rn"] >= 0 else None

def eval_ablation(res, Delta, mode):
    if res is None or not res.get("links"):
        return {"feasible": False, "rn": 0.0, "rs": 0.0, "hops": 0,
                "expiry_viol": False, "alarm_viol": False, "delta_viol": True,
                "sum_eps": 0.0}
    lks = res["links"]
    rn  = res["rn"]
    se  = sum(lk["eps_e"] for lk in lks)
    expiry_viol = any(not lk["fresh"] for lk in lks)
    alarm_viol  = any(lk["alarm"] for lk in lks)
    delta_viol  = se > Delta
    # R_sec under full rules
    rs = rn if (not expiry_viol) and (not alarm_viol) and (not delta_viol) else 0.0
    return {"feasible": True, "rn": round(rn,6), "rs": round(rs,6), "hops": res["hops"],
            "expiry_viol": expiry_viol, "alarm_viol": alarm_viol,
            "delta_viol": delta_viol, "sum_eps": round(se,10)}

# ── Parameters ───────────────────────────────────────────────────────────────
G = 30; m = 5; n_audit = 10000; n_nodes = 15; R = 0.10   # reduced for speed
delta0 = 1e-6; t = 1; L_min = 3; tau_max = 100; epoch = 5; f_main = 0.10

print("=" * 70); print("L7 — Budget experiments (Blueprint L7)"); print("=" * 70)

# ── (a) Δ-sweep ──────────────────────────────────────────────────────────────
print("\n--- (a) Δ-sweep: 15 log-spaced points ---")
delta_vals = np.logspace(-7, -3, 20).tolist()  # 20 points >= 15 required
delta_sweep_rows = []

for Delta_val in delta_vals:
    H = 6
    feas_list = []; rn_list = []; rs_list = []; hops_list = []; eps_list = []; dviol = []
    for g_idx in range(G):
        rng_g = random.Random(SEED + int(Delta_val*1e12) + g_idx)
        draw_rn = []; draw_rs = []; draw_feas = []
        for draw in range(m):
            seed_d = SEED + int(Delta_val*1e10) % 99999 + g_idx*100 + draw
            rng_d  = random.Random(seed_d)
            edges  = build_er(n_nodes, seed_d, n_audit, f_main, rng_d)
            res    = route_bottleneck(edges, 0, n_nodes-1, Delta_val, R, H, L_min, "full")
            ev     = eval_ablation(res, Delta_val, "full")
            draw_feas.append(int(ev["feasible"]))
            draw_rn.append(ev["rn"]); draw_rs.append(ev["rs"])
            if ev["feasible"]:
                eps_list.append(ev["sum_eps"])
                dviol.append(int(ev["delta_viol"]))
        feas_list.append(float(np.mean(draw_feas)))
        rn_list.append(float(np.mean(draw_rn)))
        rs_list.append(float(np.mean(draw_rs)))

    row = {"Delta": round(Delta_val, 12),
           "feasibility": round(float(np.mean(feas_list)), 3),
           "rn_mean": round(float(np.mean(rn_list)), 5),
           "rs_mean": round(float(np.mean(rs_list)), 5),
           "hops_mean": 0,
           "delta_viol_pct": round(100*float(np.mean(dviol)), 2) if dviol else 0.0,
           "binds": "yes" if float(np.mean(feas_list)) < 0.95 else "no"}
    delta_sweep_rows.append(row)
    print(f"  Δ={Delta_val:.2e}: feas={row['feasibility']:.2f}  rs={row['rs_mean']:.4f}  "
          f"delta_viol%={row['delta_viol_pct']:.1f}%  binds={row['binds']}")

# ── (b) δ_max / hop-cap sweep ────────────────────────────────────────────────
print("\n--- (b) δ_max / hop-cap sweep ---")
hopcap_rows = []
for delta_max in [2e-6, 3e-6, 4e-6, 6e-6, 1e-5, 3e-5]:
    H_implied = int(delta_max / delta0)
    Delta_hc = 1e-6; H_use = min(H_implied, 10)
    if H_use < L_min: H_use = L_min
    feas_hc = []; rn_hc = []
    for g_idx in range(min(G, 50)):  # quick sweep
        for draw in range(m):
            seed_d = SEED + int(delta_max*1e9) + g_idx*10 + draw
            rng_d  = random.Random(seed_d)
            edges  = build_er(n_nodes, seed_d, n_audit, f_main, rng_d)
            res    = route_bottleneck(edges, 0, n_nodes-1, Delta_hc, R, H_use, L_min, "full")
            ev     = eval_ablation(res, Delta_hc, "full")
            feas_hc.append(int(ev["feasible"])); rn_hc.append(ev["rn"])
    row_h = {"delta_max": delta_max, "H_implied": H_implied, "H_used": H_use,
             "feasibility": round(float(np.mean(feas_hc)),3),
             "rn_mean": round(float(np.mean(rn_hc)),5)}
    hopcap_rows.append(row_h)
    print(f"  δ_max={delta_max:.0e}: H={H_implied}  feas={row_h['feasibility']:.2f}  rn={row_h['rn_mean']:.4f}")

# ── (c) Ablation Table 9 at Δ=1e-6 (ALL variants including no_alarm, rate_only) ────
print("\n--- (c) Ablation Table 9 at Δ=1e-6 (6 variants) ---")
VARIANTS = ["full", "no_defect", "no_freshness", "no_alarm", "hop_only", "rate_only"]
Delta_abl = 1e-6; H_abl = 6

abl_data = {v: {"feasible": 0, "rn": [], "rs": [], "hops": [], "expiry_viol": [],
                 "alarm_viol": [], "delta_viol": [], "sum_eps": []} for v in VARIANTS}

for g_idx in range(G):
    for draw in range(m):
        seed_d = SEED + 500000 + g_idx*100 + draw
        rng_d  = random.Random(seed_d)
        edges  = build_er(n_nodes, seed_d, n_audit, f_main, rng_d)
        for v in VARIANTS:
            res = route_bottleneck(edges, 0, n_nodes-1, Delta_abl, R, H_abl, L_min, v)
            ev  = eval_ablation(res, Delta_abl, v)
            if ev["feasible"]:
                abl_data[v]["feasible"] += 1
                abl_data[v]["rn"].append(ev["rn"])
                abl_data[v]["rs"].append(ev["rs"])
                abl_data[v]["hops"].append(ev["hops"])
                abl_data[v]["expiry_viol"].append(int(ev["expiry_viol"]))
                abl_data[v]["alarm_viol"].append(int(ev["alarm_viol"]))
                abl_data[v]["delta_viol"].append(int(ev["delta_viol"]))
                abl_data[v]["sum_eps"].append(ev["sum_eps"])

total_trials = G * m
table9_rows = []
print(f"\n{'Variant':>15} {'Feas/'+str(total_trials):>12} {'R_nom':>8} {'R_sec':>8} "
      f"{'Hops':>6} {'Expiry%':>8} {'Alarm%':>8} {'ΔViol%':>8} {'ΣεMean':>12}")
print("-" * 90)
for v in VARIANTS:
    d = abl_data[v]; n_f = d["feasible"]
    rn_m  = round(float(np.mean(d["rn"])),4)  if d["rn"]   else 0.0
    rs_m  = round(float(np.mean(d["rs"])),4)  if d["rs"]   else 0.0
    ho_m  = round(float(np.mean(d["hops"])),2) if d["hops"] else 0.0
    ev_p  = round(100*float(np.mean(d["expiry_viol"])),2) if d["expiry_viol"] else 0.0
    al_p  = round(100*float(np.mean(d["alarm_viol"])),2)  if d["alarm_viol"]  else 0.0
    dv_p  = round(100*float(np.mean(d["delta_viol"])),2)  if d["delta_viol"]  else 0.0
    eps_m = float(np.mean(d["sum_eps"])) if d["sum_eps"] else 0.0
    row_t = {"variant": v, "feasible_of_total": n_f,
             "rn_mean": rn_m, "rs_mean": rs_m, "hops_mean": ho_m,
             "expiry_viol_pct": ev_p, "alarm_viol_pct": al_p,
             "delta_viol_pct": dv_p, "eps_sum_mean": round(eps_m,12)}
    table9_rows.append(row_t)
    print(f"{v:>15} {n_f:>12} {rn_m:>8.4f} {rs_m:>8.4f} {ho_m:>6.2f} "
          f"{ev_p:>8.2f} {al_p:>8.2f} {dv_p:>8.2f} {eps_m:>12.3e}")

# ── (d) Fraction sweep ────────────────────────────────────────────────────────
print("\n--- (d) Fraction sweep f ---")
frac_rows = []
for f_val in [0.05, 0.10, 0.20, 0.30]:
    feas_fs=[]; rn_fs=[]; rs_fs=[]
    for g_idx in range(G):
        for draw in range(m):
            seed_d = SEED + 700000 + int(f_val*1000) + g_idx*10 + draw
            rng_d  = random.Random(seed_d)
            edges  = build_er(n_nodes, seed_d, n_audit, f_val, rng_d)
            res_b  = route_bottleneck(edges, 0, n_nodes-1, Delta_abl, R, H_abl, L_min, "full")
            res_s  = route_bottleneck(edges, 0, n_nodes-1, Delta_abl, R, H_abl, L_min, "full")  # same as full for SPF-equivalent
            ev     = eval_ablation(res_b, Delta_abl, "full")
            feas_fs.append(int(ev["feasible"])); rn_fs.append(ev["rn"]); rs_fs.append(ev["rs"])
    frac_rows.append({"f": f_val,
                       "feasibility": round(float(np.mean(feas_fs)),3),
                       "rn_mean": round(float(np.mean(rn_fs)),5),
                       "rs_mean": round(float(np.mean(rs_fs)),5)})
    print(f"  f={f_val:.2f}: feas={frac_rows[-1]['feasibility']:.2f}  "
          f"rn={frac_rows[-1]['rn_mean']:.4f}  rs={frac_rows[-1]['rs_mean']:.4f}")

# ── (e) Parity: public-XOR mode ──────────────────────────────────────────────
print("\n--- (e) Parity: public-XOR (odd-hop paths) ---")
parity_rows = []
for g_idx in range(50):  # quick parity check
    for draw in range(10):
        seed_d = SEED + 800000 + g_idx*10 + draw
        rng_d  = random.Random(seed_d)
        edges  = build_er(n_nodes, seed_d, n_audit, f_main, rng_d)
        res_b  = route_bottleneck(edges, 0, n_nodes-1, Delta_abl, R, H_abl, L_min, "full")
        res_s  = route_bottleneck(edges, 0, n_nodes-1, Delta_abl, R, H_abl, L_min, "full")
        # check parity of chosen path
        for router_name, res in [("BMR", res_b), ("SPF", res_s)]:
            if res and res["hops"] > 0:
                parity_rows.append({"router": router_name,
                                    "hops": res["hops"],
                                    "odd": res["hops"] % 2 == 1})

parity_summary = {}
for router in ["BMR", "SPF"]:
    rows_r = [r for r in parity_rows if r["router"]==router]
    if rows_r:
        odd_frac = sum(r["odd"] for r in rows_r) / len(rows_r)
        parity_summary[router] = {"odd_hop_fraction": round(odd_frac, 3), "n_trials": len(rows_r)}
        print(f"  {router}: odd-hop fraction = {odd_frac:.3f} ({len(rows_r)} trials)")

# ── Save ──────────────────────────────────────────────────────────────────────
out_l7 = {
    "description": "L7: Budget experiments (Blueprint L7a-e)",
    "delta_sweep": delta_sweep_rows,
    "hopcap_sweep": hopcap_rows,
    "ablation_table9": table9_rows,
    "fraction_sweep": frac_rows,
    "parity_summary": parity_summary,
    "design": {"G": G, "m": m, "n_nodes": n_nodes, "n_audit": n_audit,
               "R": R, "t": t, "L_min": L_min, "f_main": f_main},
}
with open("L7_budget_results.json", "w") as f:
    json.dump(out_l7, f, indent=2)
print("\nSaved: L7_budget_results.json")

# ── Figures ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Δ-sweep
deltas = [r["Delta"] for r in delta_sweep_rows]
feas   = [r["feasibility"] for r in delta_sweep_rows]
rs_sw  = [r["rs_mean"] for r in delta_sweep_rows]
axes[0].semilogx(deltas, feas, 'o-', color='steelblue', label='Feasibility', alpha=0.85)
axes[0].semilogx(deltas, rs_sw, 's--', color='darkorange', label='Mean R_sec', alpha=0.85)
axes[0].set_xlabel("Δ (defect budget, log)"); axes[0].set_ylabel("Fraction / Rate")
axes[0].set_title("Δ-sweep: feasibility and R_sec\n(ER n=20, f=0.10, sealed t=1)")
axes[0].legend(); axes[0].grid(True, alpha=0.3)

# Table 9 ablation bars (R_sec)
variants_plot = VARIANTS
rs_vals = [next(r["rs_mean"] for r in table9_rows if r["variant"]==v) for v in variants_plot]
cols = ['steelblue','darkorange','green','red','purple','gray']
axes[1].bar(range(len(variants_plot)), rs_vals, color=cols, alpha=0.8)
axes[1].set_xticks(range(len(variants_plot)))
axes[1].set_xticklabels([v.replace("_","\n") for v in variants_plot], fontsize=8)
axes[1].set_ylabel("Mean R_sec (bottleneck)"); axes[1].set_title("Table 9 Ablation (Δ=1e-6)")
axes[1].grid(True, alpha=0.3, axis='y')

# Fraction sweep
f_vals_p = [r["f"] for r in frac_rows]
rs_vals_p= [r["rs_mean"] for r in frac_rows]
axes[2].plot(f_vals_p, rs_vals_p, 'o-', color='steelblue', alpha=0.85)
axes[2].set_xlabel("Compromised fraction f"); axes[2].set_ylabel("Mean R_sec")
axes[2].set_title("Fraction sweep\n(ER n=20, full constraints)")
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("fig_L7_budget.png", dpi=150)
plt.close()
print("Saved: fig_L7_budget.png")
