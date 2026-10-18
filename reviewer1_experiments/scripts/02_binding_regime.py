"""
Experiment 2: Binding-Regime Study
====================================
Addresses Reviewer 1, Comment 1 (P1 in r1_v1.md):
  "The security budgets never bind in the headline runs"
  
Goal: choose n and Δ so that Σε_e is comparable to Δ.
For ε_e ≈ delta_e (audit failure prob), we need Δ ~ H*delta_e 
where H = hop count.

With corrected pipeline: ε_e is now the security defect accumulated 
per link. For the binding regime:
  - n ∈ {1000, 1500, 2500} per setting
  - Δ swept over [1e-8, 1e-4, 1e-2, 0.05, 0.1, 0.5, 1.0]
  - H (hop count) = 3, 4, 5

Reports:
  1. Per-link ε_e at each n
  2. Path defect Σε_e for H hops
  3. Feasibility: Σε_e <= Δ  (budget binding or not)
  4. Certified path rate min(r_e) per link
  5. Delta sweep: feasibility rate vs Delta

Simulator model:
  - Links draw S ~ N(S0, sigma^2)  where sigma controls noise
  - Per-link: compute r_e, ε_e from corrected pipeline
  - Measure feasibility fraction across Monte Carlo trials
"""

import math
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# ─── Certificate functions ────────────────────────────────────────────────────
def h_pironio(S):
    if S <= 2.0: return 0.0
    if S > 2*math.sqrt(2): S = 2*math.sqrt(2)
    return 1.0 - math.log2(1 + math.sqrt(2 - S**2/4))

def tau(n, delta):
    return math.sqrt(32*math.log(8/delta)/n)

def r_e_from_S(S_hat, n, delta):
    slb = S_hat - tau(n, delta)
    return h_pironio(slb) if slb > 2.0 else 0.0

def eps_e_from_S(S_hat, n, delta):
    """Security defect = delta (the audit failure probability)."""
    return delta  # per corrected pipeline: ε_e = δ_e per link

# ─── Simulation parameters ────────────────────────────────────────────────────
S0    = 2.687   # honest CHSH
sigma = 0.05    # CHSH std dev per link (Gaussian noise model)
R     = 0.10    # minimum accepted rate

n_vals   = [1000, 1500, 2500, 5000, 10_000]
delta_e  = 1e-6
hop_vals = [3, 4, 5]
Delta_vals = np.logspace(-9, 0, 40)   # Δ sweep over 9 decades

n_trials = 10_000   # Monte Carlo trials per (n, H)
rng = np.random.default_rng(42)

results = {}

# ─── Compute per-link ε_e and path defect vs n ───────────────────────────────
per_link_eps = {}
per_link_re  = {}
for n in n_vals:
    # Nominal link (S_hat = S0)
    eps_nom = eps_e_from_S(S0, n, delta_e)
    re_nom  = r_e_from_S(S0, n, delta_e)
    per_link_eps[n] = eps_nom
    per_link_re[n]  = re_nom

print("=== Per-link ε_e and r_e (nominal honest link S₀=2.687) ===")
print(f"{'n':>8}  {'ε_e':>12}  {'r_e (bits/round)':>18}  {'4n total rounds':>16}")
print("-"*60)
for n in n_vals:
    print(f"{n:>8}  {per_link_eps[n]:>12.2e}  {per_link_re[n]:>18.6f}  {4*n:>16,}")

# ─── Monte Carlo: feasibility vs Δ sweep ─────────────────────────────────────
print("\n=== Binding-Regime Study: Feasibility vs Δ ===")

mc_results = {}
for n in [1000, 2500, 10_000]:
    mc_results[n] = {}
    for H in hop_vals:
        feasibility = []
        path_re_list = []
        for _ in range(n_trials):
            # Draw H links from Gaussian
            S_links = rng.normal(S0, sigma, size=H)
            eps_links = np.array([eps_e_from_S(s, n, delta_e) for s in S_links])
            re_links  = np.array([r_e_from_S(s, n, delta_e)  for s in S_links])
            total_eps = eps_links.sum()
            min_re    = re_links.min()
            # check r_e > R (rate constraint)
            rate_ok = min_re >= R
            # feasibility as a function of Delta
            feas_by_delta = [total_eps <= D and rate_ok for D in Delta_vals]
            feasibility.append(feas_by_delta)
            path_re_list.append(min_re)
        feas_arr = np.array(feasibility)  # (n_trials, n_delta)
        mc_results[n][H] = {
            "feasibility_rate": feas_arr.mean(axis=0).tolist(),
            "path_re_mean": float(np.mean(path_re_list)),
            "path_re_std":  float(np.std(path_re_list)),
            "path_re_CI_95": [float(np.percentile(path_re_list, 2.5)),
                              float(np.percentile(path_re_list, 97.5))],
            "median_total_eps": float(np.median([
                sum(eps_e_from_S(s, n, delta_e) for s in rng.normal(S0, sigma, H))
                for _ in range(1000)
            ])),
        }

# ─── Print binding point (Δ where feasibility crosses 50%) ───────────────────
print(f"\n{'n':>8}  {'H':>3}  {'Δ_bind(50%)':>14}  {'r_e mean':>10}  {'r_e CI 95%':>20}")
print("-"*70)
for n in [1000, 2500, 10_000]:
    for H in hop_vals:
        feas = np.array(mc_results[n][H]["feasibility_rate"])
        # find Delta where feasibility first crosses 0.5
        bind_idx = np.searchsorted(feas, 0.5)
        if bind_idx < len(Delta_vals):
            d_bind = Delta_vals[bind_idx]
        else:
            d_bind = float('inf')
        re_m = mc_results[n][H]["path_re_mean"]
        re_ci = mc_results[n][H]["path_re_CI_95"]
        print(f"{n:>8}  {H:>3}  {d_bind:>14.3e}  {re_m:>10.5f}  [{re_ci[0]:.5f}, {re_ci[1]:.5f}]")

# ─── Save JSON ────────────────────────────────────────────────────────────────
out = {
    "description": "Binding-regime study (Reviewer 1 P1 / Comment 1)",
    "S0": S0, "sigma": sigma, "R": R, "delta_e": delta_e,
    "n_trials": n_trials,
    "Delta_sweep": Delta_vals.tolist(),
    "per_link": {str(n): {"eps_e": per_link_eps[n], "r_e": per_link_re[n]} for n in n_vals},
    "mc_results": {str(n): {str(H): v for H, v in mc_results[n].items()} for n in mc_results},
}
with open("REVISION_OSN_12OCT2026/expr_r1/02_binding_regime_results.json", "w") as f:
    json.dump(out, f, indent=2)

# ─── Figures ─────────────────────────────────────────────────────────────────
# Fig A: Feasibility rate vs Δ for different n (H=4)
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

H_plot = 4
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
for i, n in enumerate([1000, 2500, 10_000]):
    feas = mc_results[n][H_plot]["feasibility_rate"]
    axes[0].semilogx(Delta_vals, feas, label=f'n={n:,}', color=colors[i])
axes[0].axhline(0.5, color='gray', linestyle='--', lw=1, label='50% feasibility')
axes[0].axhline(0.9, color='red',  linestyle=':', lw=1, label='90% feasibility')
axes[0].set_xlabel('Δ (security-defect budget per path)')
axes[0].set_ylabel('Feasibility rate')
axes[0].set_title(f'Binding-Regime Study: Feasibility vs Δ  (H={H_plot} hops)')
axes[0].legend()
axes[0].grid(True, alpha=0.3)
axes[0].set_ylim([0, 1.05])

# Fig B: certified path rate (min r_e) vs n, H=4
n_arr_plot = [1000, 2500, 10_000]
re_means = []
re_lo = []
re_hi = []
for n in n_arr_plot:
    r = mc_results[n][H_plot]
    re_means.append(r["path_re_mean"])
    re_lo.append(r["path_re_CI_95"][0])
    re_hi.append(r["path_re_CI_95"][1])
err_lo_arr = np.maximum(0, np.array(re_means) - np.array(re_lo))
err_hi_arr = np.maximum(0, np.array(re_hi) - np.array(re_means))
axes[1].errorbar(n_arr_plot, re_means,
                 yerr=[err_lo_arr, err_hi_arr],
                 fmt='o-', color='#1f77b4', capsize=4, label=f'H={H_plot}')
axes[1].axhline(R, color='red', linestyle='--', label=f'R={R}')
axes[1].set_xlabel('n (per-setting audit length)')
axes[1].set_ylabel('Certified path rate min(r_e), bits/round')
axes[1].set_title('Certified Path Rate vs n  (S₀=2.687, σ=0.05)')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_02_binding_regime.png", dpi=150)
plt.close()
print("\nSaved: fig_02_binding_regime.png")
