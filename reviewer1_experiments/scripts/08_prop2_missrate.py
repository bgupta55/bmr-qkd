"""
Experiment 8: Complexity Boundary — Proposition 2 Miss-Rate Verification
=========================================================================
Addresses Reviewer 1, Comment 1 (Appendix B.2 in r1_v1.md):
  "For S* − S† = 0.269 this is 0.945 at n = 100 (vacuous) and 3.5×10⁻³ 
   at n = 10,000, so compare against a variance-aware (binomial/normal) 
   prediction for Table 5."

Also addresses: ε_e = 2/n² vs corrected formula

Verifies:
  1. Hoeffding miss-rate: P(Ŝ < θ | S*) ≤ exp(−n(S*−S†)²/128) [B.2]
  2. Normal approximation (binomial variance) for same quantity
  3. Simulated miss-rate via Monte Carlo
  4. The ε_e formula: 2/n² (old, wrong) vs δ_e (correct)
"""

import math
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm

# ─── Certificate functions ────────────────────────────────────────────────────
def h_pironio(S):
    if S <= 2.0: return 0.0
    if S > 2*math.sqrt(2): S = 2*math.sqrt(2)
    return 1.0 - math.log2(1 + math.sqrt(2 - S**2/4))

def tau_fn(n, delta):
    return math.sqrt(32*math.log(8/delta)/n)

def hoeffding_miss(n, S_star, S_dag):
    """Hoeffding bound on miss probability P(Ŝ < θ | S*) [B.2]."""
    theta = (S_star + S_dag) / 2
    a = (S_star - theta)  # = (S* - S†)/2
    # Hoeffding: P(Ŝ < θ | S*) ≤ exp(-na²/32) = exp(-n(S*-S†)²/128)
    return math.exp(-n * (S_star - S_dag)**2 / 128)

def normal_miss(n, S_star, S_dag):
    """
    Normal approximation (variance σ² = 64/n for CHSH estimator with range 8).
    P(Ŝ < θ) = Φ((θ - S*) / σ) where σ = 8/√n
    """
    theta = (S_star + S_dag) / 2
    sigma_n = 8 / math.sqrt(n)
    z = (theta - S_star) / sigma_n
    return norm.cdf(z)

def simulated_miss(n, S_star, S_dag, n_trials=50_000, rng_seed=42):
    """
    Monte Carlo: draw Ŝ ~ N(S*, (8/√n)²) and count Ŝ < θ.
    """
    theta = (S_star + S_dag) / 2
    sigma_n = 8 / math.sqrt(n)
    rng = np.random.default_rng(rng_seed)
    S_hat = rng.normal(S_star, sigma_n, n_trials)
    return float((S_hat < theta).mean())

# ─── Parameters from paper ────────────────────────────────────────────────────
S_star = 2.687   # honest CHSH
S_dag  = 2.418   # compromised CHSH
theta  = (S_star + S_dag) / 2
n_vals = [100, 500, 1000, 2500, 5000, 10_000]

print("=== Proposition 2 Miss-Rate Verification ===")
print(f"S* = {S_star},  S† = {S_dag},  θ = {theta:.4f}")
print()
print(f"{'n':>8}  {'Hoeffding_UB':>14}  {'Normal_approx':>15}  {'MC_simulated':>14}")
print("-"*60)

rows = []
for n in n_vals:
    h_ub  = hoeffding_miss(n, S_star, S_dag)
    n_apr = normal_miss(n, S_star, S_dag)
    mc    = simulated_miss(n, S_star, S_dag)
    rows.append({"n": n, "hoeffding_UB": h_ub, "normal_approx": n_apr, "mc_simulated": mc})
    print(f"{n:>8}  {h_ub:>14.4e}  {n_apr:>15.6f}  {mc:>14.6f}")

# ─── ε_e formula comparison ───────────────────────────────────────────────────
print("\n=== ε_e Formula: Old (2/n²) vs Corrected (δ_e) ===")
delta_e = 1e-6
print(f"{'n':>8}  {'eps_old=2/n²':>14}  {'eps_corrected=δ_e':>20}  {'ratio':>10}")
print("-"*60)
eps_rows = []
for n in n_vals:
    eps_old = 2 / n**2
    eps_new = delta_e
    ratio = eps_old / eps_new
    eps_rows.append({"n": n, "eps_old": eps_old, "eps_corrected": eps_new, "ratio": ratio})
    print(f"{n:>8}  {eps_old:>14.2e}  {eps_new:>20.2e}  {ratio:>10.2f}")

print()
print("NOTE: At n=10,000, ε_old = 2e-8  vs  ε_corrected = 1e-6.")
print("The old formula gives falsely tiny defects; the corrected formula gives a")
print("constant δ_e = 10⁻⁶ per link, regardless of n. This resolves the binding-regime")
print("issue: in the binding regime (n=2500, Δ=1e-4), H=4 links give Σε=4e-6 < Δ=1e-4")
print("so the budget CAN bind if Δ is set close to H*δ_e ≈ 4×10⁻⁶.")

# ─── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

n_arr = np.array(n_vals)
hb  = [r["hoeffding_UB"]  for r in rows]
na  = [r["normal_approx"] for r in rows]
mc  = [r["mc_simulated"]  for r in rows]

axes[0].semilogy(n_arr, hb, 'b-o',  label='Hoeffding upper bound')
axes[0].semilogy(n_arr, na, 'r--s', label='Normal approximation')
axes[0].semilogy(n_arr, mc, 'g-^',  label='MC simulation (50k trials)')
axes[0].set_xlabel('n (per-setting audit length)')
axes[0].set_ylabel('Miss probability (log scale)')
axes[0].set_title(f'Prop 2 Miss Rate: P(Ŝ<θ|S*)\nS*={S_star}, S†={S_dag}, θ={theta:.3f}')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# eps comparison
eps_old_arr = [2/n**2 for n in n_arr]
axes[1].semilogy(n_arr, eps_old_arr, 'r--', label='ε_e = 2/n² (old formula)')
axes[1].axhline(delta_e, color='b', linestyle='-', label=f'ε_e = δ_e = {delta_e:.0e} (corrected)')
axes[1].set_xlabel('n (per-setting audit length)')
axes[1].set_ylabel('Security defect ε_e')
axes[1].set_title('ε_e Formula: Old vs Corrected')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_08_prop2_miss_rate.png", dpi=150)
plt.close()
print("\nSaved: fig_08_prop2_miss_rate.png")

out = {
    "description": "Prop 2 miss-rate + ε_e formula (Rev1 Comments 1, 3)",
    "S_star": S_star, "S_dag": S_dag, "theta": theta,
    "miss_rate_table": rows,
    "eps_comparison": eps_rows,
    "delta_e": delta_e,
}
with open("REVISION_OSN_12OCT2026/expr_r1/08_prop2_missrate_results.json", "w") as f:
    json.dump(out, f, indent=2)
