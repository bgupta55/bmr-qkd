"""
Experiment 1: Certificate Pipeline — Reference Values Table
============================================================
Addresses Reviewer 1, Comment 1:
  - Clarify n = per-setting vs total rounds
  - Specify delta_e numerically
  - Reproduce tau, S_lb = S_hat - tau, rₑ = h(S_lb), ε_e, γ_e

Notation (from r1_v1.md glossary):
  n      = per-setting audit length  (4n total rounds)
  S_hat  = observed CHSH (scalar input)
  tau    = Hoeffding radius = sqrt(32*ln(8/delta)/n)
  S_lb   = S_hat - tau  (lower confidence bound on S)
  h(S)   = 1 - log2(1 + sqrt(2 - S^2/4))   [Pironio min-entropy]
  rₑ     = h(S_lb)  [certified min-entropy rate, bits/round]
  R      = minimum accepted rate (routing threshold)
  ε_e    = probability defect = delta (audit failure prob)
  γ_e    = guessing excess = max(0, 0.5 - 2^{-rₑ})  (illustrative)
  θ_R    = h^{-1}(R)  [CHSH threshold for rate R]
  min_n  = smallest n such that S0-tau >= 2  (S0=honest CHSH)
"""

import math
import json
from scipy.optimize import brentq
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─── Core certificate functions ──────────────────────────────────────────────

def h_pironio(S):
    """Min-entropy per round; S in (2, 2*sqrt(2)]."""
    if S <= 2.0:
        return 0.0
    if S > 2 * math.sqrt(2):
        S = 2 * math.sqrt(2)
    return 1.0 - math.log2(1 + math.sqrt(2 - S**2 / 4))

def tau(n, delta):
    """Hoeffding radius; n = per-setting audit length."""
    return math.sqrt(32 * math.log(8 / delta) / n)

def S_lb(S_hat, n, delta):
    return S_hat - tau(n, delta)

def r_e(S_hat, n, delta):
    """Certified min-entropy rate (bits/round)."""
    slb = S_lb(S_hat, n, delta)
    return h_pironio(slb)

def epsilon_e(delta):
    """Security defect = audit-failure probability delta."""
    return delta

def gamma_e(r):
    """Guessing excess (illustrative; 0.5 - 2^{-r} for r>0)."""
    if r <= 0:
        return 0.5
    return max(0.0, 0.5 - 2**(-r))

def min_n_honest_passes(S0, delta):
    """Smallest per-setting n so that S0 - tau(n,delta) >= 2."""
    # tau^2 = 32*ln(8/delta)/n  →  n >= 32*ln(8/delta)/(S0-2)^2
    return math.ceil(32 * math.log(8 / delta) / (S0 - 2) ** 2)

def theta_R(R, eps=1e-10):
    """CHSH threshold for rate R: h^{-1}(R) via brentq."""
    if R <= 0:
        return 2.0
    if R >= 1.0:
        return 2 * math.sqrt(2)
    return brentq(lambda S: h_pironio(S) - R, 2.0 + eps, 2 * math.sqrt(2) - eps)


# ─── Build reference-values table ───────────────────────────────────────────
S0 = 2.687   # honest Werner-state CHSH (v=0.95 → S=2*sqrt(2)*v≈2.687)
S1 = 2.418   # compromised CHSH (midpoint threshold demonstration)
R  = 0.10    # minimum accepted rate (bits/round)

n_vals     = [100, 1000, 1078, 2500, 5000, 10_000]
delta_vals = [1e-2, 1e-4, 1e-6]

rows = []
for n in n_vals:
    for d in delta_vals:
        t   = tau(n, d)
        slb_honest = S0 - t
        r0  = h_pironio(slb_honest) if slb_honest > 2 else 0.0
        slb_comp   = S1 - t
        r1  = h_pironio(slb_comp)   if slb_comp > 2 else 0.0
        eps = epsilon_e(d)
        gam = gamma_e(r0)
        rows.append({
            "n_per_setting": n,
            "n_total_rounds": 4*n,
            "delta_e": d,
            "tau": round(t, 5),
            "S_lb_honest": round(slb_honest, 4),
            "r_e_honest": round(r0, 5),
            "S_lb_comp": round(slb_comp, 4),
            "r_e_comp": round(r1, 5),
            "epsilon_e": eps,
            "gamma_e_honest": round(gam, 5),
            "passes_threshold": slb_honest > 2.0,
            "comp_passes": slb_comp > 2.0,
        })

# ─── Min-n table ────────────────────────────────────────────────────────────
min_n_rows = []
for d in [1e-2, 1e-4, 1e-6]:
    mn = min_n_honest_passes(S0, d)
    min_n_rows.append({"delta_e": d, "min_n_per_setting": mn, "min_total_rounds": 4*mn})

# theta_R
thr = theta_R(R)

# ─── Save JSON ───────────────────────────────────────────────────────────────
out = {
    "description": "Certificate pipeline reference values (Reviewer 1, Comment 1)",
    "S0_honest": S0,
    "S1_compromised": S1,
    "R_min_accepted_rate": R,
    "theta_R": round(thr, 5),
    "min_n_table": min_n_rows,
    "reference_table": rows,
}
with open("REVISION_OSN_12OCT2026/expr_r1/01_cert_pipeline_results.json", "w") as f:
    json.dump(out, f, indent=2)

# ─── Print summary ────────────────────────────────────────────────────────────
print("=== Certificate Pipeline Reference Values ===")
print(f"S0 (honest CHSH)     = {S0}")
print(f"S1 (compromised CHSH)= {S1}")
print(f"R  (min accepted)    = {R}")
print(f"θ_R = h⁻¹(R)         = {thr:.5f}")
print()
print("Min-n for honest link to pass (S0=2.687):")
for r in min_n_rows:
    print(f"  δ={r['delta_e']:.0e}: n_per_setting={r['min_n_per_setting']:,}  (4n={r['min_total_rounds']:,} total rounds)")

print()
print(f"{'n':>8} {'δ':>8} {'τ':>8} {'S_lb(H)':>9} {'r_e(H)':>9} {'S_lb(C)':>9} {'r_e(C)':>9} {'passes?':>8}")
print("-"*80)
for row in rows:
    print(f"{row['n_per_setting']:>8} {row['delta_e']:>8.0e} {row['tau']:>8.4f} "
          f"{row['S_lb_honest']:>9.4f} {row['r_e_honest']:>9.5f} "
          f"{row['S_lb_comp']:>9.4f} {row['r_e_comp']:>9.5f} "
          f"{'Y' if row['passes_threshold'] else 'N':>8}")

# ─── Plot: τ vs n for various δ ──────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

n_arr = np.arange(50, 12000, 50)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
for i, d in enumerate([1e-2, 1e-4, 1e-6]):
    t_arr = np.array([tau(n, d) for n in n_arr])
    axes[0].plot(n_arr, t_arr, label=f'δ=10^{int(math.log10(d))}', color=colors[i])
axes[0].axhline(y=S0-2, color='gray', linestyle='--', label=f'S₀−2={S0-2:.3f}')
axes[0].set_xlabel('n (per-setting audit length)')
axes[0].set_ylabel('τ (Hoeffding radius)')
axes[0].set_title('Hoeffding Radius τ vs n')
axes[0].legend()
axes[0].set_xlim([0, 12000])
axes[0].grid(True, alpha=0.3)

for i, d in enumerate([1e-2, 1e-4, 1e-6]):
    r_arr = np.array([h_pironio(S0 - tau(n, d)) if S0 - tau(n, d) > 2 else 0.0 for n in n_arr])
    axes[1].plot(n_arr, r_arr, label=f'δ=10^{int(math.log10(d))}', color=colors[i])
axes[1].axhline(y=R, color='red', linestyle='--', label=f'R={R}')
axes[1].set_xlabel('n (per-setting audit length)')
axes[1].set_ylabel('r_e (certified min-entropy rate, bits/round)')
axes[1].set_title('Certified Rate r_e vs n (S₀=2.687, honest link)')
axes[1].legend()
axes[1].set_xlim([0, 12000])
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_01_cert_pipeline.png", dpi=150)
plt.close()
print("\nSaved: fig_01_cert_pipeline.png")
