"""
X2: Detection Power / ROC Table (replaces defective E8)
=========================================================
Addresses r1_v2.md §2 E8 — E8 was mathematically inconsistent:
  - Used σ = 8/√n which is 4× the physical maximum (correct: 2/√n)
  - Hoeffding UB fell BELOW the true probability (impossible)
  - Used midpoint test; correct operative test is calibrated alarm θ_α(n)

This script implements the CORRECT design per D2:
  (i)  Gate:    S_lb = Ŝ − τ ≥ 2   (certification)
  (ii) Request: r_e ≥ R = 0.10      (application floor)
  (iii) Alarm:  Ŝ < θ_α(n) = S₀,cal − z_α·σ(n),  α = 10⁻⁶  (link becomes non-routable)

Physical σ: CHSH estimator Ŝ = (1/n)Σy_i with y_i ∈ [-4,4] for n per-setting.
For Werner state v=S₀/(2√2), each term has Var(y) = 4(1-v²)·corrector ≈ 4(1-E²).
Conservative (maximum) σ = 2/√n  (verified in r1_v2.md §2 E8-a).

Per D5: use binomial sampling, σ_physical(S₀,n) = 2/√n (conservative upper bound).

ROC table as specified in r1_v2.md §2 E8, reference values table.
"""

import math
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm

# ─── Core functions ──────────────────────────────────────────────────────────
S0      = 2.687    # honest CHSH
delta_e = 1e-6     # δ_e per D1
R       = 0.10     # request floor per D2
sqrt2   = math.sqrt(2)

def h_pironio(S):
    """Min-entropy rate h(S_lb) = 1 - log2(1 + √(2 - S²/4))."""
    if S <= 2.0: return 0.0
    S = min(S, 2*sqrt2 - 1e-10)
    inner = max(0.0, 2.0 - S**2/4.0)
    return 1.0 - math.log2(1.0 + math.sqrt(inner))

def tau_paper(n, delta=delta_e):
    """Paper's τ = √(32·ln(8/δ)/n)  — kept unchanged per D11."""
    return math.sqrt(32 * math.log(8/delta) / n)

def sigma_conservative(n):
    """Conservative (maximum possible) CHSH std dev: σ = 2/√n per r1_v2.md E8-a.
    Valid for ANY quantum state; used as the Hoeffding bound argument."""
    return 2.0 / math.sqrt(n)

def sigma_werner(n, S0_val=S0):
    """Physical CHSH std dev for a Werner state with honest mean S₀.
    For a Werner state with visibility v₀ = S₀/(2√2), each CHSH correlator
    E(a,b) = ±v₀/√2.  Var(y) per round = 4·(1 − E²) = 4·(1 − v₀²/2).
    σ_Werner = √(4·(1 − v₀²/2) / n)  ≈ 1.48/√n at S₀=2.687.
    This matches r1_v2.md §2 E8 reference table exactly.
    GAP 1 FIX: use this for the 'r1_v2 reference' column; keep conservative for bound."""
    v0 = S0_val / (2 * math.sqrt(2))
    E_val = v0 / math.sqrt(2)          # per-correlator expectation
    return math.sqrt(4.0 * (1.0 - E_val**2) / n)

# Alias used throughout (conservative = default; Werner used in reference column)
def sigma_physical(n):
    return sigma_conservative(n)

def theta_alarm(n, alpha=1e-6, S0_cal=S0, use_werner=False):
    """
    Calibrated alarm threshold: θ_α(n) = S₀,cal − z_α·σ(n).
    use_werner=True  → σ_Werner (matches r1_v2.md reference table exactly).
    use_werner=False → σ_conservative=2/√n (valid upper-bound / looser alarm).
    """
    z_alpha = norm.ppf(1.0 - alpha)
    sig = sigma_werner(n, S0_cal) if use_werner else sigma_conservative(n)
    return S0_cal - z_alpha * sig

def detection_power(n, S_attack, alpha=1e-6, S0_cal=S0, use_werner=False):
    """
    Per-epoch power = P(Ŝ < θ_α(n) | true_S = S_attack).
    Ŝ ~ N(S_attack, σ²(n)).
    use_werner mirrors theta_alarm: determines which σ calibrates the threshold.
    The width of Ŝ under the attack is also set to the same σ for consistency.
    """
    theta = theta_alarm(n, alpha, S0_cal, use_werner=use_werner)
    sig   = sigma_werner(n, S0_cal) if use_werner else sigma_conservative(n)
    z = (theta - S_attack) / sig
    return float(norm.cdf(z))

def hoeffding_ub_miss(n, S_star, S_dag):
    """
    Hoeffding UB for midpoint-test miss probability P(Ŝ > θ | S† = S_dag).
    Per r1_v2.md B.2: P ≤ exp(-n(S*-S†)²/128).
    NOTE: This is valid ONLY when σ = 2/√n (range of y is 8 → range² = 64 → /2 = 32).
    Correct Hoeffding: P ≤ exp(-n·a²/32) where a = (S*-S†)/2, range = 8.
    """
    a = (S_star - S_dag) / 2.0
    return math.exp(-n * a**2 / 32.0)

def mc_miss_rate(n, S_star, S_dag, n_trials=100_000, seed=42):
    """
    MC simulation using CORRECT σ = 2/√n.
    Miss = P(Ŝ > θ | S† where θ = (S*+S†)/2 midpoint).
    """
    theta = (S_star + S_dag) / 2.0
    sigma = sigma_physical(n)
    rng = np.random.default_rng(seed)
    S_hat = rng.normal(S_dag, sigma, n_trials)
    return float((S_hat > theta).mean())

# ─── Operating points (D1) ────────────────────────────────────────────────────
n_vals = [2500, 5000, 10_000, 20_000, 40_000]

# Attack depths (S†) 
attack_depths = {
    "subtle (S†=2.62)":   2.62,
    "moderate (S†=2.55)": 2.55,
    "strong (S†=2.48)":   2.48,
    "declared (S†=2.418)":2.418,
}

z_alpha_val = norm.ppf(1 - 1e-6)
v0          = S0 / (2 * math.sqrt(2))
sig_werner_example = sigma_werner(10000)
sig_cons_example   = sigma_conservative(10000)

print("=" * 110)
print("X2: Detection Power / ROC Table  (corrected per r1_v2.md D2, D5, E8) — GAP 1 FIXED")
print("=" * 110)
print(f"\nTwo σ columns per r1_v2.md:")
print(f"  σ_conservative = 2/√n               (E8-a upper bound; valid for any state)")
print(f"  σ_Werner(S₀)   = √(4(1−E²)/n)  where E=v₀/√2  (honest-link physical; v₀=S₀/(2√2)={v0:.4f})")
print(f"  At n=10,000: σ_cons={sig_cons_example:.5f}, σ_Werner={sig_werner_example:.5f}  (ratio={sig_cons_example/sig_werner_example:.2f}×, matches r1_v2 ref)")
print(f"  z_{{10⁻⁶}} = {z_alpha_val:.4f}")
print()

# ─── Main ROC table (BOTH σ columns) ─────────────────────────────────────────
# r1_v2.md reference uses σ_Werner; we show both so referee can verify
print("CONSERVATIVE σ=2/√n  (upper bound; looser alarm):")
print(f"{'n':>8}  {'τ':>7}  {'σ_cons':>8}  {'θ_α':>8}  {'R_α':>8}", end="")
for label in attack_depths:
    short = label.split('=')[1].rstrip(')')
    print(f"  {short:>8}", end="")
print()
print("-" * 100)

roc_rows = []
for n in n_vals:
    t         = tau_paper(n)
    sig_c     = sigma_conservative(n)
    theta_c   = theta_alarm(n, use_werner=False)
    r_alpha_c = h_pironio(theta_c - t)
    powers_c  = {label: detection_power(n, S_atk, use_werner=False)
                 for label, S_atk in attack_depths.items()}
    print(f"{n:>8}  {t:>7.4f}  {sig_c:>8.5f}  {theta_c:>8.4f}  {r_alpha_c:>8.4f}", end="")
    for p in powers_c.values():
        print(f"  {p:>8.4f}", end="")
    print()
    roc_rows.append({"n": n, "tau": round(t,5),
                     "sigma_conservative": round(sig_c,6),
                     "theta_alarm_conservative": round(theta_c,5),
                     "R_alpha_conservative": round(r_alpha_c,5),
                     "powers_conservative": {lbl: round(p,6) for lbl,p in powers_c.items()}})

print()
print("WERNER σ=√(4(1−v₀²)/n)  (honest-link physical value; matches r1_v2 reference table):")
print(f"{'n':>8}  {'τ':>7}  {'σ_Werner':>9}  {'θ_α':>8}  {'R_α':>8}", end="")
for label in attack_depths:
    short = label.split('=')[1].rstrip(')')
    print(f"  {short:>8}", end="")
print()
print("-" * 100)

# r1_v2.md reference values for verification
ref_vals = {
    2500:  {"theta":2.546, "R_alpha":0.036,  "p0":0.007, "p1":0.45,  "p2":0.98, "p3":1.0},
    5000:  {"theta":2.587, "R_alpha":0.117,  "p0":0.064, "p1":0.96,  "p2":1.0,  "p3":1.0},
    10000: {"theta":2.617, "R_alpha":0.188,  "p0":0.41,  "p1":1.0,   "p2":1.0,  "p3":1.0},
    20000: {"theta":2.637, "R_alpha":0.250,  "p0":0.95,  "p1":1.0,   "p2":1.0,  "p3":1.0},
    40000: {"theta":2.652, "R_alpha":0.301,  "p0":1.0,   "p1":1.0,   "p2":1.0,  "p3":1.0},
}

for row in roc_rows:
    n         = row["n"]
    t         = row["tau"]
    sig_w     = sigma_werner(n)
    theta_w   = theta_alarm(n, use_werner=True)
    r_alpha_w = h_pironio(theta_w - t)
    powers_w  = {label: detection_power(n, S_atk, use_werner=True)
                 for label, S_atk in attack_depths.items()}
    ref       = ref_vals[n]
    theta_match = "✅" if abs(theta_w - ref["theta"]) < 0.003 else "⚠"
    p0_match    = "✅" if abs(list(powers_w.values())[0] - ref["p0"]) < 0.05 else "⚠"
    print(f"{n:>8}  {t:>7.4f}  {sig_w:>9.5f}  {theta_w:>8.4f}{theta_match}  {r_alpha_w:>8.4f}", end="")
    for p in powers_w.values():
        print(f"  {p:>8.4f}", end="")
    print(f"  {p0_match}")
    row.update({
        "sigma_werner": round(sig_w,6),
        "theta_alarm_werner": round(theta_w,5),
        "R_alpha_werner": round(r_alpha_w,5),
        "powers_werner": {lbl: round(p,6) for lbl,p in powers_w.items()},
        "ref_theta": ref["theta"], "ref_p0_subtle": ref["p0"]
    })

print()
print("NOTE: 'Conservative' column uses σ=2/√n (max possible) → wider alarm band → lower power.")
print("      'Werner' column uses σ≈1.48/√n (honest-state physical) → tighter alarm → higher power.")
print("      Both are correct; conservative is a valid bound; Werner matches r1_v2 reference.")

# ─── Mean detection latency (Werner σ — matches reference) ───────────────────
print("\n--- Mean detection latency (Werner σ; matches r1_v2 reference) ---")
print(f"{'n':>8}  {'Attack':>22}  {'Power(cons)':>13}  {'Power(Werner)':>14}  {'Latency(W,ep)':>15}  {'P35(Werner)':>14}")
print("-"*100)
latency_rows = []
for n in [5000, 10_000, 20_000]:
    for label, S_atk in attack_depths.items():
        p_c = detection_power(n, S_atk, use_werner=False)
        p_w = detection_power(n, S_atk, use_werner=True)
        mean_lat = 1.0/p_w if p_w > 1e-10 else float('inf')
        p35      = 1.0 - (1.0 - p_w)**35
        print(f"{n:>8}  {label:>22}  {p_c:>13.4f}  {p_w:>14.4f}  {mean_lat:>15.1f}  {p35:>14.4f}")
        latency_rows.append({"n": n, "attack": label, "S_atk": S_atk,
                              "power_conservative": round(p_c,6),
                              "power_werner": round(p_w,6),
                              "mean_latency_werner_epochs": round(mean_lat,2),
                              "p_detect_35_epochs_werner": round(p35,6)})

# ─── Midpoint-test correction from E8 ────────────────────────────────────────
print("\n--- Midpoint test (corrected; E8 comparison) ---")
print(f"{'n':>8}  {'Hoeffding UB':>14}  {'Normal approx':>15}  {'MC simulated':>14}  {'E8 wrong MC':>14}")
S_star, S_dag = 2.687, 2.418
print(f"S* = {S_star},  S† = {S_dag},  θ = {(S_star+S_dag)/2:.4f}")
print("-"*72)
miss_rows = []
for n in [100, 500, 1000, 2500, 5000, 10_000]:
    h_ub  = hoeffding_ub_miss(n, S_star, S_dag)
    sigma_correct = sigma_physical(n)
    sigma_wrong   = 8.0 / math.sqrt(n)   # E8's wrong value
    theta_mid = (S_star + S_dag) / 2.0
    # Normal approx with correct σ
    z_correct = (S_dag - theta_mid) / sigma_correct
    # z < 0 means S† < θ, so miss = P(Ŝ > θ | S†) = 1 - Φ(z)
    na_correct = 1.0 - norm.cdf((theta_mid - S_dag)/sigma_correct)
    # MC with correct σ
    mc_c = mc_miss_rate(n, S_star, S_dag)
    # E8's wrong approximation
    na_wrong = 1.0 - norm.cdf((theta_mid - S_dag)/sigma_wrong)
    
    print(f"{n:>8}  {h_ub:>14.4e}  {na_correct:>15.6f}  {mc_c:>14.6f}  {na_wrong:>14.6f}")
    miss_rows.append({"n": n, "hoeffding_UB": h_ub, "normal_correct": na_correct,
                      "mc_correct": mc_c, "e8_wrong": na_wrong})

print()
print("KEY FIX: E8 used σ=8/√n giving 'Normal approx' in column 'E8 wrong MC'.")
print("Correct σ=2/√n gives MUCH lower miss probability at large n.")
print("At n=10,000: correct MC ≈ 0, not 4.7% as E8 reported.")
print("The Hoeffding UB (3.5e-3) is now a valid UPPER bound on the correct value.")

# ─── Figures ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))

# Fig A: Detection power vs n — conservative AND Werner, for subtle attack only + declared
colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728']
n_plot_vals = [1000, 2500, 5000, 10000, 20000, 40000]
for i, (label, S_atk) in enumerate(attack_depths.items()):
    # Conservative (dashed) and Werner (solid)
    p_cons  = [detection_power(n, S_atk, use_werner=False) for n in n_plot_vals]
    p_wern  = [detection_power(n, S_atk, use_werner=True)  for n in n_plot_vals]
    axes[0].semilogy(n_plot_vals, p_wern, '-o', color=colors[i], label=f'{label} (Werner σ)', linewidth=2)
    axes[0].semilogy(n_plot_vals, p_cons, '--o', color=colors[i], alpha=0.4, label=f'{label} (cons σ)', linewidth=1)
axes[0].axhline(0.95, color='gray', linestyle=':', lw=1, label='95% power')
axes[0].set_xlabel('n (per-setting audit length)')
axes[0].set_ylabel('Per-epoch detection power (log)')
axes[0].set_title('Alarm Power: Werner σ (solid)\nvs Conservative σ (dashed) — GAP 1 FIXED')
axes[0].legend(fontsize=7, ncol=2)
axes[0].grid(True, alpha=0.3)

# Fig B: Detection latency (Werner σ — matches r1_v2 reference)
for i, (label, S_atk) in enumerate(attack_depths.items()):
    lat = [1.0/max(detection_power(n, S_atk, use_werner=True), 1e-10) for n in n_plot_vals]
    lat = [min(l, 1e4) for l in lat]
    axes[1].semilogy(n_plot_vals, lat, '-o', color=colors[i], label=label, linewidth=2)
axes[1].axhline(35, color='red', linestyle=':', lw=1.5, label='35-epoch window')
axes[1].set_xlabel('n (per-setting audit length)')
axes[1].set_ylabel('Mean detection latency (epochs, log)')
axes[1].set_title('Mean Detection Latency\n(Werner σ; matches r1_v2 reference)')
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.3)

# Fig C: Both alarm thresholds θ_α(n) vs n
n_plot = np.linspace(1000, 50000, 200)
theta_cons_arr  = [theta_alarm(n, use_werner=False) for n in n_plot]
theta_wern_arr  = [theta_alarm(n, use_werner=True)  for n in n_plot]
r_cons_arr = [h_pironio(theta_alarm(n, use_werner=False) - tau_paper(n)) for n in n_plot]
r_wern_arr = [h_pironio(theta_alarm(n, use_werner=True)  - tau_paper(n)) for n in n_plot]
ax2 = axes[2]
ax2.plot(n_plot, theta_wern_arr, 'b-',  label='θ_α Werner σ (matches ref)', linewidth=2)
ax2.plot(n_plot, theta_cons_arr, 'b--', label='θ_α conservative σ',         linewidth=1.5, alpha=0.6)
ax2.axhline(S0, color='blue', linestyle=':', lw=1, alpha=0.4, label=f'S₀={S0}')
ax2.set_ylabel('Alarm threshold θ_α(n)', color='blue')
ax2.set_xlabel('n (per-setting audit length)')
ax2.set_title('Alarm Threshold Comparison\nWerner vs Conservative σ')
ax3 = ax2.twinx()
ax3.plot(n_plot, r_wern_arr, 'r-',  label='r_e Werner (ref)', linewidth=2)
ax3.plot(n_plot, r_cons_arr, 'r--', label='r_e conservative', linewidth=1.5, alpha=0.6)
ax3.axhline(R, color='red', linestyle=':', lw=1, label=f'R={R}')
ax3.set_ylabel('Certified rate r_e (bits/round)', color='red')
lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax3.get_legend_handles_labels()
ax2.legend(lines1+lines2, labels1+labels2, fontsize=8, loc='lower right')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_X2_detection_power_roc.png", dpi=150)
plt.close()
print("\nSaved: fig_X2_detection_power_roc.png")

out = {
    "description": "Detection power / ROC (corrected E8 + GAP1 fix; r1_v2 X2)",
    "S0": S0, "delta_e": delta_e, "R": R,
    "sigma_model": {
        "conservative": "2/sqrt(n)  — valid upper bound for any quantum state (Hoeffding E8-a)",
        "werner": "sqrt(4*(1-v0^2)/n)  — physical honest-state value; matches r1_v2 reference table"
    },
    "gap1_status": "FIXED: roc_table now has both sigma columns; Werner column matches r1_v2 ref",
    "roc_table": roc_rows,
    "latency_table": latency_rows,
    "midpoint_miss_corrected": miss_rows
}
with open("REVISION_OSN_12OCT2026/expr_r1/X2_detection_power_results.json","w") as f:
    json.dump(out, f, indent=2)
