"""
Experiment 3: Werner State / Fidelity Formula Verification
============================================================
Addresses Reviewer 1, Comment 4:
  "The Werner fidelity F_e = (S/(2√2)+1)/2 does not equal the
   usual Werner-state Bell-state fidelity"

For Werner state: ρ = v|Φ⁺⟩⟨Φ⁺| + (1−v)I/4
  - Max CHSH:  S = 2√2 · v
  - Fidelity F = ⟨Φ⁺|ρ|Φ⁺⟩ = (1 + 3v)/4

The paper uses F = (S/(2√2) + 1)/2 = (v + 1)/2  ← NOT the Werner fidelity.
The correct expression is F = (1 + 3v)/4 = (1 + 3S/(2√2·4)) = ...

This experiment:
  1. Computes both expressions over v ∈ [0,1] and verifies the discrepancy
  2. Checks the dephased Bell-diagonal formula S = 2√2(2F−1) ↔ F=(1+S/(2√2))/2
  3. Plots both fidelity functions vs CHSH
  4. Provides the correct formula with numeric examples
"""

import math
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sqrt2 = math.sqrt(2)

def S_from_v(v):
    """Max CHSH for Werner state with visibility v."""
    return 2 * sqrt2 * v

def F_werner(v):
    """True Werner-state Bell-state fidelity F = (1+3v)/4."""
    return (1 + 3*v) / 4

def F_paper(S):
    """Fidelity formula used in paper: F = (S/(2√2) + 1)/2 = (v+1)/2."""
    return (S / (2*sqrt2) + 1) / 2

def F_bell_diagonal(S):
    """
    Dephased Bell-diagonal mixture: F|Φ⁺⟩⟨Φ⁺| + (1-F)|Φ⁻⟩⟨Φ⁻|
    → S = 2√2(2F−1)  ↔  F = (1 + S/(2√2))/2
    (same as paper's formula!)
    """
    return (1 + S / (2*sqrt2)) / 2

# ─── Numeric check ────────────────────────────────────────────────────────────
v_vals = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 1.0]

print("=== Werner State Fidelity Verification ===")
print(f"{'v':>6}  {'S=2√2·v':>10}  {'F_Werner=(1+3v)/4':>20}  {'F_paper=(v+1)/2':>18}  {'F_BellDiag':>12}  {'Match?':>8}")
print("-"*90)

rows = []
for v in v_vals:
    S   = S_from_v(v)
    Fw  = F_werner(v)
    Fp  = F_paper(S)
    Fbd = F_bell_diagonal(S)
    match = abs(Fw - Fp) < 1e-10
    rows.append({"v": v, "S": S, "F_werner": Fw, "F_paper": Fp, "F_bell_diag": Fbd, "match": match})
    print(f"{v:>6.2f}  {S:>10.5f}  {Fw:>20.6f}  {Fp:>18.6f}  {Fbd:>12.6f}  {'YES' if match else 'NO':>8}")

print()
print("Conclusion: F_paper = (v+1)/2  ≠  F_Werner = (1+3v)/4  for v≠1")
print("            F_paper = F_BellDiag  (consistent for dephased Bell-diagonal state)")
print()
print("The paper's formula is correct for a DEPHASED BELL-DIAGONAL mixture,")
print("NOT for the standard Werner state. This should be clarified in the text.")
print()

# Numeric example: S=2.687 (v=2.687/(2√2) ≈ 0.9502)
S_ex = 2.687
v_ex = S_ex / (2*sqrt2)
Fw_ex  = F_werner(v_ex)
Fp_ex  = F_paper(S_ex)
print(f"Example: S=2.687, v={v_ex:.5f}")
print(f"  F_Werner (correct)  = {Fw_ex:.6f}")
print(f"  F_paper  (formula)  = {Fp_ex:.6f}")
print(f"  Discrepancy         = {abs(Fw_ex - Fp_ex):.6f}")
print()

# Path fidelity via Werner swapping
# Swapping: v_path = prod(v_i), F_path = (1 + 3*prod(v_i))/4
H = 4   # hops
v_link = v_ex
v_path = v_link ** H
F_path_werner = F_werner(v_path)
F_path_paper  = F_paper(S_from_v(v_path))
print(f"Path fidelity (H={H} hops, v_link={v_link:.5f}):")
print(f"  v_path = {v_path:.6f}")
print(f"  F_path (Werner)  = {F_path_werner:.6f}")
print(f"  F_path (paper)   = {F_path_paper:.6f}")

# ─── Plot ─────────────────────────────────────────────────────────────────────
v_arr = np.linspace(0, 1, 500)
S_arr = 2*sqrt2*v_arr

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(v_arr, [(1+3*v)/4 for v in v_arr], 'b-',  label='F_Werner = (1+3v)/4')
axes[0].plot(v_arr, [(v+1)/2   for v in v_arr], 'r--', label='F_paper = (v+1)/2')
axes[0].set_xlabel('v (Werner visibility)')
axes[0].set_ylabel('Fidelity F')
axes[0].set_title('Werner vs Paper Fidelity Formula')
axes[0].legend()
axes[0].grid(True, alpha=0.3)
axes[0].set_xlim([0, 1])
axes[0].set_ylim([0, 1.05])

# Difference
diff = np.array([(v+1)/2 - (1+3*v)/4 for v in v_arr])
axes[1].plot(v_arr, diff, 'g-')
axes[1].axhline(0, color='gray', linestyle='--')
axes[1].set_xlabel('v (Werner visibility)')
axes[1].set_ylabel('F_paper − F_Werner')
axes[1].set_title('Fidelity Discrepancy (paper formula − Werner formula)')
axes[1].grid(True, alpha=0.3)
axes[1].set_xlim([0, 1])

plt.tight_layout()
plt.savefig("REVISION_OSN_12OCT2026/expr_r1/fig_03_werner_fidelity.png", dpi=150)
plt.close()
print("\nSaved: fig_03_werner_fidelity.png")

# ─── Save results ─────────────────────────────────────────────────────────────
out = {
    "description": "Werner fidelity formula verification (Reviewer 1, Comment 4)",
    "S0": 2.687, "v_at_S0": v_ex,
    "F_werner_at_S0": Fw_ex, "F_paper_at_S0": Fp_ex,
    "discrepancy": abs(Fw_ex - Fp_ex),
    "conclusion": (
        "Paper uses F=(S/(2√2)+1)/2 which equals the Bell-diagonal dephased-mixture "
        "formula, NOT the standard Werner-state fidelity F=(1+3v)/4. These differ "
        "for v≠1. The paper should clarify which state model is intended."
    ),
    "table": rows,
}
with open("REVISION_OSN_12OCT2026/expr_r1/03_werner_fidelity_results.json", "w") as f:
    json.dump(out, f, indent=2)
