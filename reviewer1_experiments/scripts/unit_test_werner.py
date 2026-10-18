"""
unit_test_werner.py — GAP 3 fix (r1_v2.md D8, E3-c)
======================================================
Standalone 4×4 density-matrix unit test asserting:
  (i)  S(ρ_v) = 2√2·v   (Werner state CHSH value)
  (ii) F(ρ_v) = (1+3v)/4  (Werner state fidelity — CORRECT formula per D8)
        NOT  F = (v+1)/2  (old wrong formula)

The Werner state is:
  ρ_v = (1−v)/4 · I₄  +  v · |Φ⁺⟩⟨Φ⁺|
where |Φ⁺⟩ = (|00⟩ + |11⟩)/√2.

Verified analytically in r1_v2.md §2 E3-b (fidelity) and standard QI (CHSH).

Run:  python3 unit_test_werner.py
All tests must pass (PASS / FAIL printed per assertion).
"""

import math
import numpy as np

# ─── Build Werner state 4×4 density matrix ───────────────────────────────────
def phi_plus():
    """Bell state |Φ⁺⟩⟨Φ⁺| as 4×4 matrix in {|00⟩,|01⟩,|10⟩,|11⟩} basis."""
    # |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
    phi = np.zeros(4, dtype=complex)
    phi[0] = 1/math.sqrt(2)   # |00⟩
    phi[3] = 1/math.sqrt(2)   # |11⟩
    return np.outer(phi, phi.conj())

def werner_state(v):
    """
    ρ_v = (1−v)/4 · I₄ + v · |Φ⁺⟩⟨Φ⁺|
    v ∈ [0,1]; v=0 → maximally mixed; v=1 → pure Bell state.
    """
    I4 = np.eye(4, dtype=complex)
    return (1 - v) / 4.0 * I4 + v * phi_plus()

# ─── CHSH expectation from density matrix ─────────────────────────────────────
# CHSH observable B = a·σ⊗b·σ + a·σ⊗b'·σ + a'·σ⊗b·σ − a'·σ⊗b'·σ
# Optimal settings for Werner state give S = 2√2·v.
# We compute via the CHSH operator at optimal angles:
#   Alice: a = σ_z, a' = σ_x
#   Bob:   b = (σ_z + σ_x)/√2, b' = (σ_z − σ_x)/√2

def pauli(name):
    """Pauli matrices."""
    if name == 'x':
        return np.array([[0,1],[1,0]], dtype=complex)
    if name == 'z':
        return np.array([[1,0],[0,-1]], dtype=complex)
    if name == 'y':
        return np.array([[0,-1j],[1j,0]], dtype=complex)
    if name == 'I':
        return np.eye(2, dtype=complex)

def kron(A, B):
    return np.kron(A, B)

def chsh_value(rho, v_theoretical=None):
    """
    Compute ⟨B_CHSH⟩ = Tr(B · ρ) at optimal settings.
    B = a⊗b + a⊗b' + a'⊗b − a'⊗b'
    Alice: a=σ_z, a'=σ_x
    Bob:   b=(σ_z+σ_x)/√2, b'=(σ_z−σ_x)/√2
    """
    sx = pauli('x')
    sz = pauli('z')
    b_vec  = (sz + sx) / math.sqrt(2)
    bp_vec = (sz - sx) / math.sqrt(2)
    
    # Four correlator operators
    E_ab  = kron(sz, b_vec)
    E_abp = kron(sz, bp_vec)
    E_apb = kron(sx, b_vec)
    E_apbp= kron(sx, bp_vec)
    
    # CHSH operator: a⊗b + a⊗b' + a'⊗b − a'⊗b'
    B_op = E_ab + E_abp + E_apb - E_apbp
    
    S = float(np.real(np.trace(B_op @ rho)))
    return S

def fidelity_werner(rho):
    """
    F = ⟨Φ⁺|ρ|Φ⁺⟩  = Tr(|Φ⁺⟩⟨Φ⁺| · ρ)
    For Werner state: = v + (1−v)/4
    Analytically: F = (1+3v)/4  ✅ (NOT (v+1)/2 ← old wrong formula).
    """
    phi_p = phi_plus()
    return float(np.real(np.trace(phi_p @ rho)))

# ─── Helper: pass/fail ────────────────────────────────────────────────────────
n_pass = 0
n_fail = 0

def check(name, got, expected, tol=1e-6):
    global n_pass, n_fail
    ok = abs(got - expected) < tol
    status = "PASS" if ok else "FAIL"
    if ok:
        n_pass += 1
    else:
        n_fail += 1
    print(f"  [{status}] {name}: got={got:.8f}, expected={expected:.8f} (tol={tol:.0e})")

# ─── Run tests ────────────────────────────────────────────────────────────────
print("=" * 70)
print("Werner State Unit Tests (r1_v2.md D8, E3-c)")
print("=" * 70)

test_visibilities = [0.0, 0.25, 0.50, 0.70, 0.80, 0.90, 0.95, 1.0]

print("\n--- Density matrix validity ---")
for v in test_visibilities:
    rho = werner_state(v)
    # (1) Trace = 1
    tr = float(np.real(np.trace(rho)))
    check(f"Tr(ρ_v={v}) = 1", tr, 1.0)
    # (2) Hermitian
    herm_err = float(np.max(np.abs(rho - rho.conj().T)))
    check(f"Hermitian_err(v={v}) = 0", herm_err, 0.0)
    # (3) PSD: min eigenvalue ≥ 0  (check that min_eig ≥ -tol, not that it equals 0)
    eigs = np.linalg.eigvalsh(rho)
    min_eig = float(np.min(eigs))
    check(f"min_eigval(v={v}) ≥ 0", min(min_eig, 0.0), 0.0)  # passes iff min_eig ≥ 0

print("\n--- CHSH value: S(ρ_v) = 2√2·v ---")
for v in test_visibilities:
    rho = werner_state(v)
    S_got = chsh_value(rho)
    S_expected = 2.0 * math.sqrt(2) * v
    check(f"S(v={v})", S_got, S_expected, tol=1e-6)

print("\n--- Fidelity: F(ρ_v) = (1+3v)/4  [CORRECT Werner formula per D8] ---")
for v in test_visibilities:
    rho = werner_state(v)
    F_got = fidelity_werner(rho)
    F_correct   = (1 + 3*v) / 4.0      # D8 correct formula
    F_wrong     = (v + 1) / 2.0         # old wrong formula
    check(f"F_correct(v={v})", F_got, F_correct, tol=1e-8)
    # Wrong formula should NOT match (except at v=1 both give 1)
    if v not in (0.0, 1.0):
        mismatch = abs(F_got - F_wrong)
        diff_status = "PASS" if mismatch > 1e-6 else "FAIL"
        print(f"  [{diff_status}] F_wrong(v={v})={(F_wrong):.6f} ≠ F_correct={F_correct:.6f} (diff={mismatch:.6f})")
        if diff_status == "PASS": n_pass += 1
        else: n_fail += 1

print("\n--- Specific numeric assertions from r1_v2.md ---")
# From r1_v2.md §2 E3-b, verified:
# At v=0.95: F_correct=0.9625, F_wrong=0.975 (error=0.0125)
rho_095 = werner_state(0.95)
F_095 = fidelity_werner(rho_095)
S_095 = chsh_value(rho_095)
check("F(v=0.95) = 0.9625", F_095, 0.9625, tol=1e-8)
check("S(v=0.95) = 2√2·0.95", S_095, 2*math.sqrt(2)*0.95, tol=1e-6)

# At v = S0/(2√2) where S0=2.687:
S0 = 2.687
v0 = S0 / (2*math.sqrt(2))
rho_v0 = werner_state(v0)
S_v0   = chsh_value(rho_v0)
F_v0   = fidelity_werner(rho_v0)
F_v0_wrong = (v0 + 1) / 2.0
check(f"S(v=v₀={v0:.4f}) = S₀={S0}", S_v0, S0, tol=1e-5)
check(f"F_correct(v₀) = (1+3·{v0:.4f})/4", F_v0, (1+3*v0)/4, tol=1e-8)
print(f"  [INFO] F_wrong(v₀) = (v₀+1)/2 = {F_v0_wrong:.6f}  vs F_correct = {F_v0:.6f}  (error={F_v0_wrong-F_v0:.6f})")

# 4-hop path fidelity comparison (D8: v_path=∏v_hop, F_path=(1+3·v_path)/4):
print("\n--- 4-hop path fidelity: wrong vs correct Werner formula ---")
print("  (D8: F_correct = (1+3·v_path)/4  where v_path=∏v_hop; NOT product of per-hop F)")
print(f"  {'S/hop':>8}  {'F_wrong×4':>14}  {'F_correct×4':>14}  {'Ratio W/C':>12}  {'Match r1_v2':>12}")
s_vals = [2.74, 2.687, 2.64, 2.55]
ref_wrong  = [0.93892, 0.90368, 0.87327, 0.81719]
ref_correct= [0.91052, 0.86087, 0.81924, 0.74550]
sqrt2 = math.sqrt(2)
for s_val, rw, rc in zip(s_vals, ref_wrong, ref_correct):
    v_hop = s_val / (2*sqrt2)
    v_path = v_hop**4                        # product of per-hop visibilities
    F_path_correct = (1 + 3*v_path) / 4.0   # D8 correct formula applied to path
    F_hop_wrong    = (v_hop + 1) / 2.0
    F_path_wrong   = F_hop_wrong**4          # old wrong formula chained
    r_match = "✅" if (abs(F_path_wrong - rw) < 1e-4 and abs(F_path_correct - rc) < 1e-4) else "⚠️"
    print(f"  {s_val:>8.3f}  {F_path_wrong:>14.5f}  {F_path_correct:>14.5f}  {F_path_wrong/F_path_correct:>12.3f}  {r_match}")
    check(f"F_correct_4hop(S={s_val})", F_path_correct, rc, tol=5e-4)
    check(f"F_wrong_4hop(S={s_val})", F_path_wrong, rw, tol=5e-4)

# ─── Summary ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print(f"SUMMARY: {n_pass} PASS, {n_fail} FAIL")
if n_fail == 0:
    print("✅ ALL TESTS PASSED — Werner formulas verified for r1_v2.md D8")
else:
    print(f"❌ {n_fail} test(s) FAILED")
print("=" * 70)
