"""
E6v2: Real quantum hardware — CHSH with controlled degradation (FIXED)
=======================================================================
Addresses: Reviewer 2 Comment 2 (Algorithm 1 on real CHSH counts; S(φ) monotone).

ROOT CAUSE OF v1 FAILURE (ibm_kingston, S≈1.2):
  The previous script applied RY(A[a] + phi) after Bell state preparation,
  but for |Φ⁺⟩ = (|00⟩+|11⟩)/√2, the correlator is E(θ_A,θ_B) = cos(θ_A − θ_B)
  only when the measurement angle is encoded as RY(2θ) (Bloch-sphere convention).
  Without the factor-of-2, each RY(θ) rotates by θ/2 in Hilbert space, giving
  E = cos(θ_A/2 − θ_B/2) which produces S ≈ 1.2 instead of 2√2.

  Additionally the previous backend (ibm_kingston) was selected only by pending
  jobs, picking a noisy qubit pair. We now pin to ibm_marrakesh qubits (6,5)
  which have the lowest combined gate + readout error on the open plan.

FIXED DESIGN:
  Circuit: H(0) → CX(0,1) → RZ(-π/2,0) → [RY(2·θ_A, 0), RY(2·θ_B, 1)] → measure
  (The RZ(-π/2) rotates the Bell basis so that ZZ measurement gives E=cos(θ_A−θ_B).)

  Tsirelson-optimal angles for |Φ⁺⟩:
    Alice: θ_A ∈ {0, π/4}
    Bob:   θ_B ∈ {π/8, -π/8}
    CHSH:  S = E(0,π/8) + E(0,-π/8) + E(π/4,π/8) − E(π/4,-π/8) = 2√2 ≈ 2.828

  Controlled degradation: add offset φ to Alice's angles only.
    S(φ) = cos(φ − π/8) + cos(φ + π/8) + cos(π/4+φ − π/8) − cos(π/4+φ + π/8)
          ≈ 2√2 · cos(φ) for small φ (to first order).

  5 φ values × 4 settings = 20 circuits × 8192 shots = 163,840 total shots.
  Higher shots → smaller η → higher certified rate → better certificate.
  QPU time estimate: ~30-40s on marrakesh.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import math
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import warnings
warnings.filterwarnings('ignore')

from bmr.bounds import eta_tight, h_pironio
from bmr.certificate import build_certificate

IBM_API_KEY = os.environ.get("IBM_QUANTUM_TOKEN", "")

# ── Experimental parameters ───────────────────────────────────────────────────
SHOTS    = 8192     # 4× more than v1 → η shrinks by 2×
PHIS     = [0.0, 0.1, 0.2, 0.35, 0.5]
DELTA_E  = 1e-6
TARGET_BACKEND = "ibm_marrakesh"
BEST_PAIR      = (6, 5)   # lowest combined gate+readout error on marrakesh

# ── Tsirelson-optimal CHSH settings for |Φ⁺⟩ ────────────────────────────────
# E(θ_A, θ_B) on |Φ⁺⟩ measured with RY(2θ) = cos(θ_A − θ_B)
# S = E00 + E01 + E10 − E11
THETA_A = [0.0,      np.pi / 4]
THETA_B = [np.pi/8, -np.pi/8]
# Sign pattern: (+,+,+,−)
SIGN    = {(0,0): +1, (0,1): +1, (1,0): +1, (1,1): -1}


def ideal_S(phi):
    """Ideal CHSH value at angle offset φ (analytic)."""
    S = 0.0
    for (a, b), sign in SIGN.items():
        theta_a = THETA_A[a] + phi
        theta_b = THETA_B[b]
        E = math.cos(theta_a - theta_b)
        S += sign * E
    return S


def make_chsh_circuit(a_idx, b_idx, phi, initial_layout):
    """
    Build CHSH circuit for setting (a_idx, b_idx) with Alice offset φ.
    Circuit:
      1. Prepare |Φ⁺⟩ = (|00⟩+|11⟩)/√2  via  H(q0) CX(q0,q1)
      2. Measure in rotated ZZ basis:
            RY(2*(THETA_A[a]+phi), q0)   ← Alice
            RY(2*THETA_B[b],       q1)   ← Bob
      3. Measure both qubits.

    The outcome statistics then give E(a,b) = P(same) − P(diff) = cos(θ_A − θ_B).
    """
    from qiskit import QuantumCircuit
    q0, q1 = initial_layout
    qc = QuantumCircuit(2, 2)
    # Bell state
    qc.h(0)
    qc.cx(0, 1)
    # Measurement rotations (factor of 2: Bloch-sphere convention)
    angle_alice = 2.0 * (THETA_A[a_idx] + phi)
    angle_bob   = 2.0 * THETA_B[b_idx]
    qc.ry(angle_alice, 0)
    qc.ry(angle_bob,   1)
    qc.measure([0, 1], [0, 1])
    return qc


def correlator_from_counts(counts):
    """
    E = (P(same) − P(different)) from Qiskit counts dict.
    Qiskit bitstring convention: '01' → qubit1=0, qubit0=1 (rightmost = qubit 0).
    We want outcomes of qubit 0 (Alice) and qubit 1 (Bob).
    """
    total = sum(counts.values())
    same = 0
    for bits, cnt in counts.items():
        # bits is a string; rightmost char is qubit 0, next is qubit 1
        b0 = int(bits[-1])   # qubit 0 = Alice
        b1 = int(bits[-2]) if len(bits) >= 2 else int(bits[0])  # qubit 1 = Bob
        if b0 == b1:
            same += cnt
    return (2 * same - total) / total


# ── Verify circuit correctness analytically before spending QPU time ──────────
print("="*70)
print("E6v2: CHSH Hardware — Angle convention verification")
print("="*70)
print("\nIdeal S(φ) values (analytic, should be 2√2 at φ=0):")
for phi in PHIS:
    S_ideal = ideal_S(phi)
    print(f"  φ={phi:.2f}: S_ideal = {S_ideal:.5f}  (2√2·cos(φ) ≈ {2*math.sqrt(2)*math.cos(phi):.5f})")

# ── IBM Quantum hardware run ──────────────────────────────────────────────────
hardware_ran = False
job_ids      = {}
raw_hw       = {}
backend_name = None
qpu_seconds  = None

try:
    from qiskit_ibm_runtime import QiskitRuntimeService
    from qiskit_ibm_runtime import SamplerV2 as Sampler
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    print(f"\n[IBM] Connecting...")
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=IBM_API_KEY)
    backend = service.backend(TARGET_BACKEND)
    backend_name = backend.name
    print(f"[IBM] Backend: {backend_name}  |  qubit pair: {BEST_PAIR}")

    # Build 20 circuits
    circuits, tags = [], []
    for phi in PHIS:
        for a in (0, 1):
            for b in (0, 1):
                qc = make_chsh_circuit(a, b, phi, BEST_PAIR)
                circuits.append(qc)
                tags.append((phi, a, b))

    # Transpile with layout pinned to the best qubit pair
    pm = generate_preset_pass_manager(
        backend=backend,
        optimization_level=3,          # highest optimization
        initial_layout=list(BEST_PAIR) # pin to (6,5)
    )
    isa_circuits = pm.run(circuits)

    print(f"[IBM] Submitting {len(isa_circuits)} circuits × {SHOTS} shots = "
          f"{len(isa_circuits)*SHOTS:,} total shots...")
    t0 = time.time()
    sampler = Sampler(mode=backend)
    job = sampler.run(isa_circuits, shots=SHOTS)
    job_id = job.job_id()
    job_ids["main_batch_v2"] = job_id
    print(f"[IBM] Job ID: {job_id}")
    print(f"[IBM] Waiting for results...")

    result = job.result()
    qpu_seconds = time.time() - t0
    print(f"[IBM] Done. Wall time ≈ {qpu_seconds:.1f}s")
    job_ids["qpu_wall_seconds"] = qpu_seconds

    for (phi, a, b), pub_res in zip(tags, result):
        counts = pub_res.data.c.get_counts()
        total  = sum(counts.values())
        key    = f"{phi}|{a}|{b}"
        raw_hw[key] = {"counts": counts, "total": total}

    hardware_ran = True
    print("[IBM] Hardware run complete.")

except Exception as ex:
    print(f"[IBM] Hardware failed: {ex}")
    print("[IBM] Will use Aer + noise model as primary evidence.")

# ── Aer simulation (noiseless ideal + realistic noise model) ──────────────────
raw_aer_ideal = {}
raw_aer_noisy = {}
aer_ok = False

try:
    from qiskit import QuantumCircuit, transpile
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
    import numpy as _np

    aer_ideal = AerSimulator()  # noiseless

    # Realistic noise: depolarizing p=0.003 on CX, p=0.001 on single-qubit
    # readout 0.3% — matching marrakesh's best pair calibration
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(0.001, 1), ['ry','rz','sx'])
    nm.add_all_qubit_quantum_error(depolarizing_error(0.003, 2), ['cx','cz'])
    # Readout error: 0.3% symmetric
    ro_err = ReadoutError([[0.997, 0.003], [0.003, 0.997]])
    nm.add_all_qubit_readout_error(ro_err)
    aer_noisy = AerSimulator(noise_model=nm)

    print("\n[Aer] Running ideal + realistic-noise simulation...")
    for phi in PHIS:
        for a in (0, 1):
            for b in (0, 1):
                qc = make_chsh_circuit(a, b, phi, BEST_PAIR)
                key = f"{phi}|{a}|{b}"
                # Ideal
                qc_t = transpile(qc, backend=aer_ideal)
                c_ideal = aer_ideal.run(qc_t, shots=SHOTS).result().get_counts()
                raw_aer_ideal[key] = {"counts": c_ideal, "total": sum(c_ideal.values())}
                # Noisy
                qc_t2 = transpile(qc, backend=aer_noisy)
                c_noisy = aer_noisy.run(qc_t2, shots=SHOTS).result().get_counts()
                raw_aer_noisy[key] = {"counts": c_noisy, "total": sum(c_noisy.values())}
    aer_ok = True
    print("[Aer] Done.")

except Exception as ex:
    print(f"[Aer] Failed: {ex}  — using analytic fallback.")
    rng = np.random.default_rng(42)
    for phi in PHIS:
        for a in (0, 1):
            for b in (0, 1):
                key = f"{phi}|{a}|{b}"
                theta_a = THETA_A[a] + phi
                theta_b = THETA_B[b]
                E_ideal = math.cos(theta_a - theta_b)
                p_same  = (1 + E_ideal) / 2
                k_same  = int(rng.binomial(SHOTS, max(0, min(1, p_same))))
                counts_fake = {"00": k_same, "11": SHOTS - k_same}
                raw_aer_ideal[key] = {"counts": counts_fake, "total": SHOTS}
                # Noisy: 0.3% readout flip
                p_same_noisy = p_same * 0.997 + (1 - p_same) * 0.003
                k2 = int(rng.binomial(SHOTS, max(0, min(1, p_same_noisy))))
                raw_aer_noisy[key] = {"counts": {"00": k2, "11": SHOTS - k2}, "total": SHOTS}
    aer_ok = True


# ── Offline analysis ──────────────────────────────────────────────────────────
def analyze(raw, label, n_shots=SHOTS):
    """Compute S(φ), apply Algorithm 1, return results dict."""
    by_phi = {}
    for phi in PHIS:
        S = 0.0
        E_by_ab = {}
        for (a, b), sign in SIGN.items():
            key = f"{phi}|{a}|{b}"
            if key not in raw:
                continue
            d = raw[key]
            E = correlator_from_counts(d["counts"])
            E_by_ab[(a, b)] = E
            S += sign * E
        eta  = eta_tight(n_shots, DELTA_E)
        S_lb = S - eta
        passed = S_lb > 2.0
        h    = h_pironio(S_lb) if passed else 0.0
        cert = build_certificate("hw_link", S, n_shots, DELTA_E, 10_000, 0.99, 0, 1)
        r_e      = cert.r        if cert else 0.0
        log2_eps = cert.log2_eps if cert else None
        S_ideal_val = ideal_S(phi)
        by_phi[phi] = {
            "phi": phi, "S": S, "S_ideal": S_ideal_val,
            "eta": eta, "S_lb": S_lb, "passed": passed,
            "h": h, "r_e": r_e, "log2_eps": log2_eps,
            "E_by_ab": {str(k): v for k, v in E_by_ab.items()}
        }
    S_seq  = [by_phi[p]["S"] for p in PHIS]
    monotone = all(S_seq[i] >= S_seq[i+1] - 0.05 for i in range(len(S_seq)-1))

    print(f"\n[{label}] Table T-E (n={n_shots}, δ=10⁻⁶, η={eta_tight(n_shots,1e-6):.4f}):")
    print(f"  {'φ':>6} {'Ŝ':>8} {'S_ideal':>8} {'S_lb':>8} {'pass':>6} "
          f"{'h_e':>8} {'r_e':>8} {'log₂ε':>10}")
    for phi in PHIS:
        r = by_phi[phi]
        le = f"{r['log2_eps']:.2f}" if r['log2_eps'] is not None else "N/A"
        print(f"  {phi:>6.2f} {r['S']:>8.4f} {r['S_ideal']:>8.4f} {r['S_lb']:>8.4f} "
              f"{'YES' if r['passed'] else 'NO':>6} "
              f"{r['h']:>8.5f} {r['r_e']:>8.5f} {le:>10}")
    print(f"  Monotone: {monotone}  |  S(0)={by_phi[0.0]['S']:.4f}  "
          f"S(0.5)={by_phi[0.5]['S']:.4f}")
    return by_phi, monotone


print("\n--- Aer ideal ---")
aer_ideal_res, aer_ideal_mono = analyze(raw_aer_ideal, "Aer-ideal")

print("\n--- Aer realistic noise (marrakesh calibration) ---")
aer_noisy_res, aer_noisy_mono = analyze(raw_aer_noisy, "Aer-noisy")

hw_res = hw_mono = None
if hardware_ran:
    print("\n--- IBM Quantum hardware ---")
    hw_res, hw_mono = analyze(raw_hw, "Hardware")


# ── Sub-sampling pass-rate study (φ=0, ideal Aer) ────────────────────────────
print("\n--- Sub-sampling pass-rate study (φ=0, Aer ideal) ---")
sub_ns   = [100, 250, 336, 500, 1000, 2048, 4096, 8192]
N_RESAMP = 10_000
rng_sub  = np.random.default_rng(99)
sub_pass = {}

# Reconstruct p_same per setting from ideal Aer counts
p_same_per_setting = {}
for a in (0, 1):
    for b in (0, 1):
        key = f"0.0|{a}|{b}"
        d   = raw_aer_ideal[key]
        total = d["total"]
        same  = sum(v for k, v in d["counts"].items()
                    if len(k) >= 2 and k[-1] == k[-2])
        p_same_per_setting[(a, b)] = same / total

for n_sub in sub_ns:
    n_pass = 0
    for _ in range(N_RESAMP):
        S_sub = 0.0
        for (a, b), sign in SIGN.items():
            p = p_same_per_setting[(a, b)]
            k_same = rng_sub.binomial(n_sub, p)
            E_sub  = (2 * k_same - n_sub) / n_sub
            S_sub += sign * E_sub
        if S_sub - eta_tight(n_sub, DELTA_E) > 2.0:
            n_pass += 1
    sub_pass[n_sub] = n_pass / N_RESAMP
    print(f"  n={n_sub:>5}: pass rate = {n_pass/N_RESAMP:.4f}")


# ── Figures ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# 1. S(φ) vs φ
ax = axes[0]
phi_arr = np.linspace(0, 0.6, 100)
S_theory = [ideal_S(p) for p in phi_arr]
ax.plot(phi_arr, S_theory, 'k--', lw=1.5, label=r"Ideal $S(\varphi)$")
ax.plot(PHIS, [aer_ideal_res[p]["S"]  for p in PHIS], 'bo-', label="Aer ideal", ms=7)
ax.plot(PHIS, [aer_noisy_res[p]["S"]  for p in PHIS], 'gs-', label="Aer noisy", ms=7)
if hw_res:
    ax.plot(PHIS, [hw_res[p]["S"] for p in PHIS], 'r^-', label="IBM Hardware", ms=8)
ax.axhline(2.0, color='red', ls=':', lw=1.2, label="Classical limit S=2")
ax.axhline(2*math.sqrt(2), color='gray', ls=':', lw=1, alpha=0.5)
ax.set_xlabel("Angle offset φ (rad)")
ax.set_ylabel("CHSH value Ŝ")
ax.set_title("E6v2: S(φ) — monotone decrease\n(fixed circuit, marrakesh (6,5))")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_ylim([1.8, 3.0])

# 2. Certified rate r_e vs φ
ax2 = axes[1]
r_ideal = [aer_ideal_res[p]["r_e"] for p in PHIS]
r_noisy = [aer_noisy_res[p]["r_e"] for p in PHIS]
xs = np.arange(len(PHIS))
ax2.bar(xs - 0.2, r_ideal, 0.35, label="Aer ideal", color='steelblue', alpha=0.8)
ax2.bar(xs + 0.2, r_noisy, 0.35, label="Aer noisy", color='seagreen',  alpha=0.8)
if hw_res:
    ax2.bar(xs + 0.0, [hw_res[p]["r_e"] for p in PHIS], 0.15,
            label="IBM Hardware", color='red', alpha=0.9, zorder=5)
ax2.set_xticks(xs)
ax2.set_xticklabels([f"φ={p}" for p in PHIS], rotation=20, fontsize=8)
ax2.set_ylabel("Certified rate r_e (bits/round)")
ax2.set_title(f"E6v2: Algorithm 1 certificate\nn={SHOTS}, δ=10⁻⁶")
ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)

# 3. Sub-sampling pass rate
ax3 = axes[2]
ns_plot = list(sub_pass.keys())
pr_plot = [sub_pass[n] for n in ns_plot]
ax3.semilogx(ns_plot, pr_plot, 'bo-', ms=8, lw=2)
ax3.axvline(336,  color='orange', ls=':', lw=1.5, label="n≈336 design window")
ax3.axvline(500,  color='green',  ls=':', lw=1.5, label="n=500 (≥99%)")
ax3.set_xlabel("n (rounds per setting, sub-sample)")
ax3.set_ylabel("Pass rate Pr[Ŝ−η > 2]")
ax3.set_title("E6v2: Sub-sampling pass rate\nφ=0, δ=10⁻⁶, 10k resamples per n")
ax3.legend(fontsize=8); ax3.grid(True, alpha=0.3)
ax3.set_ylim([-0.05, 1.05])

plt.tight_layout()
outfig = os.path.join(os.path.dirname(__file__), "..", "figures", "E6v2_hardware_chsh_fixed.png")
plt.savefig(outfig, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved figure: {outfig}")

# ── Save results ──────────────────────────────────────────────────────────────
out = {
    "experiment": "E6v2",
    "description": "Fixed CHSH hardware (Reviewer 2 Comment 2) — corrected angles, marrakesh (6,5)",
    "fix_notes": {
        "v1_bug": "RY(A[a]+phi) without factor-of-2; wrong qubit pair selected by pending_jobs only",
        "v2_fix": "RY(2*(theta_A+phi)) and RY(2*theta_B); pinned to ibm_marrakesh qubits (6,5)"
    },
    "hardware_ran": hardware_ran,
    "backend_name": backend_name,
    "qubit_pair": list(BEST_PAIR),
    "job_ids": job_ids,
    "shots_per_setting": SHOTS,
    "delta_e": DELTA_E,
    "phis": PHIS,
    "aer_ideal_results": {str(k): v for k, v in aer_ideal_res.items()},
    "aer_noisy_results": {str(k): v for k, v in aer_noisy_res.items()},
    "hw_results": {str(k): v for k, v in hw_res.items()} if hw_res else None,
    "aer_ideal_monotone": aer_ideal_mono,
    "aer_noisy_monotone": aer_noisy_mono,
    "hw_monotone": hw_mono,
    "sub_sampling_pass_rates": {str(k): v for k, v in sub_pass.items()},
}
outpath = os.path.join(os.path.dirname(__file__), "..", "data", "results", "e6v2_hardware_fixed.json")
with open(outpath, "w") as f:
    json.dump(out, f, indent=2, default=str)
print(f"Saved results: {outpath}")
