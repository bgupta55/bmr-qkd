# Bell-Metered Routing (BMR): Code and Data Deposit

**Paper:** *Bell-Metered Routing (BMR): Security-Budget-Constrained Path Selection for QKD over Optical Networks*  
**Journal:** Optical Switching and Networking (Elsevier) — major revision, October 2026  
**Authors:** [authors] · **Proof reader for theorems:** Sanjeev Rana

---

## Overview

BMR is a routing framework for quantum key distribution (QKD) networks. Each link runs a live CHSH Bell audit, generating a finite-key certificate `(r, ε, δ, γ, T)`. A route is a *simple path* whose summed failure bound Σ(ε+δ) stays within a budget Δ and whose length satisfies a collusion floor 2t+1. A calibrated alarm detects link degradation independently of the certificate. The destination verifies a signed, session-bound ledger.

This deposit contains all code and data needed to reproduce every number in the paper. One command regenerates all results:

```bash
python run_all.py
```

---

## Repository Structure

```
.
├── bmr/                        # Core library
│   ├── bounds.py               # η_tight / η_union radii; h_pironio entropy function
│   ├── certificate.py          # Algorithm 1: build_certificate() from audit counts
│   ├── route.py                # Algorithm 2: bmr_route() (Yen); brute_force_max_rate()
│   ├── sim.py                  # LinkSim: binomial CHSH sampling
│   └── verify.py               # Ledger verification (Algorithm 3)
│
├── experiments/                # Numbered experiments (each self-contained)
│   ├── e1_coverage.py          # η_tight / η_union coverage test (Table 4)
│   ├── e2_table.py             # Reference-value table (Table 4, §8.1)
│   ├── e3_certificate_tests.py # Unit + integration tests (Algorithm 1)
│   ├── e4_paths.py             # Walk counter-examples (§8.3)
│   ├── e5_runtime.py           # Solver timing sweep (Table 11, §8.10)
│   ├── e6_hardware.py          # IBM Quantum φ-sweep (§8.13)  ← credentials via env
│   ├── e6v2_hardware_fixed.py  # IBM Quantum 10-run repeat + fez opposite-sign test
│   ├── e10_ledger.py           # Ledger overhead (Table 10, §8.11)
│   ├── q4_collusion.py         # Collusion exposure (§8.3)
│   ├── q5_evidence_source.py   # Bell vs QBER evidence (§8.4)
│   ├── q6_common.py            # Router definitions used by sweep + dynamic
│   ├── q6_sweep.py             # 200-cell sweep (Table 6, §8.5)
│   ├── q6_dynamic.py           # Dynamic compromise study (Table 8, §8.2)
│   ├── q6_dynamic_latency.py   # Alarm latency detail (§8.2)
│   ├── q6_sweep_degenerate.py  # Degenerate topologies (ring, star, line)
│   ├── q7_e0_audit.py          # Audit design window (§4.2)
│   ├── q8_hw_repeats.py        # Hardware repeat runs (§8.13)
│   ├── q8_hw_repeats_14oct.py  # 14 Oct session, 10 runs, raw counts
│   ├── q8_hardware_reconcile.py# Session drift analysis
│   ├── q9_solver_validation.py # 36,000-instance Yen vs brute-force (§8.10)
│   ├── q9b_sparse_timing.py    # Timing on sparse graphs (§8.10)
│   ├── q11_alg3_v2.py          # Ledger attack table T0–T12 (§8.11)
│   ├── q13_scope_numbers.py    # Scope / certificate byte layout (§3.2)
│   ├── q14_twobudget_v2.py     # Two-budget solver study (§8.10, supplement)
│   ├── q15_delta_regime.py     # Regime map / Δ-feasibility (Fig. 9, §8.6)
│   ├── q15_drift_robustness.py # Drift robustness (§8.12)
│   └── generate_macros.py      # Regenerates results_macros.tex from JSON
│
├── rr_results/                 # (R)-family re-runs (new, October 2026)
│   ├── RR1_decomposition_R_family.json   # Decomposition ER+BA N=30 (Table 7)
│   ├── RR2_ring_flagship_R_family.json   # Ring flagship 500 trials (Fig. 5/6)
│   ├── RR4_dynamic_alarm.json            # Dynamic study + BMR+alarm column (Table 8)
│   └── RR5_cross_validation.json         # Yen vs brute-force 4,968 instances (§8.10)
│
├── data/results/               # All pre-computed result JSON files
│   ├── q6_sweep.json           # 200-cell sweep (main results)
│   ├── q6_dynamic.json         # Dynamic compromise (100 reps)
│   ├── q9_solver_validation.json # 36,000 solver instances
│   ├── q8_hw_repeats_14oct.json  # 10 IBM Quantum runs, raw counts + job IDs
│   └── ...                     # (see SHA256SUMS.txt for full list)
│
├── tests/                      # Unit and integration tests
├── run_all.py                  # Reproduces every result in sequence
├── run_pending_fast.py         # (R)-family RR1+RR2+RR5 re-runs (200 trials)
├── run_rr4_alarm.py            # (R)-family RR4 dynamic + alarm column
├── results_macros.tex          # LaTeX macros — every number in the paper
└── SHA256SUMS.txt              # SHA-256 checksums of all source and result files
```

---

## Key Parameters

| Symbol | Value | Role |
|--------|-------|------|
| `n` | 10,000 | audit rounds per setting (4n total) |
| `δ` | 10⁻⁶ | per-link soundness |
| `φ` | 0.99 | extraction fraction (sealed) |
| `ℓ` | 10,000 | raw rounds per key block |
| `S₀` | 2.687 | honest link baseline |
| `η_tight` | √(8 ln(1/δ)/n) | finite-sample radius (primary) |
| `R_MIN` | 0.01 | minimum acceptable rate |
| `t` | 1 | collusion floor (hop ≥ 3) |
| `Δ` | 10⁻³ | security budget (main) |
| `α` | 10⁻⁶ | alarm false-alarm rate |

---

## Quick Start

```bash
# 1. Clone and install
git clone <this-repo-url>
cd bmr-qkd-routing
pip install numpy scipy networkx

# 2. Reproduce the reference-value table (Table 4)
python experiments/e2_table.py

# 3. Run the main sweep (Table 6, ~5 min)
python experiments/q6_sweep.py

# 4. Run (R)-family decomposition (Table 7, ~10 min)
python run_pending_fast.py

# 5. Regenerate all LaTeX macros
python experiments/generate_macros.py
```

**IBM Quantum hardware scripts** (`e6_hardware.py`, `e6v2_hardware_fixed.py`, `q8_hw_repeats_14oct.py`) require a valid IBM Quantum account. Set the environment variable before running — **never paste the token into code**:

```bash
export IBM_QUANTUM_TOKEN="your-token-here"
python experiments/e6v2_hardware_fixed.py
```

---

## Results Summary

| Contribution | Key number | Source |
|---|---|---|
| Minimum audit rounds | **n_min = 235** (η_tight), 1,078 (η_union) | `e2_table.py` |
| BMR compliance | **92.8%** vs ≤78.1% (all competitors) | `q6_sweep.json` |
| Collusion exposure | **0.0%** (BMR) vs 10.6–69.8% (others) | `q4_collusion.json` |
| Rate cost of floor | **−0.0028** vs floor-free (2.2%) | `q6_sweep.json` |
| Solver validation | **0 / 40,968 mismatches** (Yen vs brute-force) | `q9_solver_validation.json` + `RR5_cross_validation.json` |
| Alarm detection | **500/500** at S†=2.55 (power 1.000), mean latency **26.5 epochs** at S†=2.62 | `RR2_ring_flagship_R_family.json` |
| IBM Quantum plausibility | Ŝ = **2.704 ± 0.066** (ten runs); fez diff **−0.189 ± 0.004** | `q8_hw_repeats_14oct.json` |

---

## Experiment Families

### (R) Family — derived ε (main text)
Uses η_tight, ε derived from each audit, n=10⁴, φ=0.99, Δ∈{10⁻³, 3×10⁻⁵}. Every trial draws an independent graph. Scripts: `q6_sweep.py`, `run_pending_fast.py`, `run_rr4_alarm.py`.

### (D) Family — declared class ε (Appendix D)
Uses ε from a fixed class table, φ=1.0, Δ=9×10⁻⁶, 50 graphs × 5 draws. Results in `data/results/` prefixed `R*` and `L*`. BMR−SPF = +0.2127 [+0.194, +0.232] under these settings.

---

## Security Notice

- **No credentials are stored in this repository.**
- All IBM Quantum scripts load tokens exclusively from environment variables.
- Job IDs in `q8_hw_repeats_14oct.json` are read-only references to public IBM job records; they contain no authentication material.
- Before reusing any script that previously contained `IBM_API_KEY`, verify the token has been rotated in the IBM Quantum dashboard.

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{bmr2026,
  title   = {Bell-Metered Routing ({BMR}): Security-Budget-Constrained
             Path Selection for {QKD} over Optical Networks},
  journal = {Optical Switching and Networking},
  year    = {2026},
  note    = {Major revision submitted October 2026}
}
```

---

## License

Code: MIT. Data: CC-BY 4.0. See `LICENSE` (to be added).

---

*Proof-reading of Corollary 1, Theorems 3a–d, Theorem 4, Propositions 1–3 and Lemma 1: Sanjeev Rana.*
