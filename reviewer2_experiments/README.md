# BMR (OSN) — Reviewer 2 release package (Q1–Q17, D1–D6, C1–C5)

Everything needed to regenerate the numbers in `R2_results_v5_mrng_output.md`.
All commands are run **from this folder** (`Reviewer 2 output mrng/`).

## Folder tree

```
Reviewer 2 output mrng/
├── README.md                       this file
├── R2_results_v5_mrng_output.md    final write-up (C1–C5, D1–D6, Q1–Q17, proofs, appendices)
├── run_all.py                      one-shot classical regeneration (no QPU)
├── results_macros.tex              LaTeX macros – no hand-typed numbers
├── bmr/
│   ├── __init__.py
│   ├── bounds.py                   eta_tight, eta_union, h_pironio, n_min
│   ├── certificate.py              build_certificate → Cert dataclass (answers C4)
│   ├── route.py                    bmr_route (budget + hop-floor + bottleneck, exact)
│   ├── sim.py                      LinkSim (epoch simulation, certificate DB)
│   └── verify.py                   Ed25519 sign / verify
├── tests/
│   └── test_integration.py         5 integration tests (T1–T5)
├── experiments/
│   ├── e1_coverage.py              C2  MC radius coverage (200 000 trials)
│   ├── e2_table.py                 C2  Table T-A (finite-sample design window)
│   ├── e3_certificate_tests.py     C4  certificate unit tests (9 cases)
│   ├── e4_paths.py                 C3  Corollary 1 brute force (m ≤ 12)
│   ├── e5_runtime.py               C3  old runtime table (kept; superseded by Q9)
│   ├── e6_hardware.py              C5  QPU – DO NOT RUN (submits to IBM Quantum)
│   ├── e6v2_hardware_fixed.py      C5  QPU fixed – DO NOT RUN
│   ├── e7b_v2_corrected.py         Q1/Q2  E7b corrected (B2 = Pareto label-setting; BA now m = 2, 14 Oct)
│   ├── e10_ledger.py               C1  Algorithm 3 v1 certificate ledger
│   ├── generate_macros.py          Q16  generates results_macros.tex
│   ├── generate_q6_figures.py      Q6  sweep + dynamic figures
│   ├── q4_collusion.py             Q4  worst-case collusion experiment
│   ├── q5_evidence_source.py       Q5  Bell vs QBER evidence source
│   ├── q6_common.py                Q6  shared router family (B0–B8, BMR)
│   ├── q6_sweep.py                 Q6  full sweep (200 cells × 1000 trials × 10 routers)
│   ├── q6_dynamic.py               Q6  dynamic mid-session compromise
│   ├── q6_dynamic_latency.py       Q6  detection latency / false-alarm / false-negative tables (14 Oct)
│   ├── q6_sweep_degenerate.py      Q6  line and star, 64 cells (14 Oct)
│   ├── q9b_sparse_timing.py        Q9  sparse-graph timing with the 2t+1 > N-1 pre-check (14 Oct)
│   ├── q7_e0_audit.py              Q7  E0 audit (old pass probability vs new)
│   ├── q8_hardware_reconcile.py    Q8  hardware reconciliation (no QPU)
│   ├── q8_hw_repeats.py            Q8  QPU repeats – DO NOT RUN (submits to IBM Quantum)
│   ├── q9_solver_validation.py     Q9  exhaustive validation (36 000 instances) + timing
│   ├── q11_alg3_v2.py              Q11 Algorithm 3 v2 (T0–T12) + timing
│   ├── q13_scope_numbers.py        Q13 r_DI sweep and Prop. 2 Monte Carlo
│   ├── q14_twobudget_v2.py         Q14 two-budget extension (Theorem 4)
│   ├── q15_delta_regime.py         Q15 Δ-sweep, regime map, fraction sweep
│   └── q15_drift_robustness.py     Q15 drift / correlation robustness, ±25 % sensitivity
├── data/
│   ├── results/                    JSON and log files for every experiment
│   └── topologies/                 SNDlib: abilene, polska, nobel-germany, geant, germany50
├── figures/                        Q4, Q5, E7b_v2, Q6 ×2, Q8, Q9, Q14, Q15 ×2, E1–E6, E10
└── superseded/
    ├── README.md                   what was withdrawn and why
    ├── e7b_sr_dijkstra.py/.json    old "SR-Dijkstra" label (now B2 / Pareto-LS)
    ├── e11_twobudget.py/.json      old E11 (two-budget; replaced by Q14)
    └── E7b_sr_dijkstra_comparison.png  old figure
```

## Environment

| Item | Value |
|---|---|
| Python | 3.14.5 (macOS arm64, 12 cores) |
| numpy / scipy / networkx / matplotlib / cryptography | as installed in system Python 3.14 |
| qiskit / qiskit-aer / qiskit-ibm-runtime | 2.5.2 / 0.17.2 / 0.50.0 — in `../expr_r2/.venv` |
| **QPU scripts** | manual only — use `QPY=../expr_r2/.venv/bin/python experiments/q8_hw_repeats.py` |

All classical scripts use `python3` from the system path.
QPU scripts (`e6_hardware.py`, `e6v2_hardware_fixed.py`, `q8_hw_repeats.py`) are **not** invoked by `run_all.py` and should never be run without network egress to IBM Quantum.

## Seeds and conventions

Parameters are those of Table 2 in `R2_results_v5_mrng_output.md`; the authoritative values are the constants at the top of each script.

| Item | Value |
|---|---|
| Seeds (per script) | E1 42, E4 43, E5 44, E10 46, E7b-v2 50 (master), Q4 404, Q5 707, Q9 99 (+1,+2,+3 for sub-experiments), Q11 2026, Q13 13, Q14 14, Q15a 99, Q15b 1515, Q6 per-cell seeds in `q6_sweep.py`, Q6-dynamic 77 |
| δ | 1 × 10⁻⁶ |
| n per setting | 500 and 10⁴ (n = 100 only in the Q7 audit and the n_min discussion) |
| ℓ / λ_EC | 10⁴ / 0 |
| φ | 0.99 (sealed), 0.30 (public-XOR) |
| Δ | 0.5 at n = 500, 1 × 10⁻³ at n = 10⁴ (see D2 caveat) |
| Honest / compromised CHSH | S₀ = 2.687, S₁ = 2.418 |
| Collusion t | 1 (hop floor m ≥ 3) |
| Q6 trials | 1000 per cell (10 graph seeds × 100 audit draws) |
| Code integrity | no git repo; scripts print a SHA-256 of themselves / of `bmr/*.py` |

## Exact commands (classical, no QPU)

```bash
# all at once (fast mode ≈ 5 min, skips Q9 timing / Q14 / Q15-drift)
cd "Reviewer 2 output mrng"
python3 run_all.py --fast

# full run (~15–25 min)
python3 run_all.py
```

### Per-experiment commands

```bash
# Integration tests (5 pass)
python3 tests/test_integration.py

# Q4 — collusion (BMR 0 %, SPF up to 69.8 %)
python3 experiments/q4_collusion.py

# Q5 — evidence source (classical mimic accepted by QBER 100 %, by BMR 0 %)
python3 experiments/q5_evidence_source.py

# Q1/Q2 — corrected E7b table (56-row self-check, B2 = Pareto label-setting)
python3 experiments/e7b_v2_corrected.py

# Q6 — full sweep (200 cells, 10 routers, 1000 trials each)
python3 experiments/q6_sweep.py
Q6_VARIANT=phi09 python3 experiments/q6_sweep.py     # sensitivity sweep, phi=0.9, Delta=1e-4 / 1e-5
python3 experiments/q6_sweep_degenerate.py                 # line, star (64 cells)
python3 experiments/q6_dynamic.py
python3 experiments/q6_dynamic_latency.py                  # latency tables (14 Oct)
python3 experiments/q9b_sparse_timing.py                   # sparse timing, pre-check applied (14 Oct)
python3 experiments/generate_q6_figures.py

# Q7 — E0 audit (old pass prob at n=100 was 98.5 %, new rule: 0.7 %)
python3 experiments/q7_e0_audit.py

# Q8 — hardware reconcile from stored raw counts (no QPU)
python3 experiments/q8_hardware_reconcile.py

# Q9 — exhaustive validation (36 000 instances, 0 mismatches) + timing
python3 experiments/q9_solver_validation.py

# Q11 — Algorithm 3 v2, T0–T12 (T1–T11 rejected, T12 accepted)
python3 experiments/q11_alg3_v2.py

# Q13 — r_DI sweep + Prop. 2 Monte Carlo
python3 experiments/q13_scope_numbers.py

# Q14 — two-budget extension (Theorem 4)
python3 experiments/q14_twobudget_v2.py

# Q15 — regime map, fraction sweep, drift / ±25 % robustness
python3 experiments/q15_delta_regime.py
python3 experiments/q15_drift_robustness.py

# C2 — Radius coverage + design window (n_min, Table T-A)
python3 experiments/e1_coverage.py
python3 experiments/e2_table.py

# C4 — certificate unit tests (9 cases)
python3 experiments/e3_certificate_tests.py

# C3 — Corollary 1 brute force
python3 experiments/e4_paths.py

# LaTeX macros (no QPU, no hand-typed numbers)
python3 experiments/generate_macros.py

# QPU repeats (MANUAL ONLY — needs ../expr_r2/.venv + IBM Quantum account)
# ../expr_r2/.venv/bin/python experiments/q8_hw_repeats.py
```

## Output files (`data/results/`)

| File | Produced by |
|---|---|
| `q4_collusion.json` | `q4_collusion.py` |
| `q5_evidence_source.json` | `q5_evidence_source.py` |
| `e7b_v2_corrected.json`, `e7b_v2_log.txt` | `e7b_v2_corrected.py` |
| `q6_sweep.json`, `q6_sweep_phi09.json`, `q6_dynamic.json`, logs | `q6_sweep.py`, `q6_dynamic.py` |
| `q6_sweep_degenerate.json`, `q6_dynamic_latency.json`, `q9b_sparse_timing.json` | `q6_sweep_degenerate.py`, `q6_dynamic_latency.py`, `q9b_sparse_timing.py` |
| `budget_only_deltas_per_cell.csv` (400 rows: both configurations) | derived from the two sweep files |
| `q8_hw_refetched.json` (job IDs, calibration, raw counts of the 12 Oct run), `q8_hw_repeats_14oct.json` | QPU fetch (read-only) |
| `pre_rmin_14oct/` | pre-fix copies, kept for audit only; do not cite |
| `q7_e0_audit.json` | `q7_e0_audit.py` |
| `q8_hw_repeats.json` (raw counts, QPU) , `q8_hardware_reconcile.json` | `q8_hw_repeats.py` (manual), `q8_hardware_reconcile.py` |
| `q9_solver_validation.json` | `q9_solver_validation.py` |
| `q11_alg3_v2.json` | `q11_alg3_v2.py` |
| `q13_scope_numbers.json` | `q13_scope_numbers.py` |
| `q14_twobudget_v2.json` | `q14_twobudget_v2.py` |
| `q15_delta_regime.json`, `q15_drift_robustness.json` | `q15_*.py` |
| `e1_coverage.json`, `e2_table.json`, `e3_tests.json`, `e4_paths.json`, `e5_runtime.json`, `e10_ledger.json` | `e1…e10` |
| `e6_hardware.json`, `e6v2_hardware_fixed.json` | historical QPU runs (not re-run) |

## Notes for the response letter

1. **B2 label change:** the router previously called "SR-Dijkstra" is now called "Pareto label-setting (= B2)" and files with the old label are in `superseded/`.
2. **n = 100 retraction:** at δ = 10⁻⁶ and S = 2.687, the Hoeffding radius η_tight = 1.051 > 0.687; n = 100 cannot certify. The minimum certifying n is 235 (tight) or 1078 (union). Corrected n = 500 is used in Q2.
3. **QPU data provenance:** ten ibm_marrakesh runs: E6v2, three repeats (`db70r0kvf2bc73d205bg`, `db70rdqmb58s738cj47g`, `db70rmqmb58s738cj4sg`), five on 14 Oct (`data/results/q8_hw_repeats_14oct.json`), and the 12 Oct run (job `db659bklf4us73c4qju0`, Ŝ(0) = 2.776, counts and calibration in `data/results/q8_hw_refetched.json`). Ŝ(0) over the ten runs has mean **2.704**, SD **0.066**, SE 0.021, n_min = 223 (tight, δ = 10⁻⁶). The nine-run values (2.696 / 0.064) and the four-run snapshot (2.650 / 0.073) are retired. The only withdrawn figures are the submitted-manuscript ibm_fez values 2.574 ± 0.040 and 2.579 ± 0.013. The failed ibm_kingston run (Ŝ = 1.226, missing factor 2 in RY(φ)) is disclosed, not averaged.
4. **Bibliography:** all references in Q17 were checked against Crossref/arXiv (metadata only; full text read for Pirandola et al. only).
5. **Pirandola baseline:** B1a/B1b implement only the single-path λ-protocol (Eq. 35). Multipath and security-optimised variants are not claimed.
6. **Q-PATH and Q-CAST:** B6 and B7 use only their objective functions. Full implementations of Q-PATH (Li et al., IEEE Trans. Commun. 2022) and Q-CAST (Shi & Qian, SIGCOMM 2020) are not reproduced.
7. **Changes of 14 Oct (third pass), details in `../Reviewer 2/answers_c5_14oct.md`:**
   - `evaluate()` in `q6_common.py` applies R_MIN in the compliance test; `q6_sweep.py` checks "compliant ⇒ feasible" in every cell. The check exposed a floating-point tie in `bmr/route.py`, fixed with `_TOL = 1e-15` (0 mismatches in 36,000; 21,527 feasible).
   - `bmr_route` returns None at once when 2t+1 > N−1.
   - E7b used `barabasi_albert_graph(n, 1)` (a tree) and now uses m = 2, as the sweep always did. BMR refusal on BA is 2.0 % in the 200-cell sweep (ring 32.6 %, all 7.2 %).
   - Sweeps, dynamic and latency studies, sparse timing and the solver validation were re-run. Pre-fix copies are in `data/results/pre_rmin_14oct/`.
   - `R2_results_v5_mrng_output.md` still quotes pre-fix numbers. Read it together with the answer files; `results_macros.tex` is current.
