# R2 Results (v5-mrng, rebuilt) — Evidence and Answers for Reviewer 2
> **Superseded numbers (14 Oct, third pass).** The sweep, E7b, dynamic, latency, sparse-timing and compliance figures below predate three fixes: R_MIN in the compliance test, the floating-point tolerance in `bmr/route.py`, and BA m = 2 in E7b. For those quantities use `results_macros.tex` and the files `Reviewer 2/answers_14oct.md` and `Reviewer 2/answers_c5_14oct.md`. Hardware: ten runs, Ŝ(0) = 2.704 ± 0.066 (SE 0.021), not the nine-run 2.696. Pre-fix copies of the data are in `data/results/pre_rmin_14oct/`.



**Manuscript:** *Bell-Metered Routing (BMR): Security-Budget-Constrained Path Selection for QKD over Optical Networks* — Optical Switching and Networking, major revision.
**Specification closed here:** `Reviewer 2/pending_v7.md` (Q1–Q17, decisions D1–D6, final gate).
**Deposit:** this folder (`Reviewer 2 output mrng/`); every number below comes from a JSON in `data/results/`; commands are in `README.md`.
**Environment:** Python 3.14.5, macOS arm64, 12 cores; `cryptography` 50.0.1; qiskit 2.5.2 / ibm_marrakesh for the hardware repeats only.
**Code integrity:** no git repo exists in this folder, so every script prints a SHA-256 of itself (e.g. Q6: `q6_common.py` 8f7824c72e36; Q7 prints the hashes of `bmr/certificate.py`, `bounds.py`, `sim.py`, `route.py`).

> **Reading rule used throughout (same as the Reviewer 1 package).** The measured number is reported even where it disagrees with an earlier statement or with the task file; where BMR does not win, that is stated. Nothing was tuned to match a number. **This file supersedes the previous v5 text**, which declared C1–C5 "CLOSED" on evidence that did not support several of its claims; §0.3 lists every retraction.

---

## 0. Status board and corrections

### 0.1 Status of every item in pending_v7

| Item | Pri | Status | Artefact (script → JSON → figure) | Headline result |
|---|---|---|---|---|
| Q1 E7b metric integrity | P0 | **Done** | `e7b_v2_corrected.py` → `e7b_v2_corrected.json` → `E7b_v2_corrected.png` | Impossible conditional rates removed; zeros averaged over all trials; self-checks pass; refusal rate reported (up to 72 % on sparse BA) |
| Q2 E7b parameters / plausibility | P0 | **Done** | `e7b_v2_corrected.py`, `q15_delta_regime.py`, `tests/test_integration.py` | n=500, φ=0.99, δ=1e-6: honest link w = 1.70e-2, so Δ=1e-3 can never be feasible; Δ=0.5 used at n=500 (plateau) and 1e-3 at n=10⁴; regime map shows the knee |
| Q3 Independent baselines | P0 | **Done** | `q6_common.py`, `q6_sweep.py` → `q6_sweep.json` | B1a/B1b implement Pirandola's λ-protocol (Eq. 35) as published; the old "SR-Dijkstra" is **renamed Pareto label-setting (= B2)**; B3, B5, B6, B7, B8 all run |
| Q4 Collusion | P0 | **Done** | `q4_collusion.py` → `q4_collusion.json` | worst-case P(min relay coalition ≤ 1): SPF 19.7–69.8 %, NoFloor 10.6–56.3 %, **BMR 0.0 %** |
| Q5 Evidence source | P1 | **Done** | `q5_evidence_source.py` → `q5_evidence_source.json` | Classical-mimic link: BMR accepts 0 %, QBER-monitor accepts 100 % (true by construction, see honesty clause) |
| Q6 Sweep, dynamic, real topologies | P0 | **Done** | `q6_sweep.py`, `q6_dynamic.py` → `q6_sweep.json`, `q6_dynamic.json` | BMR ≈ strong rate-aware routers on rate (slightly below them on budget-only rate); +14.6 pt compliance; live-feed routers recover, frozen ones collapse |
| Q7 Simulator uses `build_certificate`; E0 audit | P0 | **Done** | `bmr/sim.py`, `tests/test_integration.py` (5/5 pass), `q7_e0_audit.py` | Old code used a different tolerance and rule (§Q7); its n=100 honest pass probability was 98.5 % vs 0.7 % now |
| Q8 Hardware reconciliation | P0 | **Done** (3 new QPU repeats) | `q8_hw_repeats.py`, `q8_hardware_reconcile.py` → `q8_hw_repeats.json`, `q8_hardware_reconcile.json` | Ŝ(0) over 4 runs: mean 2.650, SD 0.073; run-to-run spread ≫ shot noise ⇒ drift |
| Q9 Solver validation / runtime | P0 | **Done** | `q9_solver_validation.py` → `q9_solver_validation.json` | 36 000 instances vs exhaustive solver, **0 mismatches**; no path cap; timeouts reported |
| Q10 Proofs | P0 | **Written** (§Q10) | `e4_paths.py` (E4a brute force m ≤ 12) | Corollary 1, Thm 3a/3b/3c, Thm 4, parity, Prop. 2, Lemma 1 |
| Q11 Algorithm 3 redesign | P1 | **Done** | `q11_alg3_v2.py` → `q11_alg3_v2.json` | T1–T11 rejected 100 %, T0 accepted 100 %, T12 accepted (documented limit); Ed25519 ≈ 170–230 µs/entry |
| Q12 Literature evidence | P0 | **Done, with stated limits** | `data/lit/pirandola_prr2026.pdf/.txt`, reading notes (§Q12) | Pirandola read in full; Kumar, Meng **not accessed** (bibliographic data and Meng abstract verified); uncited work (i) identified as Meng et al., (ii) not located |
| Q13 Security scope | P0 | **Written** (§Q13) + `q13_scope_numbers.py` | `q13_scope_numbers.json` | r_e is a certificate quantity; at S=2.418 the asymptotic DI-QKD rate is ≤ 0 |
| Q14 E11 extension | P1 | **Done** | `q14_twobudget_v2.py` → `q14_twobudget_v2.json` | Pareto = DP on every random cell; 2ᵏ-type growth on diamonds; Lagrangian misses 1.0 % |
| Q15 Δ-sweep, regime map, robustness | P1 | **Done** | `q15_delta_regime.py`, `q15_drift_robustness.py` | Feasibility knee in Δ; drift/±25 % table |
| Q16 Number-lock, deposit, credentials | P0 | **Done** | `generate_macros.py` → `results_macros.tex`, `run_all.py`, scan (§Q16) | No credential found in any file; history cannot be scanned (no repo) |
| Q17 Attribution | P0 | **Done** (§Q17) | — | Q-PATH = Li et al.; Q-CAST = Shi & Qian |

### 0.2 Answers to the five comments (where to look)

| Comment | One-sentence answer | Where |
|---|---|---|
| C1 novelty | BMR is not the first security-rate-aware QKD router; the contribution is the *integration* of live Bell certificates, an additive route budget with a collusion floor, certificate-driven re-routing and destination verification; Pirandola et al. already distribute an additive ε budget, which we now acknowledge | Q12, Q3, Q4, Q5, Q11 |
| C2 η bound and n=100 | √(32 ln(8/δ)/n) is the union-bound Hoeffding slack (derived); a tighter √(8 ln(1/δ)/n) is used; n=100 cannot certify in expectation; the old code used a different, much smaller tolerance | Appendix A, Q7, Q8 |
| C3 simple path / Dijkstra | Route = simple path; hop floor on simple paths is strongly NP-complete when t is an input; no-floor problem is polynomial; the exact solver is validated and timed uncapped | Q9, Q10 |
| C4 five certificate fields | All five fields are constructed from named inputs; the simulator provably calls the same function | Q2, Q7, Q13 |
| C5 baselines | Ten routers on identical certificate data incl. Pirandola's λ-protocol; BMR ≈ rate-aware routers on rate, better on compliance, collusion, verification and dynamic recovery | Q1, Q3, Q6, App. D |

### 0.3 Retractions and corrections of the previous v5 text

| # | Previous statement | Finding | Corrected statement |
|---|---|---|---|
| 1 | B1 "Pirandola-style SR-Dijkstra" ≡ BMR in every cell | It was Pareto label-setting on (R, Σw ≤ Δ), i.e. BMR without the floor | Renamed **Pareto label-setting (= B2, BMR-NoFloor)**. The Pirandola baseline is now B1a/B1b = λ-protocol, Eq. (35) of the paper |
| 2 | Algorithm 2 uses Yen's algorithm "with a path cap" | `bmr_route` has **no cap**; the only limit in the study is a 1 s wall-clock deadline reported as "timeout" | "No path cap; exact enumeration; timeouts reported" (Q9) |
| 3 | "median < 0.7 ms for n ≤ 100, t ≤ 3" / "sub-2 ms" | True only on sparse graphs; on dense ER (p=0.3) t=3, n=100 the median is 132 ms and t=5 times out | Runtime is graph-dependent; table with timeouts in Q9 |
| 4 | "Strict monotonicity confirmed" (hardware) | Ŝ(0)=2.5405 vs Ŝ(0.1)=2.5388 differ by 0.002, shot-noise SE 0.022 | "monotone within shot noise" |
| 5 | Sub-sampling pass rates 4.1 / 93.9 / 99.9 % at n=100/250/336 | Computed from *ideal Aer*, not from the T-E hardware run | Recomputed from the hardware correlators: 0.0 / 9.2 / 34 % (run E6v2) and 0.5 / 55–63 / 90–93 % (repeats) |
| 6 | "BMR strictly guarantees 100 % compliance" | Compliance 92.8 % over the full sweep (99.5 % on clean random n=10⁴); BMR refuses the rest | "100 % of **delivered** keys comply; BMR refuses x % of requests" (refusal reported everywhere) |
| 7 | E7b "Mean Rate (All Trials)" 0.3247 at 56.9 % feasibility, SPF 0.2569 | Violates the stated metric (implied conditional rate 0.571 > any link rate; SPF upper bound 0.2035) | Recomputed with zeros included: Q1 table |
| 8 | E4b "15.8 % coalition violated" | Different denominator and cover test | Proper relay-cover test and both denominators: 40.6 % of all, 87.3 % of non-simple (ER n=10, t=1) |
| 9 | Old n=100 code "bypassed subtracting η" | It used a different, smaller tolerance and a different rule | See Q7 table |
| 10 | "SPF and Rate-Greedy frequently select 2-hop paths" | Unmeasured | Measured (Q4): 36.6 % (ER) to 69.8 % (hub) of delivered routes are 2-hop |
| 11 | E11 "fast < 1 ms up to n=30", Lagrangian "violates budgets" | Superseded by Q14; our Lagrangian is conservative (returns *no route*), it does not output infeasible routes | "Lagrangian misses feasible instances (1.0 %)" |
| 12 | E10 timing "≈149 µs/entry", earlier "37 µs" | Different designs and machines | Q11: mean ± SD, library, machine |
| 13 | Theorem 3c "weakly NP-complete" | No numeric data in the reduction | **Strongly** NP-complete when t is an input |
| 14 | "Nine-step" table with ten rows; "certificate issued with r_e = 0" | Count / conflict with Algorithm 1 | "Ten-step"; no certificate is published, link DEGRADED |
| 15 | E1 pass criterion `fail == 0` | Wrong; the guarantee is failure rate ≤ δ | Fixed and re-run: `all_pass_tight=True`, `all_pass_union=True` (200 000 trials per cell) |
| 16 | Q-PATH = "Shi et al.", Q-CAST = "Sun et al." | Swapped/wrong | Q-PATH = Li et al.; Q-CAST = Shi & Qian |

---

## 1. Design decisions (D1–D6) and parameters used

| ID | Decision taken | Consequence reported |
|---|---|---|
| D1 | **Fixed φ** (0.99 sealed, 0.30 public-XOR) with the extraction rule unchanged; regime map shown | Q15: at n=500, φ=0.99 the per-link defect is 1.7×10⁻², so budgets below ≈0.02–0.05 are infeasible; at φ=0.9 the defect is ≈10⁻⁶ and Δ never binds |
| D2 | Δ named in every table: **0.5 at n=500, 1×10⁻³ at n=10⁴** | Both lie on the feasibility plateau (Q15). **Honest caveat:** Δ=0.5 at n=500 is not a meaningful ε-security level; it makes the budget non-binding. A meaningful budget at n=500 needs φ≈0.9 (Q15 shows 99.5 % feasible at Δ=1e-5, rate 0.084 vs 0.0896 for φ=0.99). The full 200-cell sweep **was** repeated at φ=0.9 with Δ=10⁻⁴ (n=500) and 10⁻⁵ (n=10⁴): see Q6 “Sensitivity sweep at φ = 0.9” |
| D3 | Keep public-XOR, shortened, φ=0.30 | Included in the Q6 sweep as a mode |
| D4 | Q-PATH / Q-CAST kept as objective-only references (B6, B7); no superiority claim | Appendix D |
| D5 | First failed hardware run disclosed neutrally | Q8 |
| D6 | Mehic et al. kept qualitative (its metric needs key-pool and public-channel timing that the model does not contain) | Table 1 |

### Table 2 — parameters used by every experiment (unless a section says otherwise)

| Parameter | Value |
|---|---|
| Audit rounds per setting n | {500, 10⁴} (n=100 only in the E0 audit and the n_min discussion) |
| Confidence δ | 1×10⁻⁶ (tight bound η = √(8 ln(1/δ)/n) is the operational rule) |
| Block length ℓ | 10⁴; λ_EC = 0 |
| Extraction fraction φ | 0.99 (sealed), 0.30 (public-XOR) |
| Honest / compromised CHSH | S₀ = 2.687, S₁ = 2.418 (dynamic: S₁ ∈ {2.62, 2.55, 2.48, 2.418}) |
| Collusion | t = 1 ⇒ hop floor m ≥ 3 |
| Budget Δ | 0.5 (n=500), 1×10⁻³ (n=10⁴) |
| Freshness W / rate floor | W = 10 epochs; R_min = 0.01 bit/round |
| Compromised fraction | 20 % of links (E7b-v2, Q6, Q4: none) |
| Trials | 1000 per cell (10 graph seeds × 100 audit draws) |
| Audit statistic | Ŝ ~ N(S, √((4 − S²/4)/n)) (binomial CHSH estimator), drawn only inside `LinkSim.publish_certificate` |

---

# Part A — P0 items

## Q1 — E7b metric and table integrity

**Question.** Do the E7b numbers follow the metric definitions in the same file?
**Answer.** No, the old table did not (retraction 7). It was recomputed from per-trial records with the fields prescribed in pending_v7:

```
feasible[i]    = exists_route_within_budget(instance_i)      # same for all routers (Pareto-label oracle)
found[i,r]     = route_r is not None
compliant[i,r] = found and sum_w <= Δ  (and hops >= 2t+1 for BMR in this script)
secure[i,r]    = bottleneck rate if compliant else 0.0       # mean over ALL trials
```

**Self-checks executed (`e7b_v2_log.txt`: "All self-checks passed")**: mean(secure) ≤ mean(feasible)·max link rate; compliance of delivered BMR keys = 1.0; refusal = 1 − mean(found).

**Scope note (important).** In `e7b_v2_corrected.py` the *baselines'* compliance is budget-only; the hop floor is enforced and counted for BMR. Floor compliance of every router is measured in Q4 and in the Q6 "full" metric. Feasibility here is budget-only (no floor), which is why BMR's refusal can exceed 1 − feasibility.

ER n=30 (random) and BA n=30 (scale-free), 1000 trials, 20 % compromised links; rate (bit/round, zeros included) / compliance %; Δ=0.5 at n=500, 1e-3 at n=10⁴; t=1.

| Case | Feas. % | SPF | Rate-Greedy | Pareto-LS (=B2) | BMR | BMR refusal % |
|---|---|---|---|---|---|---|
| random n=500 clean | 100.0 | 0.0700 / 99.8 | 0.0821 / 100 | 0.0949 / 100 | 0.0915 / 99.1 | 0.9 |
| random n=500 compromised | 97.7 | 0.0460 / 69.5 | 0.0754 / 97.7 | 0.0863 / 97.7 | 0.0826 / 96.7 | 3.3 |
| random n=10⁴ clean | 100.0 | 0.3290 / 100 | 0.3315 / 100 | 0.3413 / 100 | 0.3387 / 99.5 | 0.5 |
| random n=10⁴ compromised | 96.4 | 0.2071 / 62.6 | 0.2909 / 87.7 | 0.3264 / 96.4 | 0.3225 / 95.6 | 4.4 |
| scale-free n=500 clean | 99.1 | 0.0666 / 99.7 | 0.0666 / 99.1 | 0.0666 / 99.1 | 0.0395 / 63.6 | 36.4 |
| scale-free n=500 compromised | 61.0 | 0.0405 / 64.0 | 0.0403 / 61.0 | 0.0403 / 61.0 | 0.0181 / 32.7 | 67.3 |
| scale-free n=10⁴ clean | 100.0 | 0.3284 / 100 | 0.3284 / 100 | 0.3284 / 100 | 0.1924 / 59.3 | 40.7 |
| scale-free n=10⁴ compromised | 56.5 | 0.1861 / 56.5 | 0.1861 / 56.5 | 0.1861 / 56.5 | 0.0912 / 28.0 | 72.0 |

(B3 Zhang-KR and B8 QBER-Mon give the same numbers as Rate-Greedy in this script; their differences appear in Q6/Q5/dynamic. All 56 rows are in the JSON.)

**What the corrected table says.**
- The previously impossible numbers are gone; e.g. SPF random n=10⁴ compromised is now 0.2071 ≤ 0.3251 × 0.626 as required.
- On random graphs BMR delivers 95.6–99.5 % compliant routes and gives up ≈2–3 % of rate against the no-floor router (0.3387 vs 0.3413 clean, 0.3225 vs 0.3264 compromised): the *cost of the collusion floor*.
- **Availability cost, stated plainly:** on sparse scale-free graphs BMR refuses 36–72 % of requests, because with t=1 a route needs ≥3 hops while many BA source–destination pairs are at distance ≤2. In those cells BMR's secure rate is *lower* than SPF's (0.1924 vs 0.3284 at n=10⁴ clean). A reader who values collusion resistance accepts this; a reader who does not will prefer B2. Both are shown.

## Q2 — E7b parameters and plausibility at n = 500

**Parameters** (now in every caption, Table 2): ℓ=10⁴, φ=0.99, δ=1e-6, t=1 (m ≥ 3), W=10, Δ=0.5 (n=500) / 1e-3 (n=10⁴).

**Certificate of one honest link, printed from `build_certificate` (S_hat = 2.687, δ=1e-6, ℓ=10⁴, λ_EC=0, epoch 3):**

| n | φ | η | h | H | k | r_e | log₂ ε | w = ε+δ | log₂ γ |
|---|---|---|---|---|---|---|---|---|---|
| 500 | 0.99 | 0.4702 | 0.0906 | 905.8 | 896 | 0.0896 | −5.88 | **1.70×10⁻²** | +890.1 |
| 500 | 0.90 | 0.4702 | 0.0906 | 905.8 | 815 | 0.0815 | −46.38 | 1.00×10⁻⁶ | +768.6 |
| 500 | 0.30 | 0.4702 | 0.0906 | 905.8 | 271 | 0.0271 | −318.4 | 1.00×10⁻⁶ | −47.4 |
| 2500 | 0.99 | 0.2103 | 0.2490 | 2490.0 | 2465 | 0.2465 | −13.48 | 8.88×10⁻⁵ | +2451.5 |
| 10⁴ | 0.99 | 0.1051 | 0.3424 | 3423.7 | 3389 | 0.3389 | −18.37 | 3.94×10⁻⁶ | +3370.6 |
| 10⁴ | 0.30 | 0.1051 | 0.3424 | 3423.7 | 1027 | 0.1027 | −1199 | 1.00×10⁻⁶ | −172.4 |

The reviewer's prediction is confirmed: at n=500, φ=0.99 an honest link has H−k = 10 bits, ε = 2⁻⁶, w ≈ 0.017; any Δ below ≈0.02 makes a multi-hop route infeasible. Earlier runs with Δ=10⁻³ at n=500 therefore **could not** have produced 99.5 % feasibility through `build_certificate`; they came from the old pipeline (Q7). Compromised link S₁ = 2.418: no certificate at n=500 (DEGRADED); at n=10⁴ w = 2.81×10⁻³, h = 0.1408. The regime map is in Q15.

**Which regime is reported.** Δ=0.5 (n=500) and Δ=10⁻³ (n=10⁴) with φ=0.99: both lie on the plateau where the budget filters nothing at n=500 and only compromised links at n=10⁴ (a compromised link alone has w = 2.8×10⁻³ > Δ). Stated honestly: at n=500 the budget is a formal guard, not a security level.

## Q3 — Independent Pirandola baseline and the missing routers

**Pirandola et al., PRR 8, 033221 (2026), read in full** (`data/lit/pirandola_prr2026.txt`). Relevant content: static trusted-node QKD network with per-link rate R_e and per-link ε_e; ETE ε composes additively along the path (Eq. 20 and following); single-path routing by Dijkstra variants — rate-optimised (widest-path, Eq. 29), security-optimised, and the hybrid **λ protocol** with update

T_v ← min{ T_v, max(T_u, R_uv⁻¹) + λ[T_u − log₁₀(1 − ε_uv)] }  (Eq. 35),

plus multipath routing with a rate–security trade-off. It has no Bell/CHSH evidence, no collusion floor, no expiry/monitoring and no NP-hardness result.

**Implemented baselines** (`q6_common.py`, same certificate data and seeds for all routers):

| ID | Router | Note |
|---|---|---|
| B0 | SPF | hop count |
| **B1a / B1b** | **Pirandola λ-protocol**, λ=0.01 / λ=1 (Eq. 35 as published, run on the same (R_e, ε_e=w_e)) | the independent Pirandola baseline; routes do not enforce a hop floor and only implicitly a budget |
| B2 | BMR-NoFloor (Pareto label-setting, t=0) | the ablation formerly called "SR-Dijkstra" |
| B3 | Zhang-style max-min key rate | re-implementation, not Zhang's RCKTA optimiser |
| B5 | Rate-greedy | |
| B6 / B7 | Q-PATH-style / Q-CAST-style | objective-only references; no superiority claim |
| B8 | QBER monitor (Hugues-Salas-style) | QBER from Q = (1 − Ŝ/2√2)/2 with the same finite-size treatment |
| BMR | Algorithm 2 | exact, t=1 |

"As-published" (frozen database) and "live-feed" modes are both run in the dynamic study (Q6). Results: Q6.

## Q4 — Collusion experiment

**Method** (`q4_collusion.py`; n_audit=10⁴, Δ=10⁻³, φ=0.99, t=1, 1500 trials per topology, 20 % no compromise; random requests): per delivered route, (a) *worst case*: smallest relay coalition that covers every link; record P(min coalition ≤ t=1) — for a route with m ≥ 2 hops this is exactly the share of 2-hop routes; (b) *random*: each relay dishonest independently with probability p ∈ {0.05, 0.10, 0.20}, record P(all links covered). Hops, rate and Σw are the cost.

| Topology | Router | Worst-case P(min coalition ≤ 1) % | Mean hops | Mean rate | Mean Σw | Random cover % (p=0.05 / 0.10 / 0.20) |
|---|---|---|---|---|---|---|
| ER | SPF | 36.6 | 2.53 | 0.3295 | 1.10e-5 | 1.9 / 3.5 / 9.5 |
| ER | Rate-Greedy | 36.6 | 2.53 | 0.3322 | 1.01e-5 | 1.7 / 3.2 / 9.4 |
| ER | BMR-NoFloor | 10.6 | 4.62 | 0.3418 | 1.43e-5 | 0.6 / 1.2 / 4.2 |
| ER | **BMR (t=1)** | **0.0** | 5.19 | 0.3405 | 1.61e-5 | 0.0 / 0.4 / 1.8 |
| BA | SPF | 43.8 | 2.37 | 0.3301 | 1.04e-5 | 2.1 / 3.7 / 10.1 |
| BA | Rate-Greedy | 43.8 | 2.37 | 0.3322 | 9.7e-6 | 2.0 / 3.3 / 10.1 |
| BA | BMR-NoFloor | 15.3 | 4.01 | 0.3409 | 1.28e-5 | 0.6 / 1.6 / 5.1 |
| BA | **BMR (t=1)** | **0.0** | 4.62 | 0.3390 | 1.49e-5 | 0.1 / 0.5 / 2.8 |
| ring | SPF | 19.7 | 3.27 | 0.3273 | 1.43e-5 | 0.9 / 2.7 / 6.1 |
| ring | Rate-Greedy | 19.7 | 3.27 | 0.3277 | 1.41e-5 | 0.9 / 2.7 / 6.1 |
| ring | BMR-NoFloor | 16.7 | 4.40 | 0.3298 | 1.78e-5 | 0.7 / 2.4 / 5.3 |
| ring | **BMR (t=1)** | **0.0** | 7.31 | 0.3231 | 3.07e-5 | 0.1 / 0.1 / 1.6 |
| hub | SPF | 69.8 | 1.70 | 0.3335 | 7.4e-6 | 4.3 / 7.9 / 13.4 |
| hub | Rate-Greedy | 69.8 | 1.70 | 0.3377 | 6.4e-6 | 3.5 / 6.9 / 14.5 |
| hub | BMR-NoFloor | 56.3 | 2.23 | 0.3421 | 7.4e-6 | 2.9 / 6.0 / 12.7 |
| hub | **BMR (t=1)** | **0.0** | 3.33 | 0.3365 | 1.1e-5 | 0.2 / 1.1 / 3.9 |

**Reading.** (i) The worst-case figure is exactly 0 for BMR on every topology — the structural guarantee of Corollary 1. (ii) The no-floor router is *not* collusion-safe: 10.6–56.3 % of its routes have ≤1 relay, in contrast to the earlier conjecture that "the floor rarely binds"; the ring and hub topologies make the point strongest. (iii) Random coalitions: BMR reduces coverage by 3–20× at p=0.05–0.20 but is **not** zero at p ≥ 0.10 — the guarantee is structural for ≤t colluders, probabilistic beyond. (iv) Cost: +0.6–4.8 hops, −0.1 to −2 % rate relative to the best no-floor router, Σw roughly +10–115 %, all far below any plausible Δ.

## Q6 — Full sweep, dynamic scenario, live-feed mode, real topologies

### Q6a — Sweep (`q6_sweep.py`)
200 cells = {ER, BA, WS, grid, ring} × {10, 20, 30, 50} × n ∈ {500, 10⁴} × {clean, compromised} × {sealed, public-XOR} (160 cells) + 5 SNDlib networks × same (40 cells); 1000 trials per cell; 10 routers; 51 s on 12 cores. Two metrics are reported because they answer different questions:
- **full** = secure rate with budget **and** hop floor required of every router (zeros included);
- **budget-only** = budget required, floor not required (isolates pure rate performance).

Average over the 200 cells:

| Router | Rate (full) | Compliance (full) % | Rate (budget-only) | Compliance (budget-only) % | Refusal % |
|---|---|---|---|---|---|
| B0 SPF | 0.0547 | 48.2 | 0.0970 | 82.1 | 0.0 |
| B1a Pirandola λ=0.01 | 0.0998 | 76.2 | 0.1259 | 93.9 | 3.9 |
| B1b Pirandola λ=1 | 0.0670 | 58.3 | 0.1115 | 91.5 | 3.9 |
| B2 BMR-NoFloor | 0.1043 | 78.1 | 0.1263 | 93.9 | 6.1 |
| B3 Zhang-style | 0.1042 | 77.9 | 0.1263 | 93.8 | 3.9 |
| B5 Rate-greedy | 0.0733 | 60.6 | 0.1175 | 92.8 | 3.9 |
| B6 Q-PATH-style | 0.0757 | 62.1 | 0.1192 | 93.2 | 0.0 |
| B7 Q-CAST-style | 0.0517 | 44.8 | 0.0940 | 78.7 | 8.4 |
| B8 QBER monitor | 0.0944 | 77.0 | 0.1156 | 93.4 | 3.9 |
| **BMR** | **0.1235** | **92.8** | 0.1235 | 92.8 | 7.2 |

Paired differences BMR − router (per-trial pooled mean, bootstrap 95 % CI over cells, Holm-corrected Wilcoxon p; `q6_sweep.json`):

| vs | Δ rate (full) [95 % CI] | Δ rate (budget-only) [95 % CI] |
|---|---|---|
| B0 SPF | +0.0688 [0.0587, 0.0797] | +0.0265 [0.0209, 0.0324] |
| B1a λ=0.01 | +0.0237 [0.0198, 0.0280] | −0.0024 [−0.0034, −0.0016] |
| B2 NoFloor | +0.0192 [0.0160, 0.0227] | **−0.0028** [−0.0038, −0.0020] |
| B3 Zhang | +0.0193 [0.0161, 0.0229] | **−0.0028** [−0.0037, −0.0020] |
| B5 Rate-greedy | +0.0502 [0.0424, 0.0586] | +0.0060 [0.0042, 0.0078] |
| B8 QBER | +0.0291 [0.0256, 0.0329] | +0.0080 [0.0062, 0.0096] |

All p_Holm < 10⁻²⁷⁷. BMR ≥ B3 (full) in 189 of 200 cells and never below it; in the budget-only metric BMR is below B3 in all 200 cells (by 0.0028 on average).

**Decision rule of pending_v7 applied.** The rate-aware routers B2/B3/B1a match BMR within 0.003 bit/round when the floor is not charged to them. The statement that goes into the paper is therefore: **"BMR matches strong rate-aware routers on rate and adds budget, collusion and verification guarantees."** The "full" gap (+0.019) is the price those routers pay for returning routes that violate the floor, not a rate advantage of BMR; BMR's rate is ≈2 % below theirs on budget-only.

Breakdown, synthetic topologies (rate / compliance % in the full metric):

| n_audit, condition, mode | SPF | B2 NoFloor | B3 Zhang | B5 | B8 QBER | BMR |
|---|---|---|---|---|---|---|
| 500, clean, sealed | 0.0346 / 58.2 | 0.0688 / 82.1 | 0.0682 / 80.8 | 0.0452 / 62.8 | 0.0439 / 77.6 | **0.0810 / 97.2** |
| 500, compromised, sealed | 0.0147 / 27.3 | 0.0544 / 69.8 | 0.0541 / 69.3 | 0.0395 / 57.9 | 0.0375 / 65.7 | **0.0633 / 81.7** |
| 10⁴, clean, sealed | 0.1954 / 60.4 | 0.2813 / 83.8 | 0.2812 / 83.8 | 0.1975 / 60.4 | 0.2690 / 83.7 | **0.3352 / 100.0** |
| 10⁴, compromised, sealed | 0.0709 / 21.8 | 0.2274 / 68.1 | 0.2273 / 68.0 | 0.1507 / 46.1 | 0.2182 / 67.7 | **0.2671 / 80.0** |
| 10⁴, clean, public-XOR | 0.0604 / 61.7 | 0.0856 / 84.3 | 0.0856 / 84.3 | 0.0611 / 61.7 | 0.0815 / 83.7 | **0.1014 / 99.9** |
| 10⁴, compromised, public-XOR | 0.0382 / 61.6 | 0.0758 / 84.3 | 0.0758 / 84.3 | 0.0548 / 67.3 | 0.0732 / 84.3 | **0.0894 / 100.0** |

By topology (full metric, rate / compliance %): BA — SPF 0.0427/36.4, B3 0.1061/77.0, BMR 0.1339/98.0 (refusal 2.0 %); ER — 0.0428/37.7, 0.1129/80.8, 0.1366/98.1 (1.9 %); WS — 0.0515/44.4, 0.1191/85.5, 0.1383/99.8 (0.2 %); grid — 0.0651/57.6, 0.1102/82.9, 0.1274/96.4 (3.6 %); **ring — 0.0662/61.2, 0.0699/60.9, 0.0770/67.4 (refusal 32.6 %)**. **The ring is BMR's weak case**: with a t=1 floor there are only two simple paths and many requests have none that is feasible (up to 92 % refusal in compromised ring cells). Rate by network size (BMR / B3 / SPF): n=10: 0.1220/0.0861/0.0314; n=50: 0.1215/0.1135/0.0686 — BMR's advantage shrinks as networks grow.

### Real topologies (sealed mode, mean over n_audit and condition; rate / compliance %)
Topology Zoo was **unreachable** from the sandbox, so Abilene and GEANT are the SNDlib instances of the same networks (stored in `data/topologies/`).

| Network | SPF | B1a | B3 Zhang | B8 | BMR |
|---|---|---|---|---|---|
| Abilene (12 nodes) | 0.0820 / 44.0 | 0.1358 / 69.4 | 0.1401 / 70.6 | 0.1347 / 70.3 | **0.1754 / 89.1** |
| Polska (12) | 0.0688 / 36.4 | 0.1453 / 71.9 | 0.1533 / 74.2 | 0.1422 / 73.3 | **0.1987 / 97.2** |
| Nobel-Germany (17) | 0.0907 / 48.4 | 0.1583 / 79.0 | 0.1635 / 80.3 | 0.1517 / 78.7 | **0.1935 / 95.4** |
| GEANT (22) | 0.0844 / 44.8 | 0.1586 / 78.3 | 0.1657 / 80.5 | 0.1491 / 78.1 | **0.2001 / 97.0** |
| Germany50 (50) | 0.1124 / 60.5 | 0.1861 / 91.4 | 0.1919 / 92.8 | 0.1713 / 89.8 | **0.2037 / 98.4** |

### Q6b — Dynamic mid-session compromise (`q6_dynamic.py`)
80 epochs; 20 % of links compromised at epoch 25, a further 10 % at epoch 60; S₁ ∈ {2.62, 2.55, 2.48, 2.418}; ring(30), grid 5×6 ("mesh"), ER(30); n ∈ {500, 10⁴}; 100 repetitions per cell (72 cells, 216 result rows with all series). Modes: **live** (router sees current certificates) and **frozen/as-published** (epoch-0 database; BMR-frozen = static-BMR ablation). Truth: each delivered route is scored on the *current* certificates. Mean secure rate / compliance % in windows pre (0–24) | post1 (25–59) | post2 (60–79), S₁ = 2.418:

| Case | BMR-live | BMR-frozen | B3-live | B3-as-published | B5-live | B8-live | SPF |
|---|---|---|---|---|---|---|---|
| grid, n=10⁴ | 0.336 → 0.320 → 0.288 (100→97→88 %) | 0.321 → 0.064 → 0.032 (100→20→10 %) | = BMR-live | 0.320 → 0.064 → 0.032 | 0.329 → 0.238 → 0.198 (100→73→61 %) | 0.320 → 0.310 → 0.282 | 0.322 → 0.087 → 0.039 |
| ER, n=10⁴ | 0.337 → 0.297 → 0.265 (100→89→80 %) | 0.322 → 0.087 → 0.032 (100→27→10 %) | = BMR-live | 0.322 → 0.087 → 0.029 | 0.329 → 0.281 → 0.231 | 0.322 → 0.286 → 0.257 | 0.324 → 0.120 → 0.049 |
| grid, n=500 | 0.084 → 0.069 → 0.057 (100→97→89 %) | 0.052 → 0.014 → 0.008 (97→27→18 %) | 0.084 → 0.069 → 0.057 | 0.052 → 0.013 → 0.007 | 0.074 → 0.062 → 0.053 | 0.052 → 0.049 → 0.043 | 0.057 → 0.015 → 0.009 |
| ring, n=10⁴ | 0.320 → 0.038 → 0.016 (100→12→5 %) | 0.316 → 0.032 → 0.013 | 0.320 → 0.038 → 0.016 | 0.316 → 0.032 → 0.013 | same | same | 0.318 → 0.038 → 0.016 |

**Reading.** (i) *Monitoring is what matters*: every router that sees fresh certificates (BMR-live, B3-live, B8-live) recovers; every router on a frozen database collapses (BMR-frozen and B3-as-published retain 10–18 % compliance after the second compromise). (ii) **B3-live matches BMR-live on rate** — the dynamic advantage belongs to *live Bell certificates + expiry*, not to BMR's search algorithm. (iii) At n=500, a link with S₁ ≲ 2.5 is rejected in the same epoch (no certificate), whereas S₁ = 2.62 and 2.55 certify at low rate (h ≈ 0.06 / 0.03) and are avoided only through the max–min preference. (iv) On the ring all routers fall to ≈0.038 / 0.016 after the compromises (the network, not the router, is the limit). (v) B8-live is lower in the clean state because the QBER rule is conservative at the same finite-size treatment.

Figures: `figures/Q6_sweep.png`, `figures/Q6_dynamic.png` (checked visually; titles match captions).

### Sensitivity sweep at φ = 0.9 (closes the former open item “full sweep at φ = 0.9”)
Same 200 cells, 1000 trials per cell, same seeds and routers; only the extraction fraction and budgets change: sealed φ = **0.90** (public-XOR stays 0.30), Δ = **10⁻⁴** at n = 500 and **10⁻⁵** at n = 10⁴ (both on the Q15 plateau for φ = 0.9, and four to five orders of magnitude below the Δ = 0.5 used in the headline sweep, so the budget is a real ε-security level). Run: `Q6_VARIANT=phi09 python3 experiments/q6_sweep.py` → `data/results/q6_sweep_phi09.json` (log `q6_phi09_log.txt`; 53 s, code hash `8f7824c72e36`). All integrity assertions of the sweep passed (BMR never delivers a non-compliant route; BMR found ⊆ feasible).

| Router | Rate (full) | Compliance (full) % | Rate (budget-only) | Compliance (budget-only) % | Refusal % |
|---|---|---|---|---|---|
| B0 SPF | 0.0543 | 50.2 | 0.0958 | 85.3 | 0.0 |
| B1a Pirandola λ=0.01 | 0.0917 | 75.6 | 0.1162 | 93.4 | 4.0 |
| B1b Pirandola λ=1 | 0.0634 | 59.3 | 0.1061 | 93.4 | 4.0 |
| B2 BMR-NoFloor | 0.0963 | 77.7 | 0.1171 | 93.8 | 6.2 |
| B3 Zhang-style | 0.0920 | 75.3 | 0.1127 | 91.3 | 4.0 |
| B5 Rate-greedy | 0.0685 | 60.9 | 0.1105 | 93.7 | 4.0 |
| B6 Q-PATH-style | 0.0703 | 61.6 | 0.1112 | 93.0 | 0.0 |
| B7 Q-CAST-style | 0.0522 | 47.1 | 0.0937 | 82.2 | 8.4 |
| B8 QBER monitor | 0.0831 | 73.9 | 0.1030 | 90.4 | 3.9 |
| **BMR** | 0.1143 | 92.5 | 0.1143 | 92.5 | 7.5 |

Paired differences BMR − router (per-trial pooled, bootstrap 95 % CI over cells; Holm p < 10⁻¹² for all rows below):

| Router | Δrate, full (floor+budget) | Δrate, budget-only |
|---|---|---|
| B0 SPF | +0.0600 [0.0516, 0.0688] | +0.0184 [0.0150, 0.0219] |
| B1a Pirandola λ=0.01 | +0.0225 [0.0190, 0.0265] | -0.0019 [-0.0027, -0.0012] |
| B2 NoFloor | +0.0180 [0.0151, 0.0213] | -0.0028 [-0.0036, -0.0022] |
| B3 Zhang-style | +0.0223 [0.0191, 0.0260] | +0.0016 [0.0000, 0.0032] |
| B8 QBER monitor | +0.0311 [0.0277, 0.0349] | +0.0113 [0.0093, 0.0133] |

**Reading.** The conclusions of the headline sweep are unchanged at the stricter operating point: (i) with the collusion floor in the objective BMR has the highest rate (0.1143 vs 0.0963 for the next best, B2) and the highest compliance (92.5 % vs ≤ 77.7 % for every baseline); (ii) in budget-only terms BMR is **not** better than the rate-aware baselines (−0.0028 vs B2, −0.0019 vs B1a, a negligible +0.0016 vs B3) — the same honest statement as at φ = 0.99; (iii) BMR refuses 7.5 % of requests (infeasible ones), 11.8 % in compromised cells. At n = 500 (φ = 0.9, Δ = 10⁻⁴) the BMR rate is 0.0660 sealed (feasible share 90.7 %) against 0.0228 for SPF; at n = 10⁴, 0.2768 against 0.1394.

## Q7 — Does the simulator use `build_certificate`? E0 audit

### Q7a — Integration test (`tests/test_integration.py`, 5/5 PASS)
```python
sim = LinkSim([("a","b"),("b","c")], n_audit=2500, delta=1e-6, ell=10_000, phi=0.99, lam_ec=0, seed=7)
published = sim.publish_certificate(("a","b"), epoch=3)
S_hat    = sim.last_S_hat[("a","b")]
expected = build_certificate(("a","b"), S_hat, 2500, 1e-6, 10_000, 0.99, 0, 3)
assert published == expected            # dataclass equality: all five fields + h, H, k, T, n
```
Further tests: 200 epochs with a compromised link (every published tuple equals `build_certificate` and the n=500 pass-rate is in the predicted band); the pass rule `Ŝ − η_tight > 2` is the *only* rule (checked at n ∈ {100, 336, 379, 500, 2500, 10⁴}); n=100 never certifies; a static check that Q1/Q4/Q5/Q6/Q9/Q15/Q11 call the library path. `LinkSim.publish_certificate` is the only certificate-creation path of the simulator; Q4, Q5, Q9, Q15 and E7b-v2 call `build_certificate` directly.

### Q7b — E0 audit: paper vs old code vs new code (`q7_e0_audit.py`, old code = `bmr-experiments/audit/bell_audit.py`)

| Quantity | Paper (revised) | Old code did | New code does |
|---|---|---|---|
| η subtracted before r_e? | yes, η_tight = √(8 ln(1/δ)/n) | yes, but with **tol = √(ln(2/δ)/(2n))** at δ=10⁻⁵ — a single-mean Hoeffding width (0.247 at n=100), not the CHSH-sum slack (1.05 tight, 2.26 union) | η_tight (library) |
| Pass rule | Ŝ − η > 2 | **Ŝ ≥ 2 + τ_param, τ_param = 0.10** (a "detection threshold"); the η-based certificate rule was not the rule for validity | Ŝ − η_tight > 2 (only rule) |
| n at the flagship run | 500 / 10⁴ | **n = 100** "low-n" run | purged; n ∈ {500, 10⁴} |
| ε_e | ½·2^(−(H−k)/2) | 2·exp(−2n·tol²) = δ_target exactly (1×10⁻⁵ for every n) | ½·2^(−(H−k)/2) |
| δ_e | per-link audit confidence | fixed per run (1×10⁻⁵) | fixed per run (1×10⁻⁶), not link-dependent |
| γ_e | 2^k ε_e (guessing excess) | 2^(−h) − 0.5 | log₂ γ = k + log₂ ε |
| T_e | audit epoch | audit epoch | audit epoch |
| n per setting or total? | per setting | n rounds per epoch ("audit rounds") with σ_Ŝ = 2/√n, not the binomial sd | n per setting; sd = √((4 − S²/4)/n) |
| Detection | certificate rule | threshold Ŝ ≥ 2.10 (neither midpoint nor certificate rule) | certificate rule; midpoint test only in Prop. 2 |

**Pass probability of an honest link (S_true = 2.687; 20 000 audits per cell):**

| n | old rule | R1 pipeline τ(union, 10⁻⁶) | **new (η_tight, library)** |
|---|---|---|---|
| 100 | 98.5 % | 0 % | **0.7 %** |
| 250 | 100 % | 0 % | 59.2 % |
| 500 | 100 % | 0 % | 99.9 % |
| 2500 | 100 % | 100 % | 100 % |
| 10⁴ | 100 % | 100 % | 100 % |

**Verdict on the prediction in pending_v7** ("old code did not subtract η and did not apply the pass rule at n=100"): **confirmed in substance, corrected in detail.** The old code did subtract a tolerance, but a far smaller one (0.247 vs the required 1.05 at n=100) and its validity rule was a fixed threshold; hence it passed 98.5 % of honest n=100 audits where the corrected rule passes 0.7 %. The statement in the paper is changed accordingly (retraction 9). The Reviewer 1 pipeline (`expr_r1`) used τ = √(32 ln(8/δ)/n) (union bound), consistent with the manuscript, with ε_e = δ as a placeholder — different from the new library's ε_e; this is why R1 numbers and R2 numbers are not interchangeable.

## Q8 — Hardware: reconcile runs, fix the sub-sampling statistics

### Q8a — Runs on `ibm_marrakesh`, qubit pair (6, 5), 8192 shots per setting, RY(2θ) convention
| Run | Job ID | Ŝ(0) | Note |
|---|---|---|---|
| E6v2 (v5's T-E) | db6nem84qg6s73c52f8g | 2.5405 | 5 φ-values; stored correlators |
| **repeat 1** | db70r0kvf2bc73d205bg | 2.6978 | UTC 2026-10-13 10:53; raw counts stored |
| **repeat 2** | db70rdqmb58s738cj47g | 2.6775 | 10:53–10:54; raw counts stored |
| **repeat 3** | db70rmqmb58s738cj4sg | 2.6824 | 10:54; raw counts stored |
| earlier run (quoted in v5) | n/a | 2.776 | raw data not retained — cannot be verified; **not used in any statistic** |
| failed `ibm_kingston` run | db6540kvf2bc73d0m1rg | 1.2256 | 2048 shots; RY(φ) without the factor 2 and a different qubit pair; Aer for the same circuit gave 2.84. Disclosed without an attributed cause (D5) |

**Mean over the four marrakesh runs: Ŝ(0) = 2.650, SD 0.073.** The shot-noise standard error per run is 0.022. The three repeats (taken within one minute of each other) agree with each other to within 0.02; the E6v2 run is 0.12–0.16 lower, i.e. ≈5–7 shot-noise SEs. **Device drift between sessions is the supported explanation**; we did not measure the cause (gate error, crosstalk, calibration cycle) and do not claim one. Calibration snapshot at the repeats: q6 T1 287 µs, T2 286 µs, readout error 0.27 %; q5 T1 243 µs, T2 56 µs, readout error 0.78 %.

### Q8b — Sub-sampling pass rate recomputed from the hardware correlators (φ = 0), δ = 10⁻⁶, η_tight
| n | E6v2 run (2.5405) | repeats 1 / 2 / 3 | old Ideal-Aer study (source of v5's numbers) |
|---|---|---|---|
| 100 | 0.000 | 0.007 / 0.005 / 0.006 | 0.041 |
| 250 | 0.092 | 0.635 / 0.546 / 0.567 | 0.939 |
| 336 | 0.340 | 0.935 / 0.899 / 0.907 | 0.999 |
| 379 | 0.497 | 0.980 / 0.961 / 0.968 | — |
| 500 | 0.849 | 0.9998 / 0.999 / 0.999 | 1.000 |
| 700 | 0.993 | 1.0 / 1.0 / 1.0 | — |

This reproduces the reviewer's recomputation for the E6v2 run (0.03 / 9.6 / 35 / 50 / 85 / 99 %) and shows that v5's 4.1 / 93.9 / 99.9 % came from the ideal-Aer data (retraction 5). The library-rule agreement check (the sub-sampling used the same `build_certificate` rule) was 500/500. n_min(Ŝ) from the library: 379 for Ŝ = 2.5405, 262 for the four-run mean 2.65, 235 for S = 2.687.

**Expected-result statement, as measured.** On this device certification at n = 500 per setting is reliable in the repeats (> 99.9 %) and 85 % in the lowest session; n ≈ 700 is needed to be reliable in all four sessions. "Finite-size certification on real devices is sensitive to drift."

### Q8c — Gap to Aer-noisy, monotonicity, interpretation
- Aer-noisy S(0) = 2.799 and Aer-ideal 2.809, versus the hardware mean 2.650 (gap 0.15). The measured pair readout errors (0.27 % / 0.78 %) are close to the modelled 0.3 %, so readout does not explain the gap; unmodelled two-qubit-gate error, crosstalk or drift remain candidates, none measured. The calibrated noise model overestimates the device in these sessions.
- **Monotone within shot noise**: all four runs decrease with φ (strictly in every run, but E6v2 steps 0→0.1 differ by 0.002 ≪ 0.022).
- Repeats at φ ∈ {0, 0.2, 0.5}: Ŝ = 2.698/2.539/1.567, 2.677/2.482/1.537, 2.682/2.484/1.556, close to S(0)·cos 2φ (2.485, 1.458 for the mean at φ=0.2, 0.5; the φ=0.5 point sits slightly above because of the noise floor).
- **State in the text:** this is a *calibration-error emulation*, not an adversarial attack; φ is half the Bloch rotation angle (RY(2φ)); S(φ) ≈ S(0)·cos 2φ.

Figure: `figures/Q8_hardware_reconcile.png`.

## Q9 — Algorithm 2 completeness, solver validation, runtime scope

1. **Cap.** `bmr_route`/`oracle` have no path cap; Yen's enumeration is uncapped and ends when the budget is exceeded. The only limit is a 1 s wall-clock deadline (reported as "timeout"). The v5 statement "with a path cap" was wrong; no "no secure route" answer is cap-induced.
2. **Validation against exhaustive search** (`q9_solver_validation.py`): random ER graphs (p = 0.5, n = 4…12), t ∈ {0, 1, 2, 3}, 1000 instances per (n, t) cell, loose and tight Δ; exhaustive DFS gives the optimal bottleneck rate. **36 000 instances, 21 527 feasible, 0 mismatches**, 0 excluded for DFS limits; every returned path is simple, has ≥ m_min hops and weight ≤ Δ. (Bug found and fixed in `brute_force_simple_paths` during this work; the validation uses the fixed enumerator `brute_force_max_rate`.)
3. **Uncapped full-pipeline timing** (rate/freshness filter + threshold binary search + Dijkstra/Yen + Δ test + hop floor; ER with p = 0.3, i.e. dense; 100 instances per cell; 1 s limit). Median ms (timeout %):

| n | t=0 | t=1 | t=2 | t=3 | t=5 |
|---|---|---|---|---|---|
| 10 | 0.07 | 0.10 | 0.19 | 0.38 | 0.66 |
| 20 | 0.21 | 0.24 | 0.52 | 1.32 | 10.8 (1 %) |
| 30 | 0.42 | 0.51 | 1.05 | 4.6 | 286 (14 %) |
| 50 | 1.04 | 1.16 | 3.2 | 22.7 | ≥ 1000 (95 %) |
| 75 | 2.07 | 2.28 | 7.6 | 69 | ≥ 1000 (100 %) |
| 100 | 3.5 | 3.9 | 12.4 | 133 | ≥ 1000 (100 %) |
| 150 | 7.5 | 8.0 | 26 | 372 (2 %) | ≥ 1000 (100 %) |

**The statement for the paper:** for t ≤ 2 the exact solver takes ≤ 26 ms up to n = 150 on dense random graphs; t = 3 takes up to ≈0.4 s at n = 150 (2 % timeouts); t = 5 (floor 11) is only practical to n ≈ 20–30 and the problem is NP-complete in t (Thm 3c). The old E5 sub-millisecond medians were measured on sparse graphs (mean degree 4–6) and are not a bound: "median < 0.7 ms" is withdrawn (retraction 3). With the default t = 1, the pipeline stays ≤ 8 ms at n = 150.
4. **E4b denominators** (walk solver vs simple path; 2200 instances per cell, proper relay-cover test): non-simple / coalition-violated, % of all instances (violated as % of non-simple):

| n | t=1 | t=2 | t=3 |
|---|---|---|---|
| 10 | 46.5 / 40.6 (87.3) | 89.5 / 88.0 (98.3) | 99.6 / 99.5 (99.9) |
| 20 | 33.5 / 28.2 (84.2) | 79.5 / 78.4 (98.6) | 98.3 / 98.2 (99.9) |
| 30 | 28.7 / 24.8 (86.4) | 76.5 / 75.3 (98.4) | 97.6 / 97.5 (99.9) |

The earlier "15.8 %" used a different cover test (and a different base: 31 of the 196 graphs where a walk was found). Both bases are now reported; the conclusion — exact simple-path enumeration is mandatory — is stronger, not weaker.
5. **Edge weights** at n = 2500, S ~ N(2.687, 0.02), φ = 0.99, δ = 10⁻⁶: w mean 8.7×10⁻⁵ (SD 4.9×10⁻⁵, min 9.6×10⁻⁶, max 4.8×10⁻⁴); r mean 0.2466 (SD 0.0153).

Figure: `figures/Q9_solver_validation.png`.

## Q10 — Proofs (Appendix C text)

Notation: route P = (v₀=s, e₁, v₁, …, e_m, v_m=d) is a **simple path** (distinct vertices); relays Relays(P) = {v₁,…,v_{m−1}}; s and d are honest. In sealed mode the end-to-end key is K = ⊕ⱼ Kⱼ.

**Corollary 1 (collusion floor).** A coalition C ⊆ Relays(P) learns K only if it knows every hop key Kⱼ, i.e. only if C covers all m links (every link has an endpoint in C). Claim: the minimum such |C| is ⌈m/2⌉ for m ≥ 2 (and no coalition exists for m = 1). *Proof.* Each relay is incident to ≤ 2 path links, so |C| ≥ m/2; the set {v₁, v₃, v₅, …} of size ⌈m/2⌉ covers all links (link e_m needs v_{m−1}, available when m is even via v_{m−1}, and when m is odd via v_{m−1} = v_{2·(m−1)/2}… choose {v₁, v₃, …} for even m and {v₁, v₃, …, v_{m−2}} ∪ {v_{m−1}} for odd m; both have ⌈m/2⌉ elements). Hence a coalition of size ≤ t cannot cover P iff ⌈m/2⌉ > t iff **m ≥ 2t+1**. ∎ The previous formula ⌊(m−1)/2⌋+1 equals ⌈m/2⌉ for every m ≥ 2. *Bouncing-walk example:* the walk [s, a, s, a, d] has m = 4 ≥ 3 but one relay {a} covers all links (cover 1) — the bound requires **simple** paths. *Check:* E4a brute force over m = 1…12: min coalition = ⌈m/2⌉ for m ≥ 2 and ∞ for m = 1; 0 mismatches; the bouncing walk has cover 1.

**Theorem 3a (no floor: polynomial).** Maximise min-rate subject to Σ w ≤ Δ over s–d paths. Feasibility at threshold ρ (keep edges with r_e ≥ ρ, ask for a path of weight ≤ Δ) is monotone non-increasing in ρ; only ≤ |E| distinct thresholds exist; binary search needs ⌈log₂|E|⌉ tests, each one Dijkstra (w_e > 0 so the shortest walk is a simple path). Total **O((|E| + |V| log |V|)·log |E|)**.

**Theorem 3b (walk relaxation).** The minimum-weight s–d *walk* with ≥ 2t+1 edges is found by Dijkstra on the layered graph with states (v, min(h, K)), K = 2t+1: K|V| states, K|E| arcs, **O(K(|E| + |V| log |V|))**. Every simple path with ≥ K hops is a walk, so the value lower-bounds the simple-path optimum; if the returned walk is simple it is optimal. The walk can repeat nodes and then gives no collusion protection (E4b).

**Theorem 3c (floor on simple paths: strongly NP-complete when t is an input).** *In NP:* verify simplicity, length ≥ K and weight ≤ Δ. *Hardness:* reduce from Hamiltonian s–d path on G with N vertices. If N is odd add one pendant vertex d′ adjacent only to d and replace d by d′, so that N is even and a Hamiltonian s–d path has N−1 (odd) edges. Set 2t+1 = N−1, all r_e equal, all w_e = 1, Δ = N. A simple path with ≥ N−1 edges has exactly N−1 edges and visits every vertex, i.e. is Hamiltonian. The instance contains only unary numbers, so the problem is **strongly** NP-complete. (For fixed t, K is constant and the problem is polynomial by Thm 3b + enumeration of ≤ n^K candidates; the exact solver's cost grows with t as measured in Q9.)

**Theorem 4 (two budgets).** Find an s–d path with W₁ ≤ B₁ and W₂ ≤ B₂ (e.g. w = ε+δ and a second additive weight). With B₂ unconstrained this is the **Restricted Shortest Path** problem, **weakly NP-hard** (Garey & Johnson, *Computers and Intractability*, 1979; problem **ND30** “Shortest Weight-Constrained Path” — **confirmed** from the Wikipedia article on Shortest path problem, which cites ND30 and notes the FPTAS due to Hassin). Weights are non-negative so a feasible walk can be shortcut to a simple path; pseudo-polynomial DP O(|E|·B₁) and an FPTAS exist (Hassin 1992; Lorenz & Raz 2001). Empirical confirmation: Q14. Complementary result: Zhang et al. (2024) prove NP-hardness of their *different* joint routing/channel/rate/time-slot problem.

**Parity rule (public-XOR).** A parity-doubled graph (states (v, h mod 2), arcs flipping parity) gives a polynomial algorithm for the minimum-weight *odd-length walk*; the walk may repeat nodes. For simple odd paths the paper makes **no polynomial-time claim**; the exact enumerator is used and its measured cost is reported (Q9).

**Proposition 2 (change detection by the midpoint test).** Let S* be the honest and S′ < S* the compromised CHSH value, Δ_S = S* − S′, θ = (S*+S′)/2, and Ŝ = (1/n)ΣX_j, X_j ∈ {−1, +1}, j ≤ 4n. Hoeffding: P(Ŝ ≤ θ | S*) ≤ exp(−2(nΔ_S/2)²/(4n·2²)) = exp(−nΔ_S²/32), and the same bound holds for P(Ŝ > θ | S′). So a window of n rounds per setting misclassifies with probability ≤ exp(−nΔ_S²/32), decaying exponentially in n; a compromised link is detected within a geometric number of audits with failure probability that decays the same way. *Numbers (S* = 2.687, S′ = 2.418, `q13_scope_numbers.json`, Monte Carlo with the simulator's binomial sd):* n = 500: false alarm 2.1 %, miss 3.0 % (Hoeffding bound 0.32, loose); n = 2500: 3×10⁻⁶ / 1.5×10⁻⁵ (bound 3.5×10⁻³); n = 10⁴: < 5×10⁻⁷ (bound 1.5×10⁻¹⁰). **Which mechanism triggers re-routing?** *Certificate issuance* (Ŝ − η > 2) and *change detection* (Prop. 2) are different mechanisms. The simulator (`LinkSim`) and Algorithm 2 use only the certificate rule: a link without a valid fresh certificate is absent from the routing graph, which is what removes S₁ = 2.418 at n = 500. The midpoint test is an analysis device and is not used in the experiments.

**Lemma 1 (XOR sanitisation).** If Y ∈ {0,1}^ℓ is ε-secure given side information S (trace distance of (Y,S) from (U,S) ≤ ε) and X = f(S) is any function of S (including adversarial choices; extra adversary randomness is part of S), then Y ⊕ X is ε-secure given S. *Proof.* For each s the map y ↦ y ⊕ f(s) is a bijection of {0,1}^ℓ applied conditionally on s; trace distance is invariant under this joint bijection, and U ⊕ f(s) is uniform given s. ∎

**Theorem 1 corollary (budget merge).** If every audit bound holds (probability ≥ 1 − Σδⱼ) the sealed-relay key is Σεⱼ-secure; on the complement event (probability ≤ Σδⱼ) the trace distance is ≤ 1. Hence the unconditional trace distance is ≤ (1 − Σδⱼ)·Σεⱼ + Σδⱼ ≤ **Σ(εⱼ + δⱼ) = Σ wⱼ ≤ Δ**.

*Review status.* These derivations were checked by computation where computable (E4a brute force m ≤ 12; Q9 exhaustive solver comparison; Q14 solver agreement). The derivations themselves still need an independent human re-derivation by a co-author before submission (pending_v7 Q10 "How to close").

## Q12 — Literature evidence

### Reading notes (what was actually accessed)
| Work | Accessed in this study? | Notes used for Table 1 |
|---|---|---|
| **Pirandola, Harney, Spedalieri**, Phys. Rev. Research 8, 033221 (2026), DOI 10.1103/7839-3xhy | **Yes — full text** (`data/lit/`) | Static trusted-node QKD networks (paper §III); per-link (R_e, ε_e); ETE ε additive along the route (composition equations around Eq. 20); single-path rate-optimised (Eq. 29), security-optimised and hybrid λ-protocol (Eq. 35); multipath; analysed on perturbed lattices; no Bell/CHSH evidence, no collusion model or hop floor, no expiry/monitoring, no hardness result (run time O(\|E\| + \|P\| log \|P\|) is stated for the Dijkstra variants) |
| Zhang et al., IEEE TNSM 21(1) 148–160 (2024), DOI 10.1109/tnsm.2023.3290920 | Full text not re-read; **title, authors, venue, pages verified via Crossref**; title confirms joint routing, channel, key-rate and time-slot assignment | Joint routing/channel/rate/time-slot assignment; NP-hardness proved for that problem |
| Hugues-Salas et al., JOCN 11(2) A209 (2019), DOI 10.1364/jocn.11.00a209 | Full text not re-read; title (“Monitoring and Physical-Layer Attack Mitigation in SDN-Controlled QKD Networks”) and authors verified via Crossref | Real-time QBER/optical-power monitoring with SDN path switching |
| Mehic et al., IEEE/ACM ToN 28(1) 168–181 (2020), DOI 10.1109/tnet.2019.2956079 | Full text not re-read; title (“A Novel Approach to QoS Provisioning in Trusted Relay QKD Networks”) and authors verified via Crossref | Key-pool/QoS link-state routing |
| Kumar, Kar, Basu, ICCE-Taiwan 2026, pp. 339–340, DOI 10.1109/icce-taiwan71481.2026.11652792 | **No full text** (2-page paper); title “Fault-Tolerant Routing Framework for Maximizing Secret Key Rates in Multi-Hop QKD Networks” verified via Crossref | cell content unverified; the title indicates key-rate maximisation and fault tolerance, not Bell-evidence security |
| Meng et al., Phys. Scr. 101(7) 075101 (2026), DOI 10.1088/1402-4896/ae4073 | **Abstract only** | security-level- and rate-aware multipath routing with XOR key fusion and SDN control; the remaining cells unverified |

**Statement for the response:** only Pirandola et al. was read in full for this revision. For Zhang, Hugues-Salas, Mehic, Kumar and Meng the bibliographic data were verified through Crossref (Meng: abstract also read); their full texts were not accessible to us, so the Table 1 cells for Kumar and Meng remain marked as unverified in the deposit. We do not claim that all six papers were reviewed (previous v5 statement withdrawn).

### Table 1 cell corrections (supported by the text of Pirandola et al.)
| Cell (Pirandola column) | Previous | Corrected |
|---|---|---|
| Complexity characterisation | "Modified Dijkstra" | "Run time O(\|E\| + \|P\| log \|P\|) for the Dijkstra variants; no hardness result" (an algorithm is not a complexity result) |
| Collusion resistance | "Relay trust assumed" | "All nodes are assumed trusted (static trusted-node network); no collusion model" — we make no claim about the authors' earlier work |
| Live monitoring and adaptation | "Static/planned" | "Static network (properties constant in time); no monitoring" |
| End-to-end accounting | "Yes (distributed budget)" | Kept and **acknowledged**: ETE ε is the sum along the route; the additive-budget *idea* is theirs |
| Security notion | ε-security with allocated budget | Kept (composable finite-size) |

### The two uncited works named by the reviewer (searched again with Crossref, which was reachable)
The reviewer cites, without references, (i) "security-aware QKD routing with real-time security/rate requirements" and (ii) "CHSH-based validation in secure path selection". Searches by Crossref bibliographic queries (`api.crossref.org`, this session):
* **(i) Found, with high confidence.** The abstract of Meng et al., *Phys. Scr.* **101**(7) 075101 (2026), DOI 10.1088/1402-4896/ae4073, reads: "…a fusion overhead-aware link cost function, which can intelligently calculate the optimal path based on the **real-time security level and rate requirements** of the service…". This is the same paper already in Table 1; we could not read the full text, so its Table 1 cells stay "unverified", but the abstract confirms the relevant feature (security-level- and rate-aware path cost, SDN-controlled, key fusion). A second candidate with "security-aware" in the title is Ding, Bi, Yang, Verma, *Chin. Phys. B* (2026), DOI 10.1088/1674-1056/ae838c (VQC link pre-screening + admission control; no Bell test in the abstract). Both are added to the related-work paragraph as *heuristic* security-aware routing without a device-independent statistical certificate.
* **(ii) Not located.** Queries combining CHSH / Bell test / link certification / path selection / routing returned no QKD-routing paper. The only Bell-based QKD network-routing-adjacent hits were device-independent QKD papers without routing (e.g. Woodhead–Acín–Pironio, *Quantum* 5, 443 (2021)). The reference list of Pirandola et al. (read in full) does not contain such a work either. **The response states that we could not locate it, asks the reviewer for the reference, and narrows the novelty claim to “Bell-evidence-verified, budget-constrained, collusion-floored path selection, with a hardness analysis”, which does not depend on the absence of any single paper.**

## Q13 — Scope of the security statements (text for §3.4, §4.2, §9, Appendix B)

1. **A1 scope.** H = ℓ·h(Ŝ − η) applies the Pironio et al. per-round min-entropy bound to all ℓ rounds, which presupposes a product structure (collective attacks, i.i.d. rounds). No smoothing term is accounted for and λ_EC = 0. **r_e is a certificate quantity, not a deployable DI-QKD key rate.**
2. **Randomness bound vs key rate (with numbers, `q13_scope_numbers.json`).** The bound used is the min-entropy rate h(S) = 1 − log₂(1 + √(2 − S²/4)) (Pironio et al., Nature 464, 1021 (2010), randomness certification). Asymptotic DI-QKD rates under collective attacks use a different expression (Acín et al., PRL 98, 230501 (2007); Pironio et al., New J. Phys. 11, 045021 (2009)): r_DI = 1 − h₂((1 + √((S/2)² − 1))/2) − h₂(Q), with Q = (1 − S/2√2)/2 for a Werner link:

| S | Q (Werner) | h (used for r_e, asymptotic) | r_DI asymptotic |
|---|---|---|---|
| 2.687 | 0.0250 | 0.4723 | 0.539 |
| 2.55 | 0.0492 | 0.3113 | 0.234 |
| 2.418 | 0.0726 | 0.2061 | **−0.010** |

   The two quantities are not ordered (h < r_DI at S₀, h ≫ r_DI at S₁). **Consequence, stated openly:** a compromised link at S₁ = 2.418 still has positive certified min-entropy (r_e > 0 at n = 10⁴) but no positive asymptotic DI-QKD rate; replacing r_e by r_DI would make the S₁ links unusable, which would strengthen, not weaken, BMR's avoidance of compromised links.
3. **Budget merge.** One line in Appendix B: key is Σε-secure with probability ≥ 1 − Σδ ⇒ unconditionally (Σε + Σδ)-secure (derivation in Q10). The sentence "Theorem 1 proves …" points to it.
4. **γ relation.** ε-secure ⇒ |P_guess − 2^(−k)| ≤ ε (guessing is a measurement), hence P_guess ≤ 2^(−k)(1 + γ) with γ = 2^k ε (not "|P_guess − 2^(−k)| ≤ γ"). With ε = ½·2^(−(H−k)/2), γ < 1 iff H > 3k − 2, i.e. φ ≲ 1/3, which is why public-XOR uses φ = 0.30 (log₂ γ = −127 at n = 2500, −172 at n = 10⁴) while sealed mode at φ ≥ 0.9 has γ ≫ 1 (non-informative, not required).
5. **Sealed-relay channel.** Each partial XOR is sent under authenticated encryption keyed by the hop key of the *next* link; the chain K = ⊕Kⱼ remains secret if at least one hop key is unknown to the adversary, so there is no circularity in the t-collusion model (a coalition must know *all* hop keys). Stated as assumption A4 (authenticated encryption under hop keys).
6. **Other limits for §9:** trusted PKI and authenticated certificate database at the destination; a lying link endpoint (Q11 T12); single hardware backend; baselines are re-implementations; Q15 assumptions.

## Q16 — Number-lock, regeneration, deposit, credentials

1. **Macros.** `experiments/generate_macros.py` reads the JSONs and writes `results_macros.tex` (105 lines, alpha-only LaTeX names verified, no hand-typed numbers).
2. **Figures regenerated / created:** `Q4_collusion`, `Q5_evidence_source`, `E7b_v2_corrected`, `Q9_solver_validation`, `Q15_delta_regime`, `Q6_sweep`, `Q6_dynamic`, `Q8_hardware_reconcile`, `Q14_twobudget` (+ unchanged E1–E6, E10). Superseded artefacts are in `superseded/` with a README.
3. **Deposit.** Code (`bmr/`, `experiments/`, `tests/`), seeds in each script, code hashes, raw hardware counts and job IDs (`q8_hw_repeats.json`, `e6v2_hardware_fixed.json`), tests, `README.md` with one command per table, `run_all.py` (never touches a QPU).
4. **Credentials — scan results.** (a) No hard-coded `IBM_API_KEY` in any `.py` file in the tree; (b) no file under `REVISION_OSN_12OCT2026` contains the *current* token from the local Qiskit account file; (c) no zip archive contains key patterns. **Limits of the scan:** the folder is not a git repository, so history, old branches and shared folders cannot be scanned from here, and earlier keys that may have been pasted into `results_r2.md` can only be checked by the authors; **the authors must confirm in the IBM Quantum dashboard that the earlier leaked key is revoked** and scrub any repository that ever contained it. Nothing in this deposit reads or writes credentials (the QPU scripts load the locally stored account only at run time).
5. **Table 2** lists every parameter (§1).

## Q17 — Attribution and citation corrections
- Appendix D: **Q-PATH = Li, Wang, Xue, Yu, Sun, Lu, IEEE Trans. Commun. (2022)** (manuscript ref. [13]); **Q-CAST = Shi & Qian, ACM SIGCOMM (2020)** (ref. [14]); previous "Shi et al." / "Sun et al." were wrong.
- References to add (blueprint §3.12): Pirandola S., Harney C., Spedalieri G., “Security-rate trade-off in QKD: From individual links to large-scale networks”, Phys. Rev. Research **8**, 033221 (2026), DOI 10.1103/7839-3xhy — **verified**; Kumar P., Kar B., Basu D., “Fault-Tolerant Routing Framework for Maximizing Secret Key Rates in Multi-Hop QKD Networks”, ICCE-Taiwan 2026, pp. 339–340, DOI 10.1109/icce-taiwan71481.2026.11652792 — **verified**; Meng F., Bi L., Verma G., Ban X., Zhang S., “Security-aware multi-path routing with XOR-based key fusion for QKD network”, Phys. Scr. **101**(7) 075101 (2026), DOI 10.1088/1402-4896/ae4073 — **verified**; Zhang Q., Ayoub O., Gatto A., Wu J., Musumeci F., Tornatore M., “Routing, Channel, Key-Rate, and Time-Slot Assignment for QKD in Optical Networks”, IEEE TNSM 21(1) 148–160 (2024), DOI 10.1109/tnsm.2023.3290920 — **author list and DOI verified via Crossref**; Hugues-Salas E., Ntavou F., Gkounis D., Kanellos G.T., Nejabati R., Simeonidou D., “Monitoring and Physical-Layer Attack Mitigation in SDN-Controlled QKD Networks”, J. Opt. Commun. Netw. **11**(2) A209 (2019), DOI 10.1364/jocn.11.00a209 — **verified**; Mehic M., Fazio P., Rass S., Maurhart O., Peev M., Poppe A., Rozhon J., Niemiec M., Voznak M., “A Novel Approach to QoS Provisioning in Trusted Relay QKD Networks”, IEEE/ACM Trans. Netw. **28**(1) 168–181 (2020), DOI 10.1109/tnet.2019.2956079 — **verified**; Krenn S., Mir O., Lorünser T., Ramacher S., Wohner F., “Topology-Hiding Path Validation for Large-Scale QKD Networks”, ACNS 2026, LNCS, pp. 371–396, DOI 10.1007/978-3-032-32578-5_14; arXiv:2604.01831 — **verified**; Cirigliano L., Brosco V., Castellano C., Conti C., Pilozzi L., “Optimal QKD networks: capacitance versus security”, npj Quantum Inf. **10**, 44 (2024), DOI 10.1038/s41534-024-00828-7 — **verified**; Chaudhary S., Calsi D.L., Choi J., Geitz M., Nötzel J., “Distance-Security Tradeoffs for Repeaterless End-to-End QKD Networks”, arXiv:2411.17547 (2024) — **verified**; Kiktenko E.O., Tayduganov A., Fedorov A.K., “Routing Algorithm Within Multiple Non-Overlapping Paths for QKD Networks”, Entropy **26**(12) 1102 (2024), DOI 10.3390/e26121102 — **verified**; Tomamichel et al. 2011; Renner 2005; Hoeffding 1963; Hassin 1992; Lorenz–Raz 2001; plus Pironio et al. 2010, Acín et al. 2007, Pironio et al. 2009, Garey–Johnson 1979.

---

# Part B — P1 items

## Q5 — Evidence-source experiment (Bell vs QBER)

`q5_evidence_source.py`: 20 000 links per row; δ = 10⁻⁶, ℓ = 10⁴, φ = 0.99; B8 accepts iff QBER̂ + radius < Q_th = 0.11 with the same finite-size treatment; BMR accepts iff `build_certificate` returns a certificate. Link types: honest (S = 2.687, QBER 0.024), depolarised (S = 2.418, QBER 0.06), **classical-mimic** (S = 2.0, QBER 0.024 — the hardest mimic for a Bell test, correlated non-entangled source).

| n | Link | BMR accepts % | B8 (QBER) accepts % | BMR mean rate | B8 believed rate |
|---|---|---|---|---|---|
| 500 | honest | 99.97 | 100.0 | 0.0912 | 0.0896 |
| 500 | depolarised | 23.0 | 4.4 | 0.0155 | 0.0896 |
| 500 | **classical-mimic** | **0.0** | **100.0** | 0.0 | 0.0896 |
| 10⁴ | honest | 100.0 | 100.0 | 0.3394 | 0.3389 |
| 10⁴ | depolarised | 100.0 | 100.0 | 0.1396 | 0.3389 |
| 10⁴ | **classical-mimic** | **0.0** | **100.0** | 0.0 | 0.3389 |

**Honesty clause (goes into the text).** The mimic result is true by construction: a source that fakes the QBER but cannot violate CHSH is accepted by a QBER monitor and rejected by a Bell test. It shows the role of the *evidence type* in a device-independent threat model, not a rate advantage; under trusted devices QBER and CHSH carry the same information in this model. Note that the QBER of the depolarised link was set to 0.06 (Werner would give 0.073); a higher value would only make B8 reject it more often, so the row is conservative against BMR. In the depolarised case at n = 500 the QBER monitor is actually stricter than the Bell test (4.4 % vs 23 % accepted), and at n = 10⁴ both accept — but BMR then delivers only the degraded rate (0.14), whereas B8 believes 0.34.

## Q11 — Algorithm 3 design, tests and timing

**Design (v2, `bmr/verify.py`).** (i) *Per-epoch certificate signatures*: for each link and epoch the designated attestor (the lower-id endpoint) signs (e, r, ε, δ, T) once when the certificate is published. (ii) *Session binding*: the route assembler signs H(ledger) ‖ nonce ‖ s ‖ d ‖ epoch, with the nonce chosen by the destination — so link endpoints need not sign per delivery. (iii) The database check covers (Tᵢ, now] for every entry, not just epoch Tᵢ. Ed25519 via `cryptography` 50.0.1.

**Tests** (2000 trials each; `q11_alg3_v2.json`):

| Test | Attack | Result (expected) | Reject reason |
|---|---|---|---|
| T0 | benign | accepted 100 % (accept) | — |
| T1, T1b | inflate r (also re-bind with forger key) | rejected 100 % | ledger-digest mismatch |
| T2 | stale certificate | rejected 100 % | stale |
| T3, T6 | drop a hop / reorder | rejected 100 % | continuity |
| T4 | loop | rejected 100 % | simplicity |
| T5 | lower ε | rejected 100 % | ledger-digest mismatch |
| T7, T7b | DEGRADED notice in (Tᵢ, now] (late / early entry) | rejected 100 % | database supersession |
| T8 | re-signed by non-attestor | rejected 100 % | certificate signature |
| **T9** | valid old ledger + session replayed for a new delivery | rejected 100 % | nonce mismatch (replay) |
| **T10, T10b** | session / ledger issued for another (s,d′) | rejected 100 % | bound to another (s,d) |
| **T11, T11b** | entry swapped for / path substituted by another valid ledger | rejected 100 % | ledger-digest mismatch |
| **T12** | lying attestor signs an inflated certificate with its *own valid key* | **ACCEPTED 100 %** (documented limit) | none — undetectable by the ledger |
| T12-DB | same, destination database stores the independent audit record | rejected 100 % | record mismatch (database cross-check) |

T12 is a *limit*, not a defect: a signature proves who said it, not that it is true. Only the audit database (or a second independent auditor) can catch a lying endpoint. This goes into §9.

**Timing** (Ed25519, `cryptography` 50.0.1, macOS arm64, Python 3.14.5; mean ± SD over repeated verifications):

| Hops | Total µs | Per entry µs |
|---|---|---|
| 2 | 461 ± 70 | 231 ± 35 |
| 4 | 777 ± 273 | 194 ± 68 |
| 8 | 1386 ± 139 | 173 ± 17 |
| 11 | 1855 ± 214 | 169 ± 19 |

HMAC fallback ≈ 3–4 µs per entry. The earlier E10 (v1) value was ≈ 149 µs per entry (median) and the older "37 µs" is unexplained; **the numbers above supersede both** (different design, machine and library).

## Q14 — E11 extension (maps to Theorem 4, not 3c)

`q14_twobudget_v2.py`: ER (mean degree 6) and BA (m = 2) graphs, n ∈ {10, 20, 30, 50, 100, 150}, uncorrelated and anti-correlated weights, 100 instances per cell (40 at n ≥ 100); diamond chains with tight budgets up to k = 14. Solvers: exact **Pareto labels**, exact **pseudo-polynomial DP**, **FPTAS-style** scaled DP (ε = 0.1), **Lagrangian heuristic** (Dijkstra on w₁ + λw₂, bisection on λ).

- **Random graphs:** Pareto and DP disagree on **0** instances in every cell; 0 Pareto timeouts; worst p95 = 1.1 ms; median ≤ 0.24 ms at n = 150.
- **FPTAS:** missed **0** feasible instances; in 12 of 182 exactly-infeasible ER-30 anti-correlated cases it answered YES, each certified by a path within (1 + ε)·B₁ (the expected relaxation).
- **Lagrangian:** missed **8 of 782** feasible instances (1.0 %), up to 20 % in a single cell (ER-150 anti-correlated), and answers "infeasible" on the diamond chain at every k. Our Lagrangian is conservative — it returns *no route* rather than a budget-violating route; the pending_v7 phrasing "violates budgets" is replaced by "misses feasible instances".
- **Diamond chain (worst case):** Pareto labels = 3·2^(k−1) − 1 (24 575 at k = 14, 3.6 s); DP cells 237 568 (123 ms); FPTAS cells 9 280 (4.7 ms). Exact and FPTAS solvers are correct; labels grow exponentially.

Message: hardness is a worst-case statement; typical random instances are solved in microseconds–milliseconds by Pareto/DP; the heuristic is unreliable. Figure: `figures/Q14_twobudget.png`.

## Q15 — Δ-sweep, regime map, robustness

### Regime map (ER n = 30, t = 1; feasibility % / mean rate; `q15_delta_regime.json`)
| n, φ | Δ = 10⁻⁶ | 3×10⁻⁶ | 10⁻⁵ | 3×10⁻⁵ | 10⁻³ | 10⁻² | 3×10⁻² | 5×10⁻² | 0.5 |
|---|---|---|---|---|---|---|---|---|---|
| 500, 0.99 | 0 | 0 | 0 | 0 | 0 | 1.0 / 0.001 | 27.4 / 0.029 | 66.0 / 0.066 | **99.8 / 0.092** |
| 500, 0.90 | 0 | 0.3 | **99.5 / 0.084** | 99.2 | 99.5 | 100 | 98.8 | 98.4 | — |
| 10⁴, 0.99 | 0 | 0 | 27.7 / 0.095 | **99.7 / 0.339** | 100 / 0.341 | 99.7 | 100 | 99.3 | — |
| 10⁴, 0.90 | 0 | 82.4 / 0.249 | 100 / 0.309 | 99.5 | 99.8 | 99.0 | 99.7 | 99.8 | — |

The **knee** answers "your Δ never binds": for φ = 0.99 feasibility is 0 below ≈0.01 at n = 500 and below ≈10⁻⁵ at n = 10⁴; BMR's operating Δ values lie right of the knee (plateau). For φ = 0.9 the defect is ≈10⁻⁶, so Δ ≥ 10⁻⁵ is enough at both n.

### Compromised-fraction sweep (φ = 0.99, f = 0 … 0.5)
- n = 10⁴, Δ = 10⁻³: feasibility 99.5 → 96.3 (f=.1) → 94.4 → 87.2 → 76.8 → **65.8 %** (f=.5); rate 0.338 → 0.326 → 0.318 → 0.291 → 0.255 → **0.217**.
- n = 500, Δ = 0.5: feasibility 99.8 → 97.6 → 94.5 → 92.7 → 89.2 → **78.9 %**; rate 0.092 → 0.085 → 0.078 → 0.071 → 0.062 → **0.050**.

### Drift / correlation robustness (`q15_drift_robustness.json`; 1000 replicates, 80 epochs, drift from epoch 20)
Model: a 5-hop route of identical links is accepted iff all certificates exist, r ≥ R_floor and 5·w ≤ Δ; **R_floor = 1st percentile and Δ = 99th percentile of the honest distribution** (commissioning on honest data); S(e) = S₀ − slope·max(0, e − 20); audit noise AR(1) with correlation ρ. Metrics: false-alarm rate (rejections before drift), detection delay (epochs after drift starts), budget-violation rate (accepted although the true certificate violates the budget/floor). Base: slope 0.004/epoch, ρ = 0.5.

| n | φ | False alarm % | Detection delay (epochs) | Violation while accepted % |
|---|---|---|---|---|
| 500 | 0.99 | 1.1 | 18.3 | 9.07 |
| 500 | 0.90 | 1.0 | 18.5 | 10.16 |
| 10⁴ | 0.99 | 1.2 | 6.8 | 2.43 |
| 10⁴ | 0.90 | 1.1 | 6.9 | 2.53 |

Fast drift (slope 0.02): delays 6.3–6.4 (n=500) and 2.2 (n=10⁴), violations ≤ 2.6 %. ρ ∈ {0, 0.9}: delay 15.6–28.3 (n=500), 6.1–8.5 (n=10⁴), false alarms 0.9–1.5 %. **±25 % sensitivity** (slope, ρ, S₀−2 gap): false alarms ≈ 1 % in every cell **except** the gap×0.75 variant at n = 500 (S₀ = 2.515), where false alarms are **25.5 %** because the honest link then certifies only marginally at n = 500 (and delay shrinks to 3.2). Interpretation: drift is detected late at n = 500 (≈18 epochs) and the commissioned budget does not protect during the detection window (≈9–10 % of post-drift accepted epochs violate it); the audit is a conservative certificate and an abrupt-change detector, not a sensitive drift detector (a CUSUM complement is not part of the model and not claimed). A first design with Δ = 1.5× nominal gave 38 % false alarms and was discarded as badly calibrated (documented in the script).

---

# Part C — Claim ledger and gate

## Claim ledger (what the manuscript may now say)
| Claim | Evidence | Allowed? |
|---|---|---|
| n = 100 cannot certify an honest link (in expectation); n_min = 235 | App. A, E1/E2, Q7 (0.7 % pass) | **Yes** (with "in expectation"; passes arise only from sampling fluctuation) |
| Both η bounds valid; failure rate ≤ δ in 200 000 trials per cell | E1 (fixed criterion) | **Yes** |
| Five certificate fields built from named inputs; 9 unit + 5 integration tests pass | E3, Q7 | **Yes** |
| Simulator and all network experiments use `build_certificate` | Q7 tests | **Yes** |
| Corollary 1 for simple paths; walks do not inherit it | proof, E4a (m ≤ 12), Q9 E4b | **Yes** |
| Walk solvers return non-simple paths in 28.7–99.6 % (all-instance basis) | Q9 | **Yes** |
| Floor on simple paths strongly NP-complete when t is an input | Q10 proof | **Yes** (after co-author re-derivation) |
| Exact solver validated (0 mismatches over 36 000 instances); runtime scope | Q9 | **Yes** (dense-graph numbers, timeouts shown; no path cap) |
| Algorithm 1 runs on real CHSH counts | Q8 | **Yes**, with run-to-run spread 2.54–2.70 and "monotone within shot noise" |
| Ledger rejects tampered/stale/looped/replayed/substituted ledgers | Q11 T1–T11 | **Yes**; lying endpoint = documented limit |
| BMR delivers keys only over compliant routes | definition + self-checks | **Yes**, "of delivered keys", with refusal rate |
| BMR matches strong rate-aware routers on rate | Q6 (budget-only −0.003) | **Yes**; no rate superiority claim |
| BMR resists t colluders; baselines do not | Cor. 1; Q4 | **Yes**: structural for ≤ t; probabilistic beyond |
| Bell evidence catches what QBER monitoring misses | Q5 | **Yes**, only as "in a device-independent threat model, by construction" |
| BMR improves compliance by ≈ 14.6 pt over the best rate-aware router (full sweep) | Q6 (92.8 vs 78.1 %) | **Yes**, with ring/sparse-BA refusal caveats |
| Dynamic recovery | Q6b | **Yes**, "driven by live certificates; B3-live is equivalent" |
| BMR "strictly guarantees 100 % compliance" | — | **No** (retracted) |
| "SR-Dijkstra = Pirandola's algorithm" | — | **No** (retracted) |
| "All six papers reviewed" | — | **No**: only Pirandola read in full |
| Budget Δ is a meaningful security level at n = 500, φ = 0.99 | Q2/Q15 | **No** for the headline Δ = 0.5 (formal guard only); **yes** at φ = 0.9, Δ = 10⁻⁴ (Q6 sensitivity sweep: same conclusions, BMR 0.1143 / 92.5 %) |
| r_e is a DI-QKD key rate | Q13 | **No** |

## Final gate
1. **All P0 items closed with an artefact:** yes (§0.1), except limits that cannot be removed from here and are stated: Kumar/Meng full texts not accessed and one of the two uncited works not located (Q12); git-history scan impossible and key revocation to be confirmed by the authors (Q16).
2. **Claim ledger has no "No" row asserted by the manuscript** provided §0.3 corrections are applied.
3. **Where is each comment addressed?** §0.2 maps C1–C5 to sections; a reader should find the answer in under a minute.
4. **Open items for the authors (cannot be closed from here):** (a) co-author re-derivation of the Q10 proofs; (b) revocation of the old IBM key in the dashboard and a scrub of any repository that ever held it (Q16); (c) full-text check of Kumar et al. and Meng et al. for the remaining Table 1 cells; (d) the reviewer's second uncited work (CHSH-based validation in secure path selection) — ask the reviewer for the reference. **Closed in this pass:** Garey–Johnson ND30 confirmed; Zhang et al. author list and all reference metadata verified via Crossref/arXiv; uncited work (i) identified; full sweep at φ = 0.9 executed.

---

## Appendix A — Radius derivations, design window, certificate construction (computed with `bmr/bounds.py`, `bmr/certificate.py`)

### A.1 The two radii
Setting: CHSH estimate Ŝ = E₁₁ + E₁₂ + E₂₁ − E₂₂ from n rounds per setting (4n rounds), outcome products in {±1}.
* **η_union = √(32 ln(8/δ)/n).** Each of the four correlators is a mean of n terms in {±1}; two-sided Hoeffding Pr[|Ê − E| ≥ t] ≤ 2exp(−n t²/2), set to δ/4 ⇒ t = √(2 ln(8/δ)/n). Ŝ is a sum of four correlators, so with a union bound η_S = 4t = √(32 ln(8/δ)/n). `bounds.eta_union`.
* **η_tight = √(8 ln(1/δ)/n).** Ŝ is an average-weighted sum of 4n independent terms, each in an interval of width 2/n; Hoeffding gives Pr[Ŝ ≤ S − η] ≤ exp(−2η²/(4n·(2/n)²)) = exp(−η² n/8), so η = √(8 ln(1/δ)/n). `bounds.eta_tight`. This is the radius used by `build_certificate`.
* The tight radius is smaller than the union radius by a factor √(4 ln(8/δ)/ln(1/δ)) ≈ 2.14 at δ = 10⁻⁶ (table below), which is why n_min differs by a factor ≈ 4.6 (235 vs 1078 at S = 2.687).

### A.2 Table T-A — design window at δ = 10⁻⁶, true S = 2.687 (computed)
| n per setting | η_tight | η_union | Ŝ − η_tight at Ŝ = 2.687 | h(Ŝ − η_tight) | certifies? |
|---:|---:|---:|---:|---:|---|
| 100 | 1.0513 | 2.2553 | 1.636 | — | **no** (≤ 2) |
| 235 | 0.6858 | 1.4712 | 2.0012 | 0.0004 | barely (n_min) |
| 500 | 0.4702 | 1.0086 | 2.2168 | 0.0906 | yes |
| 633 | 0.4179 | 0.8964 | 2.2691 | 0.1170 | yes (n_min for S = 2.418, where Ŝ − η = 2.000) |
| 1000 | 0.3325 | 0.7132 | 2.3545 | 0.1652 | yes |
| 2500 | 0.2103 | 0.4511 | 2.4767 | 0.2490 | yes |
| 10 000 | 0.1051 | 0.2255 | 2.5819 | 0.3424 | yes |

Minimum n (computed by `bounds.n_min`, δ = 10⁻⁶): **235** at S = 2.687 (tight); **379** at S = 2.5405; **633** at S = 2.418; **1078** at S = 2.687 with η_union. **n = 100 cannot certify even an ideal-device link** (Ŝ − η_tight = 1.64 < 2); the manuscript's n = 100 operating point must be replaced (D1/Q2).

### A.3 Ten-step certificate construction (answers Reviewer comment on parameter `c`)
| # | Step | Formula (code line in `build_certificate`) |
|---|---|---|
| 1 | Inputs | link, Ŝ, n, δ, ℓ, φ, λ_EC, epoch |
| 2 | Radius | η = √(8 ln(1/δ)/n) |
| 3 | Lower bound | S_lb = Ŝ − η; if S_lb ≤ 2 → **no certificate** (DEGRADED) |
| 4 | Min-entropy rate | h = h_Pironio(S_lb) |
| 5 | Certified entropy | H = ℓ·h − λ_EC; if H ≤ 0 → no certificate |
| 6 | Output length | k = ⌊φ·H⌋ |
| 7 | Rate | r = k/ℓ |
| 8 | Leftover-hash ε | log₂ ε = −1 − (H − k)/2 (stored in log domain; underflow-safe) |
| 9 | Composability | log₂ γ = k + log₂ ε |
| 10 | Budget weight | w = ε + δ; epoch tag T |

Passing criterion (the *only* one used anywhere in the simulator or the router): Ŝ − η_tight > 2. Verified by `tests/test_integration.py` (5/5) and `experiments/e3_certificate_tests.py` (9/9).

---

## Appendix D — Baseline definitions and attribution (as run in Q6)

| ID | Router | Definition actually executed | Attribution / caveat |
|---|---|---|---|
| B0 | SPF | hop-count shortest path | connectivity only |
| B1a | P_λ(0.01) | Pirandola et al. PRR 8, 033221 (2026), Eq. (35) | implemented as published, inputs (R_e, ε_e) = certified rate and ε |
| B1b | P_λ(1) | same protocol, λ = 1 | heuristic scalarisation, no hard budget |
| B2 | BMR-NoFloor / Pareto label-setting | exact optimum with budget + bottleneck, t = 0 | previously labelled "SR-Dijkstra" — renamed |
| B3 | Zhang-style | widest (max–min certified rate) path, ties → fewest hops | Zhang et al. also solve RWA/time slots; **dropped** here — a style, not a reproduction |
| B5 | Rate-greedy | Dijkstra on Σ 1/r_e | own definition |
| B6 | Q-PATH-style | Dijkstra on −log F_e | objective-only; **Q-PATH = Li, Wang, Xue, Yu, Sun, Lu, IEEE Trans. Commun. (2022)** |
| B7 | Q-CAST-style | min hops s.t. ∏F_e ≥ 0.75 | objective-only; **Q-CAST = Shi & Qian, ACM SIGCOMM (2020)** |
| B8 | QBER-monitor | prune QBER_ub ≥ 0.11, widest path on 1 − 2h₂(QBER_ub) | same Ŝ and δ → same information as BMR; QBER measured directly in Q5 |
| BMR | budget + hop floor (t = 1 ⇒ ≥ 3 hops) + bottleneck | exact simple-path search, validated in Q9 (36 000 instances, 0 mismatches) | |

The Pirandola multipath and security-optimised variants are *not* implemented (only the single-path λ-protocol, B1a/B1b) and no claim is made about them. The label B4 is unused.

## Appendix E — Table 1 (literature matrix): status
The cell corrections for the Pirandola column are listed under Q12. Kumar and Meng cells remain flagged "unverified" (papers not accessed). Hugues-Salas and Mehic cells are qualitative and taken from the reviewer's verified characterisation.

---

*End of R2 results v5-mrng (rebuilt).*
