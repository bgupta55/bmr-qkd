### T-R4-A  Detection and false alarms (500 trials × 60 epochs, n = 10⁴, binomial noise, alarm at α = 10⁻⁶)

| quantity | value |
|---|---|
| clean link-epoch decisions (8 links × 24 pre-attack epochs + 7 links × 36 post-attack epochs, × 500) | 222,000 |
| expected false alarms at α = 10⁻⁶ | 0.222 |
| **observed false alarms** (pre-attack / non-attacked links after the attack) | **0 / 0** |
| attacked link S† = 2.418: mean Ŝ ± SD over all attacked link-epochs | 2.4180 ± 0.0160 |
| S_lb mean [min, max] | 2.1925 [2.1301, 2.2579] |
| **gate (S_lb ≥ 2) pass rate** | **1.000** (the gate does not catch this attack) |
| r_e of the attacked link (mean / max) | 0.0791 / 0.1111  (< R = 0.10) |
| alarm fires in the attack epoch (epoch 25) | 500/500; undetected 0 |
| detection-epoch histogram | {'25': 500} |

### T-R4-B  Router performance, Δ=5e-6, deterministic tie-break (arc A first)

Feasibility, R_nom|feasible and R_sec (0 if infeasible or if the selected path contains an attacked link / violates Δ) — epochs 1–24 vs 26–60; attack epoch 25 shown separately. “on A” = share of feasible selections on the attacked arc 0-1-2-3-4.

| Router | pre feas | pre R_nom\|f | pre R_sec | ep25 R_sec | post feas | post R_nom\|f | post R_sec [95% CI over trials] | post on A | alarm-on-path (post) |
|---|---|---|---|---|---|---|---|---|---|
| BMR | 1.000 | 0.2304 | 0.2304 | 0.2260 | 1.000 | 0.2261 | 0.2261 [0.2260, 0.2262] | 0.00 | 0.000 |
| SPF | 1.000 | 0.2262 | 0.2262 | 0.0000 | 1.000 | 0.0790 | 0.0000 [0.0000, 0.0000] | 1.00 | 1.000 |
| SPF_valid | 1.000 | 0.2262 | 0.2262 | 0.2260 | 1.000 | 0.2261 | 0.2261 [0.2260, 0.2262] | 0.00 | 0.000 |
| SPF_budget | 1.000 | 0.2262 | 0.2262 | 0.0000 | 1.000 | 0.0790 | 0.0000 [0.0000, 0.0000] | 1.00 | 1.000 |
| SPF_valid_budget | 1.000 | 0.2262 | 0.2262 | 0.2260 | 1.000 | 0.2261 | 0.2261 [0.2260, 0.2262] | 0.00 | 0.000 |
| BMR_hop | 1.000 | 0.2304 | 0.2304 | 0.2260 | 1.000 | 0.2261 | 0.2261 [0.2260, 0.2262] | 0.00 | 0.000 |
| QPATHstyle | 1.000 | 0.2292 | 0.2292 | 0.2260 | 1.000 | 0.2261 | 0.2261 [0.2260, 0.2262] | 0.00 | 0.000 |
| QLEAPstyle | 1.000 | 0.2292 | 0.2292 | 0.2260 | 1.000 | 0.2261 | 0.2261 [0.2260, 0.2262] | 0.00 | 0.000 |
| QCASTstyle | 1.000 | 0.2292 | 0.2292 | 0.2260 | 1.000 | 0.2261 | 0.2261 [0.2260, 0.2262] | 0.00 | 0.000 |

### T-R4-C  Router performance, Δ=5e-6, random tie-break

Feasibility, R_nom|feasible and R_sec (0 if infeasible or if the selected path contains an attacked link / violates Δ) — epochs 1–24 vs 26–60; attack epoch 25 shown separately. “on A” = share of feasible selections on the attacked arc 0-1-2-3-4.

| Router | pre feas | pre R_nom\|f | pre R_sec | ep25 R_sec | post feas | post R_nom\|f | post R_sec [95% CI over trials] | post on A | alarm-on-path (post) |
|---|---|---|---|---|---|---|---|---|---|
| BMR | 1.000 | 0.2304 | 0.2304 | 0.2264 | 1.000 | 0.2260 | 0.2260 [0.2259, 0.2261] | 0.00 | 0.000 |
| SPF | 1.000 | 0.2262 | 0.2262 | 0.1096 | 1.000 | 0.1525 | 0.1129 [0.1113, 0.1145] | 0.50 | 0.500 |
| SPF_valid | 1.000 | 0.2261 | 0.2261 | 0.2264 | 1.000 | 0.2260 | 0.2260 [0.2259, 0.2261] | 0.00 | 0.000 |
| SPF_budget | 1.000 | 0.2261 | 0.2261 | 0.1152 | 1.000 | 0.1536 | 0.1146 [0.1129, 0.1164] | 0.49 | 0.493 |
| SPF_valid_budget | 1.000 | 0.2261 | 0.2261 | 0.2264 | 1.000 | 0.2260 | 0.2260 [0.2259, 0.2261] | 0.00 | 0.000 |
| BMR_hop | 1.000 | 0.2304 | 0.2304 | 0.2264 | 1.000 | 0.2260 | 0.2260 [0.2259, 0.2261] | 0.00 | 0.000 |
| QPATHstyle | 1.000 | 0.2292 | 0.2292 | 0.2264 | 1.000 | 0.2260 | 0.2260 [0.2259, 0.2261] | 0.00 | 0.000 |
| QLEAPstyle | 1.000 | 0.2292 | 0.2292 | 0.2264 | 1.000 | 0.2260 | 0.2260 [0.2259, 0.2261] | 0.00 | 0.000 |
| QCASTstyle | 1.000 | 0.2292 | 0.2292 | 0.2264 | 1.000 | 0.2260 | 0.2260 [0.2259, 0.2261] | 0.00 | 0.000 |

### T-R4-D  Router performance, Δ=1e-6 (tight), deterministic tie-break

Feasibility, R_nom|feasible and R_sec (0 if infeasible or if the selected path contains an attacked link / violates Δ) — epochs 1–24 vs 26–60; attack epoch 25 shown separately. “on A” = share of feasible selections on the attacked arc 0-1-2-3-4.

| Router | pre feas | pre R_nom\|f | pre R_sec | ep25 R_sec | post feas | post R_nom\|f | post R_sec [95% CI over trials] | post on A | alarm-on-path (post) |
|---|---|---|---|---|---|---|---|---|---|
| BMR | 0.000 | 0.0000 | 0.0000 | 0.0000 | 0.000 | 0.0000 | 0.0000 [0.0000, 0.0000] | 0.00 | 0.000 |
| SPF | 1.000 | 0.2261 | 0.0000 | 0.0000 | 1.000 | 0.0792 | 0.0000 [0.0000, 0.0000] | 1.00 | 1.000 |
| SPF_valid | 1.000 | 0.2261 | 0.0000 | 0.0000 | 1.000 | 0.2262 | 0.0000 [0.0000, 0.0000] | 0.00 | 0.000 |
| SPF_budget | 0.000 | 0.0000 | 0.0000 | 0.0000 | 0.000 | 0.0000 | 0.0000 [0.0000, 0.0000] | 0.00 | 0.000 |
| SPF_valid_budget | 0.000 | 0.0000 | 0.0000 | 0.0000 | 0.000 | 0.0000 | 0.0000 [0.0000, 0.0000] | 0.00 | 0.000 |
| BMR_hop | 0.000 | 0.0000 | 0.0000 | 0.0000 | 0.000 | 0.0000 | 0.0000 [0.0000, 0.0000] | 0.00 | 0.000 |
| QPATHstyle | 1.000 | 0.2292 | 0.0000 | 0.0000 | 1.000 | 0.2262 | 0.0000 [0.0000, 0.0000] | 0.00 | 0.000 |
| QLEAPstyle | 1.000 | 0.2292 | 0.0000 | 0.0000 | 1.000 | 0.2262 | 0.0000 [0.0000, 0.0000] | 0.00 | 0.000 |
| QCASTstyle | 1.000 | 0.2292 | 0.0000 | 0.0000 | 1.000 | 0.2262 | 0.0000 [0.0000, 0.0000] | 0.00 | 0.000 |

### T-R4-E  Recovery after the attack (first epoch ≥ 25 with R_sec > 0), Δ = 5×10⁻⁶, deterministic tie-break vs random tie-break

| Router | det: recovered | det: mean first secure epoch | random: recovered | random: mean first secure epoch |
|---|---|---|---|---|
| BMR | 500/500 | 25.00 | 500/500 | 25.00 |
| SPF | 0/500 | – | 500/500 | 25.95 |
| SPF_valid | 500/500 | 25.00 | 500/500 | 25.00 |
| SPF_budget | 0/500 | – | 500/500 | 25.99 |
| SPF_valid_budget | 500/500 | 25.00 | 500/500 | 25.00 |
| BMR_hop | 500/500 | 25.00 | 500/500 | 25.00 |
| QPATHstyle | 500/500 | 25.00 | 500/500 | 25.00 |
| QLEAPstyle | 500/500 | 25.00 | 500/500 | 25.00 |
| QCASTstyle | 500/500 | 25.00 | 500/500 | 25.00 |
