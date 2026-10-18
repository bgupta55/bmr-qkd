"""R7_figure.py -- regenerate the hardware figure from the raw counts (data/06_ibm_chsh_results.json) + R7 Aer results.
Fixed (A+.7): legend moved outside axes, y-axis from 0, right-panel title mentions the FakeFez noisy check."""
import json, math, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

hw = json.load(open("../data/R7_hardware_recomputed.json"))
aer = json.load(open("../data/R7_aer_noiseless.json"))

fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))

# ── Left panel: hardware bars ────────────────────────────────────────────────
ax = axes[0]
labs = ["native\n(ibm_fez)", "degraded\n(ibm_fez)"]
S = [hw["native"]["S"], hw["degraded"]["S"]]

ax.bar([0, 1], S,
       yerr=[hw["native"]["se_physical"], hw["degraded"]["se_physical"]],
       capsize=5, color=["C0", "C1"], alpha=.85)
ax.errorbar([0, 1], S,
            yerr=[hw["native"]["se_conservative"], hw["degraded"]["se_conservative"]],
            fmt="none", ecolor="k", capsize=9, lw=1, label="±1 SE (conservative)")
for i, v in enumerate(["corrected", "degraded"]):
    ax.plot([i - .4, i + .4], [aer[v]["S_mean"]] * 2, "k--", lw=1.2,
            label="noiseless Aer (mean of 20)" if i == 0 else None)

ax.axhline(2, color="r", lw=.8, label="classical bound S=2")
ax.axhline(2 * math.sqrt(2), color="gray", lw=.8, ls=":", label="quantum max 2√2")
ax.set_xticks([0, 1]); ax.set_xticklabels(labs)
ax.set_ylim(0, 3.0)            # y from 0 (not 1.9) to avoid exaggeration
ax.set_ylabel("CHSH Ŝ  (8 192 shots/setting)")
ax.set_title(
    f"Hardware vs noiseless Aer\n"
    f"Δ(deg−nat) = {hw['difference_degraded_minus_native']['diff']:+.3f}, "
    f"z = {hw['difference_degraded_minus_native']['z_physical']:.2f}",
    fontsize=9)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, framealpha=.85)

# ── Right panel: noiseless candidates ───────────────────────────────────────
ax = axes[1]
names = ["corrected", "degraded", "opposite_sign", "old_halfA", "old_halfboth", "old_halfB"]
m    = [aer[v]["S_mean"] for v in names]
sd   = [aer[v]["S_sd"]   for v in names]
ideal= [aer[v]["S_ideal"] for v in names]

ax.errorbar(range(6), m, yerr=sd, fmt="o", capsize=4, label="Aer noiseless mean ± SD (20×8 192)")
ax.plot(range(6), ideal, "kx", ms=8, label="exact statevector")
ax.axhline(2.574, color="C3", ls="--", lw=.9, label="old reported 2.574 / 2.579")
ax.axhline(2.579, color="C3", ls="--", lw=.9)

ax.set_xticks(range(6))
ax.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=7)
ax.set_ylabel("CHSH Ŝ  (noiseless)")
ax.set_title(
    "No noiseless or FakeFez-noisy candidate\nreproduces 2.574",
    fontsize=9)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, framealpha=.85)
ax.grid(alpha=.3)

plt.tight_layout()
plt.savefig("../figures/fig_R7_chsh.png", dpi=150, bbox_inches="tight")
print("Saved fig_R7_chsh.png")
