"""
Q8 (part 2): hardware reconciliation.
=====================================
(a) Table of S_hat(phi=0) across ALL historical hardware runs + the new repeats.
(b) Sub-sampling pass-rate study recomputed from HARDWARE correlators
    (p_same = (1+E_ab)/2 at phi=0), per run, using the SAME certificate rule as the
    library: build_certificate(...) is not None  <=>  S_hat - eta_tight(n, delta) > 2.
    The old study used ideal-Aer counts; this one uses hardware numbers.
(c) Monotonicity of S(phi) on hardware: per-run, with shot-noise SE; the claim is
    "non-increasing within shot noise", NOT strict.
(d) Aer-noisy vs hardware gap.

All S_hat numbers are for the 4-setting CHSH sum with n=8192 shots/setting
(shot-noise SE of S ~= sqrt(4/8192) = 0.022).
"""
import sys, os, json, math, hashlib
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bmr.bounds import eta_tight
from bmr.certificate import build_certificate

ROOT = os.path.dirname(os.path.dirname(__file__))
RES = os.path.join(ROOT, "data", "results")
FIG = os.path.join(ROOT, "figures")
DELTA = 1e-6
SIGN = {(0, 0): +1, (0, 1): +1, (1, 0): +1, (1, 1): -1}
print("code sha256:", hashlib.sha256(open(__file__, "rb").read()).hexdigest()[:16])


def corr(counts):
    tot = sum(counts.values()); same = 0
    for bits, c in counts.items():
        if bits[-1] == bits[-2]:
            same += c
    return (2 * same - tot) / tot, tot


def S_from_E(E):
    return sum(SIGN[ab] * E[ab] for ab in SIGN)


# ---------------- load data ----------------
old = json.load(open(os.path.join(RES, "e6v2_hardware_fixed.json")))
runs = {}   # name -> {phi: {(a,b): E}}
meta = {}
for phi, v in old["hw_results"].items():
    runs.setdefault("marrakesh-run0 (E6v2, correlators only)", {})[float(phi)] = {
        eval(k): x for k, x in v["E_by_ab"].items()}
meta["marrakesh-run0 (E6v2, correlators only)"] = {"job": old["job_ids"]["main_batch_v2"], "shots": 8192}
rep_path = os.path.join(RES, "q8_hw_repeats.json")
nshots = 8192
if os.path.exists(rep_path):
    rp = json.load(open(rep_path))
    for i, r in enumerate(rp["runs"]):
        name = f"marrakesh-rep{i+1} (new, raw counts)"
        runs[name] = {}
        for key, cnt in r["raw_counts"].items():
            phi, a, b = key.split("|")
            E, tot = corr(cnt)
            runs[name].setdefault(float(phi), {})[(int(a), int(b))] = E
            nshots = tot
        meta[name] = {"job": r["job_id"], "shots": nshots}

# day-2 session (14 Oct 2026) repeats, same circuits / pair / shots, raw counts
rep2_path = os.path.join(RES, "q8_hw_repeats_14oct.json")
if os.path.exists(rep2_path):
    rp2 = json.load(open(rep2_path))
    for i, r in enumerate(rp2["runs"]):
        name = f"marrakesh-day2-rep{i+1} (14 Oct, raw counts)"
        runs[name] = {}
        for key, cnt in r["raw_counts"].items():
            phi, a, b = key.split("|")
            E, tot = corr(cnt)
            runs[name].setdefault(float(phi), {})[(int(a), int(b))] = E
        meta[name] = {"job": r["job_id"], "shots": tot}

# ---------------- (a) S_hat(0) table ----------------
se = math.sqrt(4.0 / nshots)
print("\n(a) S_hat(phi=0) per run (SE_shot = %.3f)" % se)
tabA = []
for name, d in runs.items():
    S0 = S_from_E(d[0.0])
    tabA.append({"run": name, "S0": S0, "job": meta[name]["job"]})
    print(f"  {name:45s} S0 = {S0:.4f}  job {meta[name]['job']}")
hist = [
    {"run": "earlier run (pre-fix, quoted in v5 text)", "S0": 2.776, "job": "n/a (not reproducible, raw data not retained)"},
    {"run": "failed ibm_kingston run (RY(phi) w/o factor 2, 2048 shots)", "S0": 1.2256, "job": "db6540kvf2bc73d0m1rg"},
]
S0s = [t["S0"] for t in tabA]
print(f"  mean over marrakesh runs = {np.mean(S0s):.4f}, SD = {np.std(S0s, ddof=1) if len(S0s)>1 else float('nan'):.4f}")

# ---------------- (b) sub-sampling from hardware ----------------
ns_list = [100, 250, 336, 379, 500, 700, 1000, 2048, 8192]
NR = 20000
rng = np.random.default_rng(2026)
tabB = {}
print("\n(b) sub-sampling pass rate from HARDWARE correlators at phi=0 (pass <=> S_hat - eta > 2)")
for name, d in runs.items():
    E = d[0.0]
    p = {ab: (1 + E[ab]) / 2 for ab in SIGN}
    row = {}
    for n in ns_list:
        S = np.zeros(NR)
        for ab, sg in SIGN.items():
            k = rng.binomial(n, p[ab], NR)
            S += sg * (2 * k - n) / n
        row[n] = float(np.mean(S - eta_tight(n, DELTA) > 2.0))
    tabB[name] = row
    print(f"  {name:45s}", " ".join(f"{n}:{row[n]:.3f}" for n in ns_list))
# cross-check equivalence with library rule on one run
name0 = list(runs)[0]
E = runs[name0][0.0]; p = {ab: (1 + E[ab]) / 2 for ab in SIGN}
chk = 0
for _ in range(500):
    n = 336
    S = sum(SIGN[ab] * (2 * rng.binomial(n, p[ab]) - n) / n for ab in SIGN)
    c = build_certificate("x", S, n, DELTA, 10_000, 0.99, 0, 1)
    chk += ((c is not None) == (S - eta_tight(n, DELTA) > 2.0))
print("  library-rule agreement (500 draws, n=336):", chk, "/ 500")
# old (ideal-Aer) numbers for contrast
print("  old ideal-Aer study:", old["sub_sampling_pass_rates"])

# ---------------- (c) monotonicity ----------------
print("\n(c) S(phi) per run (monotone within shot noise: S_i >= S_{i+1} - 2*SE*sqrt2)")
tabC = {}
for name, d in runs.items():
    ph_here = sorted(d.keys())
    seq = [S_from_E(d[ph]) for ph in ph_here]
    tol = 2 * se * math.sqrt(2)
    mono_strict = all(seq[i] > seq[i + 1] for i in range(len(seq) - 1))
    mono_noise = all(seq[i] >= seq[i + 1] - tol for i in range(len(seq) - 1))
    tabC[name] = {"phis": ph_here, "S": seq, "strict": mono_strict, "within_noise": mono_noise}
    print(f"  {name:45s}", [round(s, 3) for s in seq], f"strict={mono_strict} within-noise={mono_noise}")

# ---------------- (d) Aer gap ----------------
an = old["aer_noisy_results"]["0.0"]["S"]; ai = old["aer_ideal_results"]["0.0"]["S"]
print(f"\n(d) Aer ideal S(0) = {ai:.3f}; Aer noisy S(0) = {an:.3f}; hardware mean = {np.mean(S0s):.3f}; "
      f"gap (Aer-noisy - hw) = {an-np.mean(S0s):.3f}")
print("    Aer noise model: depolarizing 1e-3 (1q), 3e-3 (2q), readout 0.3%. Pair readout error at submission "
      "(q6/q5: 13 Oct 0.27%/0.78%; 14 Oct session 1 1.39%/2.87%; session 2 1.04%/2.19%; see the 'calibration' "
      "fields of q8_hw_repeats*.json) varies by session. Day-2 S is HIGHER despite the higher readout error, so "
      "readout error does not explain the gap or the drift; the remaining gap is attributed to unmodelled "
      "two-qubit gate error/crosstalk/drift (not measured here).")

# ---------------- figure ----------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for name, d in runs.items():
    phs = sorted(d.keys()); ax[0].plot(phs, [S_from_E(d[ph]) for ph in phs], "o-", label=name.split(" (")[0])
ax[0].axhline(2, color="k", ls="--", lw=.8); ax[0].set_xlabel(r"$\varphi$ (half Bloch rotation)")
ax[0].set_ylabel(r"$\hat S$"); ax[0].set_title("Hardware S(phi), ibm_marrakesh (6,5), n=8192/setting"); ax[0].legend(fontsize=7)
for name, row in tabB.items():
    ax[1].plot(list(row), list(row.values()), "o-", label=name.split(" (")[0])
ax[1].set_xscale("log"); ax[1].set_xlabel("n per setting"); ax[1].set_ylabel("pass rate")
ax[1].set_title("Sub-sampling pass rate (DATA: hardware correlators at phi=0)"); ax[1].legend(fontsize=7)
ax[2].bar(range(len(S0s)), S0s, color="C0")
ax[2].axhline(an, color="C1", ls="--", label="Aer noisy (optimistic model)")
ax[2].axhline(ai, color="C2", ls=":", label="Aer ideal"); ax[2].axhline(2, color="k", ls="--", lw=.8)
ax[2].set_xticks(range(len(S0s))); ax[2].set_xticklabels([f"run{i}" for i in range(len(S0s))])
ax[2].set_ylim(1.8, 3.0); ax[2].set_title(r"$\hat S(0)$ per hardware run"); ax[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "Q8_hardware_reconcile.png"), dpi=140)

json.dump({"S0_table": tabA, "historical": hist, "S0_mean": float(np.mean(S0s)),
           "S0_sd": float(np.std(S0s, ddof=1)) if len(S0s) > 1 else None, "shot_noise_se": se,
           "subsampling_hw": tabB, "subsampling_old_ideal_aer": old["sub_sampling_pass_rates"],
           "monotonicity": tabC, "aer_noisy_S0": an, "aer_ideal_S0": ai, "delta": DELTA,
           "data_source_note": "hardware correlators (E6v2 run: E_ab only; new reps: raw counts)"},
          open(os.path.join(RES, "q8_hardware_reconcile.json"), "w"), indent=1)
print("saved.")
