"""
R9_robustness.py -- R9: robustness of the audit at n = 10^4.
 (a) linear drift of a single link S: 2.687 -> 2.60 over 60 epochs: first-alarm epoch distribution (binomial noise), certificate soundness
     (S_lb <= S_true and r_e <= h(S_true)), epochs run on a drifted-but-not-alarmed link.
 (b) correlated noise (Gaussian approximation with the exact binomial sigma_W): 8 links, within-epoch common shock with correlation rho,
     and AR(1) temporal correlation rho; per-link false-alarm rate and alarm clustering at alpha=1e-6 and at a looser alpha=1e-3 (statistics).
Usage: python3 R9_robustness.py
"""
import json, math, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B
from scipy.stats import norm

N = 10_000; S0 = 2.687; EP = 60; TRIALS = 20_000
rng = np.random.default_rng(20261012)
out = {}
sigW0 = float(B.sigma_w(S0, N)); thr = S0 - B.Z_ALPHA * sigW0
print("sigma_W(S0)", sigW0, "alarm threshold", thr)

# ---------------- (a) linear drift
S_true = S0 + (2.60 - S0) * (np.arange(EP) / (EP - 1))            # epoch 1 .. 60
first = np.full(TRIALS, -1); sound = 0; tot = 0; viol_r = 0
unalarmed_drifted = []; Slb_all = []
for t in range(TRIALS // 1):
    p = (1 + S_true / 4) / 2
    K = rng.binomial(4 * N, p); Sh = 2 * K / N - 4
    al = Sh < thr
    if al.any(): first[t] = int(np.argmax(al)) + 1
    Slb = np.maximum(0, Sh - B.tau(N)); sound += int((Slb <= S_true + 1e-12).sum()); tot += EP
    viol_r += int((B.h_vec(Slb) > B.h_vec(S_true) + 1e-12).sum())
    pre = np.arange(1, EP + 1) < (first[t] if first[t] > 0 else EP + 1)
    unalarmed_drifted.append(int((pre & (S_true < S0 - 0.03)).sum()))
det = first[first > 0]
S_cross = S0 - (S0 - thr)
ep_mean_cross = 1 + (S0 - thr) / (S0 - 2.60) * (EP - 1)
out["drift"] = dict(trials=TRIALS, threshold=thr, mean_crossing_epoch=float(ep_mean_cross), detected_frac=float((first > 0).mean()),
                    first_alarm_epoch_mean=float(det.mean()), first_alarm_epoch_median=float(np.median(det)),
                    first_alarm_q05_q95=[float(np.percentile(det, 5)), float(np.percentile(det, 95))],
                    hist={int(k): int(v) for k, v in zip(*np.unique(det, return_counts=True))},
                    cert_sound_frac=sound / tot, rate_overshoot_events=viol_r,
                    mean_epochs_with_drift_gt_0p03_before_alarm=float(np.mean(unalarmed_drifted)),
                    S_true_at_epoch=dict(e1=float(S_true[0]), e30=float(S_true[29]), e47=float(S_true[46]), e60=float(S_true[-1])))
print("drift", {k: v for k, v in out["drift"].items() if k != "hist"})

# ---------------- (b) correlated noise
def sim_corr(kind, rho, trials, links=8):
    sig = float(B.sigma_w(S0, N))
    res = {}
    for alpha, z in ((1e-6, B.Z_ALPHA), (1e-3, norm.isf(1e-3))):
        res[alpha] = dict(z=float(z), alarms=0, decisions=0, epochs_ge1=0, epochs_ge2=0, epochs=0)
    for _ in range(trials):
        if kind == "common":
            c = rng.standard_normal((EP, 1)); e = rng.standard_normal((EP, links))
            Z = math.sqrt(rho) * c + math.sqrt(1 - rho) * e
        elif kind == "ar1":
            Z = np.empty((EP, links)); Z[0] = rng.standard_normal(links)
            for t in range(1, EP): Z[t] = rho * Z[t - 1] + math.sqrt(1 - rho ** 2) * rng.standard_normal(links)
        else:
            Z = rng.standard_normal((EP, links))
        for alpha, r in res.items():
            a = Z < -r["z"]
            r["alarms"] += int(a.sum()); r["decisions"] += a.size
            k = a.sum(1); r["epochs_ge1"] += int((k >= 1).sum()); r["epochs_ge2"] += int((k >= 2).sum()); r["epochs"] += EP
    return {str(a): dict(z=r["z"], rate=r["alarms"] / r["decisions"], alarms=r["alarms"], decisions=r["decisions"],
                         P_ge1_per_epoch=r["epochs_ge1"] / r["epochs"], P_ge2_per_epoch=r["epochs_ge2"] / r["epochs"],
                         P_ge2_given_ge1=(r["epochs_ge2"] / r["epochs_ge1"] if r["epochs_ge1"] else None)) for a, r in res.items()}
out["corr"] = {}
for kind, rho, tr in (("indep", 0.0, 40_000), ("common", 0.3, 40_000), ("ar1", 0.3, 40_000)):
    out["corr"][f"{kind}_rho{rho}"] = sim_corr(kind, rho, tr)
    for a, v in out["corr"][f"{kind}_rho{rho}"].items():
        print(kind, rho, "alpha", a, {k: (round(x, 6) if isinstance(x, float) else x) for k, x in v.items()})
json.dump(out, open("../data/R9_robustness.json", "w"), indent=1)
