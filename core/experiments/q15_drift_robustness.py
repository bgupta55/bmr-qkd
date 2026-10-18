"""
Q15 (remainder): drift / correlation robustness of the certificate pipeline, +-25 % sensitivity.

Model (per link, 5-hop route of identical links, route accepted iff all of):
   certificate exists, r_cert >= R_floor, 5*w_cert <= Delta.
   R_floor / Delta are COMMISSIONED on the honest distribution: R_floor = 1st percentile of r, Delta = 99th
   percentile of 5*w over 20000 honest audits (None -> r=0, w=inf).  A first design with Delta=1.5x nominal
   gave 38 % false alarms at n=500 because w depends steeply on S_hat; it was discarded as badly calibrated.
Truth: w_true, r_true = build_certificate evaluated at S_hat = S_true (the noiseless statistic).
Drift: S_true(e) = S0 - slope*max(0, e-E0), E0=20, 80 epochs.  Audit noise AR(1) with correlation rho
across epochs (marginal sd = sqrt((4-S^2/4)/n)).
Metrics (per replicate, then averaged over 1000 replicates):
   false-alarm rate : fraction of PRE-drift epochs (e<E0) where the route is rejected
   detection delay  : first epoch >= E0 where the route is rejected, minus E0 (censored at 80-E0 = 60 -> reported with fraction detected)
   budget violation : fraction of POST-drift epochs where the route is ACCEPTED but truth violates
                      (r_true < R_floor or 5*w_true > Delta or no certificate at S_true)
Sensitivity (+-25 %): slope in {0.75,1,1.25} x base; rho in {0.375,0.5,0.625}; S0-2 gap in {0.75,1,1.25} x 0.687
(upper capped at the Tsirelson bound 2.828).
"""
import sys, os, json, math, hashlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from bmr.certificate import build_certificate
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print("code sha256:", hashlib.sha256(open(__file__, "rb").read()).hexdigest()[:16])
rng = np.random.default_rng(1515)
NE, E0, REPS, H = 80, 20, 1000, 5
DELTA, ELL = 1e-6, 10_000
TS = 2 * math.sqrt(2)


def eval_cert(S, n, phi, ep):
    c = build_certificate("l", min(S, TS), n, DELTA, ELL, phi, 0, ep)
    return c


def run(n, phi, S0, slope, rho, reps=REPS):
    sd0 = math.sqrt(max(0.0, 4 - S0 * S0 / 4) / n)
    rr, ww = [], []
    for x in rng.normal(S0, sd0, 20000):
        c = eval_cert(x, n, phi, 0)
        rr.append(c.r if c else 0.0); ww.append(H * c.w if c else math.inf)
    pass_rate = float(np.mean(np.array(rr) > 0))
    Rf, Dl = float(np.percentile(rr, 1)), float(np.percentile(ww, 99))
    if pass_rate < 0.5:
        return dict(fa=1 - pass_rate, delay=0.0, detected=1.0, viol=0.0, nominal_ok=False, honest_pass=pass_rate)
    # truth pre-computation per epoch
    truth_bad = np.zeros(NE, bool); S_true = np.array([S0 - slope * max(0, e - E0) for e in range(NE)])
    for e in range(NE):
        ct = eval_cert(S_true[e], n, phi, e)
        truth_bad[e] = (ct is None) or ct.r < Rf or H * ct.w > Dl
    fa, delay, det, viol = [], [], [], []
    sd = np.array([math.sqrt(max(0.0, 4 - s * s / 4) / n) for s in S_true])
    for _ in range(reps):
        z = np.empty(NE); z[0] = rng.normal()
        for e in range(1, NE):
            z[e] = rho * z[e - 1] + math.sqrt(1 - rho * rho) * rng.normal()
        rej = np.empty(NE, bool)
        for e in range(NE):
            c = eval_cert(S_true[e] + sd[e] * z[e], n, phi, e)
            rej[e] = (c is None) or c.r < Rf or H * c.w > Dl
        fa.append(rej[:E0].mean())
        post = rej[E0:]
        idx = np.argmax(post) if post.any() else None
        det.append(post.any()); delay.append(idx if post.any() else NE - E0)
        viol.append(np.mean((~rej[E0:]) & truth_bad[E0:]))
    return dict(fa=float(np.mean(fa)), delay=float(np.mean(delay)), detected=float(np.mean(det)),
                viol=float(np.mean(viol)), nominal_ok=True, honest_pass=pass_rate, Rfloor=Rf, Delta=Dl,
                first_truth_bad=int(np.argmax(truth_bad[E0:])) if truth_bad[E0:].any() else None)


rows = []
base = dict(slope=0.004, rho=0.5, gap=0.687)
print(f"{'n':>6} {'phi':>5} | {'variant':26s} | FA%    delay(ep)  detected%  viol%   [first epoch truth is bad]")
cells = [(500, .99), (500, .90), (10_000, .99), (10_000, .90)]
variants = [("base", {})]
for f in (0.75, 1.25):
    variants += [(f"slope x{f}", dict(slope=base["slope"] * f)), (f"rho {0.5*f:.3f}", dict(rho=0.5 * f)),
                 (f"S0-2 gap x{f}", dict(gap=base["gap"] * f))]
variants += [("rho=0 (iid)", dict(rho=0.0)), ("rho=0.9", dict(rho=0.9)), ("slope 0.02 (fast)", dict(slope=0.02))]
for (n, phi) in cells:
    for name, ov in variants:
        p = {**base, **ov}
        S0 = min(2 + p["gap"], TS - 1e-9)
        r = run(n, phi, S0, p["slope"], p["rho"])
        r.update(n=n, phi=phi, variant=name, slope=p["slope"], rho=p["rho"], S0=S0)
        rows.append(r)
        print(f"{n:>6} {phi:>5} | {name:26s} | {100*r['fa']:5.1f}  {r['delay']:8.1f}  {100*r['detected']:8.1f}  {100*r['viol']:5.2f}   [{r.get('first_truth_bad')}]  S0={S0:.3f}"
              + ("" if r["nominal_ok"] else "  (honest link certifies <50%: FA is structural)"), flush=True)
json.dump({"rows": rows, "E0": E0, "epochs": NE, "reps": REPS, "H": H, "base": base}, open(os.path.join(ROOT, "data", "results", "q15_drift_robustness.json"), "w"), indent=1)
