"""
Q5 — Evidence-source experiment (P1): Bell audit vs QBER monitoring
===================================================================
Link types
  honest        S=2.687, QBER 0.024
  depolarised   S=2.418, QBER 0.060 (elevated)
  classical-mimic  correlated non-entangled source: S = 2.0 (the *strongest* classical value,
                   i.e. the hardest mimic for the Bell test), reported QBER = honest 0.024.
BMR: build_certificate (same finite-size treatment: eta_tight, delta=1e-6).
B8 : QBER-monitor, accept if QBER_hat < Q_th, where QBER_hat has the SAME finite-size treatment:
     Hoeffding one-sided radius sqrt(ln(1/delta)/(2 n_q)) with n_q = 4n sifted rounds,
     accept iff QBER_hat + radius < Q_th  (Q_th = 0.11, the BB84-style tolerable error;
     identical statistical confidence 1e-6 as BMR).
Both evaluated at n = 500 and n = 10^4.
Honesty clause: true by construction. It shows the role of the evidence type in a
device-independent threat model, not a rate advantage; under trusted devices QBER and CHSH
carry the same information in this model.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bmr.certificate import build_certificate

SEED = 707
DELTA = 1e-6; PHI = 0.99; ELL = 10_000
N_LINKS = 20_000
Q_TH = 0.11
TYPES = [("honest", 2.687, 0.024), ("depolarised", 2.418, 0.060), ("classical-mimic", 2.000, 0.024)]
rng = np.random.default_rng(SEED)

res = []
print(f"{'n':>6} {'link type':>16} | {'BMR accept%':>11} {'B8 accept%':>11} {'BMR rate|acc':>13} {'B8 would-deliver rate':>22}")
for n in (500, 10_000):
    rad_q = math.sqrt(math.log(1 / DELTA) / (2 * 4 * n))
    for name, S, q in TYPES:
        sd = math.sqrt(max(0.0, 4 - S * S / 4) / n)
        S_hat = rng.normal(S, sd, N_LINKS)
        certs = [build_certificate("e", float(x), n, DELTA, ELL, PHI, 0, 0) for x in S_hat]
        bmr_ok = np.array([c is not None for c in certs])
        r_bmr = np.mean([c.r for c in certs if c is not None]) if bmr_ok.any() else 0.0
        q_hat = rng.binomial(4 * n, q, N_LINKS) / (4 * n)
        b8_ok = (q_hat + rad_q) < Q_TH
        # rate a QBER-monitor would deliver over an accepted link: it has no Bell value, so it
        # books the *honest* nominal rate of the link class (what it believes) -> risk quantity
        r_b8_believed = 0.0
        if b8_ok.any():
            c_h = build_certificate("e", 2.687, n, DELTA, ELL, PHI, 0, 0)
            r_b8_believed = c_h.r if c_h else 0.0
        row = dict(n=n, link=name, S=S, qber=q, bmr_accept_pct=100 * float(bmr_ok.mean()),
                   b8_accept_pct=100 * float(b8_ok.mean()), bmr_mean_rate=float(r_bmr),
                   b8_believed_rate=float(r_b8_believed), qber_radius=rad_q)
        res.append(row)
        print(f"{n:>6} {name:>16} | {row['bmr_accept_pct']:>10.2f}% {row['b8_accept_pct']:>10.2f}% "
              f"{r_bmr:>13.4f} {r_b8_believed * (row['b8_accept_pct'] / 100):>22.4f}")

# self-checks: by-construction expectations
g = lambda n, k: next(r for r in res if r["n"] == n and r["link"] == k)
assert g(10_000, "classical-mimic")["bmr_accept_pct"] == 0.0           # Bell audit never certifies S<=2
assert g(500, "classical-mimic")["bmr_accept_pct"] == 0.0
assert g(10_000, "classical-mimic")["b8_accept_pct"] > 99.0            # QBER sees an honest link
assert g(10_000, "honest")["bmr_accept_pct"] > 99.9 and g(10_000, "honest")["b8_accept_pct"] > 99.9
# depolarised @ n=500: mean S_lb = 1.948 < 2 but per-audit noise (sd~0.07) lets a minority pass; report, do not assert 0
assert g(500, "depolarised")["bmr_accept_pct"] < 50.0
assert g(10_000, "depolarised")["bmr_accept_pct"] > 99.0   # S=2.418 certifies at reduced rate (0.14) at n=1e4
print("[self-check] mimic: BMR rejects 100%, B8 accepts; honest: both accept; depolarised: BMR alarm-rate 77% @500, 0% @1e4 (reduced rate)")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for a, n in zip(ax, (500, 10_000)):
    rr = [r for r in res if r["n"] == n]
    x = np.arange(3)
    a.bar(x - .2, [r["bmr_accept_pct"] for r in rr], .38, label="BMR (Bell audit)", color="green", ec="k", lw=.5)
    a.bar(x + .2, [r["b8_accept_pct"] for r in rr], .38, label="B8 (QBER monitor)", color="purple", ec="k", lw=.5)
    a.set_xticks(x); a.set_xticklabels([r["link"] for r in rr]); a.set_title(f"n = {n}"); a.grid(alpha=.3, axis="y")
ax[0].set_ylabel("link accepted as healthy [%]"); ax[0].legend(fontsize=8)
plt.suptitle("Q5: evidence source (classical mimic has S=2, honest-looking QBER)")
plt.tight_layout()
fp = os.path.join(os.path.dirname(__file__), "..", "figures", "Q5_evidence_source.png")
plt.savefig(fp, dpi=150); plt.close()
json.dump(dict(experiment="Q5", params=dict(delta=DELTA, phi=PHI, ell=ELL, Q_th=Q_TH, N_links=N_LINKS), rows=res),
          open(os.path.join(os.path.dirname(__file__), "..", "data", "results", "q5_evidence_source.json"), "w"), indent=2)
print("saved", fp)
