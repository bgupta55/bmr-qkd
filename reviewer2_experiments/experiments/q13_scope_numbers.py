"""
Q13 / Prop. 2 support numbers (all analytic or Monte-Carlo, no tuning).
  (1) r_e (min-entropy certificate, lambda_EC = 0) vs the asymptotic collective-attack DI-QKD rate
      r_DI = 1 - h2((1+sqrt((S/2)^2-1))/2) - h2(Q)   [Pironio et al., New J. Phys. 11, 045021 (2009)],
      with the Werner relation Q = (1 - S/(2 sqrt 2))/2.  r_e is NOT a DI-QKD key rate: it omits EC leakage h2(Q).
  (2) Proposition 2: midpoint test theta=(S*+S')/2.  Hoeffding bound exp(-n (S*-S')^2/32) per tail vs the exact
      binomial-CHSH Monte-Carlo error rate (sd = sqrt((4-S^2/4)/n), the simulator's own audit model).
Output: data/results/q13_scope_numbers.json
"""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from bmr.bounds import h_pironio
RES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "results")
h2 = lambda p: 0.0 if p <= 0 or p >= 1 else -p*math.log2(p) - (1-p)*math.log2(1-p)

rows1 = []
for S in (2.687, 2.55, 2.418):
    Q = (1 - S/(2*math.sqrt(2)))/2
    hmin = h_pironio(S)
    s_ent = 1 - h2((1+math.sqrt((S/2)**2 - 1))/2)           # H(A|E) bound, collective attacks
    rows1.append(dict(S=S, Q_werner=Q, h_minent=hmin, H_AE_collective=s_ent, h2_Q=h2(Q),
                      r_DI_asymptotic=s_ent - h2(Q), h_minent_minus_h2Q=hmin - h2(Q)))

rng = np.random.default_rng(13)
rows2 = []
S0, S1 = 2.687, 2.418
theta = (S0+S1)/2
for n in (500, 2500, 10000):
    sd = lambda S: math.sqrt((4 - S*S/4)/n)
    fa = float((rng.normal(S0, sd(S0), 2_000_000) <= theta).mean())      # honest flagged
    miss = float((rng.normal(S1, sd(S1), 2_000_000) > theta).mean())     # compromised missed
    hoef = math.exp(-n*(S0-S1)**2/32)
    rows2.append(dict(n=n, theta=theta, false_alarm_MC=fa, miss_MC=miss, hoeffding_bound_per_tail=hoef))
    print(f"n={n:6d}  false-alarm {fa:.2e}  miss {miss:.2e}  Hoeffding bound {hoef:.2e}")
for r in rows1: print({k: round(v, 4) for k, v in r.items()})
json.dump(dict(experiment="Q13-scope", di_vs_cert=rows1, prop2=rows2), open(os.path.join(RES, "q13_scope_numbers.json"), "w"), indent=1)
