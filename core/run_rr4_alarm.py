"""
run_rr4_alarm.py
================
RR4: Dynamic study with BMR+alarm column.
Adds the alarm filter to BMR routing in the dynamic scenario.
Alarm: a link fails the alarm if S_hat < S0_link - z_alpha * sigma_W(n)
       where z_alpha = 4.7534 (alpha=1e-6) and S0_link = 2.687 (honest baseline).

Uses the R2 package (bmr/sim.py, q6_common.py).
Run from: REVISION_OSN_12OCT2026/
  python3 run_rr4_alarm.py
"""
import sys, os, json, time, math, itertools
import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "r1_r2_150ct", "rr_results")
os.makedirs(OUTDIR, exist_ok=True)

R2_PKG = os.path.join(HERE, "Reviewer 2 output mrng")
sys.path.insert(0, R2_PKG)
sys.path.insert(0, os.path.join(R2_PKG, "experiments"))

from bmr.sim import LinkSim, S_HONEST
from bmr.route import bmr_route
from q6_common import make_topology, attach_certs, evaluate, r_bmr, M_FLOOR, R_MIN

Z_ALPHA = 4.7534
S0_DEFAULT = S_HONEST  # 2.687

def sigma_w(S0, n):
    return math.sqrt(max(0.0, 4.0 - S0**2 / 4.0) / n)

def r_bmr_alarm(G, s, d, Delta, S0=S0_DEFAULT):
    """BMR + alarm filter: exclude links where S_hat < theta_alpha."""
    import heapq
    # Build alarm-filtered graph
    H = nx.Graph()
    n_val = None
    for u, v, a in G.edges(data=True):
        cert = a.get("cert")
        if cert is None:
            continue
        n_val = cert.n
        theta = S0 - Z_ALPHA * sigma_w(S0, n_val)
        s_hat = a.get("S_hat", float("inf"))
        alarm = s_hat < theta
        if not alarm and a.get("r", 0) >= R_MIN:
            H.add_edge(u, v, **a)
    try:
        res = bmr_route(H, s, d, R_MIN, Delta, 0, 10, m_min=M_FLOOR)
        return res[1] if res else None
    except Exception:
        return None


REPS = 50  # reduced for speed (100 in the original)
N_EP, E1, E2 = 80, 25, 60
F1, F2 = 0.20, 0.10
DELTA_BY_N = {500: 0.5, 10_000: 1e-3}
S1_LIST = [2.62, 2.55, 2.48, 2.418]
TOPOS = [("ring", 30), ("grid", 30), ("ER", 30)]
WINS = {"pre": (0, 25), "post1": (25, 60), "post2": (60, 80)}

print("=== RR4: Dynamic + BMR+alarm ===")
t0 = time.time()
rng = np.random.default_rng(77)
rr4_rows = []

for (topo, n), na, S1 in itertools.product(TOPOS, [500, 10_000], S1_LIST):
    Delta = DELTA_BY_N[na]
    bmr_rates = np.zeros((REPS, N_EP))
    bma_rates = np.zeros((REPS, N_EP))
    rng_cell = np.random.default_rng(int(rng.integers(1e9)))

    for rep in range(REPS):
        for attempt in range(20):
            G0 = make_topology(topo, n, int(rng_cell.integers(1e9)))
            if nx.is_connected(G0):
                break
        dist = dict(nx.all_pairs_shortest_path_length(G0))
        pairs = [(a, b) for a in G0 for b in G0 if a < b and dist[a].get(b, 0) >= 4]
        if not pairs:
            pairs = [(a, b) for a in G0 for b in G0 if a < b and a != b]
        sv, dv = pairs[int(rng_cell.integers(len(pairs)))]
        edges = list(G0.edges())
        perm = rng_cell.permutation(len(edges))
        k1, k2 = max(1, round(F1 * len(edges))), max(1, round(F2 * len(edges)))
        set1 = [edges[i] for i in perm[:k1]]
        set2 = [edges[i] for i in perm[k1:k1+k2]]
        sim = LinkSim(edges, na, delta=1e-6, ell=10_000, phi=0.99, seed=int(rng_cell.integers(1e9)))
        G_live = G0.copy()

        for ep in range(N_EP):
            if ep == E1:
                for e in set1: sim.set_true_S(e, S1)
            if ep == E2:
                for e in set2: sim.set_true_S(e, S1)
            attach_certs(G_live, sim, epoch=ep)
            # BMR (no alarm)
            path_bmr = r_bmr(G_live, sv, dv, Delta)
            _, _, rate_bmr, _, _, _, _ = evaluate(G_live, path_bmr, Delta)
            bmr_rates[rep, ep] = rate_bmr
            # BMR + alarm
            path_bma = r_bmr_alarm(G_live, sv, dv, Delta)
            _, _, rate_bma, _, _, _, _ = evaluate(G_live, path_bma, Delta)
            bma_rates[rep, ep] = rate_bma

    row = {"topo": topo, "n_audit": na, "S1": S1, "Delta": Delta, "reps": REPS}
    for wn, (a_ep, b_ep) in WINS.items():
        row[f"BMR_live_rate_{wn}"]  = float(bmr_rates[:, a_ep:b_ep].mean())
        row[f"BMR_alarm_rate_{wn}"] = float(bma_rates[:, a_ep:b_ep].mean())
    rr4_rows.append(row)
    print(f"  {topo:5s} n={na:5d} S1={S1}: BMR_pre={row['BMR_live_rate_pre']:.4f} "
          f"BMR+alarm_pre={row['BMR_alarm_rate_pre']:.4f} "
          f"post1 BMR={row['BMR_live_rate_post1']:.4f} BMR+alarm={row['BMR_alarm_rate_post1']:.4f}")

json.dump({"experiment": "RR4-alarm", "reps": REPS, "rows": rr4_rows},
          open(os.path.join(OUTDIR, "RR4_dynamic_alarm.json"), "w"), indent=2)
print(f"\nRR4 done ({time.time()-t0:.1f}s) -> RR4_dynamic_alarm.json")
