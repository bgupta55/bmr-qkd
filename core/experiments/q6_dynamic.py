"""
experiments/q6_dynamic.py
=========================
Q6 dynamic scenario (mid-session compromise) - pending_v7 Q6.

Timeline (80 epochs): a first set of links (20 % of all links) is compromised at epoch 25, a second set
(a further 10 %) at epoch 60; compromised links have true CHSH value S1 in {2.62, 2.55, 2.48, 2.418}.
Topologies: ring (n=30), grid 5x6 ("mesh"), ER (n=30).  n_audit in {500, 10^4}; Delta as in the Q6 sweep.

Information modes:
  live            the router sees the certificates re-published in the current epoch
  frozen          the router sees only the epoch-0 database (a static design with no monitoring / "as published")
Truth: every delivered route is scored on the CURRENT live certificates:
  secure_rate = min_e r_e(now) if (all links certified now) and sum w(now) <= Delta and hops >= 3, else 0.
Reported per window (pre 0-24 | post1 25-59 | post2 60-79): mean secure rate and compliance %, zeros included.
"BMR-frozen" is the static-BMR ablation (no re-audit, no expiry).
"""
import sys, os, json, time, platform
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
sys.path.insert(0, os.path.dirname(__file__))
from multiprocessing import Pool
import numpy as np
import networkx as nx

from bmr.sim import LinkSim
from q6_common import (make_topology, attach_certs, evaluate, r_spf, r_b3, r_b5, r_b8, r_bmr, M_FLOOR)

HERE = os.path.dirname(__file__)
REPS = int(os.environ.get("Q6D_REPS", 100))
N_EP, E1, E2 = 80, 25, 60
F1, F2 = 0.20, 0.10
DELTA_BY_N = {500: 0.5, 10_000: 1e-3}
S1_LIST = [2.62, 2.55, 2.48, 2.418]
TOPOS = [("ring", 30), ("grid", 30), ("ER", 30)]
WINS = {"pre": (0, 25), "post1": (25, 60), "post2": (60, 80)}
ALGOS = {"B0-SPF": r_spf, "B3-Zhang-style": r_b3, "B5-Rate-greedy": r_b5, "B8-QBER-monitor": r_b8, "BMR": r_bmr}
SERIES = [("BMR-live", "BMR", "live"), ("BMR-frozen(static ablation)", "BMR", "frozen"),
          ("B3-live", "B3-Zhang-style", "live"), ("B3-as-published", "B3-Zhang-style", "frozen"),
          ("B5-live", "B5-Rate-greedy", "live"), ("B5-as-published", "B5-Rate-greedy", "frozen"),
          ("B8-live", "B8-QBER-monitor", "live"), ("B8-as-published", "B8-QBER-monitor", "frozen"),
          ("B0-SPF", "B0-SPF", "frozen")]


def cell(spec):
    topo, n, na, S1, seed = spec
    Delta = DELTA_BY_N[na]
    rng = np.random.default_rng(seed)
    acc = {lab: np.zeros((REPS, N_EP)) for lab, _, _ in SERIES}
    cmp_ = {lab: np.zeros((REPS, N_EP)) for lab, _, _ in SERIES}
    for rep in range(REPS):
        while True:
            G0 = make_topology(topo, n, int(rng.integers(1e9)))
            if nx.is_connected(G0):
                break
        dist = dict(nx.all_pairs_shortest_path_length(G0))
        pairs = [(a, b) for a in G0 for b in G0 if a < b and dist[a][b] >= 4]
        s, d = pairs[int(rng.integers(len(pairs)))]
        edges = list(G0.edges())
        perm = rng.permutation(len(edges))
        k1, k2 = max(1, round(F1 * len(edges))), max(1, round(F2 * len(edges)))
        set1 = [edges[i] for i in perm[:k1]]
        set2 = [edges[i] for i in perm[k1:k1 + k2]]
        sim = LinkSim(edges, na, delta=1e-6, ell=10_000, phi=0.99, seed=int(rng.integers(1e9)))
        G_live = G0.copy()
        G_frozen = None
        for ep in range(N_EP):
            if ep == E1:
                for e in set1: sim.set_true_S(e, S1)
            if ep == E2:
                for e in set2: sim.set_true_S(e, S1)
            attach_certs(G_live, sim, epoch=ep)
            if ep == 0:
                G_frozen = G_live.copy()
            for lab, an, mode in SERIES:
                Gv = G_live if mode == "live" else G_frozen
                if mode == "frozen":
                    for u, v in Gv.edges():
                        Gv[u][v]["T"] = ep            # frozen DB never expires (as-published designs have no freshness rule)
                path = ALGOS[an](Gv, s, d, Delta)
                f, c, rate, h, w, cb, rb_ = evaluate(G_live, path, Delta)   # truth = live certificates
                acc[lab][rep, ep] = rate
                cmp_[lab][rep, ep] = c
    out = []
    for lab, _, _ in SERIES:
        row = dict(topo=topo, n_audit=na, S1=S1, series=lab, Delta=Delta, reps=REPS)
        for wn, (a, b) in WINS.items():
            row["rate_" + wn] = float(acc[lab][:, a:b].mean())
            row["compl_" + wn] = float(100 * cmp_[lab][:, a:b].mean())
        row["per_epoch_rate"] = acc[lab].mean(axis=0).round(5).tolist()
        out.append(row)
    return out


def main():
    t0 = time.time()
    rng = np.random.default_rng(77)
    specs = [(t, n, na, S1, int(rng.integers(1e9))) for (t, n) in TOPOS for na in (500, 10_000) for S1 in S1_LIST]
    rows = []
    with Pool(max(1, (os.cpu_count() or 4) - 1)) as p:
        for i, r in enumerate(p.imap(cell, specs, chunksize=1)):
            rows += r
            print(f"  {i+1}/{len(specs)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(dict(experiment="Q6-dynamic",
                   header=dict(reps=REPS, epochs=N_EP, compromise_epochs=[E1, E2], fractions=[F1, F2],
                               delta_by_n=DELTA_BY_N, S1=S1_LIST, hop_floor=M_FLOOR,
                               python=platform.python_version(), seconds=round(time.time() - t0, 1)),
                   rows=rows),
              open(os.path.join(HERE, "..", "data", "results", "q6_dynamic.json"), "w"))
    print("saved q6_dynamic.json", flush=True)


if __name__ == "__main__":
    main()
