"""
experiments/q6_dynamic_latency.py
=================================
Q3 follow-up: detection latency, false-alarm and false-negative rates for the Q6-dynamic scenario
(same topologies, n_audit, S1 grid, timeline, seeds-structure, Delta and reps as q6_dynamic.py).

Two levels, both reported per (topology, n_audit, S1):
 LINK level (the certificate / evidence rule that each router relies on)
   cert rule  (BMR, B3, B5, B0 share it): link is "alarmed" iff it has NO valid certificate (r_e = 0).
   B8 rule                              : link is "alarmed" iff r_q = 0 (direct QBER bound).
   false alarm  = honest link-epoch alarmed (all epochs before the link's own onset, plus never-compromised links)
   false negative = compromised link-epoch NOT alarmed (r_e>0 although true S = S1)
   latency      = epochs from onset (25 for set 1, 60 for set 2) to the first alarmed epoch (right-censored at
                  the window end: 35 / 20 epochs); 'detected' = fraction of compromised links ever alarmed.
 ROUTE level (per router series, scored on the live certificates, as in q6_dynamic)
   FN_route   = P(route delivered AND not compliant on the live data) over post-1 epochs (25-59)
   FA_route   = P(router refuses although a compliant route exists on the live data) over pre epochs (0-24)
   refuse     = P(no route delivered) over post-1 epochs
   recover    = epochs after 25 until the router is compliant in every remaining post-1 epoch (35 = never)
"""
import sys, os, json, time, platform
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
sys.path.insert(0, os.path.dirname(__file__))
from multiprocessing import Pool
import numpy as np
import networkx as nx
from bmr.sim import LinkSim
from q6_common import (make_topology, attach_certs, evaluate, feasible_oracle, M_FLOOR)
from q6_dynamic import (REPS, N_EP, E1, E2, F1, F2, DELTA_BY_N, S1_LIST, TOPOS, ALGOS, SERIES)

HERE = os.path.dirname(__file__)
NN = {"pre": (0, E1), "post1": (E1, E2), "post2": (E2, N_EP)}


def cell(spec):
    topo, n, na, S1, seed = spec
    Delta = DELTA_BY_N[na]
    rng = np.random.default_rng(seed)
    R = {lab: dict(fn=[], fa=[], ref=[], rec=[]) for lab, _, _ in SERIES}
    L = {k: dict(fa=[], fn=[], lat1=[], det1=[], lat2=[], det2=[]) for k in ("cert", "B8")}
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
        onset = {frozenset(e): E1 for e in set1}
        onset.update({frozenset(e): E2 for e in set2})
        sim = LinkSim(edges, na, delta=1e-6, ell=10_000, phi=0.99, seed=int(rng.integers(1e9)))
        G_live = G0.copy(); G_frozen = None
        alarm = {k: {frozenset(e): [] for e in edges} for k in L}
        rt = {lab: dict(found=[], comp=[], orc=[]) for lab, _, _ in SERIES}
        for ep in range(N_EP):
            if ep == E1:
                for e in set1: sim.set_true_S(e, S1)
            if ep == E2:
                for e in set2: sim.set_true_S(e, S1)
            attach_certs(G_live, sim, epoch=ep)
            if ep == 0:
                G_frozen = G_live.copy()
            for u, v in G_live.edges():
                a = G_live[u][v]; k = frozenset((u, v))
                alarm["cert"][k].append(a["r"] <= 0)
                alarm["B8"][k].append(a["r_q"] <= 0)
            orc = int(ep < E1 and feasible_oracle(G_live, s, d, Delta))
            for lab, an, mode in SERIES:
                Gv = G_live if mode == "live" else G_frozen
                if mode == "frozen":
                    for u, v in Gv.edges():
                        Gv[u][v]["T"] = ep
                path = ALGOS[an](Gv, s, d, Delta)
                f, c, *_ = evaluate(G_live, path, Delta)
                rt[lab]["found"].append(f); rt[lab]["comp"].append(c); rt[lab]["orc"].append(orc)
        for lab, _, _ in SERIES:
            f = np.array(rt[lab]["found"]); c = np.array(rt[lab]["comp"]); o = np.array(rt[lab]["orc"])
            a, b = NN["post1"]; p0, p1 = NN["pre"]
            R[lab]["fn"].append(float(((f == 1) & (c == 0))[a:b].mean()))
            R[lab]["ref"].append(float((f == 0)[a:b].mean()))
            den = o[p0:p1].sum()
            R[lab]["fa"].append(float(((f == 0) & (o == 1))[p0:p1].sum() / den) if den else float("nan"))
            ok = c[a:b]; rec = len(ok)
            for i in range(len(ok) - 1, -1, -1):
                if ok[i] == 1: rec = i
                else: break
            R[lab]["rec"].append(rec if ok[-1] == 1 else len(ok))
        for k in L:
            for e in edges:
                key = frozenset(e); al = np.array(alarm[k][key]); on = onset.get(key, N_EP + 1)
                hon = al[:min(on, N_EP)]
                L[k]["fa"].append(float(hon.mean()))
                if key in onset:
                    cmp_ = al[on:]
                    L[k]["fn"].append(float((~cmp_).mean()))
                    hit = np.flatnonzero(cmp_)
                    tag = "1" if on == E1 else "2"
                    L[k]["det" + tag].append(int(len(hit) > 0))
                    L[k]["lat" + tag].append(int(hit[0]) if len(hit) else len(cmp_))
    def m(x):
        x = [v for v in x if v == v]
        return float(np.mean(x)) if x else float("nan")
    out = dict(topo=topo, n_audit=na, S1=S1, Delta=Delta, reps=REPS,
               link={k: {kk: m(vv) for kk, vv in L[k].items()} for k in L},
               route={lab: {kk: m(vv) for kk, vv in R[lab].items()} for lab in R})
    return out


def main():
    t0 = time.time()
    rng = np.random.default_rng(77)
    specs = [(t, n, na, S1, int(rng.integers(1e9))) for (t, n) in TOPOS for na in (500, 10_000) for S1 in S1_LIST]
    rows = []
    with Pool(max(1, (os.cpu_count() or 4) - 1)) as p:
        for i, r in enumerate(p.imap(cell, specs, chunksize=1)):
            rows.append(r); print(f"  {i+1}/{len(specs)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(dict(experiment="Q6-dynamic-latency", header=dict(reps=REPS, epochs=N_EP, compromise_epochs=[E1, E2],
               delta_by_n=DELTA_BY_N, S1=S1_LIST, python=platform.python_version(), seconds=round(time.time()-t0, 1)),
               rows=rows), open(os.path.join(HERE, "..", "data", "results", "q6_dynamic_latency.json"), "w"), indent=1)
    print("saved q6_dynamic_latency.json", flush=True)

if __name__ == "__main__":
    main()
