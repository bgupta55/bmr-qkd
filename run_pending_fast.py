"""
run_pending_fast.py
===================
Faster version: 200 trials per cell (sufficient for stable estimates at 4 decimal places).
Focuses on key cells only: ER N=30, compromised (f=0.20), S† in {2.418, 2.55, 2.62}, Delta in {1e-3, 3e-5}.

RR1: decomposition
RR2: ring flagship S†=2.55 and S†=2.62 (with alarm tracking)
RR5: 5000-instance cross-validation

Run from: REVISION_OSN_12OCT2026/
  python3 run_pending_fast.py
"""
import sys, os, json, math, time
import numpy as np
import itertools

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "r1_r2_150ct", "rr_results")
os.makedirs(OUTDIR, exist_ok=True)

R2_PKG = os.path.join(HERE, "Reviewer 2 output mrng")
sys.path.insert(0, R2_PKG)
from bmr.bounds import eta_tight, h_pironio

# ─── constants ───────────────────────────────────────────────────────────────
Z_ALPHA = 4.7534
S0_REF  = 2.687
ELL     = 10_000
PHI     = 0.99
DELTA   = 1e-6
N_AUDIT = 10_000
R_MIN   = 0.01
CLASS_S = np.array([2.720, 2.687, 2.640])
CLASS_P = np.array([0.3, 0.5, 0.2])

def sigma_w(S0, n):
    return math.sqrt(max(0.0, 4.0 - S0**2 / 4.0) / n)

def certified_rate(S_hat, n=N_AUDIT, delta=DELTA, ell=ELL, phi=PHI):
    eta = eta_tight(n, delta)
    S_lb = S_hat - eta
    if S_lb <= 2.0:
        return 0.0, 1.0, 0.0, False
    h = h_pironio(S_lb)
    H = ell * h
    k = math.floor(phi * H)
    if k <= 0 or H <= k:
        return 0.0, 1.0, h, True
    log2_eps = -1.0 - (H - k) / 2.0
    eps = 0.0 if log2_eps < -1000 else 2.0**log2_eps
    r = k / ell
    return r, eps, h, True

def alarm_fires(S_hat, S0=S0_REF, n=N_AUDIT):
    theta = S0 - Z_ALPHA * sigma_w(S0, n)
    return S_hat < theta

# ─── pre-compute derived ε for the (R) parameter table ────────────────────
print("=== Derived ε bridge (n=10^4, phi=0.99) ===")
for S in [2.720, 2.687, 2.640, 2.62, 2.55, 2.48, 2.418]:
    r, eps, h, gate = certified_rate(S, n=N_AUDIT, delta=DELTA, ell=ELL, phi=PHI)
    w = eps + DELTA
    print(f"  S={S:.3f}: r={r:.4f}, eps={eps:.3e}, w={w:.3e}, gate={gate}")


# ════════════════════════════════════════════════════════════════════════════
#  RR1 — DECOMPOSITION UNDER (R) — key cells only
# ════════════════════════════════════════════════════════════════════════════
print("\n=== RR1: Decomposition under (R) ===")

N_TRIALS = 200  # per cell (fast)
N_NODES  = 30
ER_P     = 4.0 / (N_NODES - 1)
F_COMP   = 0.20
S_ATTACK_LIST = [2.418, 2.55, 2.62]
DELTA_LIST    = [1e-3, 3e-5]


def make_er_graph(n, p, rng):
    iu, ju = np.triu_indices(n, 1)
    keep = rng.random(len(iu)) < p
    return list(zip(iu[keep].tolist(), ju[keep].tolist()))

def make_ba_graph(n, m_ba, rng):
    import networkx as nx
    G = nx.barabasi_albert_graph(n, m_ba, seed=int(rng.integers(2**31)))
    return sorted((min(a,b), max(a,b)) for a,b in G.edges())

def sample_links_R(edges, n_audit, f, rng, S_att):
    E = len(edges)
    cls = rng.choice(3, size=E, p=CLASS_P)
    S_true = CLASS_S[cls].copy()
    k_comp = max(1, int(round(f * E)))
    comp_idx = rng.choice(E, size=k_comp, replace=False)
    comp = np.zeros(E, dtype=bool)
    comp[comp_idx] = True
    S_true[comp] = S_att
    p_win = (1.0 + S_true / 4.0) / 2.0
    K = rng.binomial(4 * n_audit, p_win)
    S_hat = 2.0 * K / n_audit - 4.0
    r_arr = np.zeros(E); eps_arr = np.ones(E); gate_arr = np.zeros(E, dtype=bool); alarm_arr = np.zeros(E, dtype=bool)
    for i in range(E):
        r_arr[i], eps_arr[i], _, gate_arr[i] = certified_rate(S_hat[i])
        alarm_arr[i] = alarm_fires(S_hat[i])
    w_arr = eps_arr + DELTA
    valid_arr = gate_arr & ~alarm_arr
    return dict(r=r_arr, eps=eps_arr, w=w_arr, gate=gate_arr, alarm=alarm_arr, valid=valid_arr)

def best_path_variant(edges, n_nodes, s, d, lk, Delta, mode):
    """DFS over all simple s-d paths (Lmin=3, Hmax=8), return bottleneck rate."""
    from collections import defaultdict
    adj_list = defaultdict(list)
    for (u, v) in edges:
        adj_list[u].append(v); adj_list[v].append(u)
    edge_idx = {(u,v): i for i,(u,v) in enumerate(edges)}
    for i,(u,v) in enumerate(edges):
        edge_idx[(v,u)] = i
    r, w, valid = lk["r"], lk["w"], lk["valid"]
    best = [-1.0, math.inf, math.inf]  # [rate, hops, sw]

    def dfs(node, visited, links_on_path, hops):
        if hops > 8: return
        if node == d:
            if hops < 3: return
            sw = sum(w[li] for li in links_on_path)
            rv = min(r[li] for li in links_on_path)
            av = all(valid[li] for li in links_on_path)
            bud = sw <= Delta * (1 + 1e-9)
            rte = rv >= R_MIN
            if mode == "spf":
                if hops < best[1] or (hops == best[1]):
                    if hops < best[1]: best[0] = rv; best[1] = hops; best[2] = sw
            elif mode == "spf_valid":
                if av and (hops < best[1] or (hops == best[1] and rv > best[0])):
                    best[0] = rv; best[1] = hops; best[2] = sw
            elif mode == "spf_budget":
                if bud and rte and (hops < best[1] or (hops == best[1] and rv > best[0])):
                    best[0] = rv; best[1] = hops; best[2] = sw
            elif mode == "spf_vb":
                if av and bud and rte and (hops < best[1] or (hops == best[1] and rv > best[0])):
                    best[0] = rv; best[1] = hops; best[2] = sw
            elif mode == "bmr":
                if av and bud and rte and (rv > best[0] or (rv == best[0] and hops < best[1])):
                    best[0] = rv; best[1] = hops; best[2] = sw
            return
        for nb in adj_list[node]:
            if nb in visited: continue
            ei = edge_idx.get((node, nb), -1)
            if ei < 0: continue
            visited.add(nb)
            links_on_path.append(ei)
            dfs(nb, visited, links_on_path, hops+1)
            links_on_path.pop()
            visited.discard(nb)

    dfs(s, {s}, [], 0)
    return max(0.0, best[0])

def decomp_trial(edges, n_nodes, s, d, lk, Delta):
    r_bmr = best_path_variant(edges, n_nodes, s, d, lk, Delta, "bmr")
    r_spf = best_path_variant(edges, n_nodes, s, d, lk, Delta, "spf")
    r_d2  = best_path_variant(edges, n_nodes, s, d, lk, Delta, "spf_budget")
    r_d3  = best_path_variant(edges, n_nodes, s, d, lk, Delta, "spf_vb")
    return r_bmr, r_spf, r_d2, r_d3

rr1_results = {}
t0_rr1 = time.time()
for topo in ["er", "ba"]:
    for S_att in S_ATTACK_LIST:
        for Delta in DELTA_LIST:
            key = f"{topo}_N30_f020_Satt{S_att}_D{Delta:.0e}"
            rates = {"bmr": [], "spf": [], "d2": [], "d3": []}
            for trial in range(N_TRIALS):
                trial_rng = np.random.default_rng(np.random.SeedSequence([20261012, trial, int(S_att*1000), int(Delta*1e6)]))
                if topo == "er":
                    edges = make_er_graph(N_NODES, ER_P, trial_rng)
                else:
                    edges = make_ba_graph(N_NODES, 2, trial_rng)
                if not edges: continue
                lk = sample_links_R(edges, N_AUDIT, F_COMP, trial_rng, S_att)
                s, d_node = 0, N_NODES - 1
                r_bmr, r_spf, r_d2, r_d3 = decomp_trial(edges, N_NODES, s, d_node, lk, Delta)
                rates["bmr"].append(r_bmr); rates["spf"].append(r_spf)
                rates["d2"].append(r_d2); rates["d3"].append(r_d3)

            bm = np.array(rates["bmr"]); sp = np.array(rates["spf"])
            d2 = np.array(rates["d2"]); d3 = np.array(rates["d3"])
            n = len(bm)
            ord1_budget   = float((d2 - sp).mean())
            ord1_validity = float((d3 - d2).mean())
            objective     = float((bm - d3).mean())
            total         = float((bm - sp).mean())
            feasible      = int((bm > 0).sum())

            def boot_ci(a, b, B=1000):
                diffs = a - b
                if len(diffs) == 0: return [0.0, 0.0]
                ms = [diffs[np.random.randint(0, len(diffs), len(diffs))].mean() for _ in range(B)]
                return np.percentile(ms, [2.5, 97.5]).tolist()

            ci = boot_ci(bm, sp)
            rr1_results[key] = {
                "topo": topo, "N": N_NODES, "f": F_COMP, "S_att": S_att, "Delta": Delta,
                "n_trials": n, "n_feasible": feasible,
                "BMR_mean": float(bm.mean()), "SPF_mean": float(sp.mean()),
                "D2_mean": float(d2.mean()), "D3_mean": float(d3.mean()),
                "total_BMR_SPF": total, "BMR_SPF_ci95": ci,
                "ord1_budget": ord1_budget, "ord1_validity": ord1_validity, "objective": objective,
            }
            print(f"  {key:45s}: total={total:+.4f} [budget={ord1_budget:+.4f}, valid={ord1_validity:+.4f}, obj={objective:+.4f}] feasible={feasible}/{n}")

json.dump(rr1_results, open(os.path.join(OUTDIR, "RR1_decomposition_R_family.json"), "w"), indent=2)
print(f"RR1 done ({time.time()-t0_rr1:.1f}s) -> RR1_decomposition_R_family.json")


# ════════════════════════════════════════════════════════════════════════════
#  RR2 — RING FLAGSHIP S†=2.55 and S†=2.62 with alarm tracking
# ════════════════════════════════════════════════════════════════════════════
print("\n=== RR2: Ring flagship S†=2.55 and S†=2.62 ===")

RING_N = 8; RING_EPOCHS = 60; ATTACK_EPOCH = 25; ATT_LINK = 1
RING_TRIALS = 500
S_B = S0_REF

def ring_epoch_links(rng, attacked, S_att, n=N_AUDIT, delta=DELTA):
    E = 8
    S_true = np.full(E, S_B)
    if attacked: S_true[ATT_LINK] = S_att
    p = (1.0 + S_true / 4.0) / 2.0
    K = rng.binomial(4*n, p)
    S_hat = 2.0*K/n - 4.0
    r_arr = np.zeros(E); eps_arr = np.ones(E); gate_arr = np.zeros(E, dtype=bool); alarm_arr = np.zeros(E, dtype=bool)
    for i in range(E):
        r_arr[i], eps_arr[i], _, gate_arr[i] = certified_rate(S_hat[i])
        alarm_arr[i] = alarm_fires(S_hat[i], S0=S_B, n=n)
    w_arr = eps_arr + delta
    valid_arr = gate_arr & ~alarm_arr
    return dict(r=r_arr, eps=eps_arr, w=w_arr, gate=gate_arr, alarm=alarm_arr, valid=valid_arr, S_hat=S_hat)

def ring_choose(lk, Delta):
    """Returns (on_A_bmr, r_sec_bmr, fa_count_clean) for the 8-node ring.
    Arc A = links {0,1,2,3}, Arc B = links {7,6,5,4}."""
    arc_A = [0,1,2,3]; arc_B = [7,6,5,4]
    r, w, valid = lk["r"], lk["w"], lk["valid"]
    def arc_stats(links):
        sw = sum(w[i] for i in links)
        rv = min(r[i] for i in links)
        av = all(valid[i] for i in links)
        return av and sw <= Delta*(1+1e-9) and rv >= R_MIN, rv
    ok_A, rv_A = arc_stats(arc_A); ok_B, rv_B = arc_stats(arc_B)
    if ok_A and ok_B:
        if rv_A >= rv_B: return 1, rv_A
        return 0, rv_B
    elif ok_A: return 1, rv_A
    elif ok_B: return 0, rv_B
    return -1, 0.0

rr2_results = {}
t0_rr2 = time.time()

for S_att in [2.418, 2.55, 2.62]:
    for Delta in [1e-3, 3e-5]:
        key = f"ring8_Satt{S_att}_D{Delta:.0e}"
        rng = np.random.default_rng(np.random.SeedSequence([20261015, int(S_att*1000), int(Delta*1e9)]))
        n_clean_dec = 0; fa_pre = 0; fa_post_other = 0
        detections_alarm = []  # epoch first alarm on attacked link
        s_att_hat = []; slb_att = []; r_att_list = []; gate_att = []
        pre_r_bmr = []; post_r_bmr = []; ep25_r_bmr = []

        for trial in range(RING_TRIALS):
            detect_alarm = None
            for ep in range(1, RING_EPOCHS+1):
                attacked = ep >= ATTACK_EPOCH
                lk = ring_epoch_links(rng, attacked, S_att)
                clean = np.ones(8, dtype=bool)
                if attacked: clean[ATT_LINK] = False
                n_clean_dec += int(clean.sum())
                fa = lk["alarm"] & clean
                if not attacked:
                    fa_pre += int(fa.sum())
                else:
                    fa_post_other += int(fa.sum())
                    s_att_hat.append(lk["S_hat"][ATT_LINK])
                    eta = eta_tight(N_AUDIT, DELTA)
                    slb_att.append(max(0.0, lk["S_hat"][ATT_LINK] - eta))
                    r_att_list.append(lk["r"][ATT_LINK])
                    gate_att.append(bool(lk["gate"][ATT_LINK]))
                    if detect_alarm is None and lk["alarm"][ATT_LINK]:
                        detect_alarm = ep
                on_A, r_sec = ring_choose(lk, Delta)
                if ep < ATTACK_EPOCH: pre_r_bmr.append(r_sec)
                elif ep == ATTACK_EPOCH: ep25_r_bmr.append(r_sec)
                else: post_r_bmr.append(r_sec)
            detections_alarm.append(detect_alarm)

        rr2_results[key] = {
            "S_att": S_att, "Delta": Delta, "trials": RING_TRIALS,
            "clean_decisions": n_clean_dec,
            "fa_pre_attack": fa_pre, "fa_post_attack_other_links": fa_post_other,
            "expected_fa_at_alpha_1e-6": n_clean_dec * 1e-6,
            "S_att_hat_mean": float(np.mean(s_att_hat)) if s_att_hat else None,
            "S_lb_att_mean": float(np.mean(slb_att)) if slb_att else None,
            "S_lb_att_min": float(np.min(slb_att)) if slb_att else None,
            "S_lb_att_max": float(np.max(slb_att)) if slb_att else None,
            "r_att_mean": float(np.mean(r_att_list)) if r_att_list else None,
            "r_att_max": float(np.max(r_att_list)) if r_att_list else None,
            "gate_pass_rate": float(np.mean(gate_att)) if gate_att else None,
            "alarm_detected": int(sum(x is not None for x in detections_alarm)),
            "alarm_detected_at_ep25": int(sum(x == 25 for x in detections_alarm)),
            "alarm_mean_detect_epoch": float(np.mean([x for x in detections_alarm if x is not None])) if any(x is not None for x in detections_alarm) else None,
            "BMR_pre_r_sec": float(np.mean(pre_r_bmr)),
            "BMR_ep25_r_sec": float(np.mean(ep25_r_bmr)),
            "BMR_post_r_sec": float(np.mean(post_r_bmr)),
        }
        r = rr2_results[key]
        print(f"  {key}: FA={fa_pre}/{fa_post_other} alarm_det={r['alarm_detected']}/500 "
              f"ep25={r['alarm_detected_at_ep25']}/500 mean_ep={r['alarm_mean_detect_epoch']} "
              f"gate_pass={r['gate_pass_rate']:.3f}")

json.dump(rr2_results, open(os.path.join(OUTDIR, "RR2_ring_flagship_R_family.json"), "w"), indent=2)
print(f"RR2 done ({time.time()-t0_rr2:.1f}s) -> RR2_ring_flagship_R_family.json")


# ════════════════════════════════════════════════════════════════════════════
#  RR5 — YEN VS BRUTE FORCE CROSS-VALIDATION (5000 instances)
# ════════════════════════════════════════════════════════════════════════════
print("\n=== RR5: Cross-validation Yen vs brute-force (5000 instances) ===")
t0_rr5 = time.time()

try:
    from bmr.route import bmr_route, brute_force_max_rate
    import networkx as nx

    N_INST = 5000
    N_SIZES = list(range(4, 13))
    T_VALS  = [0, 1, 2, 3]
    ER_P_CV = 0.5

    rng_cv = np.random.default_rng(20261015)
    n_match = 0; n_mismatch = 0; n_truncated = 0; n_infeasible_both = 0
    mismatch_details = []
    timing_yen = []; timing_bf = []
    inst_per_cell = max(1, N_INST // (len(N_SIZES) * len(T_VALS)))
    total_inst = 0

    for N_cv, t_cv in itertools.product(N_SIZES, T_VALS):
        m_min = 2*t_cv + 1
        for _ in range(inst_per_cell):
            total_inst += 1
            seed_inst = int(rng_cv.integers(1e9))
            G = nx.erdos_renyi_graph(N_cv, ER_P_CV, seed=seed_inst)
            if not nx.is_connected(G) or G.number_of_edges() == 0: continue
            for u, v in G.edges():
                G[u][v]["r"] = float(rng_cv.uniform(0.01, 0.35))
                G[u][v]["w"] = float(rng_cv.uniform(1e-7, 1e-5))
                G[u][v]["T"] = 0
            nodes = list(G.nodes())
            if len(nodes) < 2: continue
            idx_pair = rng_cv.choice(len(nodes), size=2, replace=False)
            s_cv, d_cv = int(nodes[idx_pair[0]]), int(nodes[idx_pair[1]])
            Delta_cv = float(rng_cv.uniform(2e-6, 3e-5))

            tY0 = time.perf_counter()
            res_yen = bmr_route(G, s_cv, d_cv, 0.01, Delta_cv, 0, 100, m_min=m_min)
            timing_yen.append((time.perf_counter()-tY0)*1000)

            tB0 = time.perf_counter()
            r_bf, path_bf, w_bf, trunc = brute_force_max_rate(G, s_cv, d_cv, Delta_cv, m_min=m_min)
            timing_bf.append((time.perf_counter()-tB0)*1000)

            if trunc: n_truncated += 1; continue

            r_yen = res_yen[0] if res_yen is not None else None
            if r_yen is None and r_bf is None:
                n_infeasible_both += 1; n_match += 1
            elif r_yen is None or r_bf is None:
                n_mismatch += 1
                if len(mismatch_details) < 5:
                    mismatch_details.append({"N": N_cv, "t": t_cv, "r_yen": r_yen, "r_bf": r_bf})
            elif abs(r_yen - r_bf) < 1e-9:
                n_match += 1
            else:
                n_mismatch += 1
                if len(mismatch_details) < 5:
                    mismatch_details.append({"N": N_cv, "t": t_cv, "r_yen": r_yen, "r_bf": r_bf})

    rr5_result = {
        "total_instances": total_inst,
        "truncated": n_truncated, "infeasible_both": n_infeasible_both,
        "matches": n_match, "mismatches": n_mismatch,
        "mismatch_details": mismatch_details,
        "yen_median_ms": float(np.median(timing_yen)) if timing_yen else None,
        "yen_p95_ms":    float(np.percentile(timing_yen, 95)) if timing_yen else None,
        "bf_median_ms":  float(np.median(timing_bf)) if timing_bf else None,
        "bf_p95_ms":     float(np.percentile(timing_bf, 95)) if timing_bf else None,
    }
    json.dump(rr5_result, open(os.path.join(OUTDIR, "RR5_cross_validation.json"), "w"), indent=2)
    print(f"  matches={n_match}, mismatches={n_mismatch}, truncated={n_truncated}, infeasible_both={n_infeasible_both}")
    print(f"  Yen median {rr5_result['yen_median_ms']:.3f} ms, p95 {rr5_result['yen_p95_ms']:.3f} ms")
    print(f"  BF  median {rr5_result['bf_median_ms']:.3f} ms, p95 {rr5_result['bf_p95_ms']:.3f} ms")
    print(f"  -> Saved RR5_cross_validation.json")

except ImportError as e:
    print(f"  WARNING: RR5 skipped ({e})")

print(f"\nAll done. RR1={time.time()-t0_rr1:.1f}s total. Results in {OUTDIR}")
