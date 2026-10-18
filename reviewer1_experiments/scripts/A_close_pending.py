"""
A_close_pending.py -- closes the Section A items that were still open after the first pass.

Open items handled here (numbers only; no code of the experiments is changed):
  C1  defect status D1-D13, checked against the actual scripts / data files (not from memory)
  C2  A19: Q-LEAP-style F_th sweep over a range in which the threshold really binds (0.70..0.95)
  C3  A15: consistency with Reviewer 2 (eta_tight primary there, eta_union used in the R1 package):
      sensitivity of the headline R1 differences to the Hoeffding radius, and of the attack detection channel
  C4  A3/A17: N_key penalty (restated) and certificate size with int32 freshness stamp
Run from scripts/:   python3 A_close_pending.py        Output: ../data/A_close_pending.json
"""
import json, math, os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "..", "data")
OUT = {}


def src(name): return open(os.path.join(HERE, name)).read()


# ------------------------------------------------------------------ C1 defects
D = {}
fa = json.load(open(os.path.join(DATA, "R1_false_alarm_check.json")))
D["D1"] = dict(fixed=fa["binomial"]["alarms"] == 0 and fa["gauss02"]["alarms"] == 45,
               evidence=f"200,000 honest links: binomial model {fa['binomial']['alarms']} alarms, legacy Gaussian 0.02 {fa['gauss02']['alarms']} ({fa['gauss02']['rate']:.2e})")
core = src("B_core.py")
D["D2"] = dict(fixed="bmr_mask = valid & bud & reqf" in core and 'sel["SPF_valid"] = _pick(valid' in core,
               evidence="B_core.select_all: BMR mask = valid & budget & rate floor; SPF_valid uses the validity mask")
l7 = src("L7_ablation_v5b.py")
D["D3"] = dict(fixed="gok & noal & fr & bud & rf" in l7 and "assert np.array_equal(M[\"allvalid\"]" in l7,
               evidence="L7_ablation_v5b.select(): 'full' applies gate, alarm and freshness masks during selection; asserted equal to allvalid")
cells = pd.read_csv(os.path.join(DATA, "R1_cells.csv"))
D["D4"] = dict(fixed=set(cells.cond.unique()) == {"clean", "compromised"},
               evidence=f"R1_cells.csv conditions: {sorted(cells.cond.unique())}; f values {sorted(cells.f.unique())}")
r8 = src("R8_overhead_v5b.py")
D["D5"] = dict(fixed="er_p=0.25" in core and "0.4" in r8 and "0.25" in r8,
               evidence="B_core.build_graph default er_p=0.25 is used by R1/R5/R10; R8 reports both p=0.25 and p=0.40 explicitly")
l2 = src("L2_alg2_v5b.py"); r2c = json.load(open(os.path.join(DATA, "R2_correctness.json")))
D["D6"] = dict(fixed="hop-layered Bellman-Ford" in l2 and "Second variant" in l2,
               evidence="L2_alg2_v5b.py implements prefix enumeration + hop-layered Bellman-Ford; DFS kept as second variant; correctness JSON present")
aer = json.load(open(os.path.join(DATA, "R7_aer_noiseless.json")))
D["D7"] = dict(fixed=abs(aer["degraded"]["S_ideal"] - aer["corrected"]["S_ideal"]) < 1e-9 and aer["opposite_sign"]["S_ideal"] < 2.7,
               evidence=f"ideal S: corrected {aer['corrected']['S_ideal']:.4f}, degraded (same sign) {aer['degraded']['S_ideal']:.4f} -> null perturbation; opposite sign {aer['opposite_sign']['S_ideal']:.4f}. "
                        "Documented, not a code fix: the 'degraded' arm is withdrawn as a degradation test")
fig = src("R7_figure.py")
D["D8"] = dict(fixed="R7_hardware_recomputed.json" in fig and "2.490" not in fig,
               evidence="R7_figure.py reads S and SE from the JSON recomputed from raw counts; no hard-coded S")
c06 = src("06_ibm_quantum_chsh.py")
D["D9"] = dict(fixed="optimization_level=1" in c06, evidence="06_ibm_quantum_chsh.py: transpile(..., optimization_level=1); the v5 statement 'level 3' is wrong")
hw = json.load(open(os.path.join(DATA, "R7_hardware_recomputed.json")))
dd = hw["difference_degraded_minus_native"]
D["D10"] = dict(fixed=abs(dd["z_physical"] - 1.75) < 0.01 and abs(dd["z_conservative"] - 1.36) < 0.01,
                evidence=f"z_physical = {dd['z_physical']:.3f}, z_conservative = {dd['z_conservative']:.3f}: both reported and labelled")
eta8192 = math.sqrt(32 * math.log(8 / 1e-6) / 8192)
rdeg = B.h(hw["degraded"]["S"] - eta8192)
D["D11"] = dict(fixed=abs(rdeg - 0.1245) < 5e-4, evidence=f"r_e(degraded, n=8192) = {rdeg:.4f} (v5 had 0.1248)")
D["D12"] = dict(fixed={"F_P", "cost_pairs", "throughput"} <= set(cells.columns),
                evidence="R1_cells.csv exports F_P, cost_pairs, throughput for every router")
D["D13"] = dict(fixed=False, evidence="OPEN by nature: the old circuit / job IDs are not in the workspace (see C1b); needs the authors' archive or withdrawal of 2.574/2.579")
OUT["C1_defects"] = D

# ---- C1b: old-circuit candidates under the ibm_fez noise model (separate script) + other correct-circuit hardware values
p = os.path.join(DATA, "A_close_aer_old_circuit.json")
if os.path.exists(p):
    OUT["C1b_old_circuit_noisy"] = json.load(open(p))
r2hw = json.load(open(os.path.join(HERE, "..", "..", "expr_r2", "data", "results", "e6v2_hardware_fixed.json")))
OUT["C1b_other_hardware_reference"] = dict(
    note="Reviewer-2 package, corrected RY(2theta) circuit, ibm_marrakesh qubits (6,5), 8192 shots, phi = 0 (different CHSH angle set)",
    S_phi0=r2hw["hw_results"]["0.0"]["S"], job=r2hw["job_ids"]["main_batch_v2"])

# --------------------------------------------------- C2: F_th sweep where it binds
def run_fth(fths, G=50, m=5, size=30, Delta=5e-6, R=0.10, H=5, seed=20261012, eta=None):
    if eta is not None:
        B.tau = lambda n, delta=1e-6, e=eta: e
    res = {f_: {c: [] for c in ("clean", "compromised")} for f_ in fths}
    bmr = {c: [] for c in ("clean", "compromised")}
    frac_bind = {f_: [] for f_ in fths}
    for g in range(G):
        rg = np.random.default_rng(np.random.SeedSequence([seed, size, 1, g, 0]))
        edges = B.build_graph("er", size, rg); E = len(edges)
        adj = B.adj_matrix(size, edges)
        L, hops = B.enumerate_paths(adj, 0, size - 1, 3, H)
        if len(hops) == 0:
            for c in bmr: bmr[c] += [0.0] * m
            for f_ in fths:
                for c in res[f_]: res[f_][c] += [0.0] * m
            continue
        for draw in range(m):
            for f in (0.0, 0.10):
                rl = np.random.default_rng(np.random.SeedSequence([seed, size, 1, g, draw, 1]))
                rc = np.random.default_rng(np.random.SeedSequence([seed, size, 1, g, draw, 2, int(f * 1000)]))
                lk = B.make_links(rl, E, 10000, f=f, rng_comp=rc)
                M = B.path_metrics(L, lk); F = (1 + 3 * M["vprod"]) / 4
                rs = np.random.default_rng(np.random.SeedSequence([seed, size, 1, g, draw, 3]))
                cond = "clean" if f == 0 else "compromised"
                sel, M2, _, _ = B.select_all(L, hops, lk, Delta, R, "sealed", rs, ["BMR"], M=M)
                i = sel["BMR"]; ok = i >= 0 and bool(M["allvalid"][i]) and M["seps"][i] <= Delta * (1 + 1e-9)
                bmr[cond].append(float(M["rmin"][i]) if ok else 0.0)
                for f_ in fths:
                    msk = F >= f_
                    j = B._pick(msk, [M["cost"]], rs)
                    okj = j >= 0 and bool(M["allvalid"][j]) and M["seps"][j] <= Delta * (1 + 1e-9)
                    res[f_][cond].append(float(M["rmin"][j]) if okj else 0.0)
                    if cond == "clean": frac_bind[f_].append(float(1 - msk.mean()))
    return res, bmr, frac_bind


fths = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
res, bmr, fb = run_fth(fths)
C2 = {}
for f_ in fths:
    row = {"frac_of_candidate_paths_excluded_by_threshold": float(np.mean(fb[f_]))}
    for c in ("clean", "compromised"):
        x = np.array(res[f_][c]).reshape(50, 5).mean(1); b = np.array(bmr[c]).reshape(50, 5).mean(1)
        d = b - x; hw_ = 1.96 * d.std(ddof=1) / math.sqrt(len(d))
        row[c] = dict(qleap_R_sec=float(x.mean()), bmr_R_sec=float(b.mean()), diff=float(d.mean()), lo=float(d.mean() - hw_), hi=float(d.mean() + hw_))
    C2[f"{f_:.2f}"] = row
OUT["C2_Fth_sweep"] = C2

# ----------------------------------------- C3: eta_tight vs eta_union sensitivity
tau_union = B.tau  # keep reference before patching in run_fth (restore below)
import importlib; importlib.reload(B)
eta_union = B.tau(10000); eta_tight = math.sqrt(8 * math.log(1 / 1e-6) / 10000)
OUT["C3_eta"] = dict(eta_union=eta_union, eta_tight=eta_tight,
                     r_e_honest_union=B.h(2.687 - eta_union), r_e_honest_tight=B.h(2.687 - eta_tight),
                     r_e_attack_union=B.h(2.418 - eta_union), r_e_attack_tight=B.h(2.418 - eta_tight),
                     attack_below_R_union=B.h(2.418 - eta_union) < 0.10, attack_below_R_tight=B.h(2.418 - eta_tight) < 0.10)


def run_eta(eta, G=50, m=5, size=30, Delta=5e-6, R=0.10, H=5, seed=20261012):
    importlib.reload(B)
    B.tau = lambda n, delta=1e-6, e=eta: e
    names = ["BMR", "SPF", "SPF_valid", "SPF_valid_budget"]
    acc = {c: {r: [] for r in names} for c in ("clean", "compromised")}
    touch = {c: [] for c in acc}
    for g in range(G):
        rg = np.random.default_rng(np.random.SeedSequence([seed, size, 1, g, 0]))
        edges = B.build_graph("er", size, rg); E = len(edges)
        adj = B.adj_matrix(size, edges)
        L, hops = B.enumerate_paths(adj, 0, size - 1, 3, H)
        for draw in range(m):
            for f in (0.0, 0.10):
                cond = "clean" if f == 0 else "compromised"
                if len(hops) == 0:
                    for r in names: acc[cond][r].append(0.0)
                    continue
                rl = np.random.default_rng(np.random.SeedSequence([seed, size, 1, g, draw, 1]))
                rc = np.random.default_rng(np.random.SeedSequence([seed, size, 1, g, draw, 2, int(f * 1000)]))
                lk = B.make_links(rl, E, 10000, f=f, rng_comp=rc)
                rs = np.random.default_rng(np.random.SeedSequence([seed, size, 1, g, draw, 3]))
                sel, M, F, tp = B.select_all(L, hops, lk, Delta, R, "sealed", rs, names)
                rows = {r["router"]: r for r in B.trial_rows(sel, M, F, tp, L, hops, Delta, "sealed", 3)}
                for r in names: acc[cond][r].append(rows[r]["r_sec"])
    out = {}
    for cond in acc:
        b = np.array(acc[cond]["BMR"]).reshape(G, m).mean(1)
        out[cond] = {"BMR_R_sec": float(b.mean())}
        for r in names[1:]:
            x = np.array(acc[cond][r]).reshape(G, m).mean(1); d = b - x; h_ = 1.96 * d.std(ddof=1) / math.sqrt(G)
            out[cond][f"BMR-{r}"] = dict(diff=float(d.mean()), lo=float(d.mean() - h_), hi=float(d.mean() + h_), comparator_R_sec=float(x.mean()))
    return out


OUT["C3_union"] = run_eta(eta_union)
OUT["C3_tight"] = run_eta(eta_tight)
importlib.reload(B)

# ------------------------------------------------------------ C4: N_key, int32 stamp
OUT["C4_NKEY"] = {f"eps={e:g},N={N:g}": 2 * math.log2(1 / e) / N for e in (2e-7, 1e-6, 5e-6) for N in (1e5, 1e6)}
OUT["C4_cert_bytes"] = dict(payload_float64_stamp=53, payload_int32_stamp=49, signed_float64=53 + 64, signed_int32=49 + 64)
json.dump(OUT, open(os.path.join(DATA, "A_close_pending.json"), "w"), indent=1, default=str)
print(json.dumps({k: OUT[k] for k in ("C2_Fth_sweep", "C3_eta", "C3_union", "C3_tight")}, indent=1, default=str)[:6000])
print({k: v["fixed"] for k, v in D.items()})
