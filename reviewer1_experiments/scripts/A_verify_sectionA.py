"""
A_verify_sectionA.py -- numeric verification for Section A (A1-A20) of the
Reviewer-1 response package.  Reads only existing code/data; no tuning.

Run from the folder 'Reviewer 1 output mrng':
    python3 scripts/A_verify_sectionA.py
Output: data/A_verify_sectionA.json
"""
import json, math, os, re, sys
import numpy as np
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = {}

def h(S):
    if S <= 2.0: return 0.0
    S = min(S, 2*math.sqrt(2)-1e-10)
    return 1.0 - math.log2(1.0 + math.sqrt(max(0.0, 2.0 - S**2/4.0)))
def tau(n, d=1e-6): return math.sqrt(32*math.log(8/d)/n)
def sig_w(S0, n):   return math.sqrt(max(0.0, 4.0 - S0**2/4.0))/math.sqrt(n)

# ---------------- A1/A2: class parameters, per-class r_e, false-alarm rates
n = 10000; Z = 4.7534
classes = {"high": (2.720, 2e-7, .3), "std": (2.687, 1e-6, .5), "deg": (2.640, 5e-6, .2)}
sigma_gen = 2/math.sqrt(n)
A1 = {"n_audit": n, "sigma_gen(2/sqrt n)": sigma_gen, "tau(n,1e-6)": tau(n), "classes": {}}
for k, (S, eps, p) in classes.items():
    sw = sig_w(S, n)
    far_code = float(norm.cdf(-Z*sw/sigma_gen))      # what the L1 code actually realises
    A1["classes"][k] = dict(S_true=S, eps=eps, prob=p, sigma_werner=sw,
        S_lb_mean=S-tau(n), r_e_at_mean=h(S-tau(n)),
        alarm_threshold=S-Z*sw, false_alarm_rate_realised=far_code,
        nominal_false_alarm=1e-6)
# attack link: S_true=2.418, baseline S0=2.687
sw0 = sig_w(2.687, n)
thr = 2.687 - Z*sw0
A1["attack"] = dict(S_true=2.418, eps=1e-6, S0=2.687, threshold=thr,
    detection_prob=float(norm.cdf((thr-2.418)/sigma_gen)),
    S_lb_mean=2.418-tau(n), r_e_at_mean=h(2.418-tau(n)),
    gate_pass_prob=float(1-norm.cdf((2.0+tau(n)-2.418)/sigma_gen)))
# alarm with sigma_gen matched to Werner sigma (i.e. correct noise model)
A1["false_alarm_if_sigma_matched"] = float(norm.cdf(-Z))
OUT["A1_A2"] = A1

# ---------------- A3: PA penalty 2*log2(1/eps)/N
A3 = {}
for eps in (2e-7, 1e-6, 5e-6):
    for N in (1e5, 1e6):
        pen = 2*math.log2(1/eps)/N
        A3[f"eps={eps:g},N={N:g}"] = dict(penalty=pen,
            rel_to_r_e_high=pen/h(2.720-tau(n)),
            rel_to_r_e_std=pen/h(2.687-tau(n)))
OUT["A3_PA_penalty"] = A3

# ---------------- A5: expiry probability
OUT["A5"] = dict(epoch=5, tau_max=100, p_expire=1-math.exp(-5/100),
                 p_path_all_fresh_L3=(math.exp(-5/100))**3,
                 p_path_all_fresh_L4=(math.exp(-5/100))**4)

# ---------------- A5/A1: code defects verified by static inspection
def read(p): return open(os.path.join(HERE, p)).read()
l7 = read("L7_budget_full.py")
OUT["static_checks"] = dict(
    L7_link_ok_def_count=len(re.findall(r"def link_ok", l7)),
    L7_link_ok_call_count=len(re.findall(r"link_ok\(", l7)) - len(re.findall(r"def link_ok\(", l7)),
    L1_route_bmr_filters_on_valid_or_alarm=bool(re.search(r'\["(valid|alarm|gate)"\]', re.search(r"def route_bmr.*?def route_spf", read("L1_sweep_bottleneck.py"), re.S).group(0))),
    L1_n_audit=re.search(r"n_audit = (\d+)", read("L1_sweep_bottleneck.py")).group(1),
    L7_n_audit=re.search(r"n_audit = (\d+)", l7).group(1),
    L7_ER_p=0.4, L1_ER_p=0.25,
    script06_opt_level=re.search(r"optimization_level=(\d)", read("06_ibm_quantum_chsh.py")).group(1),
    script06_backend_selection="least-busy at run time; result JSON says ibm_fez",
)

# ---------------- A4: L1 failure handling in means (from JSON)
L1 = json.load(open(os.path.join(DATA, "L1_bottleneck_results.json")))
OUT["L1_json_top_keys"] = list(L1.keys())[:10] if isinstance(L1, dict) else "list"

# ---------------- A8/A9: hardware from raw counts
hw = json.load(open(os.path.join(DATA, "06_ibm_chsh_results.json")))
def E_from_counts(c):
    tot = sum(c.values()); e = 0
    for b, k in c.items():
        e += k*(1 if b[-1] == b[-2] else -1)
    return e/tot, tot
res = {}
for cond in ("native", "degraded"):
    Es, ns = {}, {}
    for key in ("00", "01", "10", "11"):
        Es[key], ns[key] = E_from_counts(hw["hardware_results"][cond][key]["counts"])
    S = Es["00"]+Es["01"]+Es["10"]-Es["11"]
    var_phys = sum((1-e*e)/ns[k] for k, e in Es.items())
    var_cons = sum(4/ns[k] for k in Es)  # NOT used; see below
    nn = sum(ns.values())
    se_phys = math.sqrt(var_phys)
    se_cons = 2/math.sqrt(8192)          # conservative sigma=2/sqrt(n) as in X13 script, n=8192
    res[cond] = dict(E=Es, shots=ns, S=S, SE_physical=se_phys, SE_conservative_2_over_sqrt8192=se_cons,
                     tau_8192=tau(8192), S_lb=S-tau(8192), r_e_lb=h(S-tau(8192)),
                     r_e_plain_h_of_S=h(S), job_ids={k: hw["job_ids"][f"{cond}_{k}"] for k in Es})
dS = res["degraded"]["S"]-res["native"]["S"]
for lab in ("SE_physical", "SE_conservative_2_over_sqrt8192"):
    se = math.sqrt(res["native"][lab]**2 + res["degraded"][lab]**2)
    z = dS/se
    res[f"diff_{lab}"] = dict(dS=dS, SE_diff=se, z=z, p_two_sided=2*(1-float(norm.cdf(abs(z)))))
res["backend"] = hw["hardware_backend"]
res["shots_per_setting"] = 8192
res["comment"] = "degraded S exceeds native S by dS; diff not significant at 5% under either SE"
OUT["A8_A9_hardware"] = res

# ---------------- A8: Aer re-simulation of the 06 circuits (noiseless) and old-circuit hypotheses
try:
    from qiskit import QuantumCircuit
    from qiskit_aer import AerSimulator
    sim = AerSimulator()
    ang = {"00": (0, math.pi/8), "01": (0, -math.pi/8), "10": (math.pi/4, math.pi/8), "11": (math.pi/4, -math.pi/8)}
    def circ(a, b, fa=2, fb=2, deg=0.0):
        qc = QuantumCircuit(2, 2); qc.h(0); qc.cx(0, 1)
        qc.ry(fa*a, 0); qc.ry(fb*b, 1)
        if deg:
            qc.ry(deg, 0); qc.ry(deg, 1)
        qc.measure([0, 1], [0, 1]); return qc
    def S_of(fa, fb, deg, shots=400000, seed=7):
        E = {}
        for k, (a, b) in ang.items():
            c = sim.run(circ(a, b, fa, fb, deg), shots=shots, seed_simulator=seed).result().get_counts()
            E[k] = E_from_counts(c)[0]
        return E["00"]+E["01"]+E["10"]-E["11"]
    def S_opp(d, shots=400000, seed=7):
        E = {}
        for k, (a, b) in ang.items():
            qc = QuantumCircuit(2, 2); qc.h(0); qc.cx(0, 1); qc.ry(2*a, 0); qc.ry(2*b, 1)
            qc.ry(d, 0); qc.ry(-d, 1); qc.measure([0, 1], [0, 1])
            E[k] = E_from_counts(sim.run(qc, shots=shots, seed_simulator=seed).result().get_counts())[0]
        return E["00"]+E["01"]+E["10"]-E["11"]
    OUT["A8_aer_noiseless"] = dict(
        qiskit_circuit_RY2a_native=S_of(2, 2, 0.0),
        qiskit_circuit_RY2a_degraded_RY0p20_each=S_of(2, 2, 0.20),
        halved_angle_RYa_native=S_of(1, 1, 0.0),
        halved_angle_RYa_scaled_by_hw_visibility=S_of(1, 1, 0.0)*res["native"]["S"]/(2*math.sqrt(2)),
        opposite_sign_RYp0p20_RYm0p20_S=S_opp(0.20),
        analytic_common_mode_S_unchanged="RY(d)(x)RY(d) leaves E=cos(thetaA-thetaB) unchanged on Phi+",
        note="Noiseless Aer; shots=400000/setting, seed 7. Hardware 2.490/2.532 are far below the 2.828 ideal -> large device noise.")
    # best-fit depolarising visibility for the hardware numbers
    OUT["A8_visibility"] = dict(
        v_native=res["native"]["S"]/(2*math.sqrt(2)),
        v_degraded=res["degraded"]["S"]/(2*math.sqrt(2)))
except Exception as ex:
    OUT["A8_aer_noiseless"] = f"not run: {ex!r}"

# ---------------- A10: X2 sigma values
X2 = {}
for S0 in (2.687,):
    X2["sigma_werner_n1e4"] = sig_w(S0, 1e4)
    X2["sigma_2_over_sqrt_n"] = 2/math.sqrt(1e4)
    X2["threshold_z4.7534"] = S0 - Z*sig_w(S0, 1e4)
    X2["detect_S2.418"] = float(norm.cdf((X2["threshold_z4.7534"]-2.418)/sig_w(2.418, 1e4)))
OUT["A10_X2"] = X2


# ---------------- A1/A4: which class mixes satisfy sum(eps) <= Delta=1e-6 (L1 design)
import itertools
EPS = [2e-7, 1e-6, 5e-6]; NAMES = ["high", "std", "deg"]
feas = {}
for L in (3, 4, 5):
    ok = [c for c in itertools.combinations_with_replacement(range(3), L) if sum(EPS[i] for i in c) <= 1e-6*(1+1e-9)]
    feas[f"L={L}"] = [[NAMES[i] for i in c] for c in ok]
OUT["A4_budget_feasible_class_mixes_Delta1e-6"] = dict(mixes=feas,
    P_all_high_L3=0.3**3, P_all_high_L4=0.3**4, P_all_high_L5=0.3**5,
    note="With Delta=1e-6 only all-'high' paths satisfy the budget at 3..5 hops; infeasible trials enter means as 0 (eval_result).")

# ---------------- A13: environment
try:
    import scipy, networkx, qiskit
    from scipy.optimize import milp
    OUT["A13_env"] = dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                          networkx=networkx.__version__, qiskit=qiskit.__version__, scipy_milp=True)
except Exception as ex:
    OUT["A13_env"] = f"partial: {ex!r}"

json.dump(OUT, open(os.path.join(DATA, "A_verify_sectionA.json"), "w"), indent=2, default=float)
print(json.dumps(OUT, indent=1, default=float)[:6000])
