"""
R7_aer_chsh.py -- R7: noiseless-Aer CHSH check of the circuit variants + per-setting correlators / SE from the raw IBM counts.
Run with:  ../../expr_r2/.venv/bin/python R7_aer_chsh.py     (needs qiskit + qiskit_aer)

Circuit family (same as 06_ibm_quantum_chsh.py): H(0); CX(0,1); RY(ca*2a + da) on q0; RY(cb*2b + db) on q1; measure.
Settings (a,b): (0,pi/8), (0,-pi/8), (pi/4,pi/8), (pi/4,-pi/8);  S = E11 + E12 + E21 - E22.
Variants:
  corrected      ca=cb=1, d=0                    ideal 2*sqrt(2)=2.828
  degraded       corrected + RY(+0.20) on BOTH qubits (as in 06_ibm_quantum_chsh.py)   -> a common shift = null perturbation
  opposite_sign  corrected + RY(+0.20) on q0, RY(-0.20) on q1                           -> relative shift 0.4 rad, ideal 2.828*cos(0.4)
  old_halfboth / old_halfA / old_halfB : angles used as Bloch angles (RY(a) instead of RY(2a)) on both / Alice / Bob
    (tested ONLY as candidate explanations of the old 2.574/2.579; the old circuit itself is not recoverable)
"""
import json, math, sys, numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

SHOTS = 8192; REPS = 20
ANG = {"11": (0, math.pi / 8), "12": (0, -math.pi / 8), "21": (math.pi / 4, math.pi / 8), "22": (math.pi / 4, -math.pi / 8)}
VAR = {"corrected": (2, 2, 0, 0), "degraded": (2, 2, .2, .2), "opposite_sign": (2, 2, .2, -.2),
       "old_halfboth": (1, 1, 0, 0), "old_halfA": (1, 2, 0, 0), "old_halfB": (2, 1, 0, 0)}


def circ(a, b, v):
    ca, cb, da, db = VAR[v]
    qc = QuantumCircuit(2, 2); qc.h(0); qc.cx(0, 1)
    qc.ry(ca * a + da, 0); qc.ry(cb * b + db, 1); qc.measure([0, 1], [0, 1]); return qc


def corr(counts):
    tot = sum(counts.values()); return sum(c * (1 if k[-1] == k[-2] else -1) for k, c in counts.items()) / tot


def main():
    sim = AerSimulator(); out = {}
    for v in VAR:
        Ss = []; Es = {k: [] for k in ANG}
        for rep in range(REPS):
            E = {}
            for k, (a, b) in ANG.items():
                qc = transpile(circ(a, b, v), sim)
                cnt = sim.run(qc, shots=SHOTS, seed_simulator=10_000 * rep + int(k) + 1000 * list(VAR).index(v)).result().get_counts()
                E[k] = corr(cnt); Es[k].append(E[k])
            Ss.append(E["11"] + E["12"] + E["21"] - E["22"])
        Ss = np.array(Ss)
        out[v] = dict(S_mean=float(Ss.mean()), S_sd=float(Ss.std(ddof=1)), S_min=float(Ss.min()), S_max=float(Ss.max()),
                      E_mean={k: float(np.mean(x)) for k, x in Es.items()}, shots=SHOTS, reps=REPS, S_all=[float(x) for x in Ss])
        print(f"{v:14s} S = {Ss.mean():.4f} ± {Ss.std(ddof=1):.4f}  [{Ss.min():.3f}, {Ss.max():.3f}]  E={ {k: round(float(np.mean(x)),3) for k,x in Es.items()} }", flush=True)
    # ideal (analytic) values from the statevector, for reference
    from qiskit.quantum_info import Statevector
    for v in VAR:
        S = 0
        for sgn, (k, (a, b)) in zip((1, 1, 1, -1), ANG.items()):
            qc = circ(a, b, v); qc.remove_final_measurements()
            p = Statevector(qc).probabilities()
            S += sgn * (p[0] + p[3] - p[1] - p[2])
        out[v]["S_ideal"] = float(S); print(f"ideal {v:14s} {S:.4f}")
    json.dump(out, open("../data/R7_aer_noiseless.json", "w"), indent=1)

    # ---- per-setting correlators / SE from the raw IBM counts
    d = json.load(open("../data/06_ibm_chsh_results.json"))
    hw = {}
    for lab in ("native", "degraded"):
        rows = {}
        for k in ("00", "01", "10", "11"):
            c = d["hardware_results"][lab][k]["counts"]; n = sum(c.values()); E = corr(c)
            p_eq = (c.get("00", 0) + c.get("11", 0)) / n
            se = math.sqrt(max(1e-12, 1 - E * E) / n)                # physical (binomial) SE of a +-1 correlator
            rows[k] = dict(E=E, n=n, se_binom=se, counts=c)
        sign = {"00": 1, "01": 1, "10": 1, "11": -1}
        S = sum(sign[k] * rows[k]["E"] for k in rows)
        se_phys = math.sqrt(sum(rows[k]["se_binom"] ** 2 for k in rows))
        se_cons = math.sqrt(4 / 8192)                                  # worst case |E|=0 for each of the four correlators
        hw[lab] = dict(settings=rows, S=S, se_physical=se_phys, se_conservative=se_cons)
        print(lab, "S =", round(S, 4), "SE physical", round(se_phys, 4), "conservative", round(se_cons, 4))
    diff = hw["degraded"]["S"] - hw["native"]["S"]
    sd_phys = math.sqrt(hw["degraded"]["se_physical"] ** 2 + hw["native"]["se_physical"] ** 2)
    sd_cons = math.sqrt(hw["degraded"]["se_conservative"] ** 2 + hw["native"]["se_conservative"] ** 2)
    hw["difference_degraded_minus_native"] = dict(diff=diff, sd_physical=sd_phys, z_physical=diff / sd_phys, sd_conservative=sd_cons, z_conservative=diff / sd_cons)
    print("diff", round(diff, 4), "z_phys", round(diff / sd_phys, 2), "z_cons", round(diff / sd_cons, 2))
    json.dump(hw, open("../data/R7_hardware_recomputed.json", "w"), indent=1)


if __name__ == "__main__":
    main()
