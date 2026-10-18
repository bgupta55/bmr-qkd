"""R9_aer_noise.py -- Aer with the ibm_fez backend noise model (FakeFez snapshot from qiskit-ibm-runtime), n_shots = n = 10^4 per setting, 20 repeats.
Run with ../../expr_r2/.venv/bin/python.  Applies Algorithm 1 (S_lb = S_hat - tau, tau(1e4, 1e-6)): gate, r_e, and the alarm test vs S0 = 2.687."""
import json, math, numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit_ibm_runtime.fake_provider import FakeFez
from R7_aer_chsh import ANG, VAR, circ, corr
N = 10_000; REPS = 20
def h(S): return 0.0 if S <= 2 else 1 - math.log2(1 + math.sqrt(max(0, 2 - S * S / 4)))
tau = math.sqrt(32 * math.log(8 / 1e-6) / N); sigW = math.sqrt(4 - 2.687 ** 2 / 4) / math.sqrt(N); thr = 2.687 - 4.7534 * sigW
fb = FakeFez(); nm = NoiseModel.from_backend(fb); sim = AerSimulator(noise_model=nm, basis_gates=nm.basis_gates)
out = {"backend_noise_model": "FakeFez (ibm_fez calibration snapshot shipped with qiskit-ibm-runtime 0.50)", "n_shots": N, "reps": REPS, "tau": tau, "alarm_threshold": thr}
for v in ("corrected", "degraded"):
    Ss = []
    for rep in range(REPS):
        E = {}
        for k, (a, b) in ANG.items():
            qc = transpile(circ(a, b, v), backend=fb, optimization_level=1, initial_layout=[0, 1])
            cnt = sim.run(qc, shots=N, seed_simulator=777 * rep + int(k)).result().get_counts()
            E[k] = corr(cnt)
        Ss.append(E["11"] + E["12"] + E["21"] - E["22"])
    Ss = np.array(Ss); Slb = Ss - tau
    out[v] = dict(S_mean=float(Ss.mean()), S_sd=float(Ss.std(ddof=1)), S_all=[float(x) for x in Ss], S_lb_mean=float(Slb.mean()),
                  gate_pass=int((Slb >= 2).sum()), r_e_mean=float(np.mean([h(max(0, s)) for s in Slb])), alarm_vs_S0=int((Ss < thr).sum()))
    print(v, {k: (round(x, 4) if isinstance(x, float) else x) for k, x in out[v].items() if k != "S_all"}, flush=True)
d = np.array(out["degraded"]["S_all"]) - np.array(out["corrected"]["S_all"])
out["degraded_minus_corrected"] = dict(mean=float(d.mean()), sd=float(d.std(ddof=1)), se=float(d.std(ddof=1) / math.sqrt(REPS)))
print(out["degraded_minus_corrected"])
json.dump(out, open("../data/R9_aer_noisemodel.json", "w"), indent=1)
