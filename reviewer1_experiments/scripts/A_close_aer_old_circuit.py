"""A_close_aer_old_circuit.py -- Section A closure, item A8 (provenance of the old 2.574 / 2.579).

Question: R7 showed that no NOISELESS candidate reaches 2.574/2.579, but two candidates (opposite_sign, old_halfB) sit at
2.606 noiselessly, only 0.03 above the old values.  Hardware noise lowers S, so the noiseless argument alone is not a proof.
This script runs every candidate under the ibm_fez noise model (FakeFez snapshot) at the OLD shot configuration
(2,000 shots x 3 repetitions, as described in the old manuscript Sec. 8.13) and with 8,192 shots, and also builds the old
'degraded' circuit as described there (two extra CNOT pairs).  Nothing is tuned: the result is reported whatever it is.

Run:  ../../expr_r2/.venv/bin/python A_close_aer_old_circuit.py     (from scripts/)
Output: ../data/A_close_aer_old_circuit.json
"""
import json, math, numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit_ibm_runtime.fake_provider import FakeFez
from R7_aer_chsh import ANG, VAR, corr

fb = FakeFez(); nm = NoiseModel.from_backend(fb)
sim = AerSimulator(noise_model=nm, basis_gates=nm.basis_gates)
CFG = [(2000, 3, 20), (8192, 1, 20)]


def circ(a, b, v):
    if v == "native_2cx_pairs":                 # old 'degraded' as described: corrected + two CNOT pairs (identity in ideal case)
        qc = QuantumCircuit(2, 2); qc.h(0); qc.cx(0, 1)
        for _ in range(2):
            qc.barrier(); qc.cx(0, 1); qc.cx(0, 1)
        qc.barrier(); qc.ry(2 * a, 0); qc.ry(2 * b, 1); qc.measure([0, 1], [0, 1]); return qc
    ca, cb, da, db = VAR[v]
    qc = QuantumCircuit(2, 2); qc.h(0); qc.cx(0, 1)
    qc.ry(ca * a + da, 0); qc.ry(cb * b + db, 1); qc.measure([0, 1], [0, 1]); return qc


variants = ["corrected", "degraded", "opposite_sign", "old_halfB", "native_2cx_pairs"]   # old_halfA/halfboth are 2.39-2.42 noiselessly (cannot reach 2.57)
out = {"noise_model": "FakeFez (ibm_fez calibration snapshot, qiskit-ibm-runtime 0.50), qubits (0,1), optimisation_level=1",
       "old_reported": {"native": [2.574, 0.040], "degraded": [2.579, 0.013], "aer_matched": [2.598, 0.025]},
       "variants_run": variants, "excluded_noiseless_below_2.45": ["old_halfA", "old_halfboth"]}
tr = {v: [transpile(circ(a, b, v), backend=fb, optimization_level=1, initial_layout=[0, 1]) for k, (a, b) in ANG.items()] for v in variants}
keys = list(ANG)
for shots, reps, batches in CFG:
    key = f"shots{shots}"; out[key] = {}
    for vi, v in enumerate(variants):
        circs = tr[v] * (reps * batches)                       # one batched run: Aer parallelises over circuits
        res = sim.run(circs, shots=shots, seed_simulator=12345 + vi).result()
        S = []
        for j in range(reps * batches):
            E = {k: corr(res.get_counts(4 * j + q)) for q, k in enumerate(keys)}
            S.append(E["11"] + E["12"] + E["21"] - E["22"])
        meas = np.array(S).reshape(batches, reps).mean(1)   # each measurement = mean of `reps` repetitions
        out[key][v] = dict(S_mean=float(meas.mean()), S_sd_between_measurements=float(meas.std(ddof=1)),
                           S_min=float(meas.min()), S_max=float(meas.max()), n_measurements=batches, reps_per_measurement=reps,
                           frac_within_0p04_of_2p574=float(np.mean(np.abs(meas - 2.574) < 0.04)),
                           frac_ge_2p574=float(np.mean(meas >= 2.574)))
        print(key, f"{v:18s} S = {meas.mean():.4f} +- {meas.std(ddof=1):.4f}  P(|S-2.574|<0.04) = {out[key][v]['frac_within_0p04_of_2p574']:.2f}", flush=True)
json.dump(out, open("../data/A_close_aer_old_circuit.json", "w"), indent=1)
