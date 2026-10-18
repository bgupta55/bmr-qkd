"""
Q8 (part 1): repeated IBM QPU CHSH runs with RAW COUNTS saved.
================================================================
Addresses Reviewer 2 / pending_v7 Q8: the earlier hardware JSON stored only
correlators E_ab; raw counts were lost, and a single run gave S(0)=2.54.
Here we run R independent jobs of the SAME fixed circuits (E6v2 design,
RY(2*theta), ibm_marrakesh qubits (6,5), 8192 shots per setting) and store the
raw counts, the backend name, job ids, timestamp and the calibration numbers of
the two qubits reported by the backend at submission time.

Credentials: QiskitRuntimeService() reads the account saved by the user in
~/.qiskit/qiskit-ibm.json (or env IBM_QUANTUM_TOKEN). Nothing is written here.

Usage: <qiskit venv python> experiments/q8_hw_repeats.py [R=3]
Open-plan QPU minutes are limited: each job is 12 circuits x 8192 shots.
"""
import sys, os, json, time, datetime, math
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
from qiskit import QuantumCircuit
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler

R = int(sys.argv[1]) if len(sys.argv) > 1 else 3
SHOTS = 8192
PHIS = [0.0, 0.2, 0.5]
PAIR = (6, 5)
BACKEND = "ibm_marrakesh"
THETA_A = [0.0, math.pi / 4]
THETA_B = [math.pi / 8, -math.pi / 8]
OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "results", "q8_hw_repeats.json")


def circ(a, b, phi):
    qc = QuantumCircuit(2, 2)
    qc.h(0); qc.cx(0, 1)
    qc.ry(2 * (THETA_A[a] + phi), 0)
    qc.ry(2 * THETA_B[b], 1)
    qc.measure([0, 1], [0, 1])
    return qc


svc = QiskitRuntimeService()
be = svc.backend(BACKEND)
calib = {}
try:
    props = be.properties()
    for q in PAIR:
        calib[f"q{q}"] = {
            "T1_us": (props.t1(q) or 0) * 1e6, "T2_us": (props.t2(q) or 0) * 1e6,
            "readout_error": props.readout_error(q)}
    calib["note"] = "backend.properties() at submission time"
except Exception as ex:  # calibration data is optional
    calib["error"] = str(ex)

tags, circuits = [], []
for phi in PHIS:
    for a in (0, 1):
        for b in (0, 1):
            tags.append((phi, a, b)); circuits.append(circ(a, b, phi))
pm = generate_preset_pass_manager(backend=be, optimization_level=3, initial_layout=list(PAIR))
isa = pm.run(circuits)

runs = []
if os.path.exists(OUT):
    runs = json.load(open(OUT)).get("runs", [])
for rep in range(R):
    t0 = time.time()
    job = Sampler(mode=be).run(isa, shots=SHOTS)
    jid = job.job_id()
    print(f"[rep {rep}] job {jid} submitted", flush=True)
    res = job.result()
    raw = {}
    for (phi, a, b), pub in zip(tags, res):
        c = pub.data.c.get_counts()
        raw[f"{phi}|{a}|{b}"] = {k: int(v) for k, v in c.items()}
    runs.append({"job_id": jid, "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 "wall_s": time.time() - t0, "raw_counts": raw})
    print(f"[rep {rep}] done wall={time.time()-t0:.1f}s", flush=True)
    json.dump({"backend": BACKEND, "qubit_pair": list(PAIR), "shots_per_setting": SHOTS,
               "phis": PHIS, "calibration": calib, "runs": runs}, open(OUT, "w"), indent=1)
print("saved", OUT)
