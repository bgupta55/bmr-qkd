"""Q4 add-on: the Q6 sweep protocol on the degenerate topologies line and star (submitted Table 6 had both).
Same one_cell(), same parameters/seeds-structure as q6_sweep.py; the 200-cell main sweep is NOT touched.
Cells: {line, star} x {10,20,30,50} x n_audit {500,10^4} x {clean, compromised} x {sealed, public-XOR}."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from multiprocessing import Pool
import numpy as np
import q6_sweep as Q
from q6_common import ROUTER_NAMES

def main():
    t0 = time.time(); rng = np.random.default_rng(962)
    specs = [(t, n, na, c, m, int(rng.integers(1e9))) for t in ("line", "star") for n in Q.SIZES
             for na in (500, 10_000) for c in ("clean", "compromised") for m in ("sealed", "public-XOR")]
    rows = []
    with Pool(max(1, (os.cpu_count() or 4) - 1)) as p:
        for i, (r, ra, fe) in enumerate(p.imap(Q.one_cell, specs, chunksize=1)):
            rows += r
            if i % 8 == 0: print(f"  {i+1}/{len(specs)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(dict(experiment="Q6-sweep-degenerate", trials=Q.TRIALS, rows=rows),
              open(os.path.join(Q.HERE, "..", "data", "results", "q6_sweep_degenerate.json"), "w"))
    print("saved", flush=True)
if __name__ == "__main__":
    main()
