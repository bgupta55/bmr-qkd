#!/usr/bin/env python3
"""
run_all.py — Execute all expr_r1 experiments in order.
Run from the REVISION_OSN_12OCT2026 directory or workspace root.
"""
import subprocess, sys, os, time

EXPR_DIR = os.path.join(os.path.dirname(__file__), "REVISION_OSN_12OCT2026", "expr_r1")

scripts = [
    "01_certificate_pipeline.py",
    "02_binding_regime.py",
    "03_werner_fidelity.py",
    "04_cluster_stats.py",
    "05_runtime_scaling.py",
    "06_ibm_quantum_chsh.py",
    "07_ablation_binding.py",
    "08_prop2_missrate.py",
]

results = {}
for script in scripts:
    path = os.path.join(EXPR_DIR, script)
    print(f"\n{'='*60}")
    print(f"Running: {script}")
    print('='*60)
    t0 = time.time()
    proc = subprocess.run([sys.executable, path], capture_output=True, text=True)
    elapsed = time.time() - t0
    status = "OK" if proc.returncode == 0 else "FAILED"
    results[script] = {"status": status, "elapsed_s": round(elapsed,2)}
    print(proc.stdout)
    if proc.stderr:
        print("STDERR:", proc.stderr[:1000])
    print(f"[{status}] {script} in {elapsed:.1f}s")

print("\n\n=== Summary ===")
for s, r in results.items():
    print(f"  {r['status']:6s}  {r['elapsed_s']:6.1f}s  {s}")
