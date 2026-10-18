"""
run_all.py -- regenerate every Reviewer-2 number from scratch WITHOUT touching a QPU.

  python3 run_all.py            # all classical experiments (~15-25 min on 12 cores)
  python3 run_all.py --fast     # skip the slowest scripts (Q9 timing, Q14, Q15 drift)

NOT run here (they submit jobs to IBM Quantum and need qiskit + a locally stored account):
  experiments/e6_hardware.py, experiments/e6v2_hardware_fixed.py, experiments/q8_hw_repeats.py
  -> run manually with  ../expr_r2/.venv/bin/python experiments/q8_hw_repeats.py
  The stored raw counts (data/results/q8_hw_repeats.json) are enough to regenerate every table:
  `python3 experiments/q8_hardware_reconcile.py` needs no QPU time.
Credentials are never read or written by this file.
"""
import sys, os, subprocess, json, time, hashlib

BASE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(BASE, "experiments")
RES = os.path.join(BASE, "data", "results")
FAST = "--fast" in sys.argv

# (id, script, question/comment, slow?)
SCRIPTS = [
    ("T-integ", "../tests/test_integration.py",  "Q7  simulator == build_certificate (5 tests)",      False),
    ("E3",   "e3_certificate_tests.py",   "C4  certificate unit tests (9)",                    False),
    ("E1",   "e1_coverage.py",            "C2  eta coverage (200,000 MC trials)",              False),
    ("E2",   "e2_table.py",               "C2  finite-sample table T-A",                       False),
    ("E4",   "e4_paths.py",               "C3  Corollary 1 brute force m<=12, E4b (old denominators)", False),
    ("E5",   "e5_runtime.py",             "C3  old runtime table (superseded by Q9, kept)",    False),
    ("E10",  "e10_ledger.py",             "C1  Algorithm 3 v1 (superseded by Q11, kept)",      False),
    ("Q4",   "q4_collusion.py",           "Q4  collusion experiment",                          False),
    ("Q5",   "q5_evidence_source.py",     "Q5  Bell vs QBER evidence source",                  False),
    ("Q1",   "e7b_v2_corrected.py",       "Q1/Q2  corrected E7b + self-checks",                False),
    ("Q7",   "q7_e0_audit.py",            "Q7  E0 audit (old code vs new)",                    False),
    ("Q8",   "q8_hardware_reconcile.py",  "Q8  hardware reconciliation from stored raw counts", False),
    ("Q9",   "q9_solver_validation.py",   "Q9  exhaustive validation + uncapped timing + E4b", True),
    ("Q11",  "q11_alg3_v2.py",            "Q11 Algorithm 3 v2, T0-T12, timing",                False),
    ("Q13",  "q13_scope_numbers.py",      "Q13 r_DI sweep + Proposition 2 Monte Carlo",        False),
    ("Q14",  "q14_twobudget_v2.py",       "Q14 two-budget solvers (Theorem 4)",                True),
    ("Q15a", "q15_delta_regime.py",       "Q15 Delta/phi regime map + compromised-fraction sweep", False),
    ("Q15b", "q15_drift_robustness.py",   "Q15 drift / correlation robustness, +-25 % sensitivity", True),
    ("Q6a",  "q6_sweep.py",               "Q6  full sweep (200 cells x 1000 trials x 10 routers)", False),
    ("Q6c",  "q6_sweep.py",               "Q6  sensitivity sweep phi=0.9 (env Q6_VARIANT=phi09, set below)", False),
    ("Q6b",  "q6_dynamic.py",             "Q6  dynamic mid-session compromise",                False),
    ("Fig",  "generate_q6_figures.py",    "Q6  figures",                                       False),
    ("Q16",  "generate_macros.py",        "Q16 results_macros.tex (no hand-typed numbers)",    False),
]

summary = []
t_all = time.time()
for sid, script, desc, slow in SCRIPTS:
    if FAST and slow:
        summary.append((sid, "SKIPPED(--fast)", 0.0, desc))
        continue
    path = os.path.normpath(os.path.join(EXP, script))
    print(f"\n=== {sid}: {desc}\n    {path}", flush=True)
    t0 = time.time()
    env = dict(os.environ)
    if sid == "Q6c":
        env["Q6_VARIANT"] = "phi09"
    r = subprocess.run([sys.executable, path], cwd=BASE, env=env)
    status = "PASSED" if r.returncode == 0 else f"FAILED({r.returncode})"
    summary.append((sid, status, time.time() - t0, desc))
    print(f"[{sid}] {status} in {time.time()-t0:.1f}s", flush=True)

print("\n" + "=" * 78)
for sid, st, el, desc in summary:
    print(f"{sid:>8} {st:>16} {el:>8.1f}s  {desc}")
ok = sum(1 for s in summary if s[1] == "PASSED")
ran = sum(1 for s in summary if not s[1].startswith("SKIPPED"))
print(f"\nPassed {ok}/{ran}   total {time.time()-t_all:.0f}s")
with open(os.path.join(RES, "run_all_summary.json"), "w") as f:
    json.dump({"summary": [dict(id=a, status=b, seconds=c, desc=d) for a, b, c, d in summary],
               "passed": ok, "ran": ran}, f, indent=2)
sys.exit(0 if ok == ran else 1)
