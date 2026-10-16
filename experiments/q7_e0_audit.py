"""
Q7 / E0 audit: what did the OLD pipeline (bmr-experiments/, the pre-revision simulator behind the first
submission) do, what did the R1 scripts (expr_r1/) do, and what does the library do now?
All 'old' numbers are computed by IMPORTING the old module, not re-typed.
"""
import sys, os, json, math, hashlib, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from bmr.bounds import eta_tight, eta_union
from bmr.certificate import build_certificate
print("code sha256:", hashlib.sha256(open(__file__, "rb").read()).hexdigest()[:16])
for f in ("bmr/certificate.py", "bmr/bounds.py", "bmr/sim.py", "bmr/route.py"):
    print("  sha256", f, hashlib.sha256(open(os.path.join(ROOT, f), "rb").read()).hexdigest()[:16])

old_path = os.path.abspath(os.path.join(ROOT, "..", "..", "bmr-experiments", "audit", "bell_audit.py"))
spec = importlib.util.spec_from_file_location("old_bell_audit", old_path)
old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)

rows = []
rng = np.random.default_rng(0)
print("\nPass probability of an HONEST link (S_true=2.687) per audit, 20000 audits, delta_old=1e-5, delta_new=1e-6")
print(f"{'n':>6} | old tol  old pass | R1 tau(union,1e-6) R1 pass | new eta(tight)  new pass(library)")
for n in (100, 250, 500, 2500, 10_000):
    N = 20000
    old_pass = 0; new_pass = 0; r1_pass = 0
    sd = 2.0 / math.sqrt(n)
    for i in range(N):
        c = old.bell_audit_certify(2.687, n, tau=0.10, delta_target=1e-5, rng=rng)
        old_pass += (c["h_bound"] > 0)           # old rule: rate>0 after subtracting its own tolerance
        s = float(rng.normal(2.687, math.sqrt(4 - 2.687 ** 2 / 4) / n ** .5))
        new_pass += build_certificate("e", s, n, 1e-6, 10_000, 0.99, 0, 0) is not None
        r1_pass += (s - eta_union(n, 1e-6) > 2.0)
    r = dict(n=n, old_tol=old.hoeffding_confidence_bound(n, 1e-5), old_pass=old_pass / N,
             r1_tau=eta_union(n, 1e-6), r1_pass=r1_pass / N,
             new_eta=eta_tight(n, 1e-6), new_pass=new_pass / N)
    rows.append(r)
    print(f"{n:>6} | {r['old_tol']:.4f}  {r['old_pass']:.3f}   | {r['r1_tau']:.4f}   {r['r1_pass']:.3f}      | {r['new_eta']:.4f}     {r['new_pass']:.3f}")

# detection rule
print("\nOld 'valid' rule: S_hat >= 2 + tau_param (tau_param=0.10 -> 2.10), NOT the midpoint 2.5525 and NOT the certificate rule.")
# Old epsilon / gamma definitions
S = 2.687; n = 10_000
oc = old.bell_audit_certify(S, n, 0.10, 1e-5, np.random.default_rng(1))
print("old eps_e (n=1e4, delta=1e-5):", oc["eps_e"], "  old gamma_e:", oc["gamma_e"], "  old h:", oc["h_bound"])
nc = build_certificate("e", S, n, 1e-6, 10_000, 0.99, 0, 1)
print("new eps_e:", nc.eps, " new log2_gamma:", nc.log2_gamma, " new h:", nc.h, " w=eps+delta:", nc.w)
json.dump({"rows": rows}, open(os.path.join(ROOT, "data", "results", "q7_e0_audit.json"), "w"), indent=1)
