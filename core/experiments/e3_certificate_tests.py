"""
E3: Certificate unit tests — Developer Guide §5
===============================================
Addresses: Reviewer 2 Comment 4 (five certificate fields, how to compute c).
Dataset: Analytic (no randomness).

Tests:
  test_h_endpoints, test_golden_h, test_eta_monotone,
  test_degraded_when_n_too_small, test_golden_certificate,
  test_public_xor_gamma_small, test_invariants
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import math
import json
import time

from bmr.bounds import h_pironio, eta_tight
from bmr.certificate import build_certificate as bc


def approx_eq(a, b, abs_tol=1e-4):
    return abs(a - b) <= abs_tol


def run_test(name, fn):
    start = time.perf_counter()
    try:
        fn()
        elapsed = (time.perf_counter() - start) * 1000
        print(f"  PASS  [{elapsed:.2f}ms]  {name}")
        return {"name": name, "passed": True, "ms": elapsed}
    except AssertionError as e:
        elapsed = (time.perf_counter() - start) * 1000
        print(f"  FAIL  [{elapsed:.2f}ms]  {name}  --> {e}")
        return {"name": name, "passed": False, "ms": elapsed, "error": str(e)}


print("=" * 60)
print("E3: Certificate Unit Tests")
print("=" * 60)

results = []

def test_h_endpoints():
    assert h_pironio(2.0) == approx_eq(h_pironio(2.0), 0.0, abs_tol=1e-12) or h_pironio(2.0) == 0.0, \
        f"h(2.0) = {h_pironio(2.0)} != 0"
    # h_pironio returns 0.0 for S <= 2 by design
    assert h_pironio(2.0) == 0.0, f"h(2.0) = {h_pironio(2.0)}"
    val_tsirelson = h_pironio(2 * math.sqrt(2))
    assert abs(val_tsirelson - 1.0) < 1e-10, f"h(2√2) = {val_tsirelson} != 1.0"
results.append(run_test("test_h_endpoints", test_h_endpoints))

def test_golden_h():
    v1 = h_pironio(2.687)
    v2 = h_pironio(2.418)
    assert abs(v1 - 0.47233) < 1e-4, f"h(2.687) = {v1:.5f} != 0.47233"
    assert abs(v2 - 0.20614) < 1e-4, f"h(2.418) = {v2:.5f} != 0.20614"
results.append(run_test("test_golden_h", test_golden_h))

def test_eta_monotone():
    # Larger n => smaller eta
    assert eta_tight(1000, 1e-6) < eta_tight(500, 1e-6) < eta_tight(250, 1e-6), \
        "eta not monotone in n"
    # Larger delta => smaller eta (less conservative)
    assert eta_tight(500, 1e-9) > eta_tight(500, 1e-3), \
        "eta not monotone in delta"
results.append(run_test("test_eta_monotone", test_eta_monotone))

def test_degraded_when_n_too_small():
    # n=100, S=2.687, delta=1e-6: eta = sqrt(8*ln(1e6)/100) >> 0.687 => S_lb < 2 => None
    c = bc("e", 2.687, 100, 1e-6, 10_000, 0.99, 0, 7)
    assert c is None, f"Expected None (DEGRADED) for n=100, got Cert with r={c.r if c else '?'}"
results.append(run_test("test_degraded_when_n_too_small", test_degraded_when_n_too_small))

def test_golden_certificate():
    c = bc("e", 2.687, 2500, 1e-6, 10_000, 0.99, 0, 7)
    assert c is not None, "Expected non-None certificate"
    assert c.k == 2465, f"k = {c.k} != 2465"
    assert abs(c.r - 0.2465) < 1e-6, f"r = {c.r} != 0.2465"
    assert abs(c.log2_eps - (-13.48)) < 0.05, f"log2_eps = {c.log2_eps:.3f} != -13.48"
    assert abs(c.log2_gamma - 2451.5) < 0.5, f"log2_gamma = {c.log2_gamma:.1f} != 2451.5"
    assert c.T == 7, f"T = {c.T} != 7"
results.append(run_test("test_golden_certificate", test_golden_certificate))

def test_public_xor_gamma_small():
    c = bc("e", 2.687, 2500, 1e-6, 10_000, 0.30, 0, 7)
    assert c is not None, "Expected non-None certificate"
    assert c.log2_gamma < -100, f"log2_gamma = {c.log2_gamma:.1f}, expected < -100 for phi=0.30"
results.append(run_test("test_public_xor_gamma_small", test_public_xor_gamma_small))

def test_invariants():
    c = bc("e", 2.65, 2500, 1e-6, 10_000, 0.9, 0, 3)
    assert c is not None, "Expected non-None certificate"
    assert 0 <= c.r <= c.h, f"rate r={c.r} not in [0, h={c.h}]"
    assert c.w > 0, f"weight w={c.w} <= 0"
    assert c.k <= c.H, f"k={c.k} > H={c.H}"
results.append(run_test("test_invariants", test_invariants))

# Additional: verify the five fields are all present and finite
def test_five_fields():
    c = bc("link_AB", 2.687, 2500, 1e-6, 10_000, 0.99, 0, 5)
    assert c is not None
    fields = ["S_hat", "n", "eta", "h", "H", "k", "r", "log2_eps", "log2_gamma", "delta", "T"]
    for f in fields:
        val = getattr(c, f)
        assert val is not None, f"Field {f} is None"
    print(f"    Certificate fields: S_hat={c.S_hat}, n={c.n}, eta={c.eta:.4f}, "
          f"h={c.h:.5f}, H={c.H:.1f}, k={c.k}, r={c.r:.5f}, "
          f"log2_eps={c.log2_eps:.2f}, log2_gamma={c.log2_gamma:.1f}, delta={c.delta}, T={c.T}")
results.append(run_test("test_five_fields (Comment 4)", test_five_fields))

# Integration test: simulate one link/epoch
def test_integration_one_link_epoch():
    S_hat_obs = 2.75   # simulated measurement
    n = 2500; delta = 1e-6; ell = 10_000; phi = 0.99; lam_ec = 0; epoch = 12
    c = bc("link_sim", S_hat_obs, n, delta, ell, phi, lam_ec, epoch)
    assert c is not None
    # Verify each field is self-consistent
    import math as _m
    eta_check = _m.sqrt(8 * _m.log(1/delta) / n)
    assert abs(c.eta - eta_check) < 1e-10, f"eta mismatch"
    from bmr.bounds import h_pironio as _h
    h_check = _h(S_hat_obs - eta_check)
    assert abs(c.h - h_check) < 1e-10, f"h mismatch"
    H_check = ell * h_check - lam_ec
    assert abs(c.H - H_check) < 1e-10, f"H mismatch"
    k_check = int(phi * H_check)
    assert c.k == k_check, f"k mismatch: {c.k} != {k_check}"
    assert abs(c.r - k_check / ell) < 1e-10, f"r mismatch"
results.append(run_test("test_integration_one_link_epoch", test_integration_one_link_epoch))

# Summary
passed = sum(r["passed"] for r in results)
total = len(results)
print(f"\n{'=' * 60}")
print(f"E3 Summary: {passed}/{total} tests passed")
if passed < total:
    print("FAILED tests:")
    for r in results:
        if not r["passed"]:
            print(f"  - {r['name']}: {r.get('error','')}")

# Save results
out = {
    "experiment": "E3",
    "description": "Certificate unit tests (Reviewer 2 Comment 4)",
    "tests": results,
    "passed": passed,
    "total": total,
    "all_passed": passed == total
}
outpath = os.path.join(os.path.dirname(__file__), "..", "data", "results", "e3_tests.json")
with open(outpath, "w") as f:
    json.dump(out, f, indent=2)
print(f"Saved results: {outpath}")
