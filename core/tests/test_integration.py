"""
Q7: integration tests tying the SIMULATOR to build_certificate.
Run:  python3 tests/test_integration.py     (or pytest tests/)
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
from bmr.certificate import build_certificate
from bmr.sim import LinkSim


def test_integration_one_link_epoch():
    """One simulated link-epoch: published tuple == build_certificate(...) on the SAME inputs."""
    sim = LinkSim([("a", "b"), ("b", "c")], n_audit=2500, delta=1e-6, ell=10_000, phi=0.99, lam_ec=0, seed=7)
    published = sim.publish_certificate(("a", "b"), epoch=3)
    S_hat = sim.last_S_hat[("a", "b")]                       # the statistic the simulator drew
    expected = build_certificate(("a", "b"), S_hat, 2500, 1e-6, 10_000, 0.99, 0, 3)
    assert published == expected, (published, expected)      # dataclass equality: every field
    assert published is not None and published.n == 2500 and published.T == 3 and published.delta == 1e-6


def test_integration_many_epochs_and_degraded():
    sim = LinkSim([(0, 1)], n_audit=500, seed=11)
    sim.set_true_S((0, 1), 2.418)                              # compromised: DEGRADED at n=500
    nones = 0
    for ep in range(200):
        pub = sim.publish_certificate((0, 1), ep)
        exp = build_certificate((0, 1), sim.last_S_hat[(0, 1)], 500, 1e-6, 10_000, 0.99, 0, ep)
        assert pub == exp
        nones += pub is None
    # S_true=2.418, n=500: S_hat ~ N(2.418, 0.071); pass needs S_hat > 2+eta = 2.470 -> P(pass) ~ 23 %
    assert 0.6 < nones / 200 < 0.95, nones


def test_pass_rule_is_eta_subtracted():
    """Certificate exists iff S_hat - eta_tight(n, delta) > 2 (no midpoint shortcut)."""
    from bmr.bounds import eta_tight
    rng = np.random.default_rng(3)
    for n in (100, 336, 379, 500, 2500, 10_000):
        for _ in range(300):
            s = float(rng.uniform(2.0, 2.83))
            c = build_certificate("x", s, n, 1e-6, 10_000, 0.99, 0, 1)
            assert (c is not None) == (s - eta_tight(n, 1e-6) > 2.0)


def test_n100_never_certifies():
    # eta_tight(100,1e-6)=1.0513 > 2*sqrt2-2 => no S in [2,2sqrt2] passes
    for s in np.linspace(2.0, 2 * math.sqrt(2), 50):
        assert build_certificate("x", float(s), 100, 1e-6, 10_000, 0.99, 0, 1) is None


def test_experiments_use_library_path():
    """Static check: every network experiment reaches build_certificate (directly or via LinkSim)."""
    root = os.path.dirname(os.path.dirname(__file__))
    for f in ("q4_collusion", "q5_evidence_source", "e7b_v2_corrected", "q9_solver_validation",
              "q15_delta_regime", "q6_common", "q6_dynamic", "q6_sweep"):
        src = open(os.path.join(root, "experiments", f + ".py")).read()
        assert ("build_certificate" in src) or ("LinkSim" in src) or ("q6_common" in src), f


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("PASS", name)
