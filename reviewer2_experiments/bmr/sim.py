"""
bmr/sim.py
==========
Link-level simulator used by every network experiment (E8-style sweep, dynamic E9, Q15, Q7 test).

The ONLY place where a certificate is created is `LinkSim.publish_certificate`, which draws
the audit statistic S_hat for the link's TRUE CHSH value (binomial/Werner standard deviation
sqrt((4 - S^2/4)/n)) and then calls `bmr.certificate.build_certificate` with the simulator's
parameters.  No experiment re-implements the certificate rule (pending_v7 Q7).
"""
import math
import numpy as np
from .certificate import build_certificate

S_HONEST = 2.687      # Werner link with visibility ~0.95
S_COMPROMISED = 2.418


class LinkSim:
    def __init__(self, links, n_audit, delta=1e-6, ell=10_000, phi=0.99, lam_ec=0, seed=0):
        self.n, self.delta, self.ell, self.phi, self.lam_ec = n_audit, delta, ell, phi, lam_ec
        self.rng = np.random.default_rng(seed)
        self.true_S = {self._key(l): S_HONEST for l in links}
        self.last_S_hat = {}
        self.published = {}

    @staticmethod
    def _key(link):
        u, v = link
        return (u, v) if u <= v else (v, u)

    def set_true_S(self, link, S):
        self.true_S[self._key(link)] = float(S)

    def observe(self, link):
        """One audit: S_hat ~ N(S_true, sqrt((4 - S^2/4)/n)) (binomial CHSH estimator)."""
        S = self.true_S[self._key(link)]
        sd = math.sqrt(max(0.0, 4.0 - S * S / 4.0) / self.n)
        s_hat = float(self.rng.normal(S, sd))
        self.last_S_hat[self._key(link)] = s_hat
        return s_hat

    def publish_certificate(self, link, epoch):
        """Audit the link and publish its certificate (Cert or None = DEGRADED)."""
        k = self._key(link)
        s_hat = self.observe(k)
        cert = build_certificate(k, s_hat, self.n, self.delta, self.ell, self.phi, self.lam_ec, epoch)
        self.published[k] = cert
        return cert
