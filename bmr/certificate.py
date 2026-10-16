"""
bmr/certificate.py
==================
Certificate construction (Developer Guide §1.2, Blueprint §3.4).
Addresses Reviewer 2, Comment 4: "does not specify how to get c from existing parameters".

The five certificate fields are:
  link  - the link identifier
  S_hat - the observed CHSH value
  n     - rounds per setting
  eta   - Hoeffding radius (tight bound)
  h     - h(S_hat - eta): min-entropy rate at the lower-confidence bound
  H     - certified entropy in bits: ell * h - lam_ec
  k     - output key length: floor(phi * H)
  r     - key rate: k / ell
  log2_eps - leftover-hash epsilon (log domain): -1 - (H - k) / 2
  log2_gamma - gamma = 2^k * eps (log domain): k + log2_eps
  delta - audit failure probability
  T     - epoch of the audit

The "c" parameter (Comment 4) is the full certificate dataclass, constructed here.
"""
from dataclasses import dataclass
import math
from .bounds import eta_tight, h_pironio


@dataclass(frozen=True)
class Cert:
    link: tuple
    S_hat: float
    n: int
    eta: float
    h: float
    H: float          # certified entropy in bits (ell*h - lam_ec)
    k: int            # output key length
    r: float          # rate = k / ell
    log2_eps: float   # eps = 0.5 * 2**(-(H-k)/2); stored in log domain (can underflow)
    log2_gamma: float # gamma = 2^k * eps (log domain)
    delta: float
    T: int            # epoch of the audit

    @property
    def eps(self):
        return 0.0 if self.log2_eps < -1000 else 2.0 ** self.log2_eps

    @property
    def gamma(self):
        return math.inf if self.log2_gamma > 60 else 2.0 ** self.log2_gamma

    @property
    def w(self):
        """Single budget weight = eps + delta (blueprint G4)."""
        return self.eps + self.delta


def build_certificate(link, S_hat, n, delta, ell, phi, lam_ec, epoch):
    """Build a certificate or return None (DEGRADED) if S_hat - eta <= 2.

    Parameters
    ----------
    link    : identifier (tuple or string) for the link
    S_hat   : observed CHSH value this epoch
    n       : rounds per setting (total 4n Bell rounds)
    delta   : audit failure probability
    ell     : raw key rounds per epoch (block length)
    phi     : compression fraction (k = floor(phi * H))
    lam_ec  : error-correction leakage (bits)
    epoch   : integer epoch index (timestamp)

    The certificate c is defined (answering Reviewer 2 Comment 4) as:
      eta  = eta_tight(n, delta)       [Hoeffding radius]
      S_lb = S_hat - eta               [lower confidence bound on true S]
      h    = h_pironio(S_lb)           [Pironio min-entropy at S_lb]
      H    = ell * h - lam_ec          [certified extractable entropy]
      k    = floor(phi * H)            [output key length]
      r    = k / ell                   [key rate per round]
      eps  = 2^{-1 - (H-k)/2}         [leftover-hash security; Tomamichel et al. 2011]
      gamma = 2^k * eps                [composability parameter]
    """
    eta = eta_tight(n, delta)
    S_lb = S_hat - eta
    if S_lb <= 2.0:
        return None   # DEGRADED: link does not certify quantum advantage
    h = h_pironio(S_lb)
    H = ell * h - lam_ec
    if H <= 0:
        return None
    k = math.floor(phi * H)
    # Leftover-hash lemma (Tomamichel et al. 2011, Eq. 1):
    # eps = 2^{-1} * 2^{-(H-k)/2}  =>  log2_eps = -1 - (H-k)/2
    log2_eps = -1.0 - (H - k) / 2.0
    # gamma = 2^k * eps  =>  log2_gamma = k + log2_eps
    log2_gamma = k + log2_eps
    return Cert(
        link=link, S_hat=S_hat, n=n, eta=eta, h=h,
        H=H, k=k, r=k / ell,
        log2_eps=log2_eps, log2_gamma=log2_gamma,
        delta=delta, T=epoch
    )
