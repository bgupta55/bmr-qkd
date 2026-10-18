"""
bmr/bounds.py
=============
Hoeffding bounds for the CHSH estimator and related quantities.
From Developer Guide §1.1.

Golden values (computed):
  h_pironio(2.687) = 0.47233
  h_pironio(2.418) = 0.20614
  eta_union(100, 1.0) = 0.8157   (but the formula uses delta not 1.0; see n_min)
  n_min(2.687, 1e-6, "tight") = 235
  n_min(2.687, 1e-6, "union") = 1078
"""
import math


def eta_tight(n: int, delta: float) -> float:
    """One-sided Hoeffding on the 4n-term CHSH sum: Pr[S_hat <= S - eta] <= delta.
    Derived from the 4n-sample Hoeffding inequality: range of S ∈ [-2√2, 2√2],
    each term ±1 scaled by 1/n, so 4n terms each in range 2/n => eta = sqrt(8 ln(1/δ)/n).
    """
    return math.sqrt(8.0 * math.log(1.0 / delta) / n)


def eta_union(n: int, delta: float) -> float:
    """Previous draft: union bound over four correlators.
    Each correlator E_ab estimated from n samples, each ±1 => Hoeffding with range 2.
    Union bound over 4 correlators at level delta/4 each: eta = sqrt(2 ln(8/δ)/n) per correlator,
    then S = sum of 4 => eta_S = 4 * sqrt(2 ln(8/δ)/n) = sqrt(32 ln(8/δ)/n).
    This is the formula in the original paper (sqrt(32 ln(8/δ)/n)).
    """
    return math.sqrt(32.0 * math.log(8.0 / delta) / n)


def h_pironio(S: float) -> float:
    """Min-entropy from Pironio et al. 2010 for CHSH value S.
    h(S) = 1 - log2(1 + sqrt(2 - S²/4)), valid for S in [2, 2√2].
    h(2) = 0, h(2√2) = 1.
    """
    if S <= 2.0:
        return 0.0
    S = min(S, 2.0 * math.sqrt(2.0))
    return 1.0 - math.log2(1.0 + math.sqrt(max(0.0, 2.0 - S * S / 4.0)))


def n_min(S: float, delta: float, kind: str = "tight") -> int:
    """Smallest n (rounds per setting) such that S - eta(n, delta) > 2.
    From S - eta > 2 => eta < S - 2 => (coefficient) * log(factor/delta) / n < (S-2)^2
    => n > (coefficient) * log(factor/delta) / (S-2)^2.
    """
    gap = S - 2.0
    if gap <= 0:
        return 10 ** 18
    if kind == "tight":
        num = 8.0 * math.log(1.0 / delta)
    else:  # "union"
        num = 32.0 * math.log(8.0 / delta)
    return math.ceil(num / gap ** 2)
