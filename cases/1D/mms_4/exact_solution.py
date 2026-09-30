"""Exact (manufactured) angular flux for the mms_3 verification case.

psi(x, mu) = {
    (1+x)^2 * cos^2(x) * g(mu), -1 <= x <= 0
    (2x+1) * (1-x^2)^3 * g(mu),  0 <= x <= 1
 }
"""

import math


def angular_function(mu: float) -> float:
    return mu**6 + 5 * mu**4 - 3 * mu**2 + 1


def exact_psi(x: float, mu: float) -> float:
    """Exact angular flux psi(x, mu) for the mms_1 MMS case, x, mu in [-1, 1]."""
    angular_contrib = angular_function(mu)
    if x <= 0:
        return (1 + x) ** 2 * math.cos(x) ** 2 * angular_contrib
    else:
        return (1 - x**2) ** 3 * (2 * x + 1) * math.exp(-3 * x**2) * angular_contrib
