"""Exact (manufactured) angular flux for the mms_3 verification case.

psi(x, mu) = {
    (sin(13x)+sin(52)) * g(mu),                -4 <= x <= 0
    sin(52)*(1-(x/4)^30)*(2-exp(-13x)) * g(mu), 0 <= x <= 4
 }
"""

import math


def angular_function(mu: float) -> float:
    return mu**6 + 5 * mu**4 - 3 * mu**2 + 1


def exact_psi(x: float, mu: float) -> float:
    """Exact angular flux psi(x, mu) for the mms_1 MMS case, x, mu in [-1, 1]."""
    angular_contrib = angular_function(mu)
    if x <= 0:
        return (math.sin(13 * x) + math.sin(52)) * angular_contrib
    else:
        return (
            math.sin(52)
            * (1 - (x / 4) ** 30)
            * (2 - math.exp(-13 * x))
            * angular_contrib
        )
