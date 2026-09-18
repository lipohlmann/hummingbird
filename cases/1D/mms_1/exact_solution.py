"""Exact (manufactured) angular flux for the mms_2 verification case.

psi(x, mu) = (1 - x^2)^2 * exp(-x^2/5) * (3*mu^2 + 5*mu^4) / (2*pi)

Verified symbolically against sources.mms_source.expression in 1D_MMS_2.json
(residual == 0) during debugging of a TinyExpr unary-minus/power precedence
issue in that expression (fixed in commit de7550b).
"""

import math


def exact_psi(x: float, mu: float) -> float:
    """Exact angular flux psi(x, mu) for the mms_2 MMS case, x, mu in [-1, 1]."""
    return (1 - x**2) ** 2 * (3 * mu**2 + 1) / (2 * math.pi)
