"""Exact (manufactured) angular flux for the mms_5 verification case.

psi(x, mu) = g(mu) * {
    cos(pi x / 8)        + (mu < 0) * (b - a) * exp(-Sigma_t1 x / mu),  -4 <= x <= 0
    (1 - (x/4)^2) / 10   + (mu > 0) * (a - b) * exp(-Sigma_t2 x / mu),   0 <= x <= 4
 }
with a = 1, b = 0.1 the interface values of the two smooth parts.

The total cross sections are read from 1D_MMS_5.json so the layer thickness
follows the case file.
"""

import json
import math
from pathlib import Path

_materials = json.loads((Path(__file__).parent / "1D_MMS_5.json").read_text())[
    "materials"
]
SIGMA_T1 = _materials["mms_material_1"]["total_xs"]
SIGMA_T2 = _materials["mms_material_2"]["total_xs"]

A = 1.0
B = 0.1


def angular_function(mu: float) -> float:
    return mu**6 + 5 * mu**4 - 3 * mu**2 + 1


def exact_psi(x: float, mu: float) -> float:
    """Exact angular flux psi(x, mu) for the mms_5 MMS case, x in [-4, 4]."""
    if x <= 0:
        spatial = math.cos(math.pi * x / 8)
        if mu < 0:
            spatial += (B - A) * math.exp(-SIGMA_T1 * x / mu)
    else:
        spatial = (1 - (x / 4) ** 2) / 10
        if mu > 0:
            spatial += (A - B) * math.exp(-SIGMA_T2 * x / mu)
    return spatial * angular_function(mu)
