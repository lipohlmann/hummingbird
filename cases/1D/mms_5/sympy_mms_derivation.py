"""Symbolic MMS source for the mms_5 (interface boundary-layer) verification case.

Thin region 1 on [-4, 0] against a thick absorber, region 2, on [0, 4].
Manufactured solution (not separable in x and mu):

    x < 0:  psi = g(mu) [ f_1(x) + (mu < 0) (b - a) exp(-Sigma_t1 x / mu) ]
    x > 0:  psi = g(mu) [ f_2(x) + (mu > 0) (a - b) exp(-Sigma_t2 x / mu) ]

    g(mu)  = mu^6 + 5 mu^4 - 3 mu^2 + 1
    f_1(x) = cos(pi x / 8)              a = f_1(0) = 1,   f_1(-4) = 0
    f_2(x) = (1 - (x/4)^2) / 10         b = f_2(0) = 0.1, f_2(4)  = 0

psi is continuous in x for every mu != 0 (so the current is continuous), but
each ordinate carries a boundary layer of thickness |mu| / Sigma_t downstream
of the interface, and psi is discontinuous in mu at x = 0.

Source that makes psi an exact solution of the S_N equations:

    Q = mu dpsi/dx + Sigma_t psi - (Sigma_s / 2) sum_n w_n psi(x, mu_n)

The layer terms are homogeneous solutions of the streaming + collision
operator (asserted below), so they enter Q only through the scattering sum.
Their angular integral has no closed form (exponential integrals), so the
scattering term uses the solver's own N_POLAR-point Gauss-Legendre rule. The
printed source is therefore exact only for n_polar = N_POLAR; regenerate it if
n_polar changes.

Printed per region: LaTeX with arbitrary cross sections, with the scattering
term left as a sum over the outgoing half of the quadrature set and written
with the paper's \\ohatx macro for mu, and a TinyExpr++ string with the
XS_VALUES cross sections substituted, ready for the "expression" field of a parsed_function source (uses
pow() rather than ^ to sidestep TinyExpr's unary-minus/power precedence).
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

DPI = 300
FIGSIZE = (6.0, 4.0)
POSITION_TOLERANCE = 1e-9
X_LABEL = "x [cm]"
SCALAR_FLUX_LABEL = "Scalar Flux [n/cm²/s]"
ANGULAR_FLUX_LABEL = "Angular Flux [n/cm²/s/str]"


def configure_style():
    """Serifed, publication-style figure defaults."""
    sns.set_theme(style="ticks")
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [
                "Nimbus Roman",
                "Liberation Serif",
                "STIXGeneral",
                "DejaVu Serif",
            ],
            "mathtext.fontset": "stix",
            "figure.dpi": DPI,
            "savefig.dpi": DPI,
            "axes.titlesize": 14,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            "legend.fontsize": 10.5,
            "legend.title_fontsize": 11,
            "legend.frameon": True,
            "legend.framealpha": 0.9,
            "legend.edgecolor": "0.8",
            "axes.linewidth": 1.0,
            "lines.linewidth": 1.8,
            "lines.markersize": 5,
            "lines.markeredgewidth": 0,
        }
    )


import math
import random

from numpy.polynomial.legendre import leggauss
from sympy import (
    Float,
    Rational,
    ccode,
    cos,
    diff,
    exp,
    integrate,
    latex,
    pi,
    simplify,
    symbols,
)

x, mu = symbols("x mu", real=True)

Sigma_t1, Sigma_s1, Sigma_t2, Sigma_s2 = symbols(
    "Sigma_t1 Sigma_s1 Sigma_t2 Sigma_s2", positive=True
)
XS_VALUES = {
    Sigma_t1: Rational(1),
    Sigma_s1: Rational(1, 2),
    Sigma_t2: Rational(50),
    Sigma_s2: Rational(10),
}

# Must match angular_treatment.n_polar in the case JSON.
N_POLAR = 12

# Stands in for g(mu) in the factored (LaTeX) form of the source.
G = symbols("G")

LATEX_NAMES = {
    mu: r"\ohatx",
    G: r"g\left(\ohatx\right)",
    Sigma_t1: r"\Sigma_{t,1}",
    Sigma_s1: r"\Sigma_{s,1}",
    Sigma_t2: r"\Sigma_{t,2}",
    Sigma_s2: r"\Sigma_{s,2}",
}

g = mu**6 + 5 * mu**4 - 3 * mu**2 + 1

regions = {
    1: {
        "f": cos(pi * x / 8),
        "Sigma_t": Sigma_t1,
        "Sigma_s": Sigma_s1,
        "x_range": (-4.0, 0.0),
    },
    2: {
        "f": (1 - (x / 4) ** 2) / 10,
        "Sigma_t": Sigma_t2,
        "Sigma_s": Sigma_s2,
        "x_range": (0.0, 4.0),
    },
}


def streaming_collision(psi, Sigma_t):
    return mu * diff(psi, x) + Sigma_t * psi


def check_expression(Q, expression, x_range, n_points=5):
    """Evaluate the printed string and Q at random points; they must agree."""
    namespace = {
        "pow": math.pow,
        "sin": math.sin,
        "cos": math.cos,
        "exp": math.exp,
        "pi": math.pi,
    }
    for _ in range(n_points):
        xv, muv = random.uniform(*x_range), random.uniform(-1, 1)
        from_string = eval(expression, namespace, {"x": xv, "mu": muv})
        from_sympy = float(Q.subs({x: xv, mu: muv}))
        assert math.isclose(from_string, from_sympy, rel_tol=1e-12, abs_tol=1e-12), (
            f"printed expression disagrees with sympy at x={xv}, mu={muv}: "
            f"{from_string} vs {from_sympy}"
        )


a = regions[1]["f"].subs(x, 0)
b = regions[2]["f"].subs(x, 0)
assert regions[1]["f"].subs(x, -4) == 0, "f_1 must vanish at the west boundary"
assert regions[2]["f"].subs(x, 4) == 0, "f_2 must vanish at the east boundary"

ordinates, weights = leggauss(N_POLAR)

for region, r in regions.items():
    # Layer carried by the ordinates leaving the interface into this region:
    # mu < 0 in region 1, mu > 0 in region 2. It must add nothing to streaming +
    # collision, and must make psi continuous at x = 0.
    jump = (b - a) if region == 1 else (a - b)
    layer = jump * exp(-r["Sigma_t"] * x / mu) * g
    assert simplify(streaming_collision(layer, r["Sigma_t"])) == 0
    assert (
        simplify((r["f"] * g + layer).subs(x, 0) - (a + b - r["f"].subs(x, 0)) * g) == 0
    )

    # S_N scattering source: (Sigma_s / 2) * sum_n w_n psi(x, mu_n).
    scalar_sum = 0
    for mu_n, w_n in zip(ordinates, weights):
        psi_n = r["f"] * g.subs(mu, Float(mu_n))
        if (mu_n < 0) == (region == 1):
            psi_n += layer.subs(mu, Float(mu_n))
        scalar_sum += Float(w_n) * psi_n

    Q = streaming_collision(r["f"] * g, r["Sigma_t"]) - r["Sigma_s"] / 2 * scalar_sum
    Q = Q.subs(XS_VALUES)
    expression = ccode(Q).replace("M_PI", "pi")
    check_expression(Q, expression, r["x_range"])
    print(f"Region {region} (x in [{r['x_range'][0]:g}, {r['x_range'][1]:g}])")

    # The N_POLAR-point rule integrates the polynomial g exactly, so the
    # separable part of the scattering sum is written with the exact integral.
    assert math.isclose(
        sum(weights * (ordinates**6 + 5 * ordinates**4 - 3 * ordinates**2 + 1)),
        float(integrate(g, (mu, -1, 1))),
    )
    first = latex(
        G * (mu * diff(r["f"], x) + r["Sigma_t"] * r["f"]), symbol_names=LATEX_NAMES
    )
    separable = latex(integrate(g, (mu, -1, 1)) * r["f"])
    half = r"\ohatxn < 0" if region == 1 else r"\ohatxn > 0"
    sigma_t, sigma_s = LATEX_NAMES[r["Sigma_t"]], LATEX_NAMES[r["Sigma_s"]]
    print("LaTeX (arbitrary cross sections):")
    print(
        "\\begin{equation}\n"
        f"  \\label{{eq:mms5-source-{region}}}\n"
        "  \\begin{split}\n"
        f"    Q_{region}\\left(x,\\ohatx\\right) &= {first} \\\\\n"
        f"    &\\quad - \\frac{{{sigma_s}}}{{2}} \\left[ {separable} "
        f"+ \\left({latex(jump)}\\right) \\sum_{{{half}}} w_n\\, "
        f"g\\left(\\ohatxn\\right) e^{{-{sigma_t} x / \\ohatxn}} \\right]\n"
        "  \\end{split}\n"
        "\\end{equation}\n"
    )
    print(f"TinyExpr++ expression:\n{expression}\n")


def angular(mu):
    return mu**6 + 5 * mu**4 / 3 * mu**2 + 1


def left_angular_flux(x, mu):
    return mu * (
        np.cos(np.pi * x / 8) - 9 / 10 * np.heaviside(-mu, 0.5) * np.exp(-1 * x / mu)
    )


def right_angular_flux(x, mu):
    return mu * (
        (1 - (x / 4) ** 2) / 10 + 9 / 10 * np.heaviside(mu, 0.5) * np.exp(-50 * x / mu)
    )


def angular_flux(xs, mu):
    fluxes = np.zeros_like(xs)
    for i in range(len(xs)):
        x = xs[i]
        if x <= 0:
            fluxes[i] = left_angular_flux(x, mu)
        else:
            fluxes[i] = right_angular_flux(x, mu)
    return fluxes


# Plot the function
fig, ax = plt.subplots(figsize=FIGSIZE)
palette = sns.color_palette("colorblind")
mus, _ = leggauss(12)
xs = np.linspace(-4, 4, 1000)

for i in range(0, len(mus), 2):
    mu = mus[i]
    ys = angular_flux(xs, mu)
    label = rf"$\mu = {mu:.3f}$"
    sns.lineplot(
        x=xs,
        y=ys,
        ax=ax,
        color=palette[int(i / 2) % len(palette)],
        label=label,
    )
ax.set_xlabel(X_LABEL)
ax.set_ylabel(ANGULAR_FLUX_LABEL)
ax.set_title("Angular Flux by Ordinate")
ax.legend(title="Ordinate", loc="upper right", ncols=2, fontsize="x-small")
fig.savefig("mms_5_analytical.png", dpi=DPI, bbox_inches="tight")
plt.close(fig)
