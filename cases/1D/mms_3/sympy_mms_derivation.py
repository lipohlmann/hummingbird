"""Symbolic MMS source for the mms_3 (two-region, jump-material) verification case.

Manufactured solution, x in [-1, 1], split at x = 0:

    psi(x, mu) = f_r(x) * g(mu)
    g(mu)      = mu^6 + 5 mu^4 - 3 mu^2 + 1
    f_1(x)     = (1 + x)^2 cos(x)^2                              x in [-1, 0]
    f_2(x)     = (2x + 1) (1 - x^2)^(3/2) exp(-3 x^2)            x in [0, 1]

Source that makes psi an exact solution of the 1D transport equation:

    Q = mu dpsi/dx + Sigma_t psi - Sigma_s / (4 pi) phi
    phi(x) = 2 pi * integral_{-1}^{1} psi dmu

Each region gets its own Q. Printed per region: LaTeX with arbitrary cross
sections, factored as g(mu) [mu f' + Sigma_t f] - (Sigma_s/2) (int g dmu) f and
written with the paper's \\ohatx macro for mu, and a TinyExpr++ string with the
XS_VALUES cross sections substituted, ready for the "expression" field of a
parsed_function source (uses pow() rather than ^ to sidestep TinyExpr's
unary-minus/power precedence).
"""

import math
import random

from sympy import (
    Rational,
    ccode,
    cos,
    diff,
    exp,
    expand,
    factor,
    integrate,
    latex,
    pi,
    sin,
    sqrt,
    symbols,
)

x, mu = symbols("x mu", real=True)

# Cross sections substituted into the TinyExpr++ string (LaTeX stays symbolic).
Sigma_t1, Sigma_s1, Sigma_t2, Sigma_s2 = symbols(
    "Sigma_t1 Sigma_s1 Sigma_t2 Sigma_s2", positive=True
)
XS_VALUES = {
    Sigma_t1: Rational(2),
    Sigma_s1: Rational(1),
    Sigma_t2: Rational(8, 5),
    Sigma_s2: Rational(11, 10),
}

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
        "f": (1 + x) ** 2 * cos(x) ** 2,
        "Sigma_t": Sigma_t1,
        "Sigma_s": Sigma_s1,
        "x_range": (-1.0, 0.0),
    },
    2: {
        "f": (2 * x + 1) * (1 - x**2) ** Rational(3, 2) * exp(-3 * x**2),
        "Sigma_t": Sigma_t2,
        "Sigma_s": Sigma_s2,
        "x_range": (0.0, 1.0),
    },
}


def mms_source(psi, Sigma_t, Sigma_s):
    phi = 2 * pi * integrate(psi, (mu, -1, 1))
    return mu * diff(psi, x) + Sigma_t * psi - Sigma_s / (4 * pi) * phi


def factored_source(f, Sigma_t, Sigma_s):
    """Same Q as mms_source for psi = f(x) g(mu), as (streaming + collision, scattering)
    with Q = first - second and g(mu) kept as the symbol G."""
    # factor() leaves sqrt(-(x - 1)*(x + 1)); expanding the argument gives sqrt(1 - x**2).
    fp = factor(diff(f, x)).replace(
        lambda e: e.is_Pow and e.exp == Rational(1, 2), lambda e: sqrt(expand(e.base))
    )
    return G * (mu * fp + Sigma_t * f), Sigma_s * integrate(g, (mu, -1, 1)) / 2 * f


def check_factored(terms, Q, x_range, n_points=5):
    """The factored terms (with G -> g) and Q must agree at random points."""
    Q_factored = terms[0] - terms[1]
    free = sorted(Q.free_symbols - {x, mu}, key=str)
    xs_check = {s: 1.5 + i for i, s in enumerate(free)}
    for _ in range(n_points):
        point = {x: random.uniform(*x_range), mu: random.uniform(-1, 1), **xs_check}
        a = float(Q_factored.subs(G, g).subs(point))
        b = float(Q.subs(point))
        assert math.isclose(
            a, b, rel_tol=1e-12, abs_tol=1e-12
        ), f"factored form disagrees with Q at {point}: {a} vs {b}"


def check_expression(Q, expression, x_range, n_points=5):
    """Evaluate the printed string and Q at random points; they must agree."""
    free = sorted(Q.free_symbols - {x, mu}, key=str)
    xs_check = {s: 1.5 + i for i, s in enumerate(free)}
    namespace = {
        "pow": math.pow,
        "sin": math.sin,
        "cos": math.cos,
        "exp": math.exp,
        "sqrt": math.sqrt,
        "pi": math.pi,
    }
    for _ in range(n_points):
        xv, muv = random.uniform(*x_range), random.uniform(-1, 1)
        local = {"x": xv, "mu": muv, **{str(s): v for s, v in xs_check.items()}}
        from_string = eval(expression, namespace, local)
        from_sympy = float(Q.subs({x: xv, mu: muv, **xs_check}))
        assert math.isclose(from_string, from_sympy, rel_tol=1e-12, abs_tol=1e-12), (
            f"printed expression disagrees with sympy at x={xv}, mu={muv}: "
            f"{from_string} vs {from_sympy}"
        )


for region, r in regions.items():
    Q_symbolic = mms_source(r["f"] * g, r["Sigma_t"], r["Sigma_s"])
    terms = factored_source(r["f"], r["Sigma_t"], r["Sigma_s"])
    check_factored(terms, Q_symbolic, r["x_range"])
    Q = Q_symbolic.subs(XS_VALUES)
    expression = ccode(Q).replace("M_PI", "pi")
    check_expression(Q, expression, r["x_range"])
    print(f"Region {region} (x in [{r['x_range'][0]:g}, {r['x_range'][1]:g}])")
    first, second = (latex(t, symbol_names=LATEX_NAMES) for t in terms)
    print("LaTeX (arbitrary cross sections):")
    print(
        "\\begin{equation}\n"
        f"  \\label{{eq:mms3-source-{region}}}\n"
        "  \\begin{split}\n"
        f"    Q_{region}\\left(x,\\ohatx\\right) &= {first} \\\\\n"
        f"    &\\quad - {second}\n"
        "  \\end{split}\n"
        "\\end{equation}\n"
    )
    print(f"TinyExpr++ expression:\n{expression}\n")
