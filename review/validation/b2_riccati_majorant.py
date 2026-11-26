"""Rational B2 majorant using the verified benchmark Riccati interpolant.

The residual verifier supplies P >= 899 I/1000 and row sums <= 12/5; the
coefficient hierarchy is then bounded as in b2-global-majorant.md. This
does not enclose the remaining numerical quantities in the grid-error budget.
"""

from fractions import Fraction as F
import json
import math
from pathlib import Path

from b2_global_majorant import add, evaluate, integrate, ladder_majorant, scale
from verify_riccati_nodes import verify

ROOT = Path(__file__).resolve().parents[2]
R = F(3, 5)
GAUSSIAN = F(899, 8000)


def coefficient_bounds():
    bounds = [{0: [F(1)]}]
    for k in (1, 2):
        source = {}
        for i in range(k):
            steps = 2 * (k - i) + 2
            for degree, poly in ladder_majorant(bounds[i], steps, R).items():
                source[degree] = add(source.get(degree, []), scale(poly, F(4, math.factorial(steps))))
        block = {}
        for degree in range(4 * k, -1, -2):
            diffusion = scale(block.get(degree + 2, []), (degree + 2) * (degree + 1))
            block[degree] = integrate(add(source.get(degree, []), diffusion))
        bounds.append(block)
    return bounds


def weight(degree):
    assert degree >= 0 and degree % 2 == 0
    if degree == 0:
        return F(1)
    # e > 1957/720, and the Gaussian exponent is at least 899/8000.
    return (F(360 * degree, 1957) / GAUSSIAN) ** (degree // 2)


def main():
    nodes = json.loads((ROOT / "review/validation/riccati-quad-nodes.json").read_text())
    certificate = verify(nodes)
    assert certificate["implied_P_min_eigenvalue_lower"] == "899/1000"
    assert certificate["implied_P_row_sum_upper"] == "12/5"
    bounds = coefficient_bounds()
    parts = []
    for k in (0, 1, 2):
        steps = 8 - 2 * k
        integrand = []
        for degree, poly in ladder_majorant(bounds[k], steps, R).items():
            integrand = add(integrand, scale(poly, weight(degree)))
        parts.append(evaluate(integrate(scale(integrand, F(4, math.factorial(steps)))), F(1)))
    total = sum(parts, F(0))
    print(json.dumps({"scope": "Normalized benchmark only; rational residual and coefficient bounds, subject to proof review.",
                      "P_error_upper_rational": certificate["residual_integral_upper_rational"],
                      "P_min_eigenvalue_lower": "899/1000", "P_row_sum_upper": "12/5",
                      "gaussian_exponent_lower": str(GAUSSIAN), "ladder_row_quarter_upper": str(R),
                      "B2_parts": [str(p) for p in parts], "B2_upper_rational": str(total),
                      "B2_upper_decimal": float(total),
                      "B2_upper_integer": (total.numerator + total.denominator - 1) // total.denominator}, indent=2))


if __name__ == "__main__":
    main()
