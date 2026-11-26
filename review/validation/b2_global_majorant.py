"""Exact rational coefficient majorant for the benchmark's normalized B_2.

The derivation and assumptions are in b2-global-majorant.md. This is a coarse
analytic bound for nu=1, rho=2/5, lambda=T=1 and the fixed terminal matrix.
It is not used by the production resource selection.

Run from the repository root:
  .venv/bin/python review/validation/b2_global_majorant.py
"""

from fractions import Fraction as F
import json
import math


def add(left, right):
    result = [F(0)] * max(len(left), len(right))
    for i, value in enumerate(left):
        result[i] += value
    for i, value in enumerate(right):
        result[i] += value
    return result


def scale(poly, factor):
    return [factor * value for value in poly]


def integrate(poly):
    return [F(0)] + [value / (i + 1) for i, value in enumerate(poly)]


def evaluate(poly, time):
    return sum(value * time**i for i, value in enumerate(poly))


def ladder_majorant(degrees, steps, row_quarter=F(9, 10)):
    """Bound homogeneous coefficient l1 norms after `steps` ladder operators."""
    current = {degree: poly[:] for degree, poly in degrees.items()}
    for _ in range(steps):
        following = {}
        for degree, poly in current.items():
            if degree:
                following[degree - 1] = add(following.get(degree - 1, []), scale(poly, degree))
            following[degree + 1] = add(following.get(degree + 1, []), scale(poly, row_quarter))
        current = following
    return current


def hierarchy_majorants():
    """Polynomial-in-time bounds U[k,d] for p_k's homogeneous l1 norms."""
    bounds = [{0: [F(1)]}]
    for k in (1, 2):
        source = {}
        for i in range(k):
            steps = 2 * (k - i) + 2
            # Two spatial directions and c_j = 2/(2j+2)!.
            factor = F(4, math.factorial(steps))
            for degree, poly in ladder_majorant(bounds[i], steps).items():
                source[degree] = add(source.get(degree, []), scale(poly, factor))
        output = {}
        for degree in range(4 * k, -1, -2):
            diffusion = scale(output.get(degree + 2, []), (degree + 2) * (degree + 1))
            output[degree] = integrate(add(source.get(degree, []), diffusion))
        bounds.append(output)
    return bounds


def gaussian_weight(degree):
    """Rational upper bound for sup r^degree exp(-(7/108) r^2)."""
    assert degree >= 0 and degree % 2 == 0
    if degree == 0:
        return F(1)
    # e > sum_{j=0}^6 1/j! = 1957/720. Every degree here is even.
    return F(38880 * degree, 13699) ** (degree // 2)


def calculate():
    bounds = hierarchy_majorants()
    parts = []
    for k in (0, 1, 2):
        steps = 2 * (2 - k) + 4
        terms = ladder_majorant(bounds[k], steps)
        integrand = []
        for degree, poly in terms.items():
            assert degree % 2 == 0
            integrand = add(integrand, scale(poly, gaussian_weight(degree)))
        # Two directions and coefficient 2/steps! in the theorem.
        parts.append(evaluate(integrate(scale(integrand, F(4, math.factorial(steps)))), F(1)))
    return bounds, parts


def main():
    bounds, parts = calculate()
    total = sum(parts, F(0))
    print(json.dumps({
        "scope": "Normalized benchmark: nu=1, rho=2/5, lambda=T=1; analytic coefficient majorant.",
        "P_lower": "14/27", "P_operator_upper": "12/5", "P_row_sum_upper": "18/5",
        "gaussian_exponent_lower": "7/108", "ladder_row_quarter_upper": "9/10",
        "coefficient_bounds": [{str(d): [str(x) for x in p] for d, p in block.items()} for block in bounds],
        "B2_parts": [str(p) for p in parts], "B2_upper_rational": str(total),
        "B2_upper_decimal": float(total), "B2_upper_integer": (total.numerator + total.denominator - 1) // total.denominator,
    }, indent=2))


if __name__ == "__main__":
    main()
