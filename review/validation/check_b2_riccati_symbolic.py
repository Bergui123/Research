"""Recompute the tightened B2 recurrence using SymPy exact polynomials."""

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
t = sp.symbols("t", nonnegative=True)
R = sp.Rational(3, 5)
zero = sp.Integer(0)


def transform(degrees, steps):
    current = dict(degrees)
    for _ in range(steps):
        following = {}
        for degree, value in current.items():
            if degree:
                following[degree - 1] = following.get(degree - 1, zero) + degree * value
            following[degree + 1] = following.get(degree + 1, zero) + R * value
        current = following
    return current


def main():
    polys = [{0: sp.Integer(1)}]
    for k in (1, 2):
        source = {}
        for i in range(k):
            steps = 2 * (k - i) + 2
            for degree, value in transform(polys[i], steps).items():
                source[degree] = source.get(degree, zero) + 4 * value / sp.factorial(steps)
        output = {}
        for degree in range(4 * k, -1, -2):
            rhs = source.get(degree, zero) + (degree + 2) * (degree + 1) * output.get(degree + 2, zero)
            output[degree] = sp.integrate(rhs, (t, 0, t)).expand()
        polys.append(output)
    parts = []
    for k in (0, 1, 2):
        steps = 8 - 2 * k
        terms = transform(polys[k], steps)
        total = zero
        for degree, value in terms.items():
            assert degree % 2 == 0
            weight = 1 if degree == 0 else (sp.Rational(360 * degree, 1957) / sp.Rational(899, 8000)) ** (degree // 2)
            total += value * weight
        parts.append(sp.integrate(4 * total / sp.factorial(steps), (t, 0, 1)))
    expected = json.loads((ROOT / "review/validation/b2-riccati-majorant.json").read_text())
    result = Fraction(int(sum(parts).p), int(sum(parts).q))
    assert result == Fraction(expected["B2_upper_rational"])
    print(json.dumps({"symbolic_recurrence_matches_rational_implementation": True,
                      "B2_upper_rational": str(result)}, indent=2))


if __name__ == "__main__":
    main()
