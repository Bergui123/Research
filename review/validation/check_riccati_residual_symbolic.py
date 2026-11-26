"""Check representative residual polynomials against a separate SymPy construction."""

from fractions import Fraction as F
import json
from pathlib import Path

import sympy as sp

from verify_riccati_nodes import bernstein, quadratic

ROOT = Path(__file__).resolve().parents[2]
z = sp.symbols("z")


def main():
    data = json.loads((ROOT / "review/validation/riccati-quad-nodes.json").read_text())
    n, den = data["segments"], data["denominator"]
    points = sorted({0, 1, n // 4, n // 2, 3 * n // 4, n - 2, n - 1})
    covariance = sp.Matrix([[1, sp.Rational(2, 5)], [sp.Rational(2, 5), 1]])
    l0, lm, l1 = 2 * (z - sp.Rational(1, 2)) * (z - 1), 4 * z * (1 - z), 2 * z * (z - sp.Rational(1, 2))
    for segment in points:
        rows = data["nodes"][2 * segment:2 * segment + 3]
        interp = [sum(sp.Rational(row[i], den) * factor for row, factor in zip(rows, (l0, lm, l1)))
                  for i in range(3)]
        p = sp.Matrix([[interp[0], interp[1]], [interp[1], interp[2]]])
        residual = n * p.diff(z) - 2 * covariance + p * p / 2
        coefficients = [quadratic(*(F(row[i], den) for row in rows)) for i in range(3)]
        matrix = [[coefficients[0], coefficients[1]], [coefficients[1], coefficients[2]]]
        for i, j in ((0, 0), (0, 1), (1, 1)):
            code = [F(0)] * 5
            code[0] = n * matrix[i][j][1] - 2 * F(int(covariance[i, j].p), int(covariance[i, j].q))
            code[1] = 2 * n * matrix[i][j][2]
            for k in (0, 1):
                for a in range(3):
                    for b in range(3):
                        code[a + b] += matrix[i][k][a] * matrix[k][j][b] / 2
            symbolic = sp.Poly(residual[i, j], z)
            for degree in range(5):
                term = symbolic.nth(degree)
                assert F(int(term.p), int(term.q)) == code[degree]
            upper = max(map(abs, bernstein(code)))
            for numerator in range(8):
                x = F(numerator, 7)
                assert abs(sum(coef * x**degree for degree, coef in enumerate(code))) <= upper
    print(json.dumps({"segments_checked": points, "matrices_per_segment": 3,
                      "scope": "Exact symbolic comparison on selected segments; the rational verifier covers all segments."}, indent=2))


if __name__ == "__main__":
    main()
