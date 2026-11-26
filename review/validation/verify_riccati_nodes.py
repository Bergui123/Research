"""Verify a uniform Riccati enclosure using only exact rational arithmetic.

Input is riccati-quad-nodes.json. See riccati-residual-bound.md for the
residual comparison argument and the Bernstein coefficient inequalities.
"""

from fractions import Fraction as F
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def bernstein(power):
    degree = len(power) - 1
    return [sum((power[j] * F(math.comb(k, j), math.comb(degree, j))
                 for j in range(k + 1)), F(0)) for k in range(degree + 1)]


def quadratic(left, middle, right):
    return [left, 4 * middle - 3 * left - right, 2 * right + 2 * left - 4 * middle]


def verify(data):
    n, denominator = data["segments"], data["denominator"]
    assert n > 0 and denominator > 0 and len(data["nodes"]) == 2 * n + 1
    assert data["components"] == ["P11", "P12", "P22"]
    nodes = [[F(value, denominator) for value in row] for row in data["nodes"]]
    assert nodes[0] == [F(2), F(3, 10), F(1)]
    c = [[F(1), F(2, 5)], [F(2, 5), F(1)]]
    integral = F(0)
    max_residual = F(0)
    max_row_sum = F(0)
    min_dominance = None
    for i in range(n):
        left, middle, right = nodes[2 * i:2 * i + 3]
        entries = [quadratic(*values) for values in zip(left, middle, right)]
        b = [bernstein(entry) for entry in entries]
        for k in range(3):
            a, q, d = b[0][k], b[1][k], b[2][k]
            # Each Bernstein coefficient matrix is above 0.9 I, so their
            # convex combination is too; q >= 0 also simplifies row sums.
            assert q >= 0 and a > F(9, 10) and d > F(9, 10)
            assert (a - F(9, 10)) * (d - F(9, 10)) > q * q
            assert a > q and d > q
            max_row_sum = max(max_row_sum, a + q, d + q)
            margin = min(a - q, d - q)
            min_dominance = margin if min_dominance is None else min(min_dominance, margin)
        p = [[entries[0], entries[1]], [entries[1], entries[2]]]
        residuals = []
        for a, b_index in ((0, 0), (0, 1), (1, 1)):
            coefficients = [F(0)] * 5
            coefficients[0] = n * p[a][b_index][1] - 2 * c[a][b_index]
            coefficients[1] = 2 * n * p[a][b_index][2]
            for k in (0, 1):
                for j in range(3):
                    for ell in range(3):
                        coefficients[j + ell] += p[a][k][j] * p[k][b_index][ell] / 2
            residuals.append(max(map(abs, bernstein(coefficients))))
        # The entrywise l1 norm bounds the symmetric matrix Frobenius norm.
        slab = residuals[0] + 2 * residuals[1] + residuals[2]
        integral += slab / n
        max_residual = max(max_residual, slab)
    assert integral < F(13, 10**9)
    # A row's l1 error is at most sqrt(2)||E||_F < 2 integral.
    assert max_row_sum + 2 * integral < F(12, 5)
    assert F(9, 10) - integral > F(899, 1000)
    return {"segments": n, "node_denominator": denominator,
            "residual_integral_upper_rational": str(integral),
            "residual_integral_upper": float(integral),
            "residual_upper_max": float(max_residual),
            "interpolant_max_row_sum": float(max_row_sum),
            "interpolant_min_diagonal_margin": float(min_dominance),
            "implied_P_min_eigenvalue_lower": "899/1000",
            "implied_P_row_sum_upper": "12/5",
            "scope": "Exact rational residual check, conditional on the stated Riccati comparison proof."}


def main():
    data = json.loads((ROOT / "review/validation/riccati-quad-nodes.json").read_text())
    print(json.dumps(verify(data), indent=2))


if __name__ == "__main__":
    main()
