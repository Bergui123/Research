"""Numerical checks of the benchmark coefficient-majorant inequalities.

These checks are corroboration, not a proof between sample times. The analytic
argument and the exact rational calculation are separate.
"""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))

import numpy as np

from b2_global_majorant import calculate, evaluate, ladder_majorant
from grid_theorem import DEG, ladder, pad, polys
from model import riccati_P


def norms_by_degree(poly):
    result = {}
    for i in range(poly.shape[0]):
        for j in range(poly.shape[1]):
            result[i + j] = result.get(i + j, 0.0) + abs(float(poly[i, j]))
    return result


def main():
    bounds, parts = calculate()
    max_p_ratio = max_ladder_ratio = 0.0
    min_eigenvalue, max_row_sum, min_diagonal_margin = np.inf, 0.0, np.inf
    for t in np.linspace(0, 1, 101):
        p = riccati_P(t)
        min_eigenvalue = min(min_eigenvalue, float(np.linalg.eigvalsh(p)[0]))
        max_row_sum = max(max_row_sum, float(np.max(np.sum(np.abs(p), axis=1))))
        min_diagonal_margin = min(min_diagonal_margin, float(p[1, 1] - abs(p[0, 1])))
        assert np.linalg.eigvalsh(p)[0] >= 14 / 27 - 1e-12
        assert np.max(np.sum(np.abs(p), axis=1)) <= 18 / 5 + 1e-12
        assert p[0, 0] >= p[1, 1] - 1e-12 and p[1, 1] >= abs(p[0, 1]) - 1e-12
        true_polys = [np.ones((1, 1)), *polys(t)]
        for k in (0, 1, 2):
            majorants = {degree: float(evaluate(poly, float(t))) for degree, poly in bounds[k].items()}
            for degree, value in norms_by_degree(true_polys[k]).items():
                upper = majorants.get(degree, 0.0)
                assert value <= upper + 1e-10
                if k > 0 and upper > 1e-12:
                    max_p_ratio = max(max_p_ratio, value / upper)
            steps = 2 * (2 - k) + 4
            transformed = ladder_majorant(bounds[k], steps)
            transformed = {degree: float(evaluate(poly, float(t))) for degree, poly in transformed.items()}
            for axis in (0, 1):
                q = ladder(pad(true_polys[k], DEG + steps + 1), p, axis, steps)
                for degree, value in norms_by_degree(q).items():
                    upper = transformed.get(degree, 0.0)
                    assert value <= upper + 1e-10
                    if upper > 1e-12:
                        max_ladder_ratio = max(max_ladder_ratio, value / upper)
    print(json.dumps({
        "sample_times": 101, "min_P_eigenvalue": min_eigenvalue,
        "max_P_row_sum": max_row_sum, "min_diagonal_dominance_margin": min_diagonal_margin,
        "max_polynomial_coefficient_ratio": max_p_ratio,
        "max_ladder_coefficient_ratio": max_ladder_ratio,
        "exact_majorant": str(sum(parts)),
        "scope": "Floating-point sample checks; they do not establish the inequalities between times.",
    }, indent=2))


if __name__ == "__main__":
    main()
