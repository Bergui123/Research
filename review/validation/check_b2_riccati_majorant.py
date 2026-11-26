"""Sample the inequalities behind the tightened coefficient majorant.

The rational residual certificate supplies all-time Riccati bounds. These
floating-point samples only corroborate the subsequent coefficient estimates.
"""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))

import numpy as np

from b2_global_majorant import evaluate, ladder_majorant
from b2_riccati_majorant import R, coefficient_bounds
from check_b2_global_majorant import norms_by_degree
from grid_theorem import DEG, ladder, pad, polys
from model import riccati_P


def main():
    bounds = coefficient_bounds()
    min_eigenvalue, max_row_sum = np.inf, 0.0
    max_p_ratio = max_ladder_ratio = 0.0
    for t in np.linspace(0, 1, 101):
        p = riccati_P(t)
        eig = float(np.linalg.eigvalsh(p)[0])
        row_sum = float(np.max(np.sum(np.abs(p), axis=1)))
        min_eigenvalue, max_row_sum = min(min_eigenvalue, eig), max(max_row_sum, row_sum)
        assert eig >= 0.899 - 1e-12 and row_sum <= 2.4 + 1e-12
        true_polys = [np.ones((1, 1)), *polys(t)]
        for k in (0, 1, 2):
            upper = {degree: float(evaluate(poly, float(t))) for degree, poly in bounds[k].items()}
            for degree, value in norms_by_degree(true_polys[k]).items():
                cap = upper.get(degree, 0.0)
                assert value <= cap + 1e-10
                if k > 0 and cap > 1e-12:
                    max_p_ratio = max(max_p_ratio, value / cap)
            steps = 8 - 2 * k
            caps = {degree: float(evaluate(poly, float(t)))
                    for degree, poly in ladder_majorant(bounds[k], steps, R).items()}
            for axis in (0, 1):
                q = ladder(pad(true_polys[k], DEG + steps + 1), p, axis, steps)
                for degree, value in norms_by_degree(q).items():
                    cap = caps.get(degree, 0.0)
                    assert value <= cap + 1e-10
                    if cap > 1e-12:
                        max_ladder_ratio = max(max_ladder_ratio, value / cap)
    print(json.dumps({"sample_times": 101, "min_P_eigenvalue": min_eigenvalue,
                      "max_P_row_sum": max_row_sum, "max_polynomial_coefficient_ratio": max_p_ratio,
                      "max_ladder_coefficient_ratio": max_ladder_ratio,
                      "scope": "Floating-point samples corroborate, but do not prove, the coefficient bounds."}, indent=2))


if __name__ == "__main__":
    main()
