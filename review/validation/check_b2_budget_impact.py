"""Illustrate how the coarse analytic B2 majorant affects two saved grids.

S and the reference rate are evaluated in floating point. The reported term is
part of the resulting *error bound*, not a lower bound on the actual error.
"""

import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))

import numpy as np

from b2_global_majorant import calculate
from grid_theorem import expansion
from model import LAM, RHO, T, X0, reference_feedback


def main():
    assert (RHO, LAM, T) == (0.4, 1.0, 1.0)
    _, parts = calculate()
    rational = sum(parts)
    bound = (rational.numerator + rational.denominator - 1) // rational.denominator
    benchmark = json.loads((ROOT / "review/validation/alternative-boundary-selection.json").read_text())["cases"][0]
    cases = [("original", benchmark["baseline"]["0.001"]),
             ("larger_box", benchmark["larger_box"]["selection"]["best"])]
    output = []
    for label, cfg in cases:
        nu, dx, d = cfg["nu"], cfg["dx"], cfg["h"] * cfg["dx"]
        nodes = np.array([X0 + d * np.eye(2)[0], X0 - d * np.eye(2)[0],
                          X0 + d * np.eye(2)[1], X0 - d * np.eye(2)[1]])
        s = expansion(nodes / math.sqrt(nu), T, dx / math.sqrt(nu), 2)
        residual = bound * (dx / math.sqrt(nu))**6
        ratios = residual / s
        assert np.all(s > 0) and np.all(ratios < 1)
        term = nu / d * float(np.max(-np.log1p(-ratios))) / np.linalg.norm(reference_feedback()[0])
        assert term > cfg["target"] / 2
        output.append({"case": label, "B2_integer_upper": bound, "delta6": (dx / math.sqrt(nu))**6,
                       "B2_delta6": residual, "minimum_S": float(s.min()),
                       "maximum_ratio": float(ratios.max()),
                       "one_term_in_relative_rate_bound": float(term),
                       "half_target": cfg["target"] / 2})
    print(json.dumps({"scope": "Floating-point illustration of the bound; no lower bound on actual error.",
                      "cases": output}, indent=2))


if __name__ == "__main__":
    main()
