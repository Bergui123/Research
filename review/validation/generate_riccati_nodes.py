"""Generate rational interpolation nodes for the benchmark Riccati ODE.

The SciPy values only choose a candidate interpolant. Correctness comes from
the independent exact-rational residual verifier, not the solver tolerance.
Run from the repository root and redirect stdout to riccati-quad-nodes.json.
"""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))

from model import LAM, RHO, T, riccati_P


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segments", type=int, default=4000)
    args = parser.parse_args()
    assert (RHO, LAM, T) == (0.4, 1.0, 1.0) and args.segments > 0
    denominator = 10**15
    nodes = []
    for i in range(2 * args.segments + 1):
        p = riccati_P(i / (2 * args.segments))
        nodes.append([round(float(p[a, b]) * denominator) for a, b in ((0, 0), (0, 1), (1, 1))])
    nodes[0] = [2 * denominator, 3 * denominator // 10, denominator]
    print(json.dumps({"segments": args.segments, "denominator": denominator,
                      "components": ["P11", "P12", "P22"], "nodes": nodes}, separators=(",", ":")))


if __name__ == "__main__":
    main()
