"""Emulate the changed benchmark candidates from the alternative search.

Run after check_alternative_boundary_search.py --all. This invokes the existing
2,000-trial emulation with its fixed seed; it does not emulate GQSP truncation
or gate synthesis, or turn the grid-error estimate into an enclosure.
"""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))

import numpy as np

from model import LAM, RHO
import step10_validation


def main():
    assert (RHO, LAM) == (0.4, 1.0)
    result = json.loads((ROOT / "review/validation/alternative-boundary-selection.json").read_text())
    benchmark = next(c for c in result["cases"] if (c["rho"], c["lambda"]) == (0.4, 1.0))
    configs = [benchmark["alternative"]["0.01"]["best"],
               benchmark["larger_box"]["selection"]["best"]]
    output = []
    for cfg, target in zip(configs, (0.01, 0.001)):
        assert cfg is not None and cfg["target"] == target
        step10_validation.RNG = np.random.default_rng(20261122)
        validation = step10_validation.emulate(cfg)
        assert validation["amp_rel_err"] <= validation["amp_budget"]
        assert validation["lower_bound"] >= 0.95
        output.append({"grid": 2**cfg["n"], "nu": cfg["nu"], "span": cfg["h"],
                       "profile": cfg["kind"], "eta_max": cfg["eta_max"],
                       "validation": validation})
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
