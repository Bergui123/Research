"""Emulate three exploratory configurations re-costed with projected box bounds.

The actual-grid emulation uses a fixed seed and excludes GQSP truncation and
gate synthesis. A passing run does not validate uncompiled circuit costs or total algorithm error.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))

import numpy as np

from model import LAM, RHO
import step10_validation


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("target", nargs="?")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--selection", type=Path, default=ROOT / "review/validation/projected-recost.json")
    args = parser.parse_args()
    if args.all:
        output = []
        for target in ("0.05", "0.01", "0.001"):
            env = os.environ.copy()
            env.update(EXEC_RHO="0.4", EXEC_LAM="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
            print(f"emulating target={target}", file=sys.stderr, flush=True)
            output.append(json.loads(subprocess.check_output([sys.executable, __file__, target, "--selection", str(args.selection)], env=env, text=True)))
        print(json.dumps(output, indent=2))
        return
    if args.target not in ("0.05", "0.01", "0.001"):
        parser.error("provide a benchmark target or --all")
    assert (RHO, LAM) == (0.4, 1.0)
    cases = json.loads(args.selection.read_text())["cases"]
    cfg = next(case["best"] for case in cases if case["target"] == args.target)
    assert cfg is not None and cfg["target"] == float(args.target)
    step10_validation.RNG = np.random.default_rng(20261122)
    result = step10_validation.emulate(cfg)
    assert result["amp_rel_err"] <= result["amp_budget"]
    assert result["lower_bound"] >= 0.95
    print(json.dumps({"grid": 2**cfg["n"], "nu": cfg["nu"], "span": cfg["h"],
                      "profile": cfg["kind"], "eta_max": cfg["eta_max"],
                      "validation": result}))


if __name__ == "__main__":
    main()
