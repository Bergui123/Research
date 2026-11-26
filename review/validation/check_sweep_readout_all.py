"""Run the existing readout-A check at all six current sweep optima.

Run from the repository root and redirect stdout to a JSON file:
  .venv/bin/python review/validation/check_sweep_readout_all.py
"""

import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
CASES = (("0", "1"), ("0.2", "1"), ("0.6", "1"),
         ("0.8", "1"), ("0.4", "0.25"), ("0.4", "4"))


def main():
    rows = []
    for rho, lam in CASES:
        env = os.environ.copy()
        env.update(EXEC_RHO=rho, EXEC_LAM=lam, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        print(f"checking rho={rho}, lambda={lam}", file=sys.stderr, flush=True)
        output = subprocess.check_output(
            [sys.executable, str(ROOT / "review/validation/check_sweep_readout.py")],
            cwd=ROOT, env=env, text=True,
        )
        row = json.loads(output)
        best_file = ROOT / f"simulation/results/rho{float(rho):g}_lam{float(lam):g}/step8.json"
        cfg = json.loads(best_file.read_text())["best"]["0.01"]
        assert (row["grid"], row["nu"], row["span"], row["profile"], row["eta_max"], row["n_p"]) == (
            2**cfg["n"], cfg["nu"], cfg["h"], cfg["kind"], cfg["eta_max"], cfg["n_p"]
        )
        assert row["passes_amplitude_budget"] and row["passes_rate_target_with_exact_amplitudes"]
        rows.append(row)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
