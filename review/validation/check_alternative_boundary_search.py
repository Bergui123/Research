"""Reselect readout-A resources with the analytic boundary bound, in temporary files.

Run from the repository root with:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
    review/validation/check_alternative_boundary_search.py --all

The saved B2, expansion coefficients and box evolution remain floating-point
estimates. This experiment does not change simulation/results or establish a
validated total error bound. It also tests one larger 256-square box at 0.1%.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))
sys.path.insert(0, str(ROOT / "review/validation"))

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

import certify_grid
from check_boundary_bound import boundary_bound
from model import LAM, RESULTS, RHO, T, grid_operator
from step5_nu_sweep import evaluate
import step8_resources as resources


def analytic_truncation(cfg, x1, x2, _b2):
    """Apply the rational boundary bound; box evolution is still floating point."""
    n, nu, dx = 2**cfg["n"], cfg["nu"], cfg["dx"]
    # Conversion to binary float must not reduce the rational boundary bound.
    ghost = math.nextafter(float(boundary_bound(cfg)[0]), math.inf)
    box, _ = grid_operator(x1, x2, dx, nu, periodic=False)
    count = np.zeros((n, n))
    count[0] += 1; count[-1] += 1; count[:, 0] += 1; count[:, -1] += 1
    source = sparse.csr_matrix((nu / dx**2 * ghost * count).reshape(-1, 1))
    augmented = sparse.bmat([[-box, source], [None, sparse.csr_matrix((1, 1))]]).tocsr()
    initial = np.zeros(n * n + 1)
    initial[-1] = 1
    u = expm_multiply(T * augmented, initial)[:-1].reshape(n, n)
    return u, ghost


def select(target):
    """Use the production candidate and cost functions on temporary inputs."""
    best = None
    candidates = accepted = 0
    for cfg in resources.candidates(target):
        candidates += 1
        result = resources.run_cost(cfg, target)
        if result is None:
            continue
        accepted += 1
        if best is None or result["T_total"] < best["T_total"]:
            best = result
    return {"candidate_profiles": candidates, "costed_candidates": accepted, "best": best}


def case(targets, extension):
    original_results = RESULTS
    original_step5 = json.loads((original_results / "step5.json").read_text())
    original_estimates = json.loads((original_results / "grid_estimates.json").read_text())
    baseline = json.loads((original_results / "step8.json").read_text())["best"]
    rows = [row for row in original_step5["grid"] if row["nu"] >= 1]
    certify_grid.truncation = analytic_truncation
    with ProcessPoolExecutor(max_workers=3) as pool:
        estimates = [item for part in pool.map(certify_grid.estimate,
                     [(row, original_estimates["B2"]) for row in rows]) for item in part]
    available = [r for r in estimates if r["numerical_error_estimate"] is not None]
    assert all(r["measured"] <= r["numerical_error_estimate"] * (1 + 1e-6) for r in available)

    with tempfile.TemporaryDirectory(prefix="research-boundary-search-") as directory:
        temp = Path(directory)
        for name in ("step5.json", "step7.json"):
            shutil.copyfile(original_results / name, temp / name)
        (temp / "grid_estimates.json").write_text(json.dumps({"B2": original_estimates["B2"], "rows": estimates}))
        resources.RESULTS = temp
        selections = {target: select(float(target)) for target in targets}
        result = {"rho": RHO, "lambda": LAM, "B2_sampled_charge": original_estimates["B2"],
                  "available_estimates": len(available), "grid_span_pairs": len(estimates),
                  "baseline": {target: baseline[target] for target in targets},
                  "alternative": selections}
        if extension:
            assert (RHO, LAM) == (0.4, 1.0)
            assert "0.001" in targets
            # Test one specified larger box. Replacing this one (nu, n) row in
            # the temporary input preserves the other search rows for this
            # experiment; it does not assert a global optimum over box sizes.
            custom = evaluate((8.0, 8, 512 / 30))
            assert custom is not None and abs(custom["dx"] - 4 / 30) < 1e-14
            custom_estimates = certify_grid.estimate((custom, original_estimates["B2"]))
            modified = dict(original_step5)
            modified["grid"] = [r for r in rows if (r["nu"], r["n"]) != (8.0, 8)] + [custom]
            modified_estimates = [r for r in estimates if (r["nu"], r["n"]) != (8.0, 8)] + custom_estimates
            (temp / "step5.json").write_text(json.dumps(modified))
            (temp / "grid_estimates.json").write_text(json.dumps({"B2": original_estimates["B2"], "rows": modified_estimates}))
            result["larger_box"] = {"row": custom, "estimates": custom_estimates,
                                    "selection": select(0.001)}
        return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs="*")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--extension", action="store_true")
    args = parser.parse_args()
    if args.all:
        cases = [("0.4", "1", ["0.05", "0.01", "0.001"], True),
                 ("0", "1", ["0.01"], False), ("0.2", "1", ["0.01"], False),
                 ("0.6", "1", ["0.01"], False), ("0.8", "1", ["0.01"], False),
                 ("0.4", "0.25", ["0.01"], False), ("0.4", "4", ["0.01"], False)]
        results = []
        for rho, lam, targets, extension in cases:
            env = os.environ.copy()
            env.update(EXEC_RHO=rho, EXEC_LAM=lam, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
            command = [sys.executable, __file__, *targets]
            if extension:
                command.append("--extension")
            print(f"checking rho={rho}, lambda={lam}, targets={targets}", file=sys.stderr, flush=True)
            results.append(json.loads(subprocess.check_output(command, env=env, text=True)))
        print(json.dumps({"scope": "Alternative selection; the boundary bound is analytic but other grid quantities remain numerical.",
                          "cases": results}, indent=2))
    else:
        if not args.targets:
            parser.error("provide targets or --all")
        print(json.dumps(case(args.targets, args.extension)))


if __name__ == "__main__":
    main()
