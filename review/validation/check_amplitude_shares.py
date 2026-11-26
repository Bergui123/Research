"""Separate the two measured amplitude errors at selected readout-A nodes.

Run from the repository root: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    .venv/bin/python review/validation/check_amplitude_shares.py --all

The calculation is floating point. It uses the existing Chebyshev recovery
and sparse exponential, and does not validate GQSP or synthesis errors. The
smallest benchmark is also checked by dense eigendecomposition.
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
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from chebyshev import p_grid, recover
from model import LAM, PHI, RESULTS, RHO, T, X0, axis_grid, grid_operator
from step5_nu_sweep import neighbours
from step7_schrodinger import g_hat_fft, profile
from step8_resources import potential_scale


def check(cfg, dense=False):
    n, nu, dx = cfg["n"], cfg["nu"], cfg["dx"]
    x1, _ = axis_grid(2**n, 2**n * dx / 2, X0[0], dx)
    x2, _ = axis_grid(2**n, 2**n * dx / 2, X0[1], dx)
    hd, potential = grid_operator(x1, x2, dx, nu)
    scale, k_bits, q5 = potential_scale(nu, x1, x2, dx)
    rounded = scale * np.floor(q5 * 2.0 ** (cfg["m"] - k_bits)) / 2**cfg["m"]
    hd_rounded = hd - diags(potential.ravel()) + diags(rounded.ravel())
    grid = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
    initial = np.exp(-np.einsum("...i,ij,...j->...", grid, PHI, grid) / (8 * nu)).ravel()
    points = [a * 2**n + b for a, b in neighbours(x1, x2, cfg["h"])]
    unrounded = expm_multiply(-T * hd, initial)[points]
    rounded_solution = expm_multiply(-T * hd_rounded, initial)[points]
    recovered, gp, window = recover(hd_rounded, initial, cfg["alpha_H"],
                                    cfg["n_p"], cfg["eta_max"], cfg["kind"],
                                    cfg["p_a"], points)
    p, p0, eta = p_grid(cfg["n_p"], cfg["eta_max"], cfg["p_a"])
    normalizer = np.sqrt(np.sum(gp**2) * np.sum(initial**2))
    window_norm = np.sqrt(np.sum(np.exp(-2 * p[window]))) / normalizer
    reference = window_norm * np.abs(unrounded)
    rounded_amp = window_norm * np.abs(rounded_solution)
    recovered_amp = np.linalg.norm(recovered, axis=0) / normalizer
    pot_error = (rounded_amp - reference) / reference
    schr_error = (recovered_amp - rounded_amp) / reference
    combined = (recovered_amp - reference) / reference
    share = 0.1 * cfg["delta_tot"]
    reference_disagreement = float(np.max(np.abs(reference / cfg["b_amp"] - 1)))
    assert reference_disagreement < 1e-9
    assert np.max(np.abs(pot_error)) <= share
    assert np.max(np.abs(schr_error)) <= share
    assert np.allclose(pot_error + schr_error, combined, rtol=1e-10, atol=1e-14)
    result = {
        "rho": RHO, "lambda": LAM, "target": cfg["target"],
        "grid": 2**n, "span": cfg["h"], "m": cfg["m"],
        "share": share, "potential_node_errors": pot_error.tolist(),
        "schrodinger_node_errors": schr_error.tolist(),
        "combined_node_errors": combined.tolist(),
        "potential_fraction_of_share": float(np.max(np.abs(pot_error)) / share),
        "schrodinger_fraction_of_share": float(np.max(np.abs(schr_error)) / share),
        "reference_disagreement": reference_disagreement,
    }
    if dense:
        assert n == 5
        eigenvalues, eigenvectors = np.linalg.eigh(hd_rounded.toarray())
        coefficients = eigenvectors.T @ initial
        kernel = np.exp(1j * np.outer(p0[window], eta)) * g_hat_fft(profile(p, cfg["kind"]), len(p))[None, :]
        shifted = np.exp(1j * np.outer(eta, eigenvalues * T))
        dense_w = kernel @ shifted @ (coefficients[:, None] * eigenvectors[points, :].T)
        dense_u = eigenvectors[points, :] @ (np.exp(-T * eigenvalues) * coefficients)
        amplitude_difference = float(np.max(np.abs(np.linalg.norm(dense_w, axis=0) /
                                                   np.linalg.norm(recovered, axis=0) - 1)))
        semigroup_difference = float(np.max(np.abs(dense_u / rounded_solution - 1)))
        assert amplitude_difference < 1e-9 and semigroup_difference < 1e-9
        result["dense_crosscheck"] = {
            "max_recovered_abs_difference": float(np.max(np.abs(dense_w - recovered))),
            "max_relative_amplitude_difference": amplitude_difference,
            "max_relative_semigroup_difference": semigroup_difference,
        }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs="*")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if args.all:
        cases = [("0.4", "1", ["0.05", "0.01", "0.001"]),
                 ("0", "1", ["0.01"]), ("0.2", "1", ["0.01"]),
                 ("0.6", "1", ["0.01"]), ("0.8", "1", ["0.01"]),
                 ("0.4", "0.25", ["0.01"]), ("0.4", "4", ["0.01"])]
        rows = []
        for rho, lam, targets in cases:
            env = os.environ.copy()
            env.update(EXEC_RHO=rho, EXEC_LAM=lam, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
            print(f"checking rho={rho}, lambda={lam}, targets={targets}", file=sys.stderr, flush=True)
            output = subprocess.check_output([sys.executable, __file__, *targets], env=env, text=True)
            rows.extend(json.loads(output))
        print(json.dumps({"scope": "Floating-point node checks; not error enclosures or circuit validation.",
                          "cases": rows}, indent=2))
    else:
        if not args.targets:
            parser.error("provide targets or --all")
        best = json.loads((RESULTS / "step8.json").read_text())["best"]
        print(json.dumps([check(best[target], dense=(RHO == 0.4 and LAM == 1.0 and target == "0.05"))
                          for target in args.targets]))


if __name__ == "__main__":
    main()
