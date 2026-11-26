"""Bound boundary values for selected readout-A grids and assess the cost of that bound.

Run from the repository root:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
    review/validation/check_boundary_bound.py --all

The boundary bound uses exact rational inequalities for exponentials and cosh.
The subsequent grid-error calculation still uses floating-point B2, polynomial
coefficients and a sparse matrix exponential; it is not a certified error bound.
"""

import argparse
from fractions import Fraction as F
import json
import math
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))

import numpy as np
from scipy import sparse
from scipy.optimize import minimize_scalar
from scipy.sparse.linalg import expm_multiply

from grid_theorem import expansion
from model import LAM, PHI, RESULTS, RHO, T, X0, axis_grid, grid_operator, reference_feedback
from step5_nu_sweep import neighbours


def ceil_fraction(value, digits=15):
    """Round a nonnegative rational upward to a fixed decimal denominator."""
    assert value >= 0
    scale = 10**digits
    return F((value.numerator * scale + value.denominator - 1) // value.denominator, scale)


def exp_negative_upper(x):
    """For x >= 0, exp(-x) <= 1 / sum_{j=0}^{32} x**j/j!."""
    assert x >= 0
    term = total = F(1)
    for j in range(1, 33):
        term *= x / j
        total += term
    return ceil_fraction(1 / total)


def cosh_minus_one_upper(x):
    """Upper bound using six terms and a geometric bound on the positive tail."""
    assert x >= 0
    term = F(1)
    total = F(0)
    for j in range(1, 7):
        term *= x * x / ((2 * j - 1) * (2 * j))
        total += term
    following = term * x * x / (13 * 14)
    ratio = x * x / (15 * 16)
    assert ratio < 1
    return ceil_fraction(total + following / (1 - ratio), digits=18)


def bound_for_a(nu, dx, period, radius, a):
    """Return an exact rational upper bound for the periodised exponential."""
    # Phi = [[2, 3/10], [3/10, 1]], and max diag(Phi^-1) = 200/191.
    # Thus u0(x) <= exp(-m x_i**2) for either coordinate.
    m = F(191, 1600) / nu
    t = a * radius - a * a / (4 * m) - 2 * nu / dx**2 * cosh_minus_one_upper(a * dx)
    if t <= 0:
        return F(1)
    e_t = exp_negative_upper(t)
    e_far = exp_negative_upper(a * (period - 2 * radius))
    e_period = exp_negative_upper(a * period)
    if e_period >= 1:
        return F(1)
    return min(F(1), ceil_fraction(e_t * (1 + e_far) / (1 - e_period), digits=12))


def periodic_distance(x, period):
    return min(x % period, (-x) % period)


def boundary_bound(cfg):
    n, nu = 2**cfg["n"], F(str(cfg["nu"]))
    # The costed family uses dx = 4/k, nu > 0 and T = 1.
    dx = F(4, round(4 / cfg["dx"]))
    assert abs(float(dx) - cfg["dx"]) < 1e-12 and T == 1 and nu > 0
    assert abs(RHO) <= 1 and LAM >= 0
    assert np.array_equal(PHI, np.array([[2.0, 0.3], [0.3, 1.0]]))
    period = n * dx
    assert abs(cfg["L"] - float(period) / 2) < 1e-12
    radii = []
    for anchor in (F(3), F(1)):
        j = round((anchor + period / 2) / dx)
        lo, hi = anchor - j * dx, anchor + (n - 1 - j) * dx
        axis, _ = axis_grid(n, cfg["L"], float(anchor), cfg["dx"])
        assert abs(axis[0] - float(lo)) < 1e-12 and abs(axis[-1] - float(hi)) < 1e-12
        radii.extend((periodic_distance(lo, period), periodic_distance(hi, period),
                      abs(lo - dx), abs(hi + dx)))
    radius = min(radii)
    assert 0 <= radius <= period / 2
    # This floating-point search only chooses a rational a; the subsequent
    # inequality is valid for any positive rational choice.
    nf, df, pf, rf = map(float, (nu, dx, period, radius))
    def approximate(a):
        t = a * rf - 400 * nf * a * a / 191 - 2 * nf / df**2 * (math.cosh(a * df) - 1)
        return math.exp(-t) * (1 + math.exp(-a * (pf - 2 * rf))) / (-math.expm1(-a * pf))
    opt = minimize_scalar(approximate, bounds=(0.005, 1.0), method="bounded")
    a = F(max(1, round(opt.x * 1000)), 1000)
    return bound_for_a(nu, dx, period, radius, a), a, radius


def check(cfg):
    n, dx, nu, h = 2**cfg["n"], cfg["dx"], cfg["nu"], cfg["h"]
    estimate_data = json.loads((RESULTS / "grid_estimates.json").read_text())
    old = next(r for r in estimate_data["rows"] if (r["nu"], r["n"], r["h"]) ==
               (nu, cfg["n"], h))
    ghost, a, radius = boundary_bound(cfg)
    x1, _ = axis_grid(n, n * dx / 2, X0[0], dx)
    x2, _ = axis_grid(n, n * dx / 2, X0[1], dx)
    box, _ = grid_operator(x1, x2, dx, nu, periodic=False)
    count = np.zeros((n, n))
    count[0] += 1; count[-1] += 1; count[:, 0] += 1; count[:, -1] += 1
    source = sparse.csr_matrix((nu / dx**2 * count).reshape(-1, 1))
    augmented = sparse.bmat([[-box, source], [None, sparse.csr_matrix((1, 1))]]).tocsr()
    initial = np.zeros(n * n + 1)
    initial[-1] = 1
    unit_u = expm_multiply(T * augmented, initial)[:-1].reshape(n, n)
    u = float(ghost) * np.array([unit_u[i, j] for i, j in neighbours(x1, x2, h)])
    d = h * dx
    nodes = np.array([X0 + d * np.eye(2)[0], X0 - d * np.eye(2)[0],
                      X0 + d * np.eye(2)[1], X0 - d * np.eye(2)[1]])
    s = expansion(nodes / math.sqrt(nu), T, dx / math.sqrt(nu), 2)
    e = (dx / math.sqrt(nu))**6 * estimate_data["B2"] + u
    ratios = e / s
    assert np.all(s > 0) and np.all(ratios < 1)
    v = reference_feedback()[0]
    ls = np.log(s)
    v_s = -nu / d * np.array([ls[0] - ls[1], ls[2] - ls[3]])
    slack = -np.log1p(-ratios)
    error = (np.linalg.norm(v_s - v) + np.linalg.norm(nu / d *
             np.array([slack[0] + slack[1], slack[2] + slack[3]]))) / np.linalg.norm(v)
    charged = max(error, old["cancellation_free"])
    return {"rho": RHO, "lambda": LAM, "target": cfg["target"], "grid": n,
            "span": h, "a": str(a), "radius": str(radius),
            "rational_boundary_upper": str(ghost), "boundary_upper": float(ghost),
            "sampled_boundary_charge": old["ghost"],
            "old_grid_charge": cfg["eps_grid"], "alternative_grid_estimate": float(error),
            "alternative_grid_charge": float(charged),
            "alternative_charge_over_target": float(charged / cfg["target"]),
            "max_E_over_S": float(np.max(ratios))}


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
            rows.extend(json.loads(subprocess.check_output([sys.executable, __file__, *targets], env=env, text=True)))
        print(json.dumps({"scope": "The boundary upper bound uses rational arithmetic; the grid-error evaluation remains floating point.",
                          "cases": rows}, indent=2))
    else:
        if not args.targets:
            parser.error("provide targets or --all")
        best = json.loads((RESULTS / "step8.json").read_text())["best"]
        print(json.dumps([check(best[target]) for target in args.targets]))


if __name__ == "__main__":
    main()
