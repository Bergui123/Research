"""Numerical grid-error sensitivity checks for the saved 1% benchmark.

Run from the repository root with:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_grid_sensitivity.py

This is a convergence diagnostic, not an interval-arithmetic certificate.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'simulation'))

import numpy as np
from scipy.sparse.linalg import expm_multiply
from grid_theorem import remainder_constant
from model import PHI, RESULTS, T, X0, axis_grid, grid_operator

cfg = json.loads((RESULTS / 'step8.json').read_text())['best']['0.01']
saved = json.loads((RESULTS / 'grid_theorem.json').read_text())['B'][2]
n = 2 ** cfg['n']
dx, nu = cfg['dx'], cfg['nu']
x1, _ = axis_grid(n, n * dx / 2, X0[0], dx)
x2, _ = axis_grid(n, n * dx / 2, X0[1], dx)
hd, _ = grid_operator(x1, x2, dx, nu)
coords = np.stack(np.meshgrid(x1, x2, indexing='ij'), axis=-1)
u0 = np.exp(-np.einsum('...i,ij,...j->...', coords, PHI, coords) / (8 * nu)).ravel()
edge = np.zeros((n, n), dtype=bool)
edge[[0, -1], :] = edge[:, [0, -1]] = True
values = expm_multiply(-hd, u0, start=0, stop=T, num=257)[:, edge.ravel()]
coarse = float(values[::8].max())
fine = float(values.max())
coarse_space = remainder_constant(2, n_s=33, step=0.08)
fine_time = remainder_constant(2, n_s=65, step=0.08)
print(json.dumps({
    'benchmark_grid': n, 'saved_B2_n33_dx004': saved,
    'B2_n33_dx008': coarse_space, 'B2_n65_dx008': fine_time,
    'wrapped_edge_max_33_times': coarse, 'wrapped_edge_max_257_times': fine,
    'edge_max_relative_increase': (fine / coarse - 1) if coarse else 0.0,
}, indent=2))
