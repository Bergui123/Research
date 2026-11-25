"""Periodic against Dirichlet boundaries: relative change of the recovered rate at N = 64, h = 1."""

import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model import PHI, RESULTS, T, X0, axis_grid, grid_operator, reference_feedback, save
from step5_nu_sweep import readout

out = []
for r in json.loads((RESULTS / "step5.json").read_text())["grid"]:
	if r["n"] != 6:
		continue
	n, dx, nu = 64, r["dx"], r["nu"]
	x1, _ = axis_grid(n, n * dx / 2, X0[0], dx)
	x2, _ = axis_grid(n, n * dx / 2, X0[1], dx)
	g = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
	u0 = np.exp(-np.einsum("...i,ij,...j->...", g, PHI, g) / (8 * nu)).ravel()
	v = [readout(x1, x2, dx, expm_multiply(-T * grid_operator(x1, x2, dx, nu, periodic=p)[0], u0).reshape(n, n), nu, 1)[0]
		for p in (True, False)]
	out.append({"nu": nu, "rel_diff": float(np.linalg.norm(v[0] - v[1]) / np.linalg.norm(reference_feedback()[0]))})
	print(out[-1])
save("boundary_check", out)
