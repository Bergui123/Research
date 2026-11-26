"""Numerical estimates of grid-error bounds for readout-A candidates.

Chain: periodic box --(truncation bound u)--> infinite lattice --(theorem, K = 2)--> S --> exact.
S = psi (1 + delta^2 p_1 + delta^4 p_2); |psi_box - S| <= delta^6 B_2 + u at each stencil node.
u solves the box problem with zero data and ghost value g. The theoretical result requires g
to bound the wrapped and infinite-lattice values at all times. The sampled floating-point
maxima used here are not validated enclosures, so these results are estimates, not certificates;
B_2 and the ghost value are charged at SAFETY times their sampled values.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import json  # noqa: E402
import math  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from grid_error import feedback_shift, predicted_error
from grid_theorem import expansion, remainder_constant
from model import PHI, RESULTS, T, X0, axis_grid, grid_operator, reference_feedback, save
from step5_nu_sweep import neighbours

TIMES = 33
SAFETY = 2.0


def truncation(cfg, x1, x2, b2):
	"""Numerically estimate the box bound and its ghost value."""
	nu, dx, n = cfg["nu"], cfg["dx"], 2 ** cfg["n"]
	delta = dx / math.sqrt(nu)
	per, pot = grid_operator(x1, x2, dx, nu)
	box, _ = grid_operator(x1, x2, dx, nu, periodic=False)
	g = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
	u0 = np.exp(-np.einsum("...i,ij,...j->...", g, PHI, g) / (8 * nu)).ravel()
	edge = np.zeros((n, n), bool)
	edge[[0, -1], :] = edge[:, [0, -1]] = True
	wrapped = expm_multiply(-per, u0, start=0, stop=T, num=TIMES)[:, edge.ravel()].max()
	# Infinite-lattice values outside the box: S + delta^6 B_2 on a surrounding band.
	pad = n // 4
	band = lambda x: np.concatenate([x[0] - dx * np.arange(pad, 0, -1), x, x[-1] + dx * np.arange(1, pad + 1)])
	e = np.stack(np.meshgrid(band(x1), band(x2), indexing="ij"), axis=-1) / math.sqrt(nu)
	outside = np.ones(e.shape[:2], bool)
	outside[pad:-pad, pad:-pad] = False
	lattice = max(expansion(e[outside], s, delta, 2).max() for s in np.linspace(0, T, TIMES)) + delta**6 * b2
	ghost = SAFETY * max(wrapped, lattice)
	count = np.zeros((n, n))
	count[0] += 1; count[-1] += 1; count[:, 0] += 1; count[:, -1] += 1
	src = sparse.csr_matrix((nu / dx**2 * ghost * count).reshape(-1, 1))
	aug = sparse.bmat([[-box, src], [None, sparse.csr_matrix((1, 1))]]).tocsr()
	start = np.zeros(n * n + 1)
	start[-1] = 1
	u = expm_multiply(T * aug, start)[:-1].reshape(n, n)
	return u, ghost


def estimate(args):
	row, b2 = args
	nu, dx, n = row["nu"], row["dx"], row["n"]
	half = 2**n * dx / 2
	x1, _ = axis_grid(2**n, half, X0[0], dx)
	x2, _ = axis_grid(2**n, half, X0[1], dx)
	delta = dx / math.sqrt(nu)
	u, ghost = truncation(row, x1, x2, b2)
	v = reference_feedback()[0]
	norm = np.linalg.norm(v)
	out = []
	for sp in row["spans"]:
		h = sp["h"]
		d = h * dx
		u_pts = np.array([u[i, j] for i, j in neighbours(x1, x2, h)])
		nodes = np.array([X0 + d * np.eye(2)[0], X0 - d * np.eye(2)[0], X0 + d * np.eye(2)[1], X0 - d * np.eye(2)[1]])
		s_val = expansion(nodes / math.sqrt(nu), T, delta, 2)
		if np.any(delta**6 * b2 + u_pts >= s_val):
			out.append({"nu": nu, "n": n, "h": h, "numerical_error_estimate": None})
			continue
		ls = np.log(s_val)
		v_s = -2 * nu * np.array([ls[0] - ls[1], ls[2] - ls[3]]) / (2 * d)
		slack = lambda e: -np.log1p(-e / s_val)
		pair = lambda e: np.linalg.norm(nu / d * np.array([slack(e)[0] + slack(e)[1], slack(e)[2] + slack(e)[3]])) / norm
		out.append({"nu": nu, "n": n, "h": h, "delta": delta, "measured": sp["rel_err_v"], "ghost": float(ghost),
			"cancellation_free": predicted_error(nu, dx, h)[1], "u_max": float(u_pts.max()),
			"higher_order": float(np.linalg.norm(v_s - v - feedback_shift(nu, dx, h)) / norm),
			"remainder_lattice": float(pair(delta**6 * b2)), "remainder": float(pair(delta**6 * b2 + u_pts)),
			"numerical_error_estimate": float(np.linalg.norm(v_s - v) / norm + pair(delta**6 * b2 + u_pts))})
	return out


if __name__ == "__main__":
	rows = [r for r in json.loads((RESULTS / "step5.json").read_text())["grid"] if r["nu"] >= 1]
	b2 = SAFETY * remainder_constant(2)
	with ProcessPoolExecutor(max_workers=8) as pool:
		out = [c for part in pool.map(estimate, [(r, b2) for r in rows]) for c in part]
	ok = [c for c in out if c["numerical_error_estimate"] is not None]
	# This regression check compares sampled numerical values; it does not certify an enclosure.
	assert all(c["measured"] <= c["numerical_error_estimate"] * (1 + 1e-6) for c in ok)
	print(RESULTS.name, f"B_2={b2:.3e}: numerical estimates available for {len(ok)} of {len(out)} (grid, span) pairs; "
		f"median estimate/measured {np.median([c['numerical_error_estimate'] / c['measured'] for c in ok]):.2f}")
	save("grid_estimates", {"B2": b2, "rows": out})
