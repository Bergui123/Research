"""Step 4 numerical exit criterion: solve the two-asset linear PDE on finite domains.

psi_s = 0.2 (psi_11 + psi_22) - 1.25 (x1^2 + 0.8 x1 x2 + x2^2) psi  (Eq. stage4-twoasset)
on [-L, L)^2 with periodic or Dirichlet boundaries, recover v* = -2 nu grad log psi
at X = (3, 1) and compare with the finite-penalty Riccati answer.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply, splu

from model import A, PHI, T, X0, axis_grid, exact_psi, grid_operator, log_det_Y, reference_feedback, save


def solve_grid(n, half_width, nu, periodic=True, dx=None):
	"""Return grid axes, psi(., T) from the exact matrix exponential, and the operator."""
	x1, dx = axis_grid(n, half_width, X0[0], dx)
	x2, _ = axis_grid(n, half_width, X0[1], dx)
	hd, _ = grid_operator(x1, x2, dx, nu, periodic)
	g = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
	psi0 = np.exp(-np.einsum("...i,ij,...j->...", g, PHI, g) / (8 * nu))
	psi = expm_multiply(-T * hd, psi0.ravel()).reshape(n, n)
	return x1, x2, dx, psi, hd, psi0


def recover(x1, x2, dx, psi, nu):
	"""Central-difference feedback and value at X (X lies on the grid)."""
	i, j = int(np.argmin(abs(x1 - X0[0]))), int(np.argmin(abs(x2 - X0[1])))
	lp = lambda a, b: np.log(psi[a, b])  # only the 5 stencil points; far tails may round below 0
	grad = np.array([lp(i + 1, j) - lp(i - 1, j), lp(i, j + 1) - lp(i, j - 1)]) / (2 * dx)
	v = -2 * nu * A @ grad
	value = -4 * nu * lp(i, j) - 2 * nu * log_det_Y(T)
	return v, value


def errors(v, value):
	v_ref, val_ref = reference_feedback()
	return float(np.linalg.norm(v - v_ref) / np.linalg.norm(v_ref)), float(abs(value - val_ref) / val_ref)


def crank_nicolson(hd, psi0, steps):
	eye = sparse.identity(hd.shape[0], format="csc")
	h = T / steps
	lu = splu((eye + h / 2 * hd).tocsc())
	rhs_op = (eye - h / 2 * hd).tocsr()
	psi = psi0.ravel().copy()
	for _ in range(steps):
		psi = lu.solve(rhs_op @ psi)
	return psi


def main(nu=0.2):
	rows = []
	for periodic in (True, False):
		for half_width in (4.0, 5.0, 6.0, 8.0):
			for n in (32, 64, 128, 256):
				x1, x2, dx, psi, hd, _ = solve_grid(n, half_width, nu, periodic)
				v, value = recover(x1, x2, dx, psi, nu)
				ev, eval_ = errors(v, value)
				g = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
				exact = exact_psi(g, nu)
				i, j = int(np.argmin(abs(x1 - X0[0]))), int(np.argmin(abs(x2 - X0[1])))
				rows.append({"bc": "periodic" if periodic else "dirichlet", "L": half_width, "N": n,
					"dx": dx, "v": v.tolist(), "rel_err_v": ev, "rel_err_V": eval_,
					"rel_err_psi_at_X": float(abs(psi[i, j] / exact[i, j] - 1)),
					"max_abs_err_psi": float(np.max(abs(psi - exact)))})
				r = rows[-1]
				print(f"{r['bc']:9s} L={half_width:4.1f} N={n:4d} dx={dx:.4f} "
					f"err_v={ev:.3e} err_V={eval_:.3e} err_psi(X)={r['rel_err_psi_at_X']:.3e}")

	# Time-grid refinement: Crank-Nicolson vs the matrix exponential on one grid.
	x1, x2, dx, psi_exp, hd, psi0 = solve_grid(64, 6.0, nu)
	v_exp, _ = recover(x1, x2, dx, psi_exp, nu)
	cn = []
	for steps in (16, 32, 64, 128, 256):
		psi_cn = crank_nicolson(hd, psi0, steps).reshape(psi_exp.shape)
		v_cn, _ = recover(x1, x2, dx, psi_cn, nu)
		cn.append({"steps": steps, "time_err_v": float(np.linalg.norm(v_cn - v_exp) / np.linalg.norm(v_exp))})
		print(f"CN steps={steps:4d} feedback time error {cn[-1]['time_err_v']:.3e}")
	assert cn[-1]["time_err_v"] < cn[0]["time_err_v"] / 100  # ~second order in time

	best = min((r for r in rows if r["bc"] == "periodic" and r["L"] >= 6), key=lambda r: r["rel_err_v"])
	save("step4", {"nu": nu, "v_ref": reference_feedback()[0].tolist(), "V_ref": reference_feedback()[1],
		"grid": rows, "crank_nicolson_L6_N64": cn, "best": best})
	return rows


if __name__ == "__main__":
	main()
