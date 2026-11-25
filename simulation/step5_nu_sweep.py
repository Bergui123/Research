"""Steps 5-6: sweep the inventory diffusion nu, grid size N = 2^n and domain L.

Stage A fixes L per nu by a truncation test at fixed dx (choosing L by minimum total
error can lock onto accidental cancellations); Stage B refines n at that L.
Recorded: Gaussian widths, grid feedback error, normalized readout probabilities at
the four central-difference neighbours of X, norm decay, and the LCU normalisation
alpha_H = 8 nu/dx^2 + V'_max (equal to the Step 6 norm bound for H = I).
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

from concurrent.futures import ProcessPoolExecutor  # noqa: E402

import numpy as np  # noqa: E402
from scipy.sparse.linalg import eigsh  # noqa: E402

from grid_error import predicted_error  # noqa: E402
from model import A, C, LAM, PHI, T, X0, circuit_dx, reference_feedback, riccati_P, save  # noqa: E402
from step4_pde import errors, recover, solve_grid  # noqa: E402

NUS = (0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
NS = (5, 6, 7, 8)


SPANS = (1, 2, 4, 8, 16)  # readout stencil half-width h in grid steps


def neighbours(x1, x2, h=1):
	i, j = int(np.argmin(abs(x1 - X0[0]))), int(np.argmin(abs(x2 - X0[1])))
	return [(i + h, j), (i - h, j), (i, j + h), (i, j - h)]


def readout(x1, x2, dx, psi, nu, h):
	"""Span-h central difference v_a = -2 nu A_aa [ln psi(X+h dx e_a) - ln psi(X-h dx e_a)]/(2 h dx)."""
	nb = neighbours(x1, x2, h)
	if min(min(ab) for ab in nb) < 1 or max(max(ab) for ab in nb) > len(x1) - 2:
		return None
	lp = [np.log(psi[a, b]) for a, b in nb]
	v = -2 * nu * np.diag(A) * np.array([lp[0] - lp[1], lp[2] - lp[3]]) / (2 * h * dx)
	return v, [float(psi[a, b] ** 2 / np.sum(psi**2)) for a, b in nb]


def evaluate(args):
	nu, n, l_min = args
	dx = circuit_dx(2**n, l_min)
	half_width = 2**n * dx / 2
	x1, x2, dx, psi, hd, psi0 = solve_grid(2**n, half_width, nu, dx=dx)
	if readout(x1, x2, dx, psi, nu, 1) is None:
		return None  # X +- dx must be interior
	v, value = recover(x1, x2, dx, psi, nu)
	ev, eval_ = errors(v, value)
	spans = []
	for h in SPANS:
		out = readout(x1, x2, dx, psi, nu, h)
		if out is not None and np.all(np.isfinite(out[0])):
			spans.append({"h": h, "rel_err_v": errors(out[0], value)[0], "p": out[1]})
	pot_max = LAM / (4 * nu) * dx**2 * (2**n) ** 2  # V'_max: comparator scale 2^{2n} (Step 8)
	return {"nu": nu, "n": n, "L": half_width, "dx": dx, "rel_err_v": ev, "rel_err_V": eval_,
		"p_neighbours": spans[0]["p"], "spans": spans,
		"norm_ratio": float(np.sum(psi**2) / np.sum(psi0**2)),
		"alpha_kin": 8 * nu / dx**2, "alpha_pot": pot_max, "alpha_H": 8 * nu / dx**2 + pot_max}


TRUNC_TOL = 1e-5  # domain-truncation share of the feedback error (10% of the tightest 1e-4 use)


def truncation(args):
	"""Feedback at fixed dx on [-L, L) for one L (truncation test, Stage A)."""
	nu, half_width, dx = args
	n_pts = int(round(2 * half_width / dx))
	x1, x2, dx, psi, _, _ = solve_grid(n_pts, half_width, nu, dx=dx)
	if max(max(ab) for ab in neighbours(x1, x2)) > n_pts - 2:
		return nu, half_width, np.full(2, np.nan)  # X + dx not interior: domain too small
	return nu, half_width, recover(x1, x2, dx, psi, nu)[0]


def choose_domain(pool):
	"""Smallest L whose truncation error vs a 1.6x larger domain is below TRUNC_TOL."""
	jobs, refs = [], {}
	for nu in NUS:
		cands = np.round(np.arange(4.0, 4.0 + 10.0 * np.sqrt(nu), 0.25 * np.sqrt(nu) + 0.25), 3)
		refs[nu] = (cands, 1.6 * cands[-1])
		dx = 2 * refs[nu][1] / 320
		jobs += [(nu, L, dx) for L in cands] + [(nu, refs[nu][1], dx)]
	out = {}
	for nu, L, v in pool.map(truncation, jobs, chunksize=1):
		out[(nu, L)] = v
	v_ref = reference_feedback()[0]
	chosen, table = {}, []
	for nu in NUS:
		cands, l_ref = refs[nu]
		for L in cands:
			e = float(np.linalg.norm(out[(nu, L)] - out[(nu, l_ref)]) / np.linalg.norm(v_ref))
			table.append({"nu": nu, "L": float(L), "trunc_err": e})
			if nu not in chosen and e < TRUNC_TOL:  # NaN compares False
				chosen[nu] = float(L)
	return chosen, table


def main():
	with ProcessPoolExecutor() as pool:
		chosen, trunc_table = choose_domain(pool)
		jobs = [(nu, n, chosen[nu]) for nu in NUS for n in NS]
		rows = [r for r in pool.map(evaluate, jobs, chunksize=1) if r is not None]
	kap = np.sqrt(LAM * np.linalg.eigvalsh(C))
	widths = {nu: {"ell": np.sqrt(2 * nu / kap).tolist(), "E0T": 0.5 * kap.sum() * T,
		"sigma_psi2_s0": np.sqrt(2 * nu / np.linalg.eigvalsh(PHI)).tolist(),
		"sigma_psi2_sT": np.sqrt(2 * nu / np.linalg.eigvalsh(riccati_P(T))).tolist()} for nu in NUS}
	for r in rows:
		print(f"nu={r['nu']:5.1f} n={r['n']} L={r['L']:6.3f} dx={r['dx']:.4f} err_v={r['rel_err_v']:.3e} "
			f"pmin={min(r['p_neighbours']):.3e} alpha={r['alpha_H']:.1f} norm={r['norm_ratio']:.3e}")
	for r in rows:  # never above the cancellation-free leading-order bound (plus truncation)
		assert r["nu"] < 1 or r["rel_err_v"] <= predicted_error(r["nu"], r["dx"], 1)[1] + 2 * TRUNC_TOL
	checks = []
	for r in rows:
		if r["nu"] == 2.0 and r["n"] in (5, 6, 7):
			_, _, _, _, hd, _ = solve_grid(2 ** r["n"], r["L"], 2.0, dx=r["dx"])
			lam_max = float(eigsh(hd, k=1, which="LA", return_eigenvectors=False)[0])
			checks.append({"n": r["n"], "lambda_max": lam_max, "alpha_H": r["alpha_H"]})
			assert lam_max <= r["alpha_H"]
	save("step5", {"domain": chosen, "truncation": trunc_table, "grid": rows, "widths": widths,
		"alpha_checks": checks, "trunc_tol": TRUNC_TOL})


if __name__ == "__main__":
	main()
