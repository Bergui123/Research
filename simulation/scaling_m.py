"""Scaling of both readouts with the number of assets M (continuum Gaussian formulas).

H = I, C = (1 - rho) I + rho 11^T, Phi = 1.5 I, X = (3, 1, 3, 1, ...), lam = T = 1, eps = 1%.
Grid spacing is held at d xi = 0.25 in xi = x / sqrt(nu); only the number of coherent runs is
modelled (no T-count). Readout A minimises over nu and the readout span d; readout B counts
runs for the uniform feedback-matrix metric. An extrapolation, not a simulation for M > 2.
"""

import numpy as np
from scipy.linalg import expm

from model import PHI, SHIFT_MARGIN, T, riccati_P, block_solution, save

EPS, DXI = 0.01, 0.25
RNG = np.random.default_rng(7)


def riccati(m, rho=0.4, lam=1.0, phi=None, c=None):
	c = (1 - rho) * np.eye(m) + rho * np.ones((m, m)) if c is None else c
	phi = 1.5 * np.eye(m) if phi is None else phi
	gen = np.block([[np.zeros((m, m)), np.eye(m) / 2], [2 * lam * c, np.zeros((m, m))]])
	yz = expm(T * gen) @ np.vstack([np.eye(m), phi])
	y, z = yz[:m], yz[m:]
	p = z @ np.linalg.inv(y)
	e0 = 0.5 * np.sum(np.sqrt(lam * np.linalg.eigvalsh(c)))
	ratio = np.sqrt(np.linalg.det(phi) / np.linalg.det((p + p.T) / 2)) / np.linalg.det(y)
	return (p + p.T) / 2, ratio, e0


def runs_a(p, x, p_succ, m):
	v = np.linalg.norm(p @ x) / 2
	best = np.inf
	for nu in np.geomspace(0.5, 64, 43):
		log_cell = m * np.log(DXI * np.sqrt(nu)) - m / 2 * np.log(4 * np.pi * nu) + 0.5 * np.linalg.slogdet(p)[1]
		for d in np.geomspace(0.25, 4, 13) * np.sqrt(nu):
			delta = 0.7 * (EPS / 2) * v * d / (2 * np.sqrt(m) * nu)
			total = 0.0
			for a in range(m):
				for sgn in (1, -1):
					y = x + sgn * d * np.eye(m)[a]
					b = np.sqrt(p_succ * np.exp(log_cell - y @ p @ y / (4 * nu)))
					total += 2 * np.pi / (delta * b)
			best = min(best, total)
	return best


def runs_b(p, p_succ, m, trials=4000):
	w, u = np.linalg.eigh(p)
	g = RNG.standard_normal((trials, m, m))
	e = (g + np.transpose(g, (0, 2, 1))) / np.sqrt(2)  # var 2 on the diagonal, 1 off it
	op = np.linalg.norm((u * np.sqrt(w) @ u.T) @ e @ (u / np.sqrt(w) @ u.T), 2, axis=(1, 2))
	return (np.quantile(op, 0.95) / (EPS / 2)) ** 2 / p_succ


if __name__ == "__main__":
	# The closed-form norm ratio must match the grid value of the benchmark (0.134).
	p2 = riccati_P(T)
	ratio2 = np.sqrt(np.linalg.det(PHI) / np.linalg.det(p2)) / np.linalg.det(block_solution(T)[0])
	assert abs(ratio2 - 0.1335) < 5e-4, ratio2
	rows = []
	for m in (1, 2, 3, 4, 5, 6, 8, 10):
		p, ratio, e0 = riccati(m)
		x = np.array([3.0, 1.0] * m)[:m]
		plain, shifted = 0.5 * ratio, 0.5 * ratio * np.exp(2 * SHIFT_MARGIN * e0 * T)
		rows.append({"M": m, "norm_ratio": ratio, "P_succ_plain": plain, "P_succ_shifted": shifted, "P_succ_exact": 0.5 * ratio * np.exp(2 * e0 * T),
			"alpha_kin_xi": 4 * m / DXI**2, "runs_A": runs_a(p, x, shifted, m), "runs_B": runs_b(p, shifted, m),
			"runs_B_plain": runs_b(p, plain, m), "grid_points_N64": 64.0**m})
		r = rows[-1]
		print(f"M={m:2d} ratio={ratio:.3e} P plain={plain:.2e} shifted={shifted:.2f} runs A={r['runs_A']:.2e} B={r['runs_B']:.2e} (B plain {r['runs_B_plain']:.2e})")
	ab = np.log([r["runs_A"] for r in rows])
	print("growth of readout A per added asset (M>=4): x%.2f" % np.exp((ab[-1] - ab[3]) / (rows[-1]["M"] - rows[3]["M"])))
	# Every mode of C with eigenvalue mu contributes a factor; M - 1 modes have mu = 1 - rho.
	k, phi = np.sqrt(0.6), 1.5
	y, z = np.cosh(k * T) + phi / (2 * k) * np.sinh(k * T), phi * np.cosh(k * T) + 2 * k * np.sinh(k * T)
	f = np.sqrt(phi * y / z) / y
	factors = {"plain": f, "margin": f * np.exp(SHIFT_MARGIN * k * T), "exact": f * np.exp(k * T)}
	for name, key in (("plain", "P_succ_plain"), ("margin", "P_succ_shifted"), ("exact", "P_succ_exact")):
		tail = (rows[-1][key] / rows[-2][key]) ** 0.5  # the one remaining mode changes slowly
		assert abs(tail / factors[name] - 1) < 0.15, (name, tail, factors[name])
	print(factors)
	save("scaling_m", {"rows": rows, "factors": factors})
