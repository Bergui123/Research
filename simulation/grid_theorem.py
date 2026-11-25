"""Numerical companion of the grid-expansion theorem (paper/Paper.tex, Theorem 3).

In xi = x / sqrt(nu) the theorem reads, with delta = dx / sqrt(nu) and psi_k = p_k psi,
  sup |psi_d - sum_{k<=K} delta^{2k} psi_k| <= B_K delta^{2K+2}.
This file (1) integrates the polynomials p_1, p_2, (2) evaluates B_1, B_2, (3) checks the closed
forms and signs of the quartic/quadratic parts, (4) checks the orders and the bounds on grids.
Sup norms and time integrals are evaluated on fine grids in floating point, not with intervals.
"""

import functools
import math

import numpy as np
from scipy.integrate import quad_vec, solve_ivp
from scipy.sparse.linalg import expm_multiply

from model import C, LAM, PHI, T, block_solution, grid_operator, log_det_Y, riccati_P, save

DEG = 8  # p_2 has degree <= 8
N1 = DEG + 1


def d(p, a):
	"""d/dxi_a of a coefficient array p[i, j] xi1^i xi2^j."""
	q = np.zeros_like(p)
	k = np.arange(1, p.shape[a])
	if a == 0:
		q[:-1] = k[:, None] * p[1:]
	else:
		q[:, :-1] = k[None, :] * p[:, 1:]
	return q


def times_linear(p, f):
	"""(f . xi) p."""
	q = np.zeros_like(p)
	q[1:] += f[0] * p[:-1]
	q[:, 1:] += f[1] * p[:, :-1]
	return q


def ladder(p, P, a, m):
	"""psi^{-1} d_a^m (p psi) for psi = const exp(-xi^T P xi / 8)."""
	for _ in range(m):
		p = d(p, a) + times_linear(p, -P[a] / 4)
	return p


def pad(p, size):
	q = np.zeros((size, size))
	n = min(len(p), size)
	q[:n, :n] = p[:n, :n]
	return q


def c(j):
	return 2 / math.factorial(2 * j + 2)


def rhs(s, y):
	P = riccati_P(s)
	p1, p2 = y.reshape(2, N1, N1)
	one = np.zeros((N1, N1))
	one[0, 0] = 1
	gen = lambda p: sum(d(d(p, a), a) - 0.5 * times_linear(d(p, a), P[a]) for a in range(2))
	f1 = sum(c(1) * ladder(one, P, a, 4) for a in range(2))
	f2 = sum(c(1) * ladder(p1, P, a, 4) + c(2) * ladder(one, P, a, 6) for a in range(2))
	return np.concatenate([(gen(p1) + f1).ravel(), (gen(p2) + f2).ravel()])


@functools.lru_cache(maxsize=None)
def solution():
	return solve_ivp(rhs, (0, T), np.zeros(2 * N1 * N1), rtol=1e-11, atol=1e-14, dense_output=True).sol


def polys(s):
	return solution()(s).reshape(2, N1, N1)


def evaluate(p, xi):
	return np.polynomial.polynomial.polyval2d(xi[..., 0], xi[..., 1], p)


def psi(xi, s):
	return np.exp(-np.einsum("...i,ij,...j->...", xi, riccati_P(s), xi) / 8 - 0.5 * log_det_Y(s))


def expansion(xi, s, delta, order):
	"""Partial sum psi (1 + delta^2 p_1 + delta^4 p_2) truncated after `order` corrections."""
	p1, p2 = polys(s)
	out = 1.0 + (order >= 1) * delta**2 * evaluate(p1, xi) + (order >= 2) * delta**4 * evaluate(p2, xi)
	return psi(xi, s) * out


@functools.lru_cache(maxsize=None)
def remainder_constant(order, n_s=33, half=14.0, step=0.04):
	"""B_K = int sum_a sum_k 2/(2(K-k)+4)! ||d_a^{2(K-k)+4} psi_k||_inf ds (Simpson in s)."""
	ax = np.arange(-half, half + step / 2, step)
	xi = np.stack(np.meshgrid(ax, ax, indexing="ij"), axis=-1)
	one = np.zeros((1, 1)) + 1
	vals = []
	for s in np.linspace(0, T, n_s):
		P, g = riccati_P(s), psi(xi, s)
		p = [one, *polys(s)]
		total = 0.0
		for k in range(order + 1):
			m = 2 * (order - k) + 4
			size = DEG + m + 1
			for a in range(2):
				total += 2 / math.factorial(m) * np.max(np.abs(evaluate(ladder(pad(p[k], size), P, a, m), xi) * g))
		vals.append(total)
	w = np.ones(n_s)
	w[1:-1:2], w[2:-1:2] = 4, 2
	return float(w @ np.array(vals) * T / (n_s - 1) / 3)


def closed_forms(xi, s=T):
	"""Quartic and quadratic parts of p_1 from the integral formulas of the theorem (nu = 1)."""
	Ys_inv = np.linalg.inv(block_solution(s)[0])

	def parts(t):
		Yt, Zt = block_solution(t)
		w = Zt @ Ys_inv @ xi
		gram = quad_vec(lambda r: (lambda Yi: Yi @ Yi.T)(np.linalg.inv(block_solution(r)[0])), t, s, epsabs=1e-12)[0] / 2
		bracket = np.diag(Zt @ np.linalg.inv(Yt) - Zt @ gram @ Zt.T)
		assert np.all(bracket > 0)
		return np.array([np.sum(w**4) / 3072, -np.sum(w**2 * bracket) / 128])
	return quad_vec(parts, 0, s, epsabs=1e-12)[0]


def positivity_check(trials=200, m=4):
	"""P_t^{-1} - int_t^s F^{-1} F^{-T}/2 dr > 0 for random H = I problems in dimension m."""
	rng = np.random.default_rng(3)
	worst = np.inf
	for _ in range(trials):
		a, b = rng.standard_normal((2, m, m))
		cov, phi = a @ a.T, b @ b.T + 0.05 * np.eye(m)
		t, s = np.sort(rng.uniform(0, 2, 2))
		y = lambda r: block_solution(r, phi, cov, 1.0, np.eye(m))
		Yt, Zt = y(t)
		gram = quad_vec(lambda r: (lambda Fi: Fi @ Fi.T)(Yt @ np.linalg.inv(y(r)[0])), t, s, epsabs=1e-10)[0] / 2
		p = Zt @ np.linalg.inv(Yt)
		worst = min(worst, np.linalg.eigvalsh(np.linalg.inv((p + p.T) / 2) - (gram + gram.T) / 2)[0])
	return float(worst)


def order_check():
	"""Sup-norm residuals against the infinite lattice (box large enough for round-off)."""
	rows = []
	for delta in (0.5, 0.4, 0.3, 0.2, 0.15):
		half = 14.0
		n = 2 * int(half / delta)
		x = (np.arange(n) - n // 2) * delta + 0.3 * delta  # arbitrary lattice shift
		xi = np.stack(np.meshgrid(x, x, indexing="ij"), axis=-1)
		hd, _ = grid_operator(x, x, delta, 1.0)
		u = expm_multiply(-T * hd, psi(xi, 0.0).ravel()).reshape(n, n)
		inner = np.max(np.abs(xi), axis=-1) < 9
		res = [float(np.max(np.abs(u - expansion(xi, T, delta, k))[inner])) for k in (0, 1, 2)]
		rows.append({"delta": delta, "res": res,
			"bound": [remainder_constant(k) * delta ** (2 * k + 2) for k in (0, 1, 2)]})
	return rows


if __name__ == "__main__":
	one = [polys(T)[0], polys(T)[1]]
	rng = np.random.default_rng(1)
	p1 = one[0]
	quart = np.array([[p1[i, j] if i + j == 4 else 0 for j in range(N1)] for i in range(N1)])
	quadr = np.array([[p1[i, j] if i + j == 2 else 0 for j in range(N1)] for i in range(N1)])
	for _ in range(5):
		xi = rng.standard_normal(2) * 2
		q4, q2 = closed_forms(xi)
		assert abs(q4 - evaluate(quart, xi)) < 1e-8 * (1 + abs(q4)) and abs(q2 - evaluate(quadr, xi)) < 1e-8 * (1 + abs(q2))
		assert q4 > 0 > q2
	print("closed forms match the ODE; q4 > 0 > q2")
	pos = positivity_check()
	print("min eigenvalue of P^-1 - J over random 4-asset problems:", pos)
	assert pos > 0
	from grid_error import chi_coefficients
	mine = [p1[i, j] for (i, j) in [(4, 0), (3, 1), (2, 2), (1, 3), (0, 4), (2, 0), (1, 1), (0, 2)]]
	assert np.allclose(mine, chi_coefficients(1.0), rtol=1e-7)  # independent earlier implementation
	print("B_0, B_1, B_2 =", [remainder_constant(k) for k in (0, 1, 2)])
	rows = order_check()
	for a, b in zip(rows, rows[1:]):
		orders = [math.log(a["res"][k] / b["res"][k]) / math.log(a["delta"] / b["delta"]) for k in (0, 1, 2)]
		print(f"delta {a['delta']}->{b['delta']}: res {['%.2e' % r for r in b['res']]} bound {['%.2e' % r for r in b['bound']]} orders {['%.2f' % o for o in orders]}")
	assert all(r["res"][k] <= r["bound"][k] for r in rows for k in (0, 1, 2) if r["res"][k] > 1e-12)
	save("grid_theorem", {"B": [remainder_constant(k) for k in (0, 1, 2)], "orders": rows, "min_eig_random": pos})
