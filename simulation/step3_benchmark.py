"""Step 3: verify the classical benchmark (WorkPlan.tex Sec. Step 3).

Checks: Eq. 16 mapping and Riccati reduction (symbolic), exact-liquidation
cost/rates/sales-only bound, direct quadrature, the discrete block-tridiagonal
convergence table, and the noncommuting finite-penalty block formula vs RK4.
"""

import sys

import numpy as np
import sympy as sp
from scipy import sparse
from scipy.integrate import quad
from scipy.sparse.linalg import spsolve

from model import C, LAM, PHI, ROOT, T, X0, block_solution, log_det_Y, riccati_P, save

sys.path.insert(0, str(ROOT.parent / "stage2"))
from BellmanPrinciple import derive_hjb  # noqa: E402


def symbolic_checks():
	"""Re-derive the two-asset HJB, the Riccati system, and the Eq. 16 reduction."""
	x1, x2, v1, v2, t, Tsym, lam = sp.symbols("x1 x2 v1 v2 t T lambda", real=True)
	e1, e2 = sp.symbols("eta1 eta2", positive=True)
	c11, c12, c22 = sp.symbols("c11 c12 c22", real=True)
	Hm, Cm = sp.diag(e1, e2), sp.Matrix([[c11, c12], [c12, c22]])
	x, v = sp.Matrix([x1, x2]), sp.Matrix([v1, v2])
	res = derive_hjb((v.T * Hm * v + lam * x.T * Cm * x)[0], [x1, x2], [v1, v2],
		[-v1, -v2], t, Tsym)
	V = res.value_function
	grad = sp.Matrix([sp.diff(V, x1), sp.diff(V, x2)])
	target = sp.diff(V, t) - (grad.T * Hm.inv() * grad)[0] / 4 + lam * (x.T * Cm * x)[0]
	assert sp.simplify(res.pde.lhs - target) == 0, "HJB differs from Step 2 Eq. pde"
	assert all(sp.simplify(res.optimal_control[vi] - (Hm.inv() * grad / 2)[i]) == 0
		for i, vi in enumerate([v1, v2])), "feedback differs from H^{-1} grad V / 2"

	# Quadratic ansatz V = x^T P x / 2 + c reduces the HJB to the Step 2 Riccati system.
	p11, p12, p22, c = (sp.Function(n)(t) for n in ("p11", "p12", "p22", "c"))
	P = sp.Matrix([[p11, p12], [p12, p22]])
	Vq = (x.T * P * x)[0] / 2 + c
	gq = sp.Matrix([sp.diff(Vq, x1), sp.diff(Vq, x2)])
	hjb_q = sp.expand(sp.diff(Vq, t) - (gq.T * Hm.inv() * gq)[0] / 4 + lam * (x.T * Cm * x)[0])
	ric = P.diff(t) - P * Hm.inv() * P / 2 + 2 * lam * Cm  # must vanish (Eq. ric)
	resid = sp.expand(hjb_q - ((x.T * ric * x)[0] / 2 + sp.diff(c, t)))
	assert resid == 0, "Riccati reduction failed"

	# BDG Eq. 16 with R = 0, theta = -V, q = x, eta = H, Sigma = C, gamma = 2 lam.
	theta = -V
	gth = sp.Matrix([sp.diff(theta, x1), sp.diff(theta, x2)])
	eq16 = sp.diff(theta, t) + (gth.T * Hm.inv() * gth)[0] / 4 - lam * (x.T * Cm * x)[0]
	assert sp.simplify(eq16 + target) == 0, "Eq. 16 mapping sign/factor error"
	return {"hjb_from_bellman": True, "riccati_reduction": True, "eq16_mapping": True}


def liquidation_checks():
	"""Exact-liquidation instance H=I, X=(3,1), x(1)=0 (Step 3 boxed equation)."""
	mu, U = np.linalg.eigh(C)
	order = np.argsort(-mu)
	mu, U = mu[order], U[:, order] * np.sign(U[0, order])  # columns (u+, u-)
	kap = np.sqrt(LAM * mu)
	y0 = U.T @ X0

	def x_of_t(t):
		return U @ (y0 * np.sinh(kap * (T - t)) / np.sinh(kap * T))

	def v_of_t(t):
		return U @ (y0 * kap * np.cosh(kap * (T - t)) / np.sinh(kap * T))

	v0 = v_of_t(0.0)
	j_closed = float(np.sum(kap / np.tanh(kap * T) * y0**2))
	j_quad = quad(lambda t: v_of_t(t) @ v_of_t(t) + LAM * x_of_t(t) @ C @ x_of_t(t),
		0, T, epsabs=1e-13, epsrel=1e-13, limit=200)[0]
	r0 = kap[1] * np.sinh(kap[0]) / (kap[0] * np.sinh(kap[1]))
	theta_max = (r0 + 1) / (r0 - 1)
	ts = np.linspace(0, T, 200001)
	vs = np.array([v_of_t(t) for t in ts])
	ratio = (kap[1] * np.cosh(kap[1] * (T - ts))) / (kap[0] * np.cosh(kap[0] * (T - ts)))
	assert np.all(np.diff(ratio) >= -1e-15), "a-/a+ not monotone in t"
	out = {
		"kappa": kap.tolist(), "v0": v0.tolist(), "J_closed": j_closed, "J_quadrature": j_quad,
		"R0": r0, "theta_max": theta_max, "min_v": vs.min(axis=0).tolist(),
		"argmin_t": ts[vs.argmin(axis=0)].tolist(),
	}
	assert abs(j_closed - 13.810512031917) < 1e-11 and abs(j_quad - j_closed) < 1e-10
	assert np.allclose(v0, [4.0488436487, 1.6639810858], atol=1e-10)
	assert np.allclose(out["min_v"], [2.5062988967, 0.6931333103], atol=1e-9)
	assert abs(theta_max - 15.989177597983) < 1e-10 and abs(r0 - 1.133429601920) < 1e-11
	return out, x_of_t


def discrete_table(x_of_t, sizes=(16, 64, 256, 1024)):
	"""Unique minimizer of the discrete objective via its block-tridiagonal system."""
	rows = []
	for n in sizes:
		tau = T / n
		main = sparse.kron(sparse.eye(n - 1), 2 * np.eye(2) + tau**2 * C)
		off = sparse.kron(sparse.diags([np.ones(n - 2), np.ones(n - 2)], [-1, 1]), np.eye(2))
		rhs = np.zeros(2 * (n - 1))
		rhs[:2] = X0
		xs = np.vstack([X0, spsolve((main - off).tocsc(), rhs).reshape(n - 1, 2), np.zeros(2)])
		dx = np.diff(xs, axis=0)
		cost = np.sum(dx**2) / tau + tau * np.einsum("ki,ij,kj->", xs[1:], C, xs[1:])
		err = np.max(np.abs(xs - np.array([x_of_t(k * tau) for k in range(n + 1)])))
		rows.append({"N": n, "max_inventory_error": err, "discrete_cost": cost})
	ref = [(6.945e-5, 13.432845715492), (4.361e-6, 13.714251771181),
		(2.726e-7, 13.786331703645), (1.704e-8, 13.804459745777)]
	for r, (e, j) in zip(rows, ref):
		assert abs(r["max_inventory_error"] / e - 1) < 1e-3 and abs(r["discrete_cost"] - j) < 1e-9
	return rows


def finite_penalty_check():
	"""Noncommuting Phi: exact block exponential vs RK4 on P' = 2 lam C - P P / 2."""
	assert not np.allclose(C @ PHI, PHI @ C)

	def rhs(p):
		return 2 * LAM * C - 0.5 * p @ p

	samples = np.arange(17) / 16
	exact = np.array([riccati_P(s) for s in samples])
	errs = {}
	for steps in (128, 256, 4096):
		h, p, path = T / steps, PHI.copy(), [PHI.copy()]
		c_int = 0.0
		for k in range(steps):
			k1 = rhs(p); k2 = rhs(p + h / 2 * k1); k3 = rhs(p + h / 2 * k2); k4 = rhs(p + h * k3)
			# c(s)/nu = int tr(A P) dr; Simpson on the RK4 stages for the scalar integral.
			c_int += h / 6 * (np.trace(p) + 2 * np.trace(p + h / 2 * k1)
				+ 2 * np.trace(p + h / 2 * k2) + np.trace(p + h * k3))
			p = p + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
			if (k + 1) % (steps // 16) == 0:
				path.append(p.copy())
		errs[steps] = {"P": float(np.max(np.abs(np.array(path) - exact))),
			"c_over_nu": abs(c_int - 2 * log_det_Y(T))}
	p1 = riccati_P(T)
	nu = 0.2
	v = 0.5 * p1 @ X0
	out = {"P1": p1.tolist(), "c1_nu0.2": 2 * nu * log_det_Y(T), "v_star": v.tolist(),
		"V": 0.5 * X0 @ p1 @ X0, "rk4_errors": errs, "detY1": np.linalg.det(block_solution(T)[0])}
	assert np.allclose(p1, [[1.9643628285, 0.4223113353], [0.4223113353, 1.7879582164]], atol=1e-10)
	assert abs(out["c1_nu0.2"] - 0.6939899962) < 1e-10
	assert np.allclose(v, [3.1576999103, 1.5274461111], atol=1e-10)
	assert abs(out["V"] - 11.0005458421) < 1e-10
	assert errs[4096]["P"] < 1e-13 and errs[256]["P"] < errs[128]["P"] / 10
	return out


if __name__ == "__main__":
	sym = symbolic_checks()
	liq, x_of_t = liquidation_checks()
	table = discrete_table(x_of_t)
	fp = finite_penalty_check()
	save("step3", {"symbolic": sym, "liquidation": liq, "discrete": table, "finite_penalty": fp})
	print("symbolic:", sym)
	print("liquidation:", {k: liq[k] for k in ("v0", "J_closed", "J_quadrature", "min_v", "theta_max")})
	for r in table:
		print(f"N={r['N']:5d}  err={r['max_inventory_error']:.4e}  J_N={r['discrete_cost']:.12f}")
	print("finite penalty:", {k: fp[k] for k in ("v_star", "V", "c1_nu0.2")})
	print("RK4 errors:", fp["rk4_errors"])
