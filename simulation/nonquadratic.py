"""Off the quadratic case: terminal penalty g(y) = y^T Phi y / 2 + KAPPA sum_a y_a^4.

The running cost stays quadratic, so the operator H_d and its circuit are those of the benchmark.
The kernel is Gaussian (Mehler), which gives the viscous rate v_nu and the inviscid rate v_0 by a
two-dimensional quadrature, respectively a Newton solve: the references are exact for every nu.
  direct:     read -2 nu grad log psi at X (continuum model of the number of runs, as scaling_m.py);
  recentred:  psi = exp(alpha.x + gamma) phi(x - xi(s), s) along the extremal of the *quadratic* model;
              phi solves the same equation and peaks near 0 (grid simulation, costed and emulated).
"""

import functools
import math
import os

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from chebyshev import p_grid, recover
from model import C, E0, LAM, P_SHIFT, PHI, T, X0, block_solution, grid_operator, save
from step10_validation import sample_qae
from step8_resources import SHARE, ae_rel_error, amplitude_budget, coherent_run, median_reps, potential_scale

KAPPA = float(os.environ.get("EXEC_KAPPA", "0.1"))
TARGETS = (0.01, 0.001)
KINDS, ETAS = ("smooth_w0.5", "smooth_w1", "smooth_w2"), (5.0, 10.0, 20.0)
MU, Q = np.linalg.eigh(LAM * C)
K = np.sqrt(MU)
COTH, CSCH = Q * (K / np.tanh(K * T)) @ Q.T, Q * (K / np.sinh(K * T)) @ Q.T
LOG_KERNEL = lambda nu: 0.5 * np.sum(np.log(K / (4 * np.pi * nu * np.sinh(K * T))))


def g(y, kappa):
	return 0.5 * np.einsum("...i,ij,...j->...", y, PHI, y) + kappa * np.sum(y**4, axis=-1)


def action(y, x, kappa):
	"""Cost of the quadratic-running-cost extremal from terminal point y to x, plus g(y)."""
	return x @ COTH @ x + np.einsum("...i,ij,...j->...", y, COTH, y) - 2 * y @ (CSCH @ x) + g(y, kappa)


def terminal_point(x, kappa):
	y = np.linalg.solve(2 * COTH + PHI, 2 * CSCH @ x)
	for _ in range(80):
		hess = 2 * COTH + PHI + np.diag(12 * kappa * y**2)
		y = y - np.linalg.solve(hess, 2 * COTH @ y - 2 * CSCH @ x + PHI @ y + 4 * kappa * y**3)
	return y, hess


def inviscid_rate(x, kappa):
	return COTH @ x - CSCH @ terminal_point(x, kappa)[0]


def viscous(x, nu, kappa):
	"""log psi(x, T) and the viscous rate, by quadrature on the whitened terminal point."""
	y, hess = terminal_point(x, kappa)
	chol = np.linalg.cholesky(np.linalg.inv(hess) * 4 * nu)
	t = np.linspace(-9, 9, 181)
	a, b = np.meshgrid(t, t, indexing="ij")
	pts = y + np.stack([a, b], -1) @ chol.T
	f0 = action(y, x, kappa)
	f = np.exp(-(action(pts, x, kappa) - f0) / (4 * nu))
	log_psi = LOG_KERNEL(nu) - f0 / (4 * nu) + math.log(f.sum() * (t[1] - t[0]) ** 2 * np.linalg.det(chol))
	return log_psi, COTH @ x - CSCH @ ((f[..., None] * pts).sum((0, 1)) / f.sum())


V0 = inviscid_rate(X0, KAPPA)
NORM = float(np.linalg.norm(V0))
rel = lambda d: float(np.linalg.norm(d) / NORM)


def bias(nu, kappa=KAPPA):
	v0 = inviscid_rate(X0, kappa)
	return float(np.linalg.norm(viscous(X0, nu, kappa)[1] - v0) / np.linalg.norm(v0))


@functools.lru_cache(maxsize=None)
def norms(nu, n=128, half=6.5):
	"""Continuum ||psi(T)||^2 and P_succ of the direct problem (grid sums in xi = x / sqrt(nu))."""
	dx = 2 * half * math.sqrt(nu) / n
	x = (np.arange(n) - n // 2) * dx
	pts = np.stack(np.meshgrid(x, x, indexing="ij"), -1)
	u0 = np.exp(-g(pts, KAPPA) / (4 * nu)).ravel()
	ut = expm_multiply(-T * grid_operator(x, x, dx, nu)[0], u0)
	return np.sum(ut**2) * dx**2, 0.5 * math.exp(-2 * P_SHIFT) * np.sum(ut**2) / np.sum(u0**2)


@functools.lru_cache(maxsize=None)
def direct_level(nu, d, dxi=0.25):
	"""Rate read with span d, the four amplitudes b_j at spacing dxi sqrt(nu), smallest p_j."""
	norm2, p_succ = norms(nu)
	nodes = [X0 + s * d * e for e in np.eye(2) for s in (1, -1)]
	lp = np.array([viscous(x, nu, KAPPA)[0] for x in nodes])
	rate = -2 * nu * np.array([lp[0] - lp[1], lp[2] - lp[3]]) / (2 * d)
	log_p = 2 * lp + 2 * math.log(dxi * math.sqrt(nu)) - math.log(norm2)
	return rate, 0.5 * (math.log(p_succ) + log_p)


def runs_for(log_b, budget, d, nu):
	delta = SHARE["ae"] * amplitude_budget(budget * NORM * d / (2 * math.sqrt(2) * nu))
	return float(np.sum(2 * np.pi / delta * np.exp(-np.minimum(log_b, 0))))


def direct(target, levels):
	"""Cheapest direct readout; `levels` > 1 extrapolates nu, 2 nu, 4 nu, ... to 0 (Richardson)."""
	best = None
	for nu in 0.0078125 * 2.0 ** (np.arange(41) / 4):  # a quarter-octave grid, so that levels share points
		nus = nu * 2.0 ** np.arange(levels)
		w = np.array([np.prod([m / (m - k) for m in nus if m != k]) for k in nus])
		for c in (0.25, 0.5, 1.0):
			data = [direct_level(k, c * math.sqrt(k)) for k in nus]
			err = rel(sum(a * r for a, (r, _) in zip(w, data)) - V0)
			if err > target / 2:
				continue
			runs = sum(runs_for(lb, (target - err) / np.sum(np.abs(w)), c * math.sqrt(k), k) for k, (_, lb) in zip(nus, data))
			if best is None or runs < best["runs"]:
				best = {"target": target, "levels": levels, "nu": float(nu), "d": c * math.sqrt(nu), "runs": runs,
					"deterministic": err, "min_b": float(np.exp(min(lb.min() for _, lb in data)))}
	return best


# --- recentred ---
XI0 = np.linalg.solve(block_solution(T)[0], X0)  # terminal point of the quadratic model
PI0, PIT = 2 * (CSCH @ X0 - COTH @ XI0), 2 * (COTH @ X0 - CSCH @ XI0)  # costate of that extremal
assert np.allclose(PI0, PHI @ XI0)


def recentred_grid(nu, n, half):
	"""Grid centred on the characteristic: x, dx, H_d, potential, initial and final phi."""
	dx = 2 * half * math.sqrt(nu) / 2**n
	x = (np.arange(2**n) - 2 ** (n - 1)) * dx
	pts = np.stack(np.meshgrid(x, x, indexing="ij"), -1)
	log0 = -(g(XI0 + pts, KAPPA) - pts @ PI0) / (4 * nu)
	u0 = np.exp(log0 - log0.max()).ravel()
	hd, pot = grid_operator(x, x, dx, nu)
	return x, dx, hd, pot, u0, expm_multiply(-T * hd, u0)


def recentred(target):
	best = None
	for nu in (0.0005, 0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5):
		v_nu = viscous(X0, nu, KAPPA)[1]
		e_bias = rel(v_nu - V0)
		if e_bias > target / 2:
			continue
		for n in (4, 5, 6, 7):
			for half in (4.5, 5.5, 6.5):
				x, dx, hd, pot, u0, ut = recentred_grid(nu, n, half)
				c, big = 2 ** (n - 1), 2**n
				rows = []
				for h in (1, 2, 4, 8):
					if h >= c // 2:
						continue
					d = h * dx
					idx = [(c + h) * big + c, (c - h) * big + c, c * big + c + h, c * big + c - h]
					lp = np.log(ut[idx])
					v_grid = PIT / 2 - 2 * nu * np.array([lp[0] - lp[1], lp[2] - lp[3]]) / (2 * d)
					le = np.array([viscous(X0 + s * d * e, nu, KAPPA)[0] for e in np.eye(2) for s in (1, -1)])
					v_span = -2 * nu * np.array([le[0] - le[1], le[2] - le[3]]) / (2 * d)
					rows.append((h, d, idx, rel(v_span - v_nu), rel(v_grid - v_span)))
				# Charged in norm, with the largest grid error over the spans: no cancellation is used.
				e_grid = max(r[4] for r in rows)
				for h, d, idx, e_span, _ in rows:
					eps = e_bias + e_span + e_grid
					if eps > target / 2:
						continue
					delta = amplitude_budget((target - eps) * NORM * d / (2 * math.sqrt(2) * nu))
					p_j = ut[idx] ** 2 / np.sum(ut**2)
					ratio = np.sum(ut**2) / np.sum(u0**2)
					proxy = 1 / (delta * math.sqrt(p_j.min() * ratio)) * (8 * nu / dx**2)
					if best is None or proxy < best["proxy"]:
						best = {"proxy": proxy, "target": target, "nu": nu, "n": n, "half": half, "h": h, "dx": dx,
							"e_bias": e_bias, "e_span": e_span, "e_grid": e_grid, "eps_det": eps, "delta_tot": delta,
							"p_min": float(p_j.min()), "reference_error": rel(PIT / 2 - V0)}
	return cost(best) if best else None  # None: the reference is too far from the optimum for this bias


def cost(b):
	"""Cheapest auxiliary profile for the chosen grid: T count, qubits and an emulation of the readout."""
	nu, n, h, dx, delta = b["nu"], b["n"], b["h"], b["dx"], b["delta_tot"]
	x, dx, hd, pot, u0, ut = recentred_grid(nu, n, b["half"])
	c, big = 2 ** (n - 1), 2**n
	idx = [(c + h) * big + c, (c - h) * big + c, c * big + c + h, c * big + c - h]
	v_scale, k_bits, q5 = potential_scale(nu, x, x, dx)
	m = min(math.ceil(math.log2(T * v_scale / math.log1p(SHARE["pot"] * delta))), k_bits)
	v_m = v_scale * np.floor(q5 * 2.0 ** (m - k_bits)) / 2**m
	hd_m = hd - diags(pot.ravel()) + diags(v_m.ravel())
	alpha_kin = 8 * nu / dx**2
	alpha = alpha_kin + v_scale
	best = None
	for kind in KINDS:
		for eta in ETAS:
			n_p = math.ceil(math.log2((alpha * T + P_SHIFT + 23.0) * eta / math.pi))
			w, gp, win = recover(hd_m, u0, alpha, n_p, eta, kind, P_SHIFT, idx)
			total = np.sum(gp**2) * np.sum(u0**2)
			a = np.sum(w**2, axis=0) / total
			ideal = np.sum(np.exp(-2 * p_grid(n_p, eta, P_SHIFT)[0][win])) * ut[idx] ** 2 / total
			amp = float(np.max(np.abs(np.sqrt(a / ideal) - 1)))
			if amp > (SHARE["schr"] + SHARE["pot"]) * delta:
				continue
			b_amp = np.sqrt(ideal)
			run = coherent_run(n, 0, alpha_kin, v_scale, m, eta, P_SHIFT, SHARE["sim"] * delta * b_amp.min())
			ts = [next(t for t in range(1, 60) if ae_rel_error(q, t) <= SHARE["ae"] * delta) for q in b_amp * (1 - 0.3 * delta)]
			runs = sum(2 * 2**t + 1 for t in ts)
			if best is None or runs * run["T_run"] < best["T_single"]:
				p_succ = float(np.sum(np.exp(-2 * p_grid(n_p, eta, P_SHIFT)[0][win])) * np.sum(ut**2) / total)
				best = {**b, **run, "kind": kind, "eta_max": eta, "P_succ": p_succ, "amp_rel_err": amp, "amp_budget": 0.2 * delta,
					"AE_log2_iterations": ts, "runs": runs, "T_single": runs * run["T_run"], "a": a.tolist(), "b_min": float(b_amp.min())}
	r = median_reps()
	best["T_total"], best["logical_qubits"] = r * best["T_single"], best["qubits_run"] + 1
	est = np.stack([np.median(sample_qae(a, 2**t, (2000, r)).reshape(2000, r), axis=1) for a, t in zip(best["a"], best["AE_log2_iterations"])], -1)
	lb = 0.5 * np.log(est)
	v = PIT / 2 - 2 * nu * np.stack([lb[:, 0] - lb[:, 1], lb[:, 2] - lb[:, 3]], -1) / (2 * h * dx)
	err = np.linalg.norm(v - V0, axis=-1) / NORM
	ok = int(np.sum(err <= best["target"]))
	best["within"], best["within_lower"] = ok / 2000, float(beta.ppf(0.001, ok, 2001 - ok))
	best["q95_err"] = float(np.quantile(err, 0.95))
	return best


if __name__ == "__main__":
	# The quadrature reproduces the quadratic model exactly, and the recentred grid the viscous rate.
	assert bias(1.0, 0.0) < 1e-12 and bias(0.05, 0.0) < 1e-12
	x, dx, hd, pot, u0, ut = recentred_grid(0.1, 7, 5.5)
	lp = np.log(ut.reshape(128, 128))
	v = PIT / 2 - 0.2 * np.array([lp[65, 64] - lp[63, 64], lp[64, 65] - lp[64, 63]]) / (2 * dx)
	assert rel(v - viscous(X0, 0.1, KAPPA)[1]) < 1e-3
	curve = [{"nu": float(nu), **{f"kappa{k:g}": bias(nu, k) for k in (0.0, 0.02, 0.1, 0.5)}} for nu in np.geomspace(0.01, 8, 25)]
	amps = []
	for nu in np.geomspace(0.02, 8, 25):  # both at spacing 0.25 sqrt(nu) and span d = 0.5 sqrt(nu)
		_, _, _, _, u0, ut = recentred_grid(nu, 6, 8.0)
		idx = [34 * 64 + 32, 30 * 64 + 32, 32 * 64 + 34, 32 * 64 + 30]
		rec = math.sqrt(0.5 * math.exp(-2 * P_SHIFT) / np.sum(u0**2)) * ut[idx].min()
		amps.append({"nu": float(nu), "direct": float(np.exp(direct_level(nu, 0.5 * math.sqrt(nu))[1].min())), "recentred": float(rec)})
	w0 = float(action(terminal_point(X0, KAPPA)[0], X0, KAPPA))
	out = {"kappa": KAPPA, "v0": V0.tolist(), "W0": w0, "reference_rate": (PIT / 2).tolist(), "reference_error": rel(PIT / 2 - V0), "bias": curve,
		"slope": bias(0.01) / 0.01, "direct_amplitudes": amps, "direct": [], "recentred": []}
	for eps in TARGETS:
		for levels in (1, 2, 3):
			d = direct(eps, levels)
			if d:
				out["direct"].append(d)
			print(d)
		r = recentred(eps)
		if r is None:
			print("no recentred configuration at", eps)
			continue
		out["recentred"].append(r)
		print({k: r[k] for k in ("target", "nu", "n", "half", "h", "e_bias", "e_span", "e_grid", "delta_tot", "p_min", "kind", "eta_max",
			"n_p", "alpha_H", "P_succ", "amp_rel_err", "amp_budget", "AE_log2_iterations", "runs", "T_run", "T_total", "logical_qubits", "within", "within_lower")})
		assert r["within_lower"] >= 0.95 and r["amp_rel_err"] <= r["amp_budget"]
	save("nonquadratic" if KAPPA == 0.1 else f"nonquadratic_kappa{KAPPA:g}", out)
