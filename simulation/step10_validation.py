"""Step 10: end-to-end classical emulation of the costed pointwise-readout pipeline (readout A),
and the lambda -> 0 limit. Readout B is validated inside step9_sampling.py.

For each Step 8 optimum: the grid operator with the comparator-rounded potential, the exact
Schroedingerised p-grid evolution with the costed (N_p, L_p, profile, window) through a
Chebyshev expansion (any grid size), the exact readout probabilities a_j, then canonical
amplitude estimation sampled from its exact outcome law with M_j = 2^t_j and median-of-r.
QSP truncation and gate-synthesis errors are budgeted in Step 8 but not emulated.
"""

import json

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from chebyshev import p_grid, recover
from model import A, PHI, RESULTS, T, X0, axis_grid, block_solution, grid_operator, reference_feedback, save
from step5_nu_sweep import neighbours
from step8_resources import potential_scale

RNG = np.random.default_rng(20261122)


def readout_probabilities(cfg):
	nu, n, dx = cfg["nu"], cfg["n"], cfg["dx"]
	x1, _ = axis_grid(2**n, 2**n * dx / 2, X0[0], dx)
	x2, _ = axis_grid(2**n, 2**n * dx / 2, X0[1], dx)
	hd, pot = grid_operator(x1, x2, dx, nu)
	v_scale, k_bits, q5 = potential_scale(nu, x1, x2, dx)
	v_m = v_scale * np.floor(q5 * 2.0 ** (cfg["m"] - k_bits)) / 2 ** cfg["m"]
	hd_m = hd - diags(pot.ravel()) + diags(v_m.ravel())  # what the comparator actually encodes
	g = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
	u0 = np.exp(-np.einsum("...i,ij,...j->...", g, PHI, g) / (8 * nu)).ravel()
	pts = [a * 2**n + b for a, b in neighbours(x1, x2, cfg["h"])]
	w, gp, win = recover(hd_m, u0, cfg["alpha_H"], cfg["n_p"], cfg["eta_max"], cfg["kind"], cfg["p_a"], pts)
	norm = np.sum(gp**2) * np.sum(u0**2)
	u_t = expm_multiply(-T * hd, u0)[pts]  # unrounded grid solution, for the amplitude-error check
	p, _, _ = p_grid(cfg["n_p"], cfg["eta_max"], cfg["p_a"])
	return np.sum(w**2, axis=0) / norm, np.sum(np.exp(-2 * p[win])) * u_t**2 / norm


def sample_qae(a, m_iter, size):
	"""Sample all M distinct outcomes of canonical QAE, including the tails."""
	theta = np.arcsin(np.sqrt(a)) / np.pi
	delta = np.arange(m_iter, dtype=float) / m_iter - theta
	probability = (np.sinc(m_iter * delta) / np.sinc(delta)) ** 2
	probability /= probability.sum()
	# The opposite eigenphase has y -> -y modulo M; sin^2(pi*y/M) is
	# unchanged, so one eigenphase gives the same distribution of estimates.
	y = RNG.choice(m_iter, size=size, p=probability)
	return np.sin(np.pi * y / m_iter) ** 2


def feedback(a, cfg):
	with np.errstate(divide="ignore"):
		lb = 0.5 * np.log(a)  # a rare failed estimate a = 0 gives an infinite (counted) error
	return -2 * cfg["nu"] * np.diag(A) * np.stack([lb[..., 0] - lb[..., 1], lb[..., 2] - lb[..., 3]],
		axis=-1) / (2 * cfg["h"] * cfg["dx"])


def emulate(cfg, trials=2000):
	a, a_ideal = readout_probabilities(cfg)
	v_ref = reference_feedback()[0]
	v_exact_a = feedback(a, cfg)
	r = cfg["median_reps"]
	# The costed exponents assume amplitudes no smaller than b_j (1 - 0.3 delta); check it on this grid.
	b_ratio = float(np.min(np.sqrt(a) / np.array(cfg["b_amp"])))
	assert b_ratio >= 1 - 0.3 * cfg["delta_tot"]
	est = np.stack([np.median(sample_qae(a[j], cfg["AE_iterations"][j], (trials, r)).reshape(trials, r),
		axis=1) for j in range(4)], axis=-1)
	err = np.linalg.norm(feedback(est, cfg) - v_ref, axis=-1) / np.linalg.norm(v_ref)
	single = np.stack([sample_qae(a[j], cfg["AE_iterations"][j], trials) for j in range(4)], axis=-1)
	err1 = np.linalg.norm(feedback(single, cfg) - v_ref, axis=-1) / np.linalg.norm(v_ref)
	three = np.stack([np.median(sample_qae(a[j], cfg["AE_iterations"][j], (trials, 3)).reshape(trials, 3),
		axis=1) for j in range(4)], axis=-1)
	err3 = np.linalg.norm(feedback(three, cfg) - v_ref, axis=-1) / np.linalg.norm(v_ref)
	ok = int(np.sum(err <= cfg["target"]))
	return {"target": cfg["target"], "b_ratio": b_ratio, "lower_bound": float(beta.ppf(0.001, ok, trials - ok + 1)), "median3_frac_within": float(np.mean(err3 <= cfg["target"])), "a": a.tolist(), "amp_rel_err": float(np.max(np.abs(np.sqrt(a / a_ideal) - 1))), "amp_budget": 0.2 * cfg["delta_tot"],
		"deterministic_rel_err": float(np.linalg.norm(v_exact_a - v_ref) / np.linalg.norm(v_ref)),
		"median_err": float(np.median(err)), "q95_err": float(np.quantile(err, 0.95)),
		"frac_within_target": float(np.mean(err <= cfg["target"])),
		"single_shot_frac_within": float(np.mean(err1 <= cfg["target"])),
		"single_shot_median_err": float(np.median(err1)), "trials": trials,
		"v_mean": feedback(est, cfg).mean(axis=0).tolist(), "err_samples": err.tolist()}


def zero_risk_limit(nu=4.0, n=7, half_width=12.0):
	"""lambda = 0: grid PDE vs Riccati; Phi = alpha I, alpha -> inf gives straight-line v = X/T."""
	out = []
	for alpha_pen in (None, 10.0, 100.0):
		phi = PHI if alpha_pen is None else alpha_pen * np.eye(2)
		x1, dx = axis_grid(2**n, half_width, X0[0])
		x2, _ = axis_grid(2**n, half_width, X0[1])
		hd, _ = grid_operator(x1, x2, dx, nu, lam=0.0)
		g = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
		psi = expm_multiply(-T * hd, np.exp(-np.einsum("...i,ij,...j->...", g, phi, g) / (8 * nu)).ravel())
		psi = psi.reshape(2**n, 2**n)
		i, j = int(np.argmin(abs(x1 - X0[0]))), int(np.argmin(abs(x2 - X0[1])))
		grad = np.array([np.log(psi[i + 1, j] / psi[i - 1, j]), np.log(psi[i, j + 1] / psi[i, j - 1])]) / (2 * dx)
		v_grid = -2 * nu * grad
		y, z = block_solution(T, phi, lam=0.0)
		v_ric = 0.5 * (z @ np.linalg.inv(y)) @ X0
		out.append({"Phi": "benchmark" if alpha_pen is None else f"{alpha_pen:g} I", "v_grid": v_grid.tolist(),
			"v_riccati": v_ric.tolist(), "rel_err": float(np.linalg.norm(v_grid - v_ric) / np.linalg.norm(v_ric)),
			"straight_line": (X0 / T).tolist()})
		print(out[-1])
	return out


def main():
	best = json.loads((RESULTS / "step8.json").read_text())["best"]
	runs = []
	for key in sorted(best, key=float, reverse=True):
		res = emulate(best[key])
		runs.append(res)
		print(f"eps={key}: amp err={res['amp_rel_err']:.1e} (budget {res['amp_budget']:.1e}) deterministic err={res['deterministic_rel_err']:.2e}, median={res['median_err']:.2e}, "
			f"q95={res['q95_err']:.2e}, within target={res['frac_within_target']:.3f} "
			f"(single-shot {res['single_shot_frac_within']:.3f}, median-of-3 {res['median3_frac_within']:.3f})")
		assert res["lower_bound"] >= 0.95
		assert res["amp_rel_err"] <= res["amp_budget"]  # Schroedingerisation + rounding shares
	save("step10", {"emulation": runs, "zero_risk": zero_risk_limit() if RESULTS.name == "results" else None})


if __name__ == "__main__":
	main()
