"""The same problem without the auxiliary variable: e^{-T H_d} as a polynomial of the same walk.

H_d is Hermitian with spectrum in [0, 2 alpha], x = H_d/alpha - 1 in [-1, 1] is what the walk encodes, and
  e^{-T H_d} = e^{-beta (1 + x)},  beta = alpha T,  = e^{-beta} [I_0(beta) + 2 sum_k (-1)^k I_k(beta) T_k(x)].
The degree-K series has uniform Bessel tail tau and norm at most 1 + tau. The resource estimate scales it by
1/(1 + tau), giving a contraction with error at most 2 tau, then uses 2K GQSP signal calls on one signal qubit.
GQSP phase synthesis is not compiled here.
Same grids, same amplitude budget and same arithmetic as the costed configurations; the share of the
Schroedingerisation error goes to the polynomial. Success amplitude psi(T)/||psi(0)||, no shift.
"""

import json
import math

import numpy as np
from qualtran.bloqs.mcmt import MultiControlZ
from scipy.sparse.linalg import expm_multiply
from scipy.special import ive

from model import PHI, RESULTS, T, X0, axis_grid, grid_operator, save
from step5_nu_sweep import neighbours
from step8_resources import SHARE, ae_rel_error, median_reps, t_rot, tcount
from step9_sampling import runs_needed


def order(beta, tol):
	"""Smallest K with 2 e^{-beta} sum_{k > K} I_k(beta) <= tol."""
	k = np.arange(0, int(40 * math.sqrt(beta) + 400))
	tail = 2 * np.cumsum(ive(k, beta)[::-1])[::-1]
	return int(np.argmax(tail <= tol)) - 1


def bessel_tail(beta, k_max):
	"""Uniform tail bound for the degree-k_max semigroup polynomial."""
	upper = max(k_max + 1, int(beta + 40 * math.sqrt(beta) + 400))
	return float(2 * np.sum(ive(np.arange(k_max + 1, upper + 1), beta)))


def series(hd, u0, alpha, k_max):
	"""Apply the truncated series (used to check the order)."""
	beta = alpha * T
	x = lambda v: (hd @ v) / alpha - v
	prev, cur = u0, x(u0)
	out = ive(0, beta) * prev - 2 * ive(1, beta) * cur
	for k in range(2, k_max + 1):
		prev, cur = cur, 2 * x(cur) - prev
		out += 2 * (-1) ** k * ive(k, beta) * cur
	return out


def run(cfg, frac, eps_state):
	"""T count and qubits of one coherent run; cfg is a costed configuration."""
	n, m, alpha = cfg["n"], cfg["m"], cfg["alpha_walk"]
	k = order(alpha * T, eps_state / 4)  # rescaling below 1 doubles the tail
	tail = bessel_tail(alpha * T, k)
	scale = 1 / (1 + tail)
	steps = 2 * k
	n_rot = 4 * steps
	eps_rot = eps_state / (2 * n_rot)
	b_rot = math.ceil(math.log2(1 / eps_rot))
	parts = cfg["T_walk_parts"]
	# Controlled SELECT: one AND for the kinetic control, one for the phase of the potential branch.
	walk = parts["select_kinetic"] + parts["potential_oracle"] + 8 + tcount(MultiControlZ((1,) * (4 + m))) + 4 * t_rot(eps_rot)
	prep = 4 * (2 * 2 ** (2 * n) + 2 * n * b_rot)
	qubits = 2 * n + (3 + m) + 1 + 1 + cfg["qubits_parts"]["arith_work"] + b_rot
	return k, steps, prep + steps * walk, qubits, tail, scale


def grid(cfg, readout):
	n = 2 ** cfg["n"]
	if readout == "A":
		nu, x = cfg["nu"], [axis_grid(n, n * cfg["dx"] / 2, X0[a], cfg["dx"])[0] for a in (0, 1)]
	else:
		nu, x = 1.0, [(np.arange(n) - n // 2) * cfg["dx"]] * 2
	hd, _ = grid_operator(x[0], x[1], cfg["dx"], nu)
	pts = np.stack(np.meshgrid(*x, indexing="ij"), -1)
	return x, hd, np.exp(-np.einsum("...i,ij,...j->...", pts, PHI, pts) / (8 * nu)).ravel()


def check(hd, u0, alpha, k, eps_state):
	"""Check the contraction-scaled polynomial against the finite-grid semigroup."""
	ref = expm_multiply(-T * hd, u0)
	tail = bessel_tail(alpha * T, k)
	scale = 1 / (1 + tail)
	approx = scale * series(hd, u0, alpha, k)
	err = np.linalg.norm(approx - ref) / np.linalg.norm(u0)
	assert err <= 2 * tail + 1e-12, (err, 2 * tail)
	assert err <= eps_state / 2 + 1e-12, (err, eps_state)
	return approx, float(err)


def cost_a(cfg):
	x, hd, u0 = grid(cfg, "A")
	ref = expm_multiply(-T * hd, u0)
	indices = [i * len(x[1]) + j for i, j in neighbours(x[0], x[1], cfg["h"])]
	b_ref = ref[indices] / np.linalg.norm(u0)
	delta = cfg["delta_tot"]
	eps = (SHARE["sim"] + SHARE["schr"]) * delta * b_ref.min()
	k, steps, t_run, qubits, tail, scale = run(cfg, 2, eps)
	approx, err = check(hd, u0, cfg["alpha_walk"], k, eps)
	b = approx[indices] / np.linalg.norm(u0)
	ts = [next(t for t in range(1, 60) if ae_rel_error(q, t) <= SHARE["ae"] * delta) for q in b * (1 - 0.3 * delta)]
	runs = sum(2 * 2**t + 1 for t in ts)
	total = median_reps() * runs * t_run
	return {"readout": "A", "target": cfg["target"], "K": k, "steps": steps, "steps_schr": cfg["walk_steps_per_run"], "T_run": t_run,
		"runs": runs, "runs_schr": cfg["calls_per_estimate_set"], "T_total": total, "T_schr": cfg["T_total"], "ratio": cfg["T_total"] / total,
		"logical_qubits": qubits + 1, "qubits_schr": cfg["logical_qubits"], "P_succ": float(np.sum(approx**2) / np.sum(u0**2)),
		"series_error": err, "tail": tail, "polynomial_scale": scale}


def cost_b(cfg):
	_, hd, u0 = grid(cfg, "B")
	ref = expm_multiply(-T * hd, u0)
	p_ref = float(np.sum(ref**2) / np.sum(u0**2))
	eps = 0.01 * cfg["target"] * p_ref * cfg["sigma_min"] / (8 * cfg["L"] ** 2)
	k, steps, t_run, qubits, tail, scale = run(cfg, 0, eps)
	approx, err = check(hd, u0, cfg["alpha_walk"], k, eps)
	p = float(np.sum(approx**2) / np.sum(u0**2))
	runs = runs_needed(cfg["accepted_samples"], p)
	return {"readout": "B", "target": cfg["target"], "K": k, "steps": steps, "steps_schr": cfg["walk_steps_per_run"], "T_run": t_run,
		"runs": runs, "runs_schr": cfg["runs"], "T_total": runs * t_run, "T_schr": cfg["T_total"], "ratio": cfg["T_total"] / (runs * t_run),
		"logical_qubits": qubits, "qubits_schr": cfg["logical_qubits"], "P_succ": p,
		"series_error": err, "tail": tail, "polynomial_scale": scale}


if __name__ == "__main__":
	assert 2 * np.sum(ive(np.arange(order(155.0, 1e-8) + 1, 3000), 155.0)) <= 1e-8
	a = json.loads((RESULTS / "step8.json").read_text())["best"].values()
	b = json.loads((RESULTS / "step9_sampling.json").read_text())["best"].values()
	from model import riccati_P
	for c in b:
		c["sigma_min"] = 2 / float(np.linalg.eigvalsh(riccati_P(T))[-1])
	out = [cost_a(c) for c in a] + [cost_b(c) for c in b]
	for r in out:
		print(r)
	save("direct_polynomial", out)
