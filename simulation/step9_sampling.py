"""Step 9, readout B: sample the recovered register instead of estimating amplitudes.

For the quadratic model |psi(x,T)|^2 is a centred Gaussian with covariance 2 nu P_T^{-1}, so the
second moments of accepted position samples give P_T and hence v*(x) = P_T x / 2 for EVERY x.
By scale invariance nu = 1 without loss of generality and X need not lie on the grid.
Error metric: ||(P_hat - P_T) P_T^{-1}||_2, which bounds the relative-to-P_T-x feedback error
uniformly; the systematic parts (grid, then rounding + Schroedingerisation) are added in norm
so that they cannot cancel. The run count is chosen so that enough samples are accepted with
probability 0.99, and the coherent-simulation error is charged to the budget.
"""

import math

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta, binom

from chebyshev import recover
from model import P_SHIFT, PHI, T, X0, grid_operator, riccati_P, save
from step8_resources import TARGETS, coherent_run, potential_scale

RNG = np.random.default_rng(20261124)
P = riccati_P(T)
NU = 1.0
Q_COV, Q_ACC = 0.96, 0.99  # covariance quantile and acceptance; the product exceeds 0.95
COHERENT = 0.02  # share of eps charged for GQSP truncation and rotation synthesis


def grid(n, half_width, dx=None):
	dx = 2 * half_width / 2**n if dx is None else dx
	x = (np.arange(round(2 * half_width / dx)) - round(half_width / dx)) * dx
	g = np.stack(np.meshgrid(x, x, indexing="ij"), axis=-1)
	hd, pot = grid_operator(x, x, dx, NU)
	u0 = np.exp(-np.einsum("...i,ij,...j->...", g, PHI, g) / (8 * NU)).ravel()
	return x, dx, hd, pot, u0, g.reshape(-1, 2)


def p_from(prob, pts):
	prob = prob / prob.sum()
	return 2 * NU * np.linalg.inv(np.einsum("k,ki,kj->ij", prob, pts, pts))


def rel(p_hat, p_ref):
	return float(np.linalg.norm((p_hat - p_ref) @ np.linalg.inv(P), 2))


def truncation_domain(eps):
	"""Smallest L whose P changes by < 0.02 eps when the box grows 1.6x at fixed dx = 0.1."""
	def p_l(half_width):
		_, _, hd, _, u0, pts = grid(None, half_width, 0.1)
		return p_from(expm_multiply(-T * hd, u0) ** 2, pts)
	for half_width in np.arange(4.0, 9.01, 0.5):
		if rel(p_l(half_width), p_l(1.6 * half_width)) < 0.02 * eps:
			return float(half_width)
	raise RuntimeError("no domain found")


def quantile_constant(trials=40000, n=4000):
	"""kappa with Q_COV quantile of ||(P_hat - P) P^-1||_2 = kappa / sqrt(N) (Gaussian sampling, large N)."""
	w, u = np.linalg.eigh(P)
	root, iroot = u * np.sqrt(w) @ u.T, u / np.sqrt(w) @ u.T
	g = RNG.standard_normal((trials, 3))
	e = np.empty((trials, 2, 2))  # GOE limit of the whitened sample covariance
	e[:, 0, 0], e[:, 1, 1] = math.sqrt(2) * g[:, 0], math.sqrt(2) * g[:, 1]
	e[:, 0, 1] = e[:, 1, 0] = g[:, 2]
	return float(np.quantile(np.linalg.norm(root @ e @ iroot, 2, axis=(1, 2)), Q_COV))


def runs_needed(shots, p):
	"""Smallest number of runs with at least `shots` accepted with probability Q_ACC."""
	r = math.ceil(shots / p)
	lo, hi = r, r + math.ceil(4 * math.sqrt(shots * (1 - p)) / p) + 10
	assert binom.sf(shots - 1, hi, p) >= Q_ACC
	while lo < hi:
		mid = (lo + hi) // 2
		lo, hi = (lo, mid) if binom.sf(shots - 1, mid, p) >= Q_ACC else (mid + 1, hi)
	return lo


def configuration(eps, n, half_width, kind, eta_max, kappa):
	x, dx, hd, pot, u0, pts = grid(n, half_width)
	v_scale, k_bits, q5 = potential_scale(NU, x, x, dx)
	m = min(math.ceil(math.log2(T * v_scale / (0.02 * eps))), k_bits)
	v_m = v_scale * np.floor(q5 * 2.0 ** (m - k_bits)) / 2**m
	hd_m = hd - diags(pot.ravel()) + diags(v_m.ravel())  # what the comparator actually encodes
	u_t = expm_multiply(-T * hd, u0)
	p_pde = p_from(u_t**2, pts)
	rho = rel(p_pde, P)
	if rho > eps / 2:
		return None
	alpha_kin = 8 * NU / dx**2
	alpha = alpha_kin + v_scale
	n_p = math.ceil(math.log2((alpha * T + P_SHIFT + 23.0) * eta_max / math.pi))
	w, gp, win = recover(hd_m, u0, alpha, n_p, eta_max, kind, P_SHIFT)
	prob = np.sum(w**2, axis=0)
	p_succ = float(prob.sum() / (np.sum(gp**2) * np.sum(u0**2)))
	rho += rel(p_from(prob, pts), p_pde)
	if rho > eps / 2:
		return None
	sigma_min = 2 * NU / np.linalg.eigvalsh(P)[-1]
	# A state error e shifts the second-moment matrix by at most 8 L^2 e / sqrt(P) to first order, and the
	# metric by that over sigma_min: at most COHERENT eps / 2 for this e; COHERENT eps is charged.
	eps_state = 0.5 * COHERENT * eps * p_succ * sigma_min / (8 * half_width**2)
	run = coherent_run(n, 0, alpha_kin, v_scale, m, eta_max, P_SHIFT, eps_state)
	assert run["n_p"] == n_p
	shots = math.ceil((kappa / (eps * (1 - COHERENT) - rho)) ** 2)
	runs = runs_needed(shots, p_succ)
	return {**run, "target": eps, "n": n, "L": half_width, "dx": dx, "kind": kind, "eta_max": eta_max,
		"P_succ": p_succ, "rho_systematic": rho, "accepted_samples": shots, "runs": runs,
		"T_total": runs * run["T_run"], "logical_qubits": run["qubits_run"], "_prob": prob, "_pts": pts}


def validate(cfg, trials=2000):
	"""Draw the costed number of runs, reject by P_succ, and compare the estimate with Riccati."""
	prob = cfg["_prob"] / cfg["_prob"].sum()
	v = 0.5 * P @ X0
	uniform, at_x = [], []
	for _ in range(trials):
		counts = RNG.multinomial(RNG.binomial(cfg["runs"], cfg["P_succ"]), prob)
		p_hat = p_from(counts.astype(float), cfg["_pts"])
		uniform.append(rel(p_hat, P))
		at_x.append(np.linalg.norm(0.5 * p_hat @ X0 - v) / np.linalg.norm(v))
	ok = int(np.sum(np.array(uniform) <= cfg["target"]))
	return {"frac_uniform_within": ok / trials, "lower_bound": float(beta.ppf(0.001, ok, trials - ok + 1)),
		"q95_uniform": float(np.quantile(uniform, 0.95)), "frac_at_X_within": float(np.mean(np.array(at_x) <= cfg["target"])),
		"median_at_X": float(np.median(at_x)), "errors_at_X": at_x}


def main():
	kappa = quantile_constant()
	out = {"kappa": kappa, "best": {}}
	for eps in TARGETS:
		half_width = truncation_domain(eps)
		best = None
		for n in range(3, 8):
			found = [c for kind in ("kink", "smooth_w0.5", "smooth_w1", "smooth_w2") for eta in (5.0, 10.0, 20.0, 40.0)
				if (c := configuration(eps, n, half_width, kind, eta, kappa))]
			for c in found:
				if best is None or c["T_total"] < best["T_total"]:
					best = c
			if found and n >= best["n"] + 1:
				break  # a finer grid is already more expensive
		best["validation"] = validate(best)
		v = best["validation"]
		print(f"eps={eps:g}: T={best['T_total']:.3e} qubits={best['logical_qubits']} N={2**best['n']} L={half_width} "
			f"{best['kind']} eta={best['eta_max']:g} n_p={best['n_p']} alpha={best['alpha_H']:.0f} steps={best['walk_steps_per_run']} "
			f"T/step={best['T_walk_step']:.0f} T/run={best['T_run']:.2e} P={best['P_succ']:.3f} rho={best['rho_systematic']:.2e} "
			f"shots={best['accepted_samples']} runs={best['runs']} | uniform within {v['frac_uniform_within']:.3f}, at X {v['frac_at_X_within']:.3f}")
		assert v["lower_bound"] >= 0.95  # one-sided Clopper-Pearson at level 0.001
		out["best"][str(eps)] = {k: v for k, v in best.items() if not k.startswith("_")}
	save("step9_sampling", out)


if __name__ == "__main__":
	main()
