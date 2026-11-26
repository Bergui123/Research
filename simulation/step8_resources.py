"""Steps 8-9: logical-qubit and T-gate estimate for the full pipeline.

Pipeline per coherent run U:
  prepare g(p) (x) psi(x,0) -> QFT_p -> prod_j C_j[exp(+-i tau_j H_d)] -> QFT_p^dag.
The offset-binary eta register gives eta = sum_j s_j 2^{j-1} d_eta (s_j = +-1), so segment j
evolves for tau_j = 2^{j-1} d_eta T with a bit-controlled sign (Clifford Z per SELECT call).
The walk encodes H_d - (alpha_H/2) I with normalisation alpha = alpha_H/2 = 4 nu/dx^2 + V'/2: the
identity part of the LCU is a translation in p, i.e. a relabelled window, and costs no gates.
Each segment is a GQSP polynomial in the qubitization walk (Berry et al. 2024, Thm 2, Eqs. 9-11):
K_j + 1 directionally controlled and 2 controlled walk steps, K_j >= alpha tau_j even. With Jacobi-Anger
tail t the |+> -> |+> block errs by at most 2t (rescaling included) and leaks at most sqrt(4t) to |->;
leak and return over n_p segments add at most 2 n_p^2 t, so t = eps_state / (4 n_p (n_p + 1)) suffices.

Readout A (this file): amplitude estimation of b_j = sqrt(P_succ p_j) at X +- h dx e_a, with
M = 2^t Grover iterates and a semiclassical QFT (one recycled control qubit).
Readout B (step9_sampling.py): sampling the recovered register.

The grid error charged to the budget is the larger of the cancellation-free leading estimate
(grid_error.py) and a numerical error-bound estimate (certify_grid.py). These values are not
validated upper enclosures; the resulting accuracy and confidence claims are conditional.
SELECT is charged at the AND count of the compiled circuit (select_circuit.py); 1 Toffoli = 1 AND = 4 T.
"""

import functools
import json
import math
import os

import numpy as np
from qualtran.bloqs.mcmt import MultiControlZ
from qualtran.resource_counting import QECGatesCost, get_cost_value
from scipy.special import jv
from scipy.stats import binom

from chebyshev import p_grid
from grid_error import predicted_error
from model import C, LAM, RESULTS, RHO, T, X0, axis_grid, reference_feedback, save
from select_circuit import count as select_count
from step7_schrodinger import in_window, profile

TARGETS = tuple(float(t) for t in os.environ.get("EXEC_TARGETS", "0.05,0.01,0.001").split(","))
SHARE = {"ae": 0.7, "schr": 0.1, "sim": 0.1, "pot": 0.1}  # of the amplitude-error budget
P_AE = 8 / np.pi**2  # Brassard et al. single-run success probability
CONF = 0.95  # joint confidence for the four amplitude estimates
W10 = round(10 * RHO)  # 5 x^T C x / dx^2 = 5 z1^2 + 5 z2^2 + W10 z1 z2 needs an integer weight
assert abs(W10 - 10 * RHO) < 1e-12


@functools.lru_cache(maxsize=None)
def tcount(bloq):
	return get_cost_value(bloq, QECGatesCost()).total_t_count()


def t_rot(eps):
	"""Ross-Selinger z-rotation synthesis, leading term 3 log2(1/eps)."""
	return 3 * math.log2(1 / eps) + 2


@functools.lru_cache(maxsize=None)
def jacobi_anger_order(tau, eps):
	"""Smallest even K >= tau (to the search step) with 2 sum_{k>K} |J_k(tau)| <= eps."""
	k = max(2, 2 * math.ceil(tau / 2))
	while True:
		tail = 2 * np.sum(np.abs(jv(np.arange(k + 1, k + 400), tau)))
		if tail <= eps:
			return k
		k += 2 * max(1, int(tau ** (1 / 3) / 4))


def median_reps():
	"""Smallest odd r with P(median of r estimates correct)^4 >= CONF."""
	r = 1
	while (1 - binom.sf((r - 1) // 2, r, 1 - P_AE)) ** 4 < CONF:
		r += 2
	return r


def ae_rel_error(b, t):
	"""Relative error bound on b from M = 2^t iterates: pi/(M b) + pi^2/(2 M^2 b^2)."""
	y = math.pi / (2.0**t * b)
	return y + y**2 / 2


def amplitude_budget(log_budget):
	"""Relative error per node whose worst log error is at most log_budget."""
	return -math.expm1(-log_budget)


def ae_pair(b_pair, budget):
	"""Powers of two (t+, t-) of least total cost with err(t+) + err(t-) <= budget."""
	best = None
	for t0 in range(1, 60):
		rest = budget - ae_rel_error(b_pair[0], t0)
		if rest <= 0:
			continue
		t1 = next(t for t in range(1, 60) if ae_rel_error(b_pair[1], t) <= rest)
		cost = 2**t0 + 2**t1
		if best is None or cost < best[0]:
			best = (cost, t0, t1)
		if t1 == 1:
			break
	return best[1:]


def potential_scale(nu, x1, x2, dx):
	"""Comparator scale V' (alpha_pot) and the integer bits of 5 x^T C x / dx^2."""
	g1, g2 = np.meshgrid(x1 / dx, x2 / dx, indexing="ij")
	z1, z2 = np.rint(4 * g1), np.rint(4 * g2)
	if np.all(np.abs(4 * g1 - z1) < 1e-9) and np.all(np.abs(4 * g2 - z2) < 1e-9):
		# The circuit uses quarter-integer coordinates and an integer 10*rho.
		# Compute its polynomial as integer arithmetic before scaling by 1/16;
		# otherwise roundoff just below a truncation boundary changes q by one.
		z1, z2 = z1.astype(np.int64), z2.astype(np.int64)
		q5 = (5 * z1**2 + W10 * z1 * z2 + 5 * z2**2) / 16.0
	else:
		# General coordinates are useful for scale diagnostics; costed grids
		# follow the exact quarter-integer branch above.
		q5 = 5 * (C[0, 0] * g1**2 + 2 * C[0, 1] * g1 * g2 + C[1, 1] * g2**2)
	# frexp(x)[1] is the least k with 2**k > x, including powers of two.
	k = math.frexp(float(q5.max()))[1]
	return LAM * dx**2 / (4 * nu) * 2**k / 5, k, q5


def auxiliary_qubits(alpha_h, eta_max, p_a):
	"""n_p from the safe rule 2 L_p >= alpha_H T + p_a + 3 + 20."""
	return math.ceil(math.log2((alpha_h * T + p_a + 3.0 + 20.0) * eta_max / math.pi))


def success_probability(norm_ratio, n_p, eta_max, p_a, kind):
	"""P_succ on the actual auxiliary grid: window sum of e^{-2p} over the norm of g."""
	p, _, _ = p_grid(n_p, eta_max, p_a)
	return float(np.sum(np.exp(-2 * p[in_window(p, p_a)])) / np.sum(profile(p, kind) ** 2) * norm_ratio)


def coherent_run(n, frac, alpha_kin, v_scale, m, eta_max, p_a, eps_state):
	"""T-count and qubits of one coherent run (preparation, evolution, p-window flag)."""
	b = n + frac + 1
	alpha_h = alpha_kin + v_scale  # spectrum of H_d in [0, alpha_H]
	alpha = alpha_h / 2  # walk normalisation
	n_p = auxiliary_qubits(alpha_h, eta_max, p_a)
	d_eta = 2 * eta_max / (2**n_p - 1)
	orders = [jacobi_anger_order(alpha * d_eta * 2 ** (j - 1) * T, eps_state / (4 * n_p * (n_p + 1))) for j in range(n_p)]
	q_run = sum(k + 3 for k in orders)
	n_rot = 4 * q_run + 2 * n_p
	eps_rot = eps_state / (2 * n_rot)
	b_rot = math.ceil(math.log2(1 / eps_rot))

	ands, kin_ands, work = select_count(n, frac, m, W10)  # compiled circuit, 4 T per AND
	sel_kin, oracle_b = 4 * kin_ands, 4 * (ands - kin_ands)
	refl = 2 * tcount(MultiControlZ((1,) * (4 + m)))
	rot = 4 * t_rot(eps_rot)
	walk = sel_kin + oracle_b + refl + rot
	prep = 4 * (2 * 2 ** (2 * n) + 2 * n * b_rot) + 4 * (2 * 2**n_p + n_p * b_rot)
	qft = 2 * 4 * n_p**2
	t_run = prep + qft + q_run * walk + 2 * n_p * t_rot(eps_rot) + 4 * 2 * n_p
	# The workspace is the peak of the compiled SELECT; reflections and state preparation reuse it.
	assert work >= b_rot + 2 * n
	parts = {"x": 2 * n, "p": n_p, "prep_lcu": 3 + m, "gqsp": 1, "arith_work": work, "phase_grad": b_rot}
	return {"n_p": n_p, "walk_steps_per_run": q_run, "walk_steps_std_qsp": sum(2 * k for k in orders),
		"T_walk_step": walk, "T_walk_parts": {"select_kinetic": sel_kin, "potential_oracle": oracle_b,
			"reflections": refl, "rotations": rot}, "T_prep": prep, "T_qft": qft, "T_run": t_run,
		"eps_rot": eps_rot, "alpha_H": alpha_h, "alpha_walk": alpha, "alpha_kin": alpha_kin, "alpha_pot": v_scale, "b": b, "m": m,
		"qubits_parts": parts, "qubits_run": sum(parts.values())}


def run_cost(cfg, target):
	"""Readout A: total T-count and logical qubits, or None if the budget cannot be met."""
	nu, n, dx, h, eps_grid, pc = cfg["nu"], cfg["n"], cfg["dx"], cfg["h"], cfg["eps_grid"], cfg["pgrid"]
	v_norm = np.linalg.norm(reference_feedback()[0])
	log_budget = (target - eps_grid) * v_norm * h * dx / (2 * math.sqrt(2) * nu)
	delta_tot = amplitude_budget(log_budget)
	if delta_tot <= 0 or pc["amp_err_env"] > SHARE["schr"] * delta_tot:
		return None
	x1, _ = axis_grid(2**n, cfg["L"], X0[0], dx)
	x2, _ = axis_grid(2**n, cfg["L"], X0[1], dx)
	v_scale, k_bits, _ = potential_scale(nu, x1, x2, dx)
	n_p = auxiliary_qubits(8 * nu / dx**2 + v_scale, pc["eta_max"], pc["p_a"])
	p_succ = success_probability(cfg["norm_ratio"], n_p, pc["eta_max"], pc["p_a"], pc["kind"])
	b_amp = np.sqrt(p_succ * np.array(cfg["p"]))
	# dx = 4/k: offsets are quarter-integers, so 5Q is exact with k_bits + 4 bits (2 fractional
	# bits per coordinate). By positivity (Feynman-Kac) rounding V by <= V' 2^-m changes every
	# amplitude by a relative factor <= exp(T V' 2^-m) - 1.
	m = min(math.ceil(math.log2(T * v_scale / math.log1p(SHARE["pot"] * delta_tot))), k_bits + 4)
	run = coherent_run(n, 2, 8 * nu / dx**2, v_scale, m, pc["eta_max"], pc["p_a"],
		SHARE["sim"] * delta_tot * b_amp.min())
	# Estimation is relative to the actual prepared amplitude, which may be
	# smaller than its reference by the three deterministic shares together.
	b_min = b_amp * (1 - (1 - SHARE["ae"]) * delta_tot)
	ts = [next(t for t in range(1, 60) if ae_rel_error(b, t) <= SHARE["ae"] * delta_tot)
		for b in b_min]
	calls = sum(2 * 2**t + 1 for t in ts)
	qubits = run["qubits_run"] + 1  # recycled control qubit of the semiclassical QFT
	marks = 4 * (2 * run["n_p"] + 2 * n + 3 + m + 1 + run["n_p"] + 2 * n + 3 + m + 1)  # S_chi and S_0
	r = median_reps()
	single = calls * run["T_run"] + sum(2**t for t in ts) * marks
	return {**run, "target": target, "nu": nu, "n": n, "L": cfg["L"], "dx": dx, "h": h, "eps_grid": eps_grid,
		"eps_grid_measured": cfg["eps_grid_measured"], "kind": pc["kind"], "eta_max": pc["eta_max"],
		"p_a": pc["p_a"], "P_succ": p_succ, "amp_err_schr": pc["amp_err_env"], "delta_tot": delta_tot,
		"b_amp": b_amp.tolist(), "AE_log2_iterations": ts, "AE_iterations": [2**t for t in ts],
		"calls_per_estimate_set": calls, "median_reps": r, "T_total": r * single, "T_total_no_boost": single,
		"logical_qubits": qubits}


def block_encoding_check(n=2, m=2, nu=2.0, dx=0.5):
	"""Circuit-level matrix check of the LCU block encoding and its qubitization walk.

	Registers: ancilla (t, k1 k0, r[m]) (x) system (x1[n], x2[n]). Builds PREP and the Hermitian
	SELECT (direction bit flipped after each shift) as explicit unitaries, then checks SELECT^2 = I,
	<0|PREP^dag SELECT PREP|0> = (K + V_m - alpha_H/2)/alpha with V_m the m-bit comparator potential,
	|V_m - V| <= V' 2^-m, and the walk eigenphases +-arccos((E - alpha)/alpha).
	"""
	from model import grid_operator
	npts = 2**n
	x1, _ = axis_grid(npts, npts * dx / 2, X0[0], dx)
	x2, _ = axis_grid(npts, npts * dx / 2, X0[1], dx)
	h_exact, pot = grid_operator(x1, x2, dx, nu)
	v_scale, k, q5_grid = potential_scale(nu, x1, x2, dx)
	q5 = q5_grid.ravel()
	q_top = np.floor(q5 * 2**m / 2**k)  # comparator compares r < top m bits of 5Q
	a_kin = 4 * nu / dx**2
	alpha = a_kin + v_scale / 2
	dim = npts * npts
	eye_n = np.eye(npts)
	shift = np.roll(eye_n, 1, axis=0)  # S|j> = |j+1 mod N>
	shifts = {(0, 0): np.kron(shift, eye_n), (0, 1): np.kron(shift.T, eye_n),
		(1, 0): np.kron(eye_n, shift), (1, 1): np.kron(eye_n, shift.T)}
	n_anc = 1 + 2 + m
	sel = np.zeros((2**n_anc * dim, 2**n_anc * dim))
	for anc in range(2**n_anc):
		t, k1, k0, r = anc >> (m + 2) & 1, anc >> (m + 1) & 1, anc >> m & 1, anc & (2**m - 1)
		if t == 0:
			out = anc ^ (1 << m)  # flip direction bit: makes SELECT Hermitian
			u = -shifts[(k1, k0)]
		else:
			out, u = anc, -np.diag(np.where(r < q_top, -1.0, 1.0))
		sel[out * dim:(out + 1) * dim, anc * dim:(anc + 1) * dim] = u
	ry = np.array([[math.sqrt(a_kin / alpha), -math.sqrt(v_scale / 2 / alpha)],
		[math.sqrt(v_scale / 2 / alpha), math.sqrt(a_kin / alpha)]])
	had = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
	prep = ry
	for _ in range(2 + m):
		prep = np.kron(prep, had)
	prep_full = np.kron(prep, np.eye(dim))
	assert np.allclose(sel @ sel, np.eye(len(sel)))
	block = (prep_full.T @ sel @ prep_full)[:dim, :dim]
	v_m = v_scale * q_top / 2**m
	h_m = h_exact.toarray() - np.diag(pot.ravel()) + np.diag(v_m)
	assert np.allclose(block * alpha, h_m - alpha * np.eye(dim), atol=1e-10)
	assert np.max(np.abs(v_m - pot.ravel())) <= v_scale * 2**-m + 1e-12
	refl = np.eye(len(sel))
	refl[dim:, dim:] *= -1  # 2|0><0|_anc - I on the prepared basis
	walk = prep_full @ refl @ prep_full.T @ sel
	phases = np.angle(np.linalg.eigvals(walk))
	for e in np.linalg.eigvalsh(h_m):
		assert np.min(np.abs(np.abs(phases) - math.acos((e - alpha) / alpha))) < 1e-7
	return {"n": n, "m": m, "alpha": alpha, "max_pot_rounding": float(np.max(np.abs(v_m - pot.ravel()))),
		"bound": v_scale * 2**-m, "dim": len(sel)}


def candidates(target):
	s5 = json.loads((RESULTS / "step5.json").read_text())
	s7 = [r for r in json.loads((RESULTS / "step7.json").read_text())["rows"] if r["safe_Lp"]]
	for r in s7:  # envelope: never trust an amplitude error that is worse again at larger eta_max
		r["amp_err_env"] = max(q["amp_err"] for q in s7 if (q["nu"], q["n"], q["kind"], q["p_a"])
			== (r["nu"], r["n"], r["kind"], r["p_a"]) and q["eta_max"] >= r["eta_max"])
	grid_estimates = {(c["nu"], c["n"], c["h"]): c["numerical_error_estimate"]
		for c in json.loads((RESULTS / "grid_estimates.json").read_text())["rows"]}
	for r in s5["grid"]:
		if r["nu"] < 1.0:
			continue  # readout probabilities below 1e-7
		for sp in r["spans"]:
			estimate = grid_estimates.get((r["nu"], r["n"], sp["h"]))
			if estimate is None:
				continue  # no numerical error-bound estimate for this grid and span
			eps = max(predicted_error(r["nu"], r["dx"], sp["h"])[1], estimate)
			if eps > target / 2:
				continue
			pgrids = sorted((p for p in s7 if p["nu"] == r["nu"] and p["n"] == min(r["n"], 6)),
				key=lambda p: p["eta_max"] / math.sqrt(p["P_succ"]))
			v_norm = np.linalg.norm(reference_feedback()[0])
			delta = amplitude_budget((target - eps) * v_norm * sp["h"] * r["dx"] / (2 * math.sqrt(2) * r["nu"]))
			ok = [p for p in pgrids if p["amp_err_env"] <= SHARE["schr"] * delta]
			for pc in ok[:3]:  # cost is proportional to eta_max / sqrt(P_succ) to leading order
				yield {"nu": r["nu"], "n": r["n"], "L": r["L"], "dx": r["dx"], "h": sp["h"], "eps_grid": eps,
					"eps_grid_measured": sp["rel_err_v"], "p": sp["p"], "norm_ratio": r["norm_ratio"], "pgrid": pc}


def main():
	best, per_nu = {}, {}
	for target in TARGETS:
		for cfg in candidates(target):
			res = run_cost(cfg, target)
			if res is None:
				continue
			if target not in best or res["T_total"] < best[target]["T_total"]:
				best[target] = res
			k2 = (target, cfg["nu"])
			if k2 not in per_nu or res["T_total"] < per_nu[k2]["T_total"]:
				per_nu[k2] = res
	for target, r in best.items():
		print(f"eps={target:g}: T={r['T_total']:.3e} (single set {r['T_total_no_boost']:.3e}) "
			f"qubits={r['logical_qubits']} nu={r['nu']} N={2**r['n']} h={r['h']} {r['kind']} eta={r['eta_max']:.1f} "
			f"n_p={r['n_p']} alpha={r['alpha_H']:.0f} steps={r['walk_steps_per_run']} T/step={r['T_walk_step']:.0f} "
			f"t={r['AE_log2_iterations']} r={r['median_reps']} eps_grid={r['eps_grid']:.1e} (measured {r['eps_grid_measured']:.1e})")
	save("step8", {"block_encoding_check": block_encoding_check() if RHO == 0.4 else None,
		"best": {str(k): v for k, v in best.items()},
		"per_nu": [v for _, v in sorted(per_nu.items())], "shares": SHARE, "confidence": CONF})


if __name__ == "__main__":
	# Self-checks of the cost helpers.
	assert jacobi_anger_order(0.0, 1e-9) == 2
	assert 2 * np.sum(np.abs(jv(np.arange(jacobi_anger_order(500.0, 1e-10) + 1, 1400), 500.0))) <= 1e-10
	t0, t1 = ae_pair(np.array([1e-3, 4e-3]), 2e-2)
	assert ae_rel_error(1e-3, t0) + ae_rel_error(4e-3, t1) <= 2e-2
	assert median_reps() == 11
	main()
