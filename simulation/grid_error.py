"""Leading-order grid error of the recovered feedback (derived for the paper).

With psi = exp(phi) and D2 = d^2 + (dx^2/12) d^4 + O(dx^4), phi_d = phi + dx^2 chi where
  chi_s = nu Lap chi - (P_s x / 2).grad chi
          + (nu/12) sum_a [phi_a^4 + 6 phi_a^2 phi_aa + 3 phi_aa^2],   chi(0) = 0,
  phi_a = -(P_s x)_a / (4 nu),  phi_aa = -(P_s)_aa / (4 nu).
chi is an even quartic polynomial, so its coefficients solve a small linear ODE. The quartic
and quadratic parts push the feedback in opposite directions and cancel near some nu; the
cancellation-free bound adds their magnitudes instead, and is what Step 8 budgets.
"""

import functools

import numpy as np
from scipy.integrate import solve_ivp

from model import PHI, T, X0, reference_feedback, riccati_P

MONO = [(4, 0), (3, 1), (2, 2), (1, 3), (0, 4), (2, 0), (1, 1), (0, 2)]  # x1^i x2^j
IDX = {m: k for k, m in enumerate(MONO)}


def _add(out, powers, value):
	if powers in IDX:
		out[IDX[powers]] += value


def _rhs(s, c, nu):
	p = riccati_P(s)
	out = np.zeros(len(MONO))
	for (i, j), v in zip(MONO, c):
		_add(out, (i - 2, j), nu * i * (i - 1) * v)
		_add(out, (i, j - 2), nu * j * (j - 1) * v)
		# drift -(P x / 2) . grad
		_add(out, (i, j), -0.5 * (p[0, 0] * i + p[1, 1] * j) * v)
		_add(out, (i - 1, j + 1), -0.5 * p[0, 1] * i * v)
		_add(out, (i + 1, j - 1), -0.5 * p[0, 1] * j * v)
	for a in range(2):
		f = -p[a] / (4 * nu)  # phi_a = f . x
		faa = -p[a, a] / (4 * nu)
		quart = np.array([f[0]**4, 4 * f[0]**3 * f[1], 6 * f[0]**2 * f[1]**2, 4 * f[0] * f[1]**3, f[1]**4])
		out[:5] += nu / 12 * quart
		out[5:] += nu / 12 * 6 * faa * np.array([f[0]**2, 2 * f[0] * f[1], f[1]**2])
	return out


@functools.lru_cache(maxsize=None)
def chi_coefficients(nu):
	return solve_ivp(_rhs, (0, T), np.zeros(len(MONO)), args=(nu,), rtol=1e-10, atol=1e-14).y[:, -1]


def feedback_shift(nu, dx, h, part=slice(None)):
	"""Predicted v_grid - v_exact at X for spacing dx and readout span h (chosen monomials)."""
	c = np.zeros(len(MONO))
	c[part] = chi_coefficients(nu)[part]

	def chi(x):
		return sum(v * x[0]**i * x[1]**j for (i, j), v in zip(MONO, c))
	d = h * dx
	return np.array([-2 * nu * dx**2 * (chi(X0 + d * e) - chi(X0 - d * e)) / (2 * d) for e in np.eye(2)])


def predicted_error(nu, dx, h):
	"""(leading-order relative error, cancellation-free relative bound)."""
	q4, q2 = feedback_shift(nu, dx, h, slice(0, 5)), feedback_shift(nu, dx, h, slice(5, 8))
	ref = np.linalg.norm(reference_feedback()[0])
	return float(np.linalg.norm(q4 + q2) / ref), float(np.linalg.norm(np.abs(q4) + np.abs(q2)) / ref)


if __name__ == "__main__":
	import json
	from model import RESULTS, save
	rows = []
	for r in json.loads((RESULTS / "step5.json").read_text())["grid"]:
		for sp in r["spans"]:
			lead, cf = predicted_error(r["nu"], r["dx"], sp["h"])
			rows.append({"nu": r["nu"], "n": r["n"], "h": sp["h"], "measured": sp["rel_err_v"],
				"leading": lead, "cancellation_free": cf})
	fine = [r for r in rows if r["n"] == 8 and r["h"] == 1 and r["nu"] >= 0.5]
	for r in fine:
		print(f"nu={r['nu']:5.1f} n=8 h=1 measured={r['measured']:.4e} leading={r['leading']:.4e} cf={r['cancellation_free']:.4e}")
		# Near the sign change the error itself almost vanishes, so compare on the uncancelled scale;
		# 1e-5 is the Step 5 domain-truncation tolerance.
		assert abs(r["leading"] - r["measured"]) < 0.02 * r["cancellation_free"] + 1e-5
	for nu, n in ((4.0, 5), (2.0, 5), (16.0, 6)):
		print(nu, n, [(r["h"], f"{r['measured']:.1e}", f"{r['leading']:.1e}", f"{r['cancellation_free']:.1e}")
			for r in rows if r["nu"] == nu and r["n"] == n])
	print("measured <= cancellation-free bound:",
		np.mean([r["measured"] <= r["cancellation_free"] for r in rows if r["nu"] >= 1]))
	save("grid_error", rows)
