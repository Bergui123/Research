"""Step 7: classical simulation of the Schroedingerised (enlarged) system.

Warped phase w(s, p, x) with w(0) = g(p) u0(x); for each eigenpair (E, phi_E) of the
grid operator H_d the exact dynamics is a shift, w_E(s, p) = g(p + E s) u0_E. On a
periodic p-grid the circuit implements this shift spectrally:
  w_hat_k(T) = exp(i eta_k H_d T) w_hat_k(0),  eta_k = (k - (N_p - 1)/2) pi / L_p,
the half-integer frequency grid used by the offset-binary circuit of Step 8.
Recovery post-selects p in [p_a, p_a + 3]; ideally w(T, p) = e^{-p} u(T) there, for any
p_a >= -E_min T. P_succ = sum_{window} e^{-2p} ||u(T)||^2 / (sum g^2 ||u0||^2).

Two p-profiles: the source's g = e^{-|p|} (kink at 0) and a C-infinity profile
that equals e^{-p} for p >= 0 and e^{p} for p <= -w (smooth blend between).
Only p >= p_a >= 0 is read out, so g need only be exact there.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")

from concurrent.futures import ProcessPoolExecutor  # noqa: E402

import json  # noqa: E402

import numpy as np  # noqa: E402

from model import P_SHIFT, PHI, RESULTS, T, X0, axis_grid, circuit_dx, grid_operator, save  # noqa: E402
from step5_nu_sweep import SPANS, neighbours  # noqa: E402

NUS = (1.0, 2.0, 4.0, 8.0, 16.0)
ETAS = (1.25, 2.5, 5.0, 10.0, 20.0, 40.0)
KINDS = ("kink", "smooth_w0.5", "smooth_w1", "smooth_w2")
PAS = (P_SHIFT, 0.0, 1.0)
WINDOW = 3.0


def smooth_step(y):
	"""C-infinity step: 0 for y <= -1, 1 for y >= 0, all derivatives vanish at both ends."""
	f = lambda z: np.where(z > 0, np.exp(-1.0 / np.where(z > 0, z, 1.0)), 0.0)
	return f(y + 1) / (f(y + 1) + f(-y))


def profile(p, kind):
	if kind == "kink":
		return np.exp(-np.abs(p))
	w = float(kind.split("_w")[1])
	s = smooth_step(np.clip(p / w, -1, 0))
	# Clip the exponents: s = 0 times exp(-p) = inf would give NaN once L_p > 709.
	grow, decay = np.exp(-np.maximum(p, -w)), np.exp(np.minimum(p, 0))
	return np.where(p >= 0, np.exp(-np.maximum(p, 0)), s * grow + (1 - s) * decay)


def g_hat_fft(gp, n_pp):
	"""g_hat_k = (1/N_p) sum_l exp(-i eta_k p_l) g_l on the half-integer frequency grid."""
	c = (n_pp - 1) / 2
	k = np.arange(n_pp)
	return np.exp(1j * np.pi * (k - c)) * np.fft.fft(np.exp(2j * np.pi * c * k / n_pp) * gp) / n_pp


def p_offset(l_p, dp, p_a):
	"""Offset of the p grid that puts a point at the window start p_a: the circuit prepares
	g(p + offset) on the plain grid, so the offset costs nothing."""
	return (p_a + l_p) % dp


def in_window(p, p_a):
	return (p >= p_a - 1e-9) & (p <= p_a + WINDOW + 1e-9)


def p_grids(e_max):
	"""(N_p, L_p) pairs: safe L_p >= (E_max T + p_b + 20)/2 so wrapped modes carry e^{-20},
	plus two small domains that show the wrap-around failure."""
	l_min = (e_max * T + max(PAS) + WINDOW + 20.0) / 2
	out = []
	for eta_max in ETAS:
		n_pp = 2 ** int(np.ceil(np.log2(2 * l_min * eta_max / np.pi)))
		out.append((n_pp, np.pi * n_pp / (2 * eta_max), True))
	for l_p in (16.0, 32.0):
		out += [(2**n_p, l_p, False) for n_p in (8, 10)]
	return out


def simulate_grid(args):
	nu, n, l_min = args
	dx = circuit_dx(2**n, l_min)
	x1, _ = axis_grid(2**n, 2**n * dx / 2, X0[0], dx)
	x2, _ = axis_grid(2**n, 2**n * dx / 2, X0[1], dx)
	hd, _ = grid_operator(x1, x2, dx, nu)
	energies, vecs = np.linalg.eigh(hd.toarray())
	g = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
	u0 = np.exp(-np.einsum("...i,ij,...j->...", g, PHI, g) / (8 * nu)).ravel()
	coef = vecs.T @ u0
	u_t = vecs @ (coef * np.exp(-energies * T))
	pts, spans = [], []
	for h in SPANS:
		nb = neighbours(x1, x2, h)
		if min(min(ab) for ab in nb) >= 1 and max(max(ab) for ab in nb) <= 2**n - 2:
			spans.append((h, list(range(len(pts), len(pts) + 4))))
			pts += [a * 2**n + b for a, b in nb]
	d = coef[:, None] * vecs[pts, :].T  # (nE, n_pts): c_E phi_E(x_j)
	u_pts = u_t[pts]
	rows = []
	for n_pp, l_p, safe in p_grids(energies[-1]):
		dp, deta = 2 * l_p / n_pp, np.pi / l_p
		p0 = -l_p + dp * np.arange(n_pp)
		eta = (np.arange(n_pp) - (n_pp - 1) / 2) * deta
		bd = np.zeros((n_pp, len(pts)), complex)
		for blk in range(0, n_pp, 1024):  # sum_E exp(i eta_k E T) c_E phi_E(x_j), chunked
			bd[blk:blk + 1024] = np.exp(1j * np.outer(eta[blk:blk + 1024], energies * T)) @ d
		for kind in KINDS:
			for p_a in PAS:
				p = p0 + p_offset(l_p, dp, p_a)
				gp = profile(p, kind)
				g_hat = g_hat_fft(gp, n_pp)
				if n_pp <= 256:
					assert np.allclose(np.exp(1j * np.outer(p0, eta)) @ g_hat, gp, atol=1e-12)
				total0 = np.sum(gp**2) * np.sum(u0**2)
				win = in_window(p, p_a)
				if p_a + WINDOW >= l_p - 1:
					continue
				w = (np.exp(1j * np.outer(p0[win], eta)) * g_hat[None, :]) @ bd  # w(T, p_l, x_j)
				b2 = np.sum(np.abs(w) ** 2, axis=0) / total0
				b2_ideal = np.sum(np.exp(-2 * p[win])) * u_pts**2 / total0
				fb = {}
				for h, idx in spans:
					lb, li = 0.5 * np.log(b2[idx]), np.log(np.abs(u_pts[idx]))
					v_q = -2 * nu * np.array([lb[0] - lb[1], lb[2] - lb[3]]) / (2 * h * dx)
					v_g = -2 * nu * np.array([li[0] - li[1], li[2] - li[3]]) / (2 * h * dx)
					fb[h] = float(np.linalg.norm(v_q - v_g) / np.linalg.norm(v_g))
				rows.append({"nu": nu, "n": n, "n_p": int(np.log2(n_pp)), "L_p": l_p, "safe_Lp": safe,
					"kind": kind, "p_a": p_a, "eta_max": float(eta.max()),
					"amp_err": float(np.max(np.abs(np.sqrt(b2 / b2_ideal) - 1))), "fb_err": fb,
					"P_succ": float(np.sum(np.exp(-2 * p[win])) * np.sum(u_t**2) / total0)})
	assert energies[0] * T >= -P_SHIFT  # the shifted window needs c <= E_min
	weights = np.abs(coef) ** 2
	spectrum = {"nu": nu, "n": n, "E_min": float(energies[0]), "E_max": float(energies[-1]),
		"weight_above": {str(e): float(weights[energies > e].sum() / weights.sum())
			for e in (5, 10, 20, 40, 80)}}
	return rows, spectrum


def main():
	domain = json.loads((RESULTS / "step5.json").read_text())["domain"]
	jobs = [(nu, n, domain[str(nu)]) for nu in NUS for n in (5, 6)]
	rows, spectra = [], []
	with ProcessPoolExecutor(max_workers=4) as pool:
		for r, s in pool.map(simulate_grid, jobs):
			rows += r
			spectra.append(s)
			print(s)
	# Smallest eta_max meeting each amplitude tolerance, per profile (the walk cost scales with it).
	summary = []
	for tol in (1e-3, 1e-4, 1e-5, 1e-6):
		for kind in KINDS:
			for nu in NUS:
				ok = [r for r in rows if r["kind"] == kind and r["nu"] == nu and r["n"] == 6
					and r["safe_Lp"] and r["amp_err"] <= tol]
				if ok:
					b = min(ok, key=lambda r: (r["eta_max"] / np.sqrt(r["P_succ"]), r["n_p"]))
					summary.append({"tol": tol, **{k: b[k] for k in
						("kind", "nu", "n_p", "L_p", "p_a", "eta_max", "amp_err", "P_succ")}})
					print(summary[-1])
	save("step7", {"rows": rows, "spectra": spectra, "summary": summary})


if __name__ == "__main__":
	main()
