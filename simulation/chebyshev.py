"""Exact Schroedingerised recovery for any grid, without diagonalising H_d.

w(T, p_l, .) = f_l(H_d) psi(0) with f_l(E) = sum_k g_hat_k exp(i eta_k (p_l + E T)), the
trigonometric interpolant of g at p_l + E T. f_l is band limited (|eta| <= eta_max), so its
Chebyshev series on [0, alpha] converges super-exponentially past degree alpha eta_max T / 2.
"""

import numpy as np
from scipy.fft import dct

from model import T
from step7_schrodinger import g_hat_fft, in_window, p_offset, profile


def p_grid(n_p, eta_max, p_a):
	"""p points (one of them at the window start p_a), their unshifted positions, and eta."""
	n_pp = 2**n_p
	d_eta = 2 * eta_max / (n_pp - 1)
	l_p = np.pi / d_eta
	p0 = -l_p + 2 * l_p / n_pp * np.arange(n_pp)
	return p0 + p_offset(l_p, 2 * l_p / n_pp, p_a), p0, (np.arange(n_pp) - (n_pp - 1) / 2) * d_eta


def window_coefficients(alpha, n_p, eta_max, kind, p_a):
	"""Chebyshev coefficients c[k, l] of f_l on [0, alpha] for every window point p_l."""
	p, p0, eta = p_grid(n_p, eta_max, p_a)
	gp = profile(p, kind)
	g_hat = g_hat_fft(gp, len(p))
	win = in_window(p, p_a)
	kernel = np.exp(1j * np.outer(eta, p0[win])) * g_hat[:, None]
	deg = int(1.05 * alpha * eta_max * T / 2 + 200)
	nodes = alpha / 2 * (1 + np.cos(np.pi * (np.arange(deg) + 0.5) / deg))
	f = np.empty((deg, win.sum()))
	for i in range(0, deg, 512):
		f[i:i + 512] = (np.exp(1j * np.outer(nodes[i:i + 512] * T, eta)) @ kernel).real
	c = dct(f, type=2, axis=0) / deg
	c[0] /= 2
	return c, gp, win


def recover(hd, u0, alpha, n_p, eta_max, kind, p_a, pts=None):
	"""w(T, p_l, x) at grid indices `pts` (all if None), with g and the window mask. Real part:
	g is real and the half-integer grid is symmetric, so the imaginary part is rounding only."""
	c, gp, win = window_coefficients(alpha, n_p, eta_max, kind, p_a)
	sel = slice(None) if pts is None else pts
	prev, cur = u0, (2 / alpha) * (hd @ u0) - u0
	w = np.outer(c[0], prev[sel]) + np.outer(c[1], cur[sel])
	for k in range(2, len(c)):
		prev, cur = cur, 2 * ((2 / alpha) * (hd @ cur) - cur) - prev
		w += np.outer(c[k], cur[sel])
	return w, gp, win


if __name__ == "__main__":  # self-check against dense diagonalisation
	from model import PHI, X0, axis_grid, circuit_dx, grid_operator
	nu, dx = 2.0, circuit_dx(32, 8.225)
	x1, _ = axis_grid(32, 16 * dx, X0[0], dx)
	x2, _ = axis_grid(32, 16 * dx, X0[1], dx)
	hd, _ = grid_operator(x1, x2, dx, nu)
	g = np.stack(np.meshgrid(x1, x2, indexing="ij"), axis=-1)
	u0 = np.exp(-np.einsum("...i,ij,...j->...", g, PHI, g) / (8 * nu)).ravel()
	e, v = np.linalg.eigh(hd.toarray())
	p, p0, eta = p_grid(9, 10.0, -0.9)
	gp = profile(p, "smooth_w0.5")
	win = in_window(p, -0.9)
	ref = (np.exp(1j * np.outer(p0[win], eta)) * g_hat_fft(gp, 512)) @ np.exp(1j * np.outer(eta, e * T)) @ ((v.T @ u0)[:, None] * v.T)
	w, _, _ = recover(hd, u0, 91.0, 9, 10.0, "smooth_w0.5", -0.9)
	print("max deviation from dense:", np.max(np.abs(w - ref)), " imag part:", np.max(np.abs(ref.imag)))
	assert np.max(np.abs(w - ref)) < 1e-10
