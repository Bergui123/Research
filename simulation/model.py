"""Shared benchmark data and exact references for Steps 3-10.

Notation follows WorkPlan.tex: H impact (here I), C covariance, lam risk
aversion, Phi finite terminal penalty, s = T - t time remaining, nu the
inventory diffusion of Step 4, A = H^{-1}.
"""

import json
import math
import os
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy import sparse

ROOT = Path(__file__).resolve().parent
RHO = float(os.environ.get("EXEC_RHO", "0.4"))
LAM = float(os.environ.get("EXEC_LAM", "1.0"))
# Sweep points write next to, never over, the benchmark results.
RESULTS = ROOT / ("results" if (RHO, LAM) == (0.4, 1.0) else f"results/rho{RHO:g}_lam{LAM:g}")
FIGURES = ROOT / "figures"

H = np.eye(2)
A = np.linalg.inv(H)
C = np.array([[1.0, RHO], [RHO, 1.0]])
T = 1.0
X0 = np.array([3.0, 1.0])
# Recovery window start: modes shift by E s with E >= E_min, so w = e^{-p} psi(T) already holds on
# p >= -E_min T. E_0 = sum sqrt(lam mu_i)/2 is the continuum ground level; 2% margin for the grid.
E0 = 0.5 * float(np.sum(np.sqrt(LAM * np.linalg.eigvalsh(C))))
SHIFT_MARGIN = 0.98
P_SHIFT = -SHIFT_MARGIN * E0 * T
PHI = np.array([[2.0, 0.3], [0.3, 1.0]])  # Step 4 finite penalty, C@PHI != PHI@C


def block_solution(s, phi=PHI, c=C, lam=LAM, a=A):
	"""Y(s), Z(s) of Y' = A Z/2, Z' = 2 lam C Y, Y(0)=I, Z(0)=Phi (exact, via expm)."""
	m = len(phi)
	gen = np.block([[np.zeros((m, m)), a / 2], [2 * lam * c, np.zeros((m, m))]])
	yz = expm(s * gen) @ np.vstack([np.eye(m), phi])
	return yz[:m], yz[m:]


def riccati_P(s, phi=PHI):
	"""P_s = Z Y^{-1}: Step 2 Riccati matrix at time remaining s."""
	y, z = block_solution(s, phi)
	p = z @ np.linalg.inv(y)
	return (p + p.T) / 2


def log_det_Y(s, phi=PHI):
	"""log det Y(s); the value correction is c(s) = 2 nu log det Y(s)."""
	return np.log(np.linalg.det(block_solution(s, phi)[0]))


def reference_feedback(x=X0, s=T, phi=PHI):
	"""Exact signed selling rate v* = P_s x / 2 (H = I) and value x^T P_s x / 2."""
	p = riccati_P(s, phi)
	return 0.5 * A @ p @ x, 0.5 * x @ p @ x


def exact_psi(points, nu, s=T, phi=PHI):
	"""Full-space Gaussian solution of the linear PDE (Step 4, Eq. stage4-gaussian)."""
	p = riccati_P(s, phi)
	quad = np.einsum("...i,ij,...j->...", points, p, points)
	return np.exp(-quad / (8 * nu) - 0.5 * log_det_Y(s, phi))


def axis_grid(n_points, half_width, anchor, dx=None):
	"""N uniform points spanning ~[-L, L) that contain `anchor` exactly.

	Placing the readout point X on the grid avoids interpolation error.
	"""
	dx = 2 * half_width / n_points if dx is None else dx
	j0 = int(round((anchor + half_width) / dx))
	return anchor + (np.arange(n_points) - j0) * dx, dx


def circuit_dx(n_points, l_min):
	"""Smallest dx = 4/k with N dx/2 >= l_min.

	With X on the grid, the offsets o_a = J_a - X_a/dx are then quarter-integers, so the
	circuit's coordinate arithmetic is exact with 2 fractional bits (Step 8).
	"""
	return 4.0 / math.floor(2 * n_points / l_min)


def second_difference(n_points, dx, periodic=True):
	"""Central second difference D2 (periodic = cyclic shifts, as the circuit uses)."""
	e = np.ones(n_points)
	d2 = sparse.diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1], format="lil")
	if periodic:
		d2[0, -1] = d2[-1, 0] = 1.0
	return (d2 / dx**2).tocsr()


def grid_operator(x1, x2, dx, nu, periodic=True, c=C, lam=LAM, a=A):
	"""Discrete H_d = -nu sum_a A_aa D2_a + lam/(4nu) x^T C x (WorkPlan Eq. disc-operator).

	State ordering: index = i1 * N2 + i2 (x1 slow, x2 fast).
	"""
	n1, n2 = len(x1), len(x2)
	kin = -nu * (a[0, 0] * sparse.kron(second_difference(n1, dx, periodic), sparse.eye(n2))
		+ a[1, 1] * sparse.kron(sparse.eye(n1), second_difference(n2, dx, periodic)))
	g1, g2 = np.meshgrid(x1, x2, indexing="ij")
	pot = lam / (4 * nu) * (c[0, 0] * g1**2 + 2 * c[0, 1] * g1 * g2 + c[1, 1] * g2**2)
	return (kin + sparse.diags(pot.ravel())).tocsr(), pot


def save(name, data):
	RESULTS.mkdir(parents=True, exist_ok=True)
	(RESULTS / f"{name}.json").write_text(json.dumps(data, indent=1, default=float))
