"""Finite-difference check of the affine gauge/translation identity in Recentre.lean."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


SEED = 8042026
NU = 0.37
A = np.array([[1.4, 0.25], [0.25, 0.9]])
Q = np.array([[0.7, 0.15], [0.15, 0.5]])
TIME_STEP = 2e-5
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "review" / "validation" / "recentre-numeric-results.json"


def potential(x: np.ndarray) -> float:
	x1, x2 = x
	return float(0.3 * (1 - np.cos(x1)) + 0.2 * (1 - np.cos(x2))
		+ 0.15 * np.sin(x1 + 2 * x2) + 0.03 * np.cos(2 * x1 - x2))


def grad_potential(x: np.ndarray) -> np.ndarray:
	x1, x2 = x
	phase = x1 + 2 * x2
	other = 2 * x1 - x2
	return np.array([0.3 * np.sin(x1) + 0.15 * np.cos(phase) - 0.06 * np.sin(other),
		0.2 * np.sin(x2) + 0.30 * np.cos(phase) + 0.03 * np.sin(other)])


def gaussian_jet(y: np.ndarray, c: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
	phi = float(np.exp(-0.5 * y @ Q @ y + c @ y))
	g = -Q @ y + c
	return phi, phi * g, phi * (np.outer(g, g) - Q)


def run() -> dict:
	rng = np.random.default_rng(SEED)
	nu = NU
	xi0 = rng.uniform(-0.45, 0.45, size=2)
	alpha0 = rng.uniform(-0.4, 0.4, size=2)
	gamma0 = float(rng.uniform(-0.2, 0.2))
	c = rng.uniform(-0.2, 0.2, size=2)
	z0 = np.concatenate([xi0, alpha0, [gamma0]])

	def state_rhs(_s: float, z: np.ndarray) -> np.ndarray:
		xi, alpha, gamma = z[:2], z[2:4], z[4]
		g = grad_potential(xi)
		gamma_dot = nu * alpha @ A @ alpha + g @ xi - potential(xi)
		return np.concatenate([-2 * nu * A @ alpha, -g, [gamma_dot]])

	def state_at(s: float) -> np.ndarray:
		if s == 0:
			return z0.copy()
		return solve_ivp(state_rhs, (0.0, s), z0, rtol=2e-13, atol=2e-14,
			method="DOP853").y[:, -1]

	def transformed_rhs(y: np.ndarray) -> float:
		phi, _grad, hess = gaussian_jet(y, c)
		vtilde = potential(xi0 + y) - potential(xi0) - grad_potential(xi0) @ y
		return float(nu * np.sum(A * hess) - vtilde * phi)

	def psi(x: np.ndarray, s: float) -> float:
		z = state_at(s)
		xi, alpha, gamma = z[:2], z[2:4], z[4]
		y = x - xi
		phi, _grad, _hess = gaussian_jet(y, c)
		# This local time jet satisfies the transformed PDE at s=0 by construction.
		return float(np.exp(alpha @ x + gamma) * (phi + s * transformed_rhs(y)))

	# Confirm A is strictly positive definite and V is genuinely non-quadratic.
	eigenvalues = np.linalg.eigvalsh(A)
	assert np.all(eigenvalues > 0)
	points = xi0 + rng.uniform(-0.5, 0.5, size=(8, 2))
	space_steps = (0.08, 0.04, 0.02, 0.01)
	rows = []
	for dx in space_steps:
		residuals = []
		for x in points:
			f0 = psi(x, 0.0)
			fip = [psi(x + dx * np.eye(2)[i], 0.0) for i in range(2)]
			fim = [psi(x - dx * np.eye(2)[i], 0.0) for i in range(2)]
			dxx = (fip[0] - 2 * f0 + fim[0]) / dx**2
			dyy = (fip[1] - 2 * f0 + fim[1]) / dx**2
			fpp = psi(x + dx * np.array([1.0, 1.0]), 0.0)
			fpm = psi(x + dx * np.array([1.0, -1.0]), 0.0)
			fmp = psi(x + dx * np.array([-1.0, 1.0]), 0.0)
			fmm = psi(x - dx * np.array([1.0, 1.0]), 0.0)
			dxy = (fpp - fpm - fmp + fmm) / (4 * dx**2)
			sp = psi(x, TIME_STEP)
			sm = psi(x, -TIME_STEP)
			lhs = (sp - sm) / (2 * TIME_STEP)
			rhs = nu * (A[0, 0] * dxx + 2 * A[0, 1] * dxy + A[1, 1] * dyy) - potential(x) * f0
			residuals.append(float(abs(lhs - rhs) / max(1.0, abs(lhs), abs(rhs))))
		rows.append({"dx": dx, "max_relative_residual": max(residuals), "residuals": residuals})
	return {
		"seed": SEED,
		"nu": nu,
		"A": A.tolist(),
		"A_eigenvalues": eigenvalues.tolist(),
		"potential": "0.3(1-cos x1)+0.2(1-cos x2)+0.15 sin(x1+2x2)+0.03 cos(2x1-x2)",
		"xi0": xi0.tolist(),
		"alpha0": alpha0.tolist(),
		"gamma0": gamma0,
		"phi0": "exp(-0.5*y^T Q y + c^T y), with its transformed-PDE time derivative imposed at s=0",
		"time_difference_step": TIME_STEP,
		"space_checks": rows,
		"scope": "Finite-difference residual check of the identity at s=0 for random points. It tests the formula with SPD A and non-quadratic V; it is not a global PDE solve or a proof of the Lean theorem.",
	}


if __name__ == "__main__":
	data = run()
	OUT.write_text(json.dumps(data, indent=2) + "\n")
	print(json.dumps({"xi0": data["xi0"], "alpha0": data["alpha0"],
		"max_residual_by_dx": {str(r["dx"]): r["max_relative_residual"] for r in data["space_checks"]},
		"output": str(OUT)}, indent=2))
