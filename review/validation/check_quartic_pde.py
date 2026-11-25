"""Independent quartic HJB check: shooting for v_0 and spectral heat PDE for v_nu."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


KAPPA = 0.1
NU = 1.0
LAMBDA = 1.0
T = 1.0
C = np.array([[1.0, 0.4], [0.4, 1.0]])
PHI = np.array([[2.0, 0.3], [0.3, 1.0]])
X = np.array([3.0, 1.0])
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "review" / "validation" / "quartic-pde-results.json"


def terminal_cost(y: np.ndarray) -> float:
	return float(0.5 * y @ PHI @ y + KAPPA * np.sum(y**4))


def grad_terminal_cost(y: np.ndarray) -> np.ndarray:
	return PHI @ y + 4 * KAPPA * y**3


def shoot(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
	z0 = np.concatenate([y, 0.5 * grad_terminal_cost(y)])
	def rhs(_s: float, z: np.ndarray) -> np.ndarray:
		return np.concatenate([z[2:], LAMBDA * C @ z[:2]])
	result = solve_ivp(rhs, (0.0, T), z0, method="DOP853", rtol=2e-12, atol=2e-13)
	return result.y[:2, -1], result.y[2:, -1]


def solve_inviscid() -> dict:
	def endpoint_residual(y: np.ndarray) -> np.ndarray:
		return shoot(y)[0] - X

	opt = root(endpoint_residual, np.array([0.5, 0.2]), method="hybr", options={"xtol": 1e-12})
	if not opt.success and np.linalg.norm(endpoint_residual(opt.x)) > 1e-10:
		raise RuntimeError(f"shooting root failed: {opt.message}")
	xT, vT = shoot(opt.x)
	return {"terminal_point": opt.x.tolist(), "endpoint": xT.tolist(), "endpoint_residual": float(np.linalg.norm(xT - X)),
		"initial_velocity": (0.5 * grad_terminal_cost(opt.x)).tolist(), "terminal_velocity": vT.tolist(),
		"v0": vT.tolist(), "terminal_cost": terminal_cost(opt.x)}


def spectral_heat(n: int, dt: float) -> dict:
	# Periodic Fourier pseudospectral box [-8,8)^2. The quartic Gaussian is
	# negligible at the boundary; X falls exactly on the grid for n multiple of 16.
	length = 8.0
	dx = 2 * length / n
	xgrid = -length + dx * np.arange(n)
	ix = int(round((X[0] + length) / dx))
	iy = int(round((X[1] + length) / dx))
	if not (np.isclose(xgrid[ix], X[0]) and np.isclose(xgrid[iy], X[1])):
		raise ValueError("grid does not contain the readout point")
	x1, x2 = np.meshgrid(xgrid, xgrid, indexing="ij")
	g = 0.5 * (PHI[0, 0] * x1**2 + 2 * PHI[0, 1] * x1 * x2 + PHI[1, 1] * x2**2) + KAPPA * (x1**4 + x2**4)
	u = np.exp(-g / (4 * NU))
	v = LAMBDA / (4 * NU) * (C[0, 0] * x1**2 + 2 * C[0, 1] * x1 * x2 + C[1, 1] * x2**2)
	k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
	k2 = k[:, None] ** 2 + k[None, :] ** 2
	diffusion = np.exp(-NU * k2 * dt)
	potential_half = np.exp(-0.5 * dt * v)
	steps = round(T / dt)
	if not np.isclose(steps * dt, T):
		raise ValueError("dt must divide T")
	for _ in range(steps):
		u *= potential_half
		u = np.fft.ifft2(np.fft.fft2(u) * diffusion).real
		u *= potential_half
	if u[ix, iy] <= 0:
		raise FloatingPointError("non-positive heat solution at readout point")
	fft_u = np.fft.fft2(u)
	ux = np.fft.ifft2(1j * k[:, None] * fft_u).real
	uy = np.fft.ifft2(1j * k[None, :] * fft_u).real
	grad_log = np.array([ux[ix, iy], uy[ix, iy]]) / u[ix, iy]
	vnu = -2 * NU * grad_log
	return {"n": n, "dt": dt, "dx": dx, "steps": steps, "psi_at_X": float(u[ix, iy]),
		"grad_log_psi": grad_log.tolist(), "vnu": vnu.tolist()}


def run() -> dict:
	inv = solve_inviscid()
	cases = [spectral_heat(128, 0.002), spectral_heat(256, 0.002),
		spectral_heat(256, 0.001), spectral_heat(384, 0.002), spectral_heat(384, 0.001)]
	for case in cases:
		case["relative_bias_vs_shooting"] = float(np.linalg.norm(np.array(case["vnu"]) - np.array(inv["v0"])) / np.linalg.norm(inv["v0"]))
	return {
		"method": "Independent Hamiltonian shooting for the inviscid characteristic; Fourier pseudospectral Strang splitting of the 2D heat equation for nu=1.",
		"parameters": {"kappa": KAPPA, "nu": NU, "lambda": LAMBDA, "T": T, "C": C.tolist(), "Phi": PHI.tolist(), "X": X.tolist()},
		"inviscid_shooting": inv,
		"viscous_pde_cases": cases,
		"scope": "Finite periodic box and second-order time splitting are numerical approximations. Resolution/timestep convergence is reported; this is independent of the Mehler-kernel quadrature in nonquadratic.py.",
	}


if __name__ == "__main__":
	data = run()
	OUT.write_text(json.dumps(data, indent=2) + "\n")
	print(json.dumps({"v0": data["inviscid_shooting"]["v0"],
		"cases": [{"n": c["n"], "dt": c["dt"], "vnu": c["vnu"], "relative_bias": c["relative_bias_vs_shooting"]} for c in data["viscous_pde_cases"]],
		"output": str(OUT)}, indent=2))
