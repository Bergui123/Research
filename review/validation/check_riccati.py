"""Compare the block-exponential Riccati solution with an independent ODE solve.

Run: OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_riccati.py
The sample is numerical evidence, not a proof for all positive matrices.
"""
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

rng = np.random.default_rng(420613)
worst_rel = 0.0
worst_sym = 0.0
min_y_det = float('inf')
min_p_eig = float('inf')
cases = 0
for dim in (1, 2, 4):
    for _ in range(20):
        matrices = rng.normal(size=(3, dim, dim))
        h = matrices[0].T @ matrices[0] + .5*np.eye(dim)
        c = matrices[1].T @ matrices[1]
        phi = matrices[2].T @ matrices[2] + .1*np.eye(dim)
        lam = float(rng.uniform(0, 2))
        duration = float(rng.uniform(.1, 2))
        hinv = np.linalg.inv(h)
        def ode(_, y):
            p = y.reshape(dim, dim)
            return (2*lam*c - .5*p@hinv@p).ravel()
        sol = solve_ivp(ode, (0, duration), phi.ravel(), rtol=1e-11, atol=1e-13)
        if not sol.success:
            raise RuntimeError(sol.message)
        p_ode = sol.y[:, -1].reshape(dim, dim)
        generator = np.block([[np.zeros((dim, dim)), .5*hinv],
                              [2*lam*c, np.zeros((dim, dim))]])
        yz = expm(duration*generator) @ np.vstack((np.eye(dim), phi))
        y, z = yz[:dim], yz[dim:]
        p_block = np.linalg.solve(y.T, z.T).T
        worst_rel = max(worst_rel, float(np.linalg.norm(p_ode-p_block) / np.linalg.norm(p_block)))
        worst_sym = max(worst_sym, float(np.linalg.norm(p_block-p_block.T)))
        min_y_det = min(min_y_det, float(np.linalg.det(y)))
        min_p_eig = min(min_p_eig, float(np.linalg.eigvalsh((p_block+p_block.T)/2).min()))
        cases += 1
assert worst_rel < 1e-9 and worst_sym < 1e-9 and min_y_det > 0 and min_p_eig > 0
print(json.dumps({'seed': 420613, 'cases': cases, 'dimensions': [1,2,4],
                  'max_relative_difference': worst_rel, 'max_symmetry_residual': worst_sym,
                  'min_det_y': min_y_det, 'min_eigenvalue_p': min_p_eig}, indent=2))
