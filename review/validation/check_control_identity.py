"""Numerically check the completion-of-squares identity for random controls.

Run: OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_control_identity.py
These finite cases do not replace the analytic proof.
"""
import json

import numpy as np
from scipy.integrate import quad_vec, solve_ivp

rng = np.random.default_rng(15021)
max_identity_error = 0.0
min_gap = float('inf')
for _ in range(20):
    dim = 2
    z = rng.normal(size=(3, dim, dim))
    h = z[0].T @ z[0] + .5*np.eye(dim)
    c = z[1].T @ z[1]
    phi = z[2].T @ z[2] + .2*np.eye(dim)
    hinv = np.linalg.inv(h)
    lam = float(rng.uniform(0, 2))
    horizon = float(rng.uniform(.3, 2))
    x0, a, b = rng.normal(size=(3, dim))
    def ode(_, y):
        p = y.reshape(dim, dim)
        return (2*lam*c - .5*p@hinv@p).ravel()
    sol = solve_ivp(ode, (0, horizon), phi.ravel(), dense_output=True,
                    rtol=1e-11, atol=1e-13)
    if not sol.success:
        raise RuntimeError(sol.message)
    def quantities(t):
        v = a + b*np.sin(t)
        x = x0 - a*t + b*(np.cos(t)-1)
        p = sol.sol(horizon-t).reshape(dim, dim)
        difference = v - .5*hinv@p@x
        return np.array([v@h@v + lam*x@c@x, difference@h@difference])
    integrals = quad_vec(quantities, 0, horizon, epsabs=1e-10, epsrel=1e-10)[0]
    xend = x0 - a*horizon + b*(np.cos(horizon)-1)
    initial = .5*x0@sol.sol(horizon).reshape(dim,dim)@x0
    gap = float(integrals[0] + .5*xend@phi@xend - initial)
    max_identity_error = max(max_identity_error, abs(gap-integrals[1]))
    min_gap = min(min_gap, gap)
assert max_identity_error < 1e-8 and min_gap >= 0
print(json.dumps({'seed': 15021, 'cases': 20,
                  'max_absolute_identity_residual': max_identity_error,
                  'smallest_cost_gap': min_gap}, indent=2))
