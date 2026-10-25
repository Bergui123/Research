"""Reproduce the Steps 1--3 audit: python3 verify_research.py.

Requires numpy, scipy, sympy. All costs below exclude permanent/spread offsets.
The discrete convergence check sets permanent impact and spread to zero.
"""

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import spsolve

from BellmanPrinciple import derive_hjb


def check_hjb():
    x1, x2, v1, v2, t, T = sp.symbols("x1 x2 v1 v2 t T", real=True)
    h1, h2, lam = sp.symbols("h1 h2 lam", positive=True)
    a, b, d = sp.symbols("a b d", real=True)
    x, v = sp.Matrix([x1, x2]), sp.Matrix([v1, v2])
    H, C = sp.diag(h1, h2), sp.Matrix([[a, b], [b, d]])
    result = derive_hjb((v.T * H * v + lam * x.T * C * x)[0],
                        [x1, x2], [v1, v2], [-v1, -v2], t, T)
    grad = sp.Matrix([sp.diff(result.value_function, z) for z in x])
    expected = -(grad.T * H.inv() * grad)[0] / 4 + lam * (x.T * C * x)[0]
    assert sp.simplify(result.minimized_hamiltonian - expected) == 0
    p, q, r = sp.symbols("p q r", real=True)
    P = sp.Matrix([[p, q], [q, r]])
    Pdot = P * H.inv() * P / 2 - 2 * lam * C
    residual = (x.T * Pdot * x)[0] / 2 - (x.T * P * H.inv() * P * x)[0] / 4
    assert sp.expand(residual + lam * (x.T * C * x)[0]) == 0
    print("HJB and Riccati symbolic residuals: zero; constant equation: c'=0")


def check_instance():
    C = np.array([[1., .4], [.4, 1.]])
    X = np.array([3., 1.])
    eigenvalues, U = np.linalg.eigh(C)
    k = np.sqrt(eigenvalues)

    def trajectory(t):
        return U @ (np.sinh(k * (1 - t)) / np.sinh(k) * (U.T @ X))

    def velocity(t):
        return U @ (k * np.cosh(k * (1 - t)) / np.sinh(k) * (U.T @ X))

    ratio = (k[0] / np.sinh(k[0])) / (k[1] / np.sinh(k[1]))
    theta_max = (ratio + 1) / (ratio - 1)
    assert 3 < theta_max
    assert np.allclose(trajectory(0), X) and np.allclose(trajectory(1), 0)
    # Dense sampling supplements, but does not replace, the analytic sign proof.
    sampled_min = np.array([velocity(t) for t in np.linspace(0, 1, 10001)]).min(axis=0)
    assert np.all(sampled_min > 0)
    P0 = 2 * U @ np.diag(k / np.tanh(k)) @ U.T
    value = X @ P0 @ X / 2
    integral = quad(lambda t: velocity(t) @ velocity(t) + trajectory(t) @ C @ trajectory(t),
                    0, 1, epsabs=1e-11, epsrel=1e-11)[0]
    assert abs(value - integral) < 1e-10
    print(f"Stage 3: eigenvalues={eigenvalues}; kappa={k}")
    print(f"ratio={ratio:.12f}; theta_max={theta_max:.12f}")
    print(f"v(0)={velocity(0)}; sampled min v={sampled_min}")
    print(f"J={value:.12f}; quadrature error={abs(value-integral):.3e}")
    previous = None
    for N in (16, 64, 256, 1024):
        dt = 1 / N
        # Stationarity of sum |x[k]-x[k-1]|^2/dt + dt*x[k]'C*x[k].
        lap = diags([-np.ones(N-2), 2*np.ones(N-1), -np.ones(N-2)], [-1, 0, 1])
        A = kron(lap, eye(2)) + dt**2 * kron(eye(N-1), C)
        rhs = np.zeros((N-1, 2))
        rhs[0] = X
        interior = spsolve(A.tocsc(), rhs.ravel()).reshape(-1, 2)
        path = np.vstack((X, interior, np.zeros(2)))
        velocity_discrete = -np.diff(path, axis=0) / dt
        assert velocity_discrete.min() > 0
        grid = np.linspace(0, 1, N+1)
        error = np.max(np.abs(path - np.array([trajectory(t) for t in grid])))
        Jd = dt * (np.sum(velocity_discrete**2) + np.einsum("ni,ij,nj->", path[1:], C, path[1:]))
        assert previous is None or error < previous
        previous = error
        print(f"N={N}: max inventory error={error:.3e}; discrete J={Jd:.12f}")
    assert previous < 1e-6


def check_finite_penalty():
    # Unequal impacts and a penalty that does not commute with the risk matrix.
    H = np.diag([1., 2.])
    C = np.array([[1., .4], [.4, 1.]])
    Phi = np.array([[3., .3], [.3, 2.]])
    hs = np.diag(np.sqrt(np.diag(H)))
    hi = np.diag(1 / np.sqrt(np.diag(H)))
    mu, U = np.linalg.eigh(hi @ C @ hi)
    K = U @ np.diag(np.sqrt(mu)) @ U.T
    G = hi @ Phi @ hi / 2
    assert np.linalg.norm(K @ G - G @ K) > .01

    def closed(s):
        E = expm(-K * s)
        B = np.linalg.solve(K, np.eye(2) - E @ E) / 2 + E @ np.linalg.inv(G + K) @ E
        return 2 * hs @ (np.linalg.inv(B) - K) @ hs

    def ode(s, flat):
        P = flat.reshape(2, 2)
        return (2*C - P @ np.linalg.solve(H, P) / 2).ravel()

    sol = solve_ivp(ode, [0, 1], Phi.ravel(), rtol=1e-11, atol=1e-12, dense_output=True)
    assert sol.success
    boundary_error = np.max(np.abs(closed(0)-Phi))
    ode_error = max(np.max(np.abs(closed(s)-sol.sol(s).reshape(2, 2)))
                    for s in np.linspace(0, 1, 101))
    assert boundary_error < 1e-12 and ode_error < 1e-9
    old_boundary = 2 * hs @ (np.linalg.inv(G + np.linalg.inv(K)) - K) @ hs
    print(f"Finite Phi: terminal error={boundary_error:.3e}; ODE error={ode_error:.3e}")
    print(f"Previous formula terminal error={np.max(np.abs(old_boundary-Phi)):.6f}")


def check_cole_hopf():
    x1, x2 = sp.symbols("x1 x2", real=True)
    nu, h1, h2, lam = sp.symbols("nu h1 h2 lam", positive=True)
    p, q, r, a, b, d = sp.symbols("p q r a b d", real=True)
    x = sp.Matrix([x1, x2])
    P, C = sp.Matrix([[p, q], [q, r]]), sp.Matrix([[a, b], [b, d]])
    Hi = sp.diag(1/h1, 1/h2)
    W = (x.T * P * x)[0] / 2
    psi = sp.exp(-W / (4*nu))
    Ps = 2*lam*C - P*Hi*P/2
    Ws = (x.T*Ps*x)[0]/2 + nu*sp.trace(Hi*P)
    rhs = nu*sum(Hi[j,j]*sp.diff(psi, x[j], 2) for j in range(2))
    rhs -= lam*(x.T*C*x)[0]*psi/(4*nu)
    assert sp.simplify(-Ws/(4*nu) - rhs/psi) == 0
    for j in range(2):
        assert sp.simplify(-4*nu*sp.diff(psi, x[j])/psi - (P*x)[j]) == 0
    print("Cole--Hopf PDE and gradient symbolic residuals: zero")


if __name__ == "__main__":
    check_hjb()
    check_instance()
    check_finite_penalty()
    check_cole_hopf()
    print("All audit checks passed.")
