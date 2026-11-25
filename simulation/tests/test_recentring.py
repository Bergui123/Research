"""Independent checks of the recentring identity, the quartic references and the GQSP leakage bound."""

import unittest

import numpy as np
from scipy import sparse
from scipy.integrate import solve_bvp, solve_ivp
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

import nonquadratic as nq
from model import C, LAM, PHI, T, X0, axis_grid, grid_operator


class QuarticReferenceTest(unittest.TestCase):
    def test_inviscid_rate_by_pontryagin_bvp(self):
        """x(s), pi(s) in time remaining: x' = pi/2, pi' = 2 lam C x, pi(0) = grad g(x(0)), x(T) = X."""
        rhs = lambda s, y: np.vstack([y[2:] / 2, 2 * LAM * C @ y[:2]])
        bc = lambda a, b: np.concatenate([a[2:] - PHI @ a[:2] - 4 * nq.KAPPA * a[:2] ** 3, b[:2] - X0])
        s = np.linspace(0, T, 40)
        guess = np.vstack([np.outer(X0, s), np.ones((2, 40))])
        sol = solve_bvp(rhs, bc, s, guess, tol=1e-10, max_nodes=20000)
        self.assertTrue(sol.success)
        np.testing.assert_allclose(sol.sol(T)[2:] / 2, nq.V0, rtol=1e-7)

    def test_viscous_rate_by_grid(self):
        """Second-order grid values of the direct equation, extrapolated in the spacing."""
        nu, rates = 1.0, []
        for k in (4, 8):
            n, dx = 24 * k // 2 * 2, 1.0 / k
            x1, _ = axis_grid(n, n * dx / 2, X0[0], dx)
            x2, _ = axis_grid(n, n * dx / 2, X0[1], dx)
            hd, _ = grid_operator(x1, x2, dx, nu)
            pts = np.stack(np.meshgrid(x1, x2, indexing="ij"), -1)
            psi = expm_multiply(-T * hd, np.exp(-nq.g(pts, nq.KAPPA) / (4 * nu)).ravel()).reshape(n, n)
            i, j = int(np.argmin(abs(x1 - X0[0]))), int(np.argmin(abs(x2 - X0[1])))
            lp = np.log(psi)
            rates.append(-2 * nu * np.array([lp[i + 1, j] - lp[i - 1, j], lp[i, j + 1] - lp[i, j - 1]]) / (2 * dx))
        rate = (4 * rates[1] - rates[0]) / 3
        np.testing.assert_allclose(rate, nq.viscous(X0, nu, nq.KAPPA)[1], rtol=2e-4)


class RecentringIdentityTest(unittest.TestCase):
    def test_one_asset_with_quartic_running_cost(self):
        """H = 1/a, V = v2 x^2 + v4 x^4, an arbitrary characteristic (not the optimal one)."""
        nu, a, v2, v4, tf = 0.3, 1.7, 0.4, 0.15, 0.4
        v, dv = lambda x: v2 * x**2 + v4 * x**4, lambda x: 2 * v2 * x + 4 * v4 * x**3
        ode = lambda s, y: [-2 * nu * a * y[1], -dv(y[0]), nu * a * y[1] ** 2 + dv(y[0]) * y[0] - v(y[0])]
        path = solve_ivp(ode, (0, tf), [0.5, -0.8, 0.1], rtol=1e-11, atol=1e-12, dense_output=True).sol
        n, half = 1200, 6.0
        dx = 2 * half / n
        x = (np.arange(n) - n // 2) * dx
        # fourth-order second difference; the solutions vanish at the ends
        e = np.ones(n)
        lap = sparse.diags([-e[:-2] / 12, 4 * e[:-1] / 3, -2.5 * e, 4 * e[:-1] / 3, -e[:-2] / 12], [-2, -1, 0, 1, 2], format="csr") / dx**2
        psi0 = lambda z: np.exp(-0.9 * z**2 + 0.2 * z)
        psi = expm_multiply(tf * (nu * a * lap - sparse.diags(v(x))), psi0(x))
        xi0, al0, ga0 = path(0)
        phi = np.exp(-(al0 * (x + xi0) + ga0)) * psi0(x + xi0)
        steps = 400
        for k in range(steps):  # midpoint exponential: second order in the step
            xi = path((k + 0.5) * tf / steps)[0]
            pot = v(xi + x) - v(xi) - dv(xi) * x
            phi = expm_multiply(tf / steps * (nu * a * lap - sparse.diags(pot)), phi)
        xi, al, ga = path(tf)
        y = np.linspace(-1.5, 1.5, 31)
        back = np.exp(al * (y + xi) + ga) * np.interp(y, x, phi)
        np.testing.assert_allclose(back, np.interp(y + xi, x, psi), rtol=2e-4)


class LeakageTest(unittest.TestCase):
    def test_leak_and_return_bound(self):
        """Segments with success block A = (1 - a) U: the product errs by at most n a + n^2 a (1 + small)."""
        rng = np.random.default_rng(5)
        dim, n, a = 6, 12, 1e-5
        for _ in range(50):
            total, ideal = np.eye(2 * dim), np.eye(dim)
            for _ in range(n):
                q, _ = np.linalg.qr(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
                w1, w2 = (np.linalg.qr(rng.normal(size=(dim, dim)))[0] for _ in range(2))
                c, s = 1 - a, np.sqrt(1 - (1 - a) ** 2)
                seg = np.block([[c * q, -s * q @ w1], [s * w2, c * w2 @ w1]])
                self.assertLess(np.linalg.norm(seg.conj().T @ seg - np.eye(2 * dim)), 1e-10)
                total, ideal = seg @ total, q @ ideal
            self.assertLessEqual(np.linalg.norm(total[:dim, :dim] - ideal, 2), (n * a + n**2 * a) * 1.01)


if __name__ == "__main__":
    unittest.main()
