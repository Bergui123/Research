import unittest

import numpy as np
from step7_schrodinger import g_hat_fft


class AuxiliaryFourierTest(unittest.TestCase):
    def test_fft_matches_half_integer_basis_and_inverse(self):
        for n in (4, 8, 16, 32):
            p = -1 + 2 * np.arange(n) / n
            eta = np.pi * (np.arange(n) - (n - 1) / 2)
            values = np.random.default_rng(n).normal(size=n)
            basis = np.exp(1j * np.outer(p, eta))
            coeff = g_hat_fft(values, n)
            with self.subTest(n=n):
                np.testing.assert_allclose(coeff, basis.conj().T @ values / n, atol=2e-14)
                np.testing.assert_allclose(basis @ coeff, values, atol=2e-14)
                np.testing.assert_allclose(np.exp(1j * np.outer(p + 2, eta)) @ coeff,
                                           -values, atol=2e-14)


if __name__ == '__main__':
    unittest.main()
