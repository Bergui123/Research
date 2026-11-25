import unittest

import numpy as np
from scipy.special import ive

import direct_polynomial


class DirectPolynomialTest(unittest.TestCase):
    def test_finite_grid_check_returns_contraction_scaled_output(self):
        h = np.diag([0.0, 1.0, 2.0])
        u0 = np.array([1.0, 0.0, 0.0])

        output, _ = direct_polynomial.check(h, u0, alpha=1.0, k=1, eps_state=4.0)

        tail = 2.0 * np.sum(ive(np.arange(2, 100), 1.0))
        expected = (ive(0, 1.0) + 2.0 * ive(1, 1.0)) / (1.0 + tail)
        self.assertAlmostEqual(output[0], expected, places=12)
        self.assertLess(output[0], 1.0)


if __name__ == "__main__":
    unittest.main()
