"""Checks for the omitted Jacobi--Anger tail bound."""
from fractions import Fraction as F
import unittest

import numpy as np
from scipy.special import jv

from check_bessel_tail import omitted_upper


class BesselTailTests(unittest.TestCase):
    def test_bounds_independent_finite_sums(self):
        for x, n in [(F(1), 12), (F(5), 25), (F(20), 50)]:
            cap = omitted_upper(x, n)
            measured = 2 * np.sum(np.abs(jv(np.arange(n, n + 200), float(x))))
            self.assertGreater(float(cap), measured)

    def test_enforces_domain_and_finite_decay(self):
        self.assertGreater(omitted_upper(F(20), 50), 0)
        with self.assertRaises(AssertionError):
            omitted_upper(F(-1), 50)
        with self.assertRaises(AssertionError):
            omitted_upper(F(20), 20)


if __name__ == '__main__':
    unittest.main()
