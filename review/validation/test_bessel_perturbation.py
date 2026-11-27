"""Checks for the argument-perturbation bound."""
from fractions import Fraction as F
import unittest

import numpy as np
from scipy.special import jv

from check_bessel_perturbation import derivative_upper, input_interval


class BesselPerturbationTests(unittest.TestCase):
    def test_derivative_bound_covers_independent_finite_differences(self):
        lo, hi, n = F(1), F(101, 100), 12
        left = 2*np.sum(np.abs(jv(np.arange(n, n+200), float(lo))))
        right = 2*np.sum(np.abs(jv(np.arange(n, n+200), float(hi))))
        self.assertLess(abs(right-left), float((hi-lo)*derivative_upper(hi, n)))

    def test_interval_contains_center_and_exact_formula(self):
        cfg = {'alpha_walk': 10.5, 'eta_max': 7.0, 'n_p': 7}
        center = 10.5 * (2*7.0/(2**7-1)) * 2**(3-1)
        lo, hi = input_interval(cfg, 3, center)
        exact = F.from_float(10.5)*F.from_float(7.0)*8/127
        self.assertLessEqual(lo, min(F.from_float(center), exact))
        self.assertGreaterEqual(hi, max(F.from_float(center), exact))


if __name__ == '__main__':
    unittest.main()
