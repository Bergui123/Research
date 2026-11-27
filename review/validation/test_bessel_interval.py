"""Cross-check the ball-arithmetic Jacobi--Anger tail evaluator."""
import unittest
from fractions import Fraction as F

import numpy as np
from flint import arb
from scipy.special import jv

from check_bessel_interval import enclose_tail


class BesselIntervalTests(unittest.TestCase):
    def test_small_tail_exceeds_independent_finite_sum(self):
        result = enclose_tail(5.0, 20, 1e-8)
        observed = 2 * np.sum(np.abs(jv(np.arange(21, 221), 5.0)))
        self.assertTrue(result['upper'] > arb(observed))
        self.assertTrue(result['upper'] < arb(1e-8))

    def test_fails_when_threshold_is_too_small(self):
        with self.assertRaises(AssertionError):
            enclose_tail(5.0, 20, 1e-30)

    def test_small_argument_does_not_amplify_recurrence_uncertainty(self):
        result = enclose_tail(1.3360799846266507, 12, 5.3947700215042114e-11)
        self.assertTrue(result['upper'] < arb(5.3947700215042114e-11))

    def test_exact_rational_argument_and_threshold(self):
        result = enclose_tail(F(4, 3), 12, F(1, 10**10))
        self.assertTrue(result['upper'] < arb('1e-10'))

    def test_recurrence_branch_exceeds_independent_finite_sum(self):
        result = enclose_tail(200.0, 250, 1e-5)
        observed = 2 * np.sum(np.abs(jv(np.arange(251, 451), 200.0)))
        self.assertTrue(result['upper'] > arb(observed))
        self.assertTrue(result['upper'] < arb(1e-5))


if __name__ == '__main__':
    unittest.main()
