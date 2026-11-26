"""Checks for the parameterized polynomial-error comparison."""
from fractions import Fraction as F
import unittest

from sweep_polynomial import error_bounds


class SweepPolynomialTest(unittest.TestCase):
    def test_zero_error_and_positive_residual_propagation(self):
        a=[{d:F(0) for d in (0,2,4)},{d:F(0) for d in (0,2,4,6,8)}]
        r=[dict(b) for b in a]
        bounds=error_bounds(a,r,F(1),F(0))
        self.assertTrue(all(not any(poly) for block in bounds for poly in block.values()))
        r[0][4]=F(1)
        bounds=error_bounds(a,r,F(1),F(0))
        self.assertGreater(sum(bounds[0][0]),0)
        self.assertTrue(all(v>=0 for block in bounds for p in block.values() for v in p))


if __name__=='__main__':unittest.main()
