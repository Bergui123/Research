import unittest
from fractions import Fraction as F

from check_expansion_enclosures import exp_negative
from verify_polynomial_nodes import error_bounds, ladder_error, residuals, round_up


class PolynomialCertificateTests(unittest.TestCase):
    def test_rational_exp_encloses_independent_series(self):
        # e is between these positive partial sums with an elementary tail bound.
        lo, hi = exp_negative(F(1))
        partial = sum((F(1, __import__('math').factorial(k)) for k in range(41)), F(0))
        tail = F(1, __import__('math').factorial(41)) * F(42, 41)
        self.assertLess(lo, 1 / (partial + tail))
        self.assertGreater(hi, 1 / partial)
        self.assertEqual(exp_negative(F(0)), (F(1), F(1)))

    def test_outward_rounding(self):
        self.assertEqual(round_up(F(1, 3), 100), F(34, 100))
        self.assertEqual(round_up(F(7, 10), 100), F(7, 10))

    def test_residual_sees_nonzero_potential_for_zero_polynomials(self):
        zero = [F(0)] * 3
        one = [F(1), F(0), F(0)]
        matrix = [[one, zero], [zero, one]]
        block1 = {(i, d-i): zero for d in (0, 2, 4) for i in range(d+1)}
        block2 = {(i, d-i): zero for d in (0, 2, 4, 6, 8) for i in range(d+1)}
        r = residuals(matrix, [block1, block2], 1)
        self.assertNotEqual(r[0][(0, 0)][0], 0)
        self.assertNotEqual(r[1][(0, 0)][0], 0)

    def test_error_propagation_is_nonnegative_and_sensitive(self):
        approx = [{d: F(1) for d in (0, 2, 4)}, {d: F(1) for d in (0, 2, 4, 6, 8)}]
        res = [{d: F(0) for d in block} for block in approx]
        a = error_bounds(approx, res)
        res[0][4] = F(1)
        b = error_bounds(approx, res)
        self.assertGreater(sum(b[0][0]), sum(a[0][0]))
        self.assertGreater(sum(b[1][0]), sum(a[1][0]))
        self.assertTrue(all(v >= 0 for block in b for poly in block.values() for v in poly))
        self.assertTrue(ladder_error({0: [F(1)]}, {}, 4))


if __name__ == '__main__':
    unittest.main()
