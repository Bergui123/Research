"""Focused exact checks for the rational Riccati residual certificate."""

from fractions import Fraction as F
import json
import math
from pathlib import Path
import unittest

from verify_riccati_nodes import bernstein, quadratic, verify

ROOT = Path(__file__).resolve().parents[2]


class RiccatiCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "review/validation/riccati-quad-nodes.json").read_text())

    def test_quadratic_matches_all_three_nodes(self):
        left, middle, right = F(7, 5), F(11, 8), F(13, 9)
        p = quadratic(left, middle, right)
        value = lambda x: sum(c * x**i for i, c in enumerate(p))
        self.assertEqual(value(F(0)), left)
        self.assertEqual(value(F(1, 2)), middle)
        self.assertEqual(value(F(1)), right)

    def test_bernstein_conversion_matches_polynomial_at_rational_points(self):
        coeffs = [F(2, 5), F(-3, 7), F(5, 11), F(7, 13), F(-1, 3)]
        b = bernstein(coeffs)
        for k in range(8):
            x = F(k, 7)
            power_value = sum(c * x**i for i, c in enumerate(coeffs))
            bernstein_value = sum(b[i] * math.comb(4, i) * x**i * (1 - x)**(4 - i) for i in range(5))
            self.assertEqual(power_value, bernstein_value)
            self.assertLessEqual(abs(power_value), max(map(abs, b)))

    def test_certificate_and_invalid_nodes(self):
        result = verify(self.data)
        self.assertLess(F(result["residual_integral_upper_rational"]), F(13, 10**9))
        corrupted = dict(self.data)
        corrupted["nodes"] = list(self.data["nodes"])
        corrupted["nodes"][0] = list(corrupted["nodes"][0])
        corrupted["nodes"][0][0] += 1
        with self.assertRaises(AssertionError):
            verify(corrupted)
        corrupted["nodes"][0] = self.data["nodes"][0]
        corrupted["nodes"][1] = list(corrupted["nodes"][1])
        corrupted["nodes"][1][2] = 0
        with self.assertRaises(AssertionError):
            verify(corrupted)


if __name__ == "__main__":
    unittest.main()
