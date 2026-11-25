import math
import unittest

import numpy as np
from step8_resources import potential_scale


class ComparatorScaleTest(unittest.TestCase):
    def test_quarter_grid_uses_exact_comparator_value(self):
        _, k, q = potential_scale(1, np.array([-6.0]), np.array([3.5]), 1)
        m = k + 4
        self.assertEqual(q[0, 0], 157.25)
        self.assertEqual(int(np.floor(q[0, 0] * 2 ** (m - k))), 2516)

    def test_power_of_two_and_neighbours_fit_in_m_bits(self):
        for multiplier in (1 - 1e-10, 1, 1 + 1e-10):
            x = np.array([math.sqrt(8 / 5) * multiplier])
            _, k, q = potential_scale(1, x, np.array([0.0]), 1)
            with self.subTest(maximum=float(q.max())):
                self.assertGreater(2 ** k, q.max())
                m = 3
                top = np.floor(q * 2 ** (m - k))
                self.assertLess(top.max(), 2 ** m)
                self.assertTrue(np.all(q - top * 2 ** (k - m) >= 0))
                self.assertTrue(np.all(q - top * 2 ** (k - m) < 2 ** (k - m)))


if __name__ == '__main__':
    unittest.main()
