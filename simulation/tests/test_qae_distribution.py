import unittest

import numpy as np
import step10_validation


class QAEDistributionTest(unittest.TestCase):
    def test_four_outcomes_match_exact_distribution(self):
        theta = 0.11
        m = 4
        a = np.sin(np.pi * theta) ** 2
        y = np.arange(m)
        delta = y / m - theta
        probability = (np.sinc(m * delta) / np.sinc(delta)) ** 2
        probability /= probability.sum()
        step10_validation.RNG = np.random.default_rng(521)
        samples = step10_validation.sample_qae(a, m, 200000)
        expected = (probability[0], probability[1] + probability[3], probability[2])
        for value, p in zip((0.0, 0.5, 1.0), expected):
            with self.subTest(value=value):
                self.assertLess(abs(np.mean(np.isclose(samples, value)) - p), 0.004)

    def test_large_grid_retains_tail_outcomes(self):
        m = 65536
        theta = 7208.5 / m
        a = np.sin(np.pi * theta) ** 2
        step10_validation.RNG = np.random.default_rng(813)
        samples = step10_validation.sample_qae(a, m, 1000000)
        # The old sampler discarded |k| > 20000 around either eigenphase.
        folded_y = np.arcsin(np.sqrt(samples)) * m / np.pi
        self.assertGreater(np.count_nonzero(folded_y > 27209), 0)


if __name__ == '__main__':
    unittest.main()
