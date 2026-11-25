import math
import unittest

import step8_resources


class AmplitudeBudgetTest(unittest.TestCase):
    def test_negative_point_one_is_counterexample_to_linear_log_bound(self):
        relative_error = -0.1
        exact_log_error = abs(math.log1p(relative_error))
        self.assertAlmostEqual(exact_log_error, -math.log(0.9))
        self.assertGreater(exact_log_error, abs(relative_error))

    def test_finite_budget_covers_composed_errors_at_threshold(self):
        log_budget = 0.1
        delta = step8_resources.amplitude_budget(log_budget)
        self.assertAlmostEqual(-math.log1p(-delta), log_budget)
        for relative_error in (-delta, -delta * (1 - 1e-12),
                               delta * (1 - 1e-12), delta):
            self.assertLessEqual(abs(math.log1p(relative_error)),
                                 log_budget + 1e-15)
        self.assertLess(delta, log_budget)

    def test_three_deterministic_sources_and_estimation_compose_at_four_nodes(self):
        log_budget = 0.1
        delta = step8_resources.amplitude_budget(log_budget)
        lower, upper = 1 - delta, 1 / (1 - delta)
        deterministic_parts = (
            (0.10, 0.10, 0.10),
            (-0.15, -0.10, -0.05),
            (0.20, -0.05, 0.05),
            (-0.05, 0.20, 0.05),
        )
        estimation_parts = (0.70, -0.70, 0.65, -0.55)

        for parts, estimation_share in zip(deterministic_parts,
                                           estimation_parts):
            # Each deterministic source is bounded against the same reference
            # amplitude. The triangle inequality bounds their combined error.
            deterministic = sum(parts) * delta
            estimation = estimation_share * delta
            ratio = (1 + deterministic) * (1 + estimation)
            self.assertGreaterEqual(ratio, lower - 1e-15)
            self.assertLessEqual(ratio, upper + 1e-15)

    def test_unequal_errors_at_all_four_nodes_bound_two_component_rate(self):
        log_budget = 0.1
        delta = step8_resources.amplitude_budget(log_budget)
        nu, d = 1.7, 2
        node_errors = (0.55 * delta, -0.85 * delta,
                       -0.35 * delta, 0.71 * delta)
        component_log_errors = (
            math.log1p(node_errors[0]) - math.log1p(node_errors[1]),
            math.log1p(node_errors[2]) - math.log1p(node_errors[3]),
        )
        component_bounds = tuple(
            -math.log1p(-abs(node_errors[2 * axis]))
            -math.log1p(-abs(node_errors[2 * axis + 1]))
            for axis in range(2)
        )

        for error, bound in zip(component_log_errors, component_bounds):
            self.assertLessEqual(abs(error), bound + 1e-15)
        rate_error = nu / d * math.hypot(*component_log_errors)
        rate_bound = nu / d * math.hypot(*component_bounds)
        self.assertLessEqual(rate_error, rate_bound + 1e-15)
        self.assertLessEqual(rate_bound,
                             2 * math.sqrt(2) * nu * log_budget / d + 1e-15)


if __name__ == '__main__':
    unittest.main()
