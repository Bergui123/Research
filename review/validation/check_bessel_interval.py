"""Enclose reconstructed Jacobi--Anger tails using FLINT ball arithmetic.

Arguments and thresholds are reconstructed from the saved floating-point
values. See bessel-interval.md for the remaining limitations.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys

import flint
from flint import arb, ctx, fmpq

from check_bessel_tail import omitted_upper, saved_cases

ROOT = Path(__file__).resolve().parents[2]


def enclose_tail(x, k, threshold):
    """Return a proven ball upper bound below threshold, or raise."""
    x = F(x)
    threshold = F(threshold)
    assert math.isfinite(x) and 0 < x < k+1
    assert math.isfinite(threshold) and threshold > 0
    # Large-order direct evaluation requires more precision to overcome
    # cancellation. Once enclosed, the recurrence uses 256-bit balls.
    ctx.prec = max(256, math.ceil(float(x)*1.4)+512)
    ax = arb(fmpq(x.numerator, x.denominator))
    prev = ax.bessel_j(k+1)
    curr = ax.bessel_j(k+2)
    ctx.prec = 256
    ax, prev, curr = +ax, +prev, +curr
    limit = arb(fmpq(threshold.numerator, threshold.denominator))
    total = 2*abs(prev)
    used = 1
    for n in range(k+2, k+401):
        if used % 8 == 1:
            remaining = omitted_upper(x, n)
            upper = (total + arb(fmpq(remaining.numerator, remaining.denominator))).upper()
            if upper < limit:
                return {'upper': upper, 'terms': used}
        assert total < limit
        total += 2*abs(curr)
        if x < 128:
            # Upward recurrence amplifies ball widths when order >> argument.
            prev, curr = curr, ax.bessel_j(n+1)
        else:
            prev, curr = curr, (2*n/ax)*curr-prev
        used += 1
    raise AssertionError('Could not enclose tail below the threshold')


def main():
    input_path = ROOT/'review/validation/bessel-tail.json'
    data = json.loads(input_path.read_text())
    configurations = dict(saved_cases())
    rows = []
    for case in data['cases']:
        cfg = configurations[case['case']]
        n_p, q_run = cfg['n_p'], cfg['walk_steps_per_run']
        n_rot = 4*q_run + 2*n_p
        assert cfg['eps_rot'] > sys.float_info.min
        # eps_rot is a rounded binary division of eps_state by 2*n_rot.
        # Step down one float, then allow one relative rounding unit for the
        # final positive normal division that produced the tail threshold.
        threshold = (F.from_float(math.nextafter(cfg['eps_rot'], -math.inf)) *
                     (2*n_rot) / (4*n_p*(n_p+1)) * (1-F(1, 2**52)))
        assert 0 < threshold <= F.from_float(case['budget_float'])
        segments = []
        for s in case['segments']:
            # Also check the exact expression formed from the saved binary
            # alpha and eta, before the floating-point operations in the run.
            x_exact = (F.from_float(cfg['alpha_walk']) * F.from_float(cfg['eta_max']) *
                       2**s['segment'] / (2**n_p-1))
            checks = []
            for x in (F.from_float(s['x_float']), x_exact):
                result = enclose_tail(x, s['K'], threshold)
                ratio = (result['upper']/arb(fmpq(threshold.numerator, threshold.denominator))).upper()
                assert ratio < arb(fmpq(999821, 1000000))
                assert result['terms'] <= 113
                checks.append({'terms_enclosed': result['terms'],
                               'tail_over_threshold_upper': str(ratio.fmpq()),
                               'ratio_decimal_approx': float(ratio)})
            segments.append({'segment': s['segment'], 'binary_argument': checks[0],
                             'rational_argument': checks[1]})
        rows.append({'case': case['case'], 'segments': segments})
    assert len(rows) == 27 and sum(len(r['segments']) for r in rows) == 317
    print(json.dumps({'scope': 'FLINT ball enclosure at reconstructed binary and exact-formula arguments and a lower threshold from saved eps_rot. Does not enclose underlying physical-input uncertainty, historical individual orders, or compile the circuit.',
                      'python_flint_version': flint.__version__, 'cases': rows}, indent=2))


if __name__ == '__main__':
    main()
