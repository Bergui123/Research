"""Enclose each reconstructed Bessel tail on a small argument interval.

See bessel-perturbation.md for the derivation and numerical assumptions.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys

import flint
from flint import arb, fmpq

from check_bessel_interval import enclose_tail
from check_bessel_tail import omitted_upper, saved_cases

ROOT = Path(__file__).resolve().parents[2]


def input_interval(cfg, j, binary):
    """Contain binary result and formula with inputs one adjacent float apart."""
    alpha, eta, n_p = cfg['alpha_walk'], cfg['eta_max'], cfg['n_p']
    assert alpha > 0 and eta > 0 and 0 <= j < n_p
    scale = F(2**j, 2**n_p-1)
    lower = F.from_float(math.nextafter(alpha, -math.inf))*F.from_float(math.nextafter(eta, -math.inf))*scale
    upper = F.from_float(math.nextafter(alpha, math.inf))*F.from_float(math.nextafter(eta, math.inf))*scale
    center = F.from_float(binary)
    return min(lower, center), max(upper, center)


def derivative_upper(xmax, first):
    """Bound the Lipschitz constant of 2 sum_{k>=first}|J_k(x)| on [0,xmax]."""
    assert first >= 2 and 0 < xmax < first-1
    return (omitted_upper(xmax, first-1) + omitted_upper(xmax, first+1))/2


def threshold_lower(cfg):
    n_p, q_run = cfg['n_p'], cfg['walk_steps_per_run']
    assert cfg['eps_rot'] > sys.float_info.min
    return (F.from_float(math.nextafter(cfg['eps_rot'], -math.inf)) *
            (2*(4*q_run+2*n_p)) / (4*n_p*(n_p+1)) * (1-F(1, 2**52)))


def main():
    cases = dict(saved_cases())
    data = json.loads((ROOT/'review/validation/bessel-tail.json').read_text())
    output = []
    for case in data['cases']:
        cfg = cases[case['case']]
        threshold = threshold_lower(cfg)
        assert 0 < threshold <= F.from_float(case['budget_float'])
        segments = []
        assert len(case['segments']) == cfg['n_p']
        assert sum(s['K']+3 for s in case['segments']) == cfg['walk_steps_per_run']
        for j, segment in enumerate(case['segments']):
            assert segment['segment'] == j
            center = F.from_float(segment['x_float'])
            # Reproduce the original floating-point expression before using
            # the saved center and the reconstructed K as checker inputs.
            actual = cfg['alpha_walk']*(2*cfg['eta_max']/(2**cfg['n_p']-1))*2**(j-1)
            assert center == F.from_float(actual)
            lo, hi = input_interval(cfg, j, segment['x_float'])
            radius = max(center-lo, hi-center)
            correction = radius*derivative_upper(hi, segment['K']+1)
            assert 0 <= correction < threshold
            result = enclose_tail(center, segment['K'], threshold-correction)
            bound = (result['upper']+arb(fmpq(correction.numerator, correction.denominator))).upper()
            ratio = (bound/arb(fmpq(threshold.numerator, threshold.denominator))).upper()
            assert ratio < arb(fmpq(999821, 1000000))
            segments.append({'segment': j, 'interval_lower': str(lo), 'interval_upper': str(hi),
                             'threshold_lower': str(threshold), 'correction_upper': str(correction),
                             'ratio_upper': str(ratio.fmpq()), 'ratio_decimal_approx': float(ratio),
                             'terms_enclosed': result['terms']})
        output.append({'case': case['case'], 'segments': segments})
    assert len(output) == 27 and sum(len(r['segments']) for r in output) == 317
    print(json.dumps({'scope': 'Entire interval between computed binary argument and formula with saved alpha and eta varied by one adjacent float; conservative reconstructed threshold. Other input and historical-order uncertainties and the circuit are not validated.',
                      'python_flint_version': flint.__version__, 'cases': output}, indent=2))


if __name__ == '__main__':
    main()
