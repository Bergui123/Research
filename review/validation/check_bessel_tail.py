"""Bound the uncomputed Jacobi--Anger tail for saved resource configurations.

The bound uses rational arithmetic for the binary floating-point argument.
The selected order and the first 399 terms still use SciPy, so the whole
tail is not a validated numerical enclosure. See bessel-tail.md.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.special import jv

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'simulation'))
from step8_resources import jacobi_anger_order


def omitted_upper(x, first):
    """Rational bound for 2 sum_{k>=first}|J_k(x)|, for 0<x<first."""
    assert 0 < x < first and isinstance(first, int)
    # Any rational r>1 is valid; the floating-point choice need not be exact.
    opt = (first + math.sqrt(first**2 - float(x)**2)) / float(x)
    r = F(max(4097, min(8192, math.floor(4096 * opt))), 4096)
    y = (r - 1) / (r + 1)
    log_lower = 2 * sum((y**(2*k+1) / (2*k+1) for k in range(32)), F(0))
    d = first * log_lower - x * (r - 1/r) / 2
    assert d > 0
    # e>2 gives exp(-d) <= 2^-floor(d); no transcendental evaluation.
    return 2 * F(1, 2**(d.numerator // d.denominator)) / (1 - 1/r)


def saved_cases():
    for pattern in ('simulation/results/**/step8.json',
                    'simulation/results/**/step9_sampling.json'):
        for path in sorted(ROOT.glob(pattern)):
            for target, cfg in sorted(json.loads(path.read_text())['best'].items()):
                if cfg is not None:
                    yield f'{path.relative_to(ROOT)}:{target}', cfg
    for path in sorted((ROOT/'review/validation').glob('sweep-recost-*.json')):
        yield str(path.relative_to(ROOT)), json.loads(path.read_text())['best']
    path = ROOT/'review/validation/projected-recost-56.json'
    for case in json.loads(path.read_text())['cases']:
        yield f'{path.relative_to(ROOT)}:{case["target"]}', case['best']


def main():
    records = []
    for name, cfg in saved_cases():
        n_p = cfg['n_p']
        q_run = cfg['walk_steps_per_run']
        # coherent_run stores eps_rot = eps_state/(2*(4*q_run+2*n_p)).
        eps = (cfg['eps_rot'] * (2 * (4*q_run + 2*n_p)) /
               (4*n_p*(n_p+1)))
        d_eta = 2 * cfg['eta_max'] / (2**n_p - 1)
        segments = []
        for j in range(n_p):
            x = cfg['alpha_walk'] * d_eta * 2**(j-1)
            k = jacobi_anger_order(x, eps)
            measured = float(2*np.sum(np.abs(jv(np.arange(k+1, k+400), x))))
            omitted = omitted_upper(F.from_float(x), k+400)
            assert measured < eps and F.from_float(measured) + omitted < F.from_float(eps)
            assert omitted < F(341, 10**18)*F.from_float(eps)
            assert F.from_float(measured) < F(99977, 100000)*F.from_float(eps)
            segments.append({'segment': j, 'x_float': x, 'K': k,
                             'computed_tail_float': measured,
                             'omitted_upper': str(omitted),
                             'omitted_over_budget_upper': float(omitted/F.from_float(eps))})
        assert sum(s['K']+3 for s in segments) == q_run
        records.append({'case': name, 'budget_float': eps, 'segments': segments})
    assert len(records) == 27
    print(json.dumps({'scope': 'Omitted tail bound for the binary floating-point arguments and reconstructed selected orders; first 399 terms, parameters and complete algorithm remain numerical or uncompiled.',
                      'cases': records}, indent=2))


if __name__ == '__main__':
    main()
