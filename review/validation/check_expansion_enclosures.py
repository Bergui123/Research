"""Bound the terminal expansion at three benchmark stencils using rational arithmetic.

Run verify_polynomial_nodes.py and compare it with polynomial-residual-bound.json first.
"""
import argparse
from fractions import Fraction as F
import gzip
import json
import math
from pathlib import Path
import sys

from verify_riccati_nodes import quadratic
from verify_polynomial_nodes import EPS

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'simulation'))


def exp_negative(x):
    """Rational lower and upper bounds on exp(-x), for 0 <= x < 30."""
    assert 0 <= x < 30
    term = total = F(1)
    n = 35
    for j in range(1, n + 1):
        term *= x / j
        total += term
    next_term = term * x / (n + 1)
    tail_upper = next_term / (1 - x / (n + 2))
    return 1 / (total + tail_upper), 1 / total


def lower(x, denom=10**15):
    return F(x.numerator * denom // x.denominator, denom)


def upper(x, denom=10**15):
    return -lower(-x, denom)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--case-set',choices=('b2-875','b2-56'),default='b2-875')
    args=parser.parse_args()
    nodes = json.loads(gzip.decompress((ROOT / 'review/validation/polynomial-quad-nodes.json.gz').read_bytes()))
    pn = json.loads((ROOT / 'review/validation/riccati-quad-nodes.json').read_text())
    cert = json.loads((ROOT / 'review/validation/polynomial-residual-bound.json').read_text())
    assert nodes['segments'] == pn['segments'] == cert['segments'] == 4000
    n, q, pq = nodes['segments'], nodes['denominator'], pn['denominator']
    integrals = [F(0)] * 3
    for i in range(n):
        for k in range(3):
            a = quadratic(*[F(pn['nodes'][2*i+j][k], pq) for j in range(3)])
            integrals[k] += (a[0] + a[1]/2 + a[2]/3) / n
    terminal_p = [F(v, pq) for v in pn['nodes'][-1]]
    indices = [[tuple(x) for x in block] for block in nodes['indices']]
    poly = [{ij: F(v, q) for ij, v in zip(indices[k], nodes['nodes'][-1][k])} for k in range(2)]
    errors = [{int(d): F(v) for d, v in block.items()} for block in cert['coefficient_error_l1_bounds_at_T']]
    result = []
    cases = [('0.05', F(1, 4), 16), ('0.01', F(4, 33), 16), ('0.001', F(2, 15), 16)]
    if args.case_set == 'b2-56':
        cases = [('0.05', F(1, 2), 8), ('0.01', F(1, 4), 16), ('0.001', F(2, 15), 16)]
    for target, dx, h in cases:
        nu, span = F(8), h * dx
        node_values = []
        for x, y in [(3 + span, F(1)), (3 - span, F(1)), (F(3), 1 + span), (F(3), 1 - span)]:
            r2, delta2 = (x*x + y*y) / nu, dx*dx / nu
            factor, factor_error = F(1), F(0)
            for k in range(2):
                delta_power = delta2**(k + 1)
                value = sum(v * x**i * y**j / nu**((i+j)//2)
                            for (i, j), v in poly[k].items())
                factor += delta_power * value
                factor_error += delta_power * sum(v * r2**(d//2) for d, v in errors[k].items())
            assert factor > factor_error > 0
            exponent = ((terminal_p[0]*x*x + 2*terminal_p[1]*x*y + terminal_p[2]*y*y) / (8*nu)
                        + (integrals[0] + integrals[2]) / 4)
            exponent_error = EPS * (r2 / 8 + F(1, 2))
            assert exponent > exponent_error
            psi_low = exp_negative(exponent + exponent_error)[0]
            psi_high = exp_negative(exponent - exponent_error)[1]
            lo, hi = lower(psi_low * (factor-factor_error)), upper(psi_high * (factor+factor_error))
            node_values.append({'x': str(x), 'y': str(y), 'lower': str(lo), 'upper': str(hi),
                                'lower_decimal': float(lo), 'upper_decimal': float(hi),
                                'width_upper': float(hi-lo), 'factor_error_upper': float(factor_error)})
        result.append({'target': target, 'nu': str(nu), 'dx': str(dx), 'h': h, 'nodes': node_values})
    print(json.dumps({'scope': 'Expansion S at the three exploratory benchmark stencils only; relies on separately verified coefficient and Riccati certificates.', 'cases': result}, indent=2))


if __name__ == '__main__':
    main()
