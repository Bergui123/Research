"""Rational one-dimensional comparison bounds for three benchmark box responses.

See projected-box-bound.md for the comparison proof and limitations. The
only nonlinear functions below are bounded with rational series; propagation
uses outward-rounded fixed-point integer operations and a uniformization tail.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

from check_boundary_bound import bound_for_a
from check_expansion_enclosures import exp_negative
from verify_polynomial_nodes import EPS

ROOT = Path(__file__).resolve().parents[2]
Q = 10**18
SUBSTEPS = 128
TERMS = 60


def ceil_scaled(x, q=Q):
    return (x.numerator * q + x.denominator - 1) // x.denominator


def ceil_div(n, d):
    return (n + d - 1) // d


def response_upper(anchor, dx, n, nu, h, potential_coeff=F(21,25)):
    """Upper-bound the unit-source 1D response for V=potential_coeff*x^2/(4 nu)."""
    assert potential_coeff >= 0
    period = n * dx
    index = round((anchor + period / 2) / dx)
    assert 0 <= index-h and index+h < n
    points = [anchor + (i-index)*dx for i in range(n)]
    # For a sweep, potential_coeff = lambda*(1-rho^2).
    potential = [potential_coeff*x*x/(4*nu) for x in points]
    rate = nu / dx**2
    lam = math.ceil(2*rate + max(potential))
    diagonal = [ceil_scaled(1 - (2*rate+v)/lam) for v in potential]
    side = ceil_scaled(rate/lam)
    assert min(diagonal) >= 0 and lam > 0
    z = F(lam, SUBSTEPS)
    _, exp_hi = exp_negative(z)
    weights, term = [], F(1)
    for j in range(TERMS+1):
        if j:
            term *= z/j
        weights.append(ceil_scaled(exp_hi*term))
    assert z < TERMS + 2
    tail = ceil_scaled(exp_hi * term*z/(TERMS+1) / (1-z/(TERMS+2)))
    values = [0]*n
    for _ in range(SUBSTEPS):
        power = values[:]
        accum = [weights[0]*v for v in power]
        for weight in weights[1:]:
            power = [min(Q, ceil_div(diagonal[i]*power[i] + side*(
                (power[i-1] if i else Q) + (power[i+1] if i+1 < n else Q)), Q))
                for i in range(n)]
            for i, value in enumerate(power):
                accum[i] += weight * value
        values = [min(Q, ceil_div(v, Q) + tail) for v in accum]
    return {"center": F(values[index], Q), "plus": F(values[index+h], Q),
            "minus": F(values[index-h], Q), "uniformization_rate": lam,
            "tail_units_per_substep": tail, "anchor_index": index,
            "left": points[0], "right": points[-1]}


def round_down(x, q=Q):
    return F(x.numerator*q // x.denominator, q)


def round_up(x, q=Q):
    return -round_down(-x, q)


def log_interval(x):
    """Outward rational log interval via the atanh series, x > 0."""
    assert x > 0
    if x < 1:
        lo, hi = log_interval(1/x)
        return -hi, -lo
    y = (x-1)/(x+1)
    assert 0 <= y < 1
    total, power = F(0), y
    count = 48
    for k in range(count):
        total += 2*power/(2*k+1)
        power *= y*y
    tail = 2*power/((2*count+1)*(1-y*y))
    return round_down(total), round_up(total+tail)


def sqrt_down(x):
    assert x >= 0
    return F(math.isqrt(x.numerator*Q**2//x.denominator), Q)


def sqrt_up(x):
    return sqrt_down(x) + F(1, Q)


def reference_interval(pnodes):
    den = pnodes['denominator']
    p11, p12, p22 = [F(x, den) for x in pnodes['nodes'][-1]]
    # ||P-Phat||_F < EPS; |(E x)_i|/2 <= EPS (|3|+|1|)/2.
    error = 2*EPS
    values = [(3*p11+p12)/2, (3*p12+p22)/2]
    return [(v-error, v+error) for v in values]


def min_abs_interval(interval):
    lo, hi = interval
    return F(0) if lo <= 0 <= hi else min(abs(lo), abs(hi))


def max_abs_interval(interval):
    return max(map(abs, interval))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--b2',type=int,choices=(56,875),default=875)
    parser.add_argument('--case-set',choices=('b2-875','b2-56'),default='b2-875')
    args=parser.parse_args()
    assert args.case_set != 'b2-56' or args.b2 == 56
    if args.b2==56:
        cert=json.loads((ROOT/'review/validation/b2-polynomial-majorant.json').read_text())
        assert F(cert['B2_upper_rational']) < 56
    expansion_file=('expansion-enclosures-56.json' if args.case_set == 'b2-56'
                    else 'expansion-enclosures.json')
    expansions = json.loads((ROOT/'review/validation'/expansion_file).read_text())['cases']
    pn = json.loads((ROOT/'review/validation/riccati-quad-nodes.json').read_text())
    reference = reference_interval(pn)
    ref_lower = sqrt_down(sum(min_abs_interval(v)**2 for v in reference))
    assert ref_lower > 0
    # Rational a and r are inputs to a bound valid for every a>0. They agree
    # with the exact grid-derived minima; no optimizer result is trusted here.
    parameters = [
        ('0.05', F(1,4), 128, 16, F(3,10)),
        ('0.01', F(4,33), 256, 16, F(3,10)),
        ('0.001', F(2,15), 256, 16, F(3,10)),
    ]
    if args.case_set == 'b2-56':
        parameters = [
            ('0.05', F(1,2), 64, 8, F(3,10)),
            ('0.01', F(1,4), 128, 16, F(3,10)),
            ('0.001', F(2,15), 256, 16, F(3,10)),
        ]
    cases = []
    for (target, dx, n, h, a), expansion in zip(parameters, expansions):
        nu, period = F(8), n*dx
        assert expansion['target'] == target and F(expansion['dx']) == dx and expansion['h'] == h
        r1 = response_upper(F(3), dx, n, nu, h)
        r2 = response_upper(F(1), dx, n, nu, h)
        # Min distance of wrapped edge points and first outside points.
        def circ(x):
            y=x % period
            return min(y,period-y)
        radius=min(circ(r['left']) for r in (r1,r2))
        radius=min(radius,*(circ(r['right']) for r in (r1,r2)),
                   *(abs(r['left']-dx) for r in (r1,r2)),
                   *(abs(r['right']+dx) for r in (r1,r2)))
        assert 0 <= radius <= period/2
        ghost=bound_for_a(nu, dx, period, radius, a)
        unit=[r1['plus']+r2['center'], r1['minus']+r2['center'],
              r1['center']+r2['plus'], r1['center']+r2['minus']]
        box=[ghost*u for u in unit]
        lattice=args.b2*(dx*dx/nu)**3
        slo=[F(v['lower']) for v in expansion['nodes']]
        shi=[F(v['upper']) for v in expansion['nodes']]
        assert all(lattice+u < s for u,s in zip(box,slo))
        slacks=[log_interval(s/(s-(lattice+u)))[1] for s,u in zip(slo,box)]
        factor=nu/(h*dx)
        differences=[]
        for component in range(2):
            plus,minus=2*component,2*component+1
            lower=log_interval(slo[plus]/shi[minus])[0]
            upper=log_interval(shi[plus]/slo[minus])[1]
            vs=(-factor*upper,-factor*lower)
            ref=reference[component]
            differences.append(max_abs_interval((vs[0]-ref[1],vs[1]-ref[0])))
        discrepancy=sqrt_up(sum(v*v for v in differences))
        remainder=factor*sqrt_up((slacks[0]+slacks[1])**2+(slacks[2]+slacks[3])**2)
        charge=round_up((discrepancy+remainder)/ref_lower, 10**15)
        cases.append({'target':target,'n':n,'nu':str(nu),'dx':str(dx),'h':h,
                      'boundary_a':str(a),'radius':str(radius),'ghost_upper':str(ghost),
                      'box_upper':[str(round_up(v,10**15)) for v in box],
                      'box_upper_decimal':[float(v) for v in box],
                      'lattice_upper':str(lattice),'S_lower':[str(v) for v in slo],
                      'E_below_S':True,'reference_norm_lower':str(round_down(ref_lower,10**15)),
                      'relative_grid_bound_upper':str(charge),
                      'relative_grid_bound_upper_decimal':float(charge),
                      'half_target':str(F(target)/2),'within_half_target':charge < F(target)/2,
                      'uniformization_rates':[r1['uniformization_rate'],r2['uniformization_rate']],
                      'tail_units_per_substep':[r1['tail_units_per_substep'],r2['tail_units_per_substep']]})
    result={'scope':'Three exploratory benchmark configurations, ideal rational grid; conditional on the preceding boundary, B2 and expansion proofs. Circuit cost and total algorithm error are not established.',
            'fixed_point_denominator':Q,'substeps':SUBSTEPS,'last_poisson_term':TERMS,
            'cases':cases}
    if args.b2 != 875:
        result['B2_upper_applied']=args.b2
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
