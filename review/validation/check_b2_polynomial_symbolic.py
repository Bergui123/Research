"""Independently evaluate the coefficient bound using SymPy polynomials."""
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
x=sp.symbols('x')


def rational(value):
    f=Fraction(value)
    return sp.Rational(f.numerator,f.denominator)


def main():
    data=json.loads((ROOT/'review/validation/polynomial-residual-bound.json').read_text())
    polys=[sp.Integer(1)]
    for approx,error in zip(data['interpolant_degree_l1_bounds'],data['coefficient_error_l1_bounds_at_T']):
        polys.append(sum((rational(a)+rational(error[d]))*x**int(d) for d,a in approx.items()))
    parts=[]
    for k,poly in enumerate(polys):
        steps=8-2*k
        for _ in range(steps):
            poly=sp.diff(poly,x)+sp.Rational(3,5)*x*poly
            poly=sp.expand(poly)
        bound=sp.Integer(0)
        for (degree,),coeff in sp.Poly(poly,x).terms():
            assert degree%2==0 and coeff>=0
            w=sp.Integer(1) if degree==0 else (sp.Rational(360*degree,1957)/sp.Rational(899,8000))**(degree//2)
            bound+=coeff*w
        parts.append(4*bound/sp.factorial(steps))
    expected=json.loads((ROOT/'review/validation/b2-polynomial-majorant.json').read_text())
    for value,saved in zip(parts,expected['parts_rational']):
        assert value==rational(saved)
    assert sum(parts)==rational(expected['B2_upper_rational'])
    print(json.dumps({'symbolic_generating_polynomial_matches_exact_recurrence':True,
                      'parts_rational':[str(p) for p in parts],
                      'B2_upper_rational':str(sum(parts))},indent=2))


if __name__=='__main__':main()
