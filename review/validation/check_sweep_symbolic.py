"""Independent SymPy generating-polynomial check for a sweep B2 certificate."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

x=sp.symbols('x')


def q(value):
    f=Fraction(value)
    return sp.Rational(f.numerator,f.denominator)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate',required=True,type=Path)
    args=parser.parse_args()
    data=json.loads(args.certificate.read_text())
    r=q(data['P_row_quarter_upper'])
    gaussian=q(data['P_min_eigen_lower'])/8
    polynomials=[sp.Integer(1)]
    for approx,error in zip(data['interpolant_degree_l1_bounds'],data['coefficient_error_l1_bounds_at_T']):
        polynomials.append(sum((q(v)+q(error[d]))*x**int(d) for d,v in approx.items()))
    parts=[]
    for k,poly in enumerate(polynomials):
        steps=8-2*k
        for _ in range(steps):poly=sp.expand(sp.diff(poly,x)+r*x*poly)
        total=sp.Integer(0)
        for (d,),coefficient in sp.Poly(poly,x).terms():
            assert d%2==0 and coefficient>=0
            weight=sp.Integer(1) if d==0 else (sp.Rational(360*d,1957)/gaussian)**(d//2)
            total+=coefficient*weight
        parts.append(4*total/sp.factorial(steps))
    assert sum(parts)==q(data['B2_upper_rational'])
    print(json.dumps({'rho':data['rho'],'lambda':data['lambda'],
                      'symbolic_majorant_agrees':True,'B2_upper_rational':str(sum(parts))}))


if __name__=='__main__':main()
