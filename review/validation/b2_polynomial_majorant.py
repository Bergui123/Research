"""Benchmark B2 bound using the rational polynomial residual certificate.

The all-time coefficient enclosures are saved in polynomial-residual-bound.json.
Use --verify-certificate to rerun the slower exact rational verifier first.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
from pathlib import Path

from b2_global_majorant import ladder_majorant
from b2_riccati_majorant import R, weight
from verify_polynomial_nodes import verify

ROOT=Path(__file__).resolve().parents[2]
VAL=ROOT/'review/validation'


def certificate():
    return json.loads((VAL/'polynomial-residual-bound.json').read_text())


def coefficient_bounds(data):
    out=[{0:[F(1)]}]
    for k in (0,1):
        approx=data['interpolant_degree_l1_bounds'][k]
        error=data['coefficient_error_l1_bounds_at_T'][k]
        expected={str(d) for d in range(0,4*(k+1)+1,2)}
        assert set(approx)==set(error)==expected
        block={}
        for d in sorted(expected,key=int):
            value=F(approx[d])+F(error[d])
            assert value>=0
            block[int(d)]=[value]
        out.append(block)
    return out


def calculate(data):
    parts=[]
    for k,block in enumerate(coefficient_bounds(data)):
        steps=8-2*k
        ladder=ladder_majorant(block,steps,R)
        assert all(d%2==0 and len(poly)==1 for d,poly in ladder.items())
        parts.append(F(4,math.factorial(steps))*sum(poly[0]*weight(d) for d,poly in ladder.items()))
    return parts


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify-certificate',action='store_true')
    args=parser.parse_args()
    data=certificate()
    if args.verify_certificate:
        nodes=json.loads(gzip.decompress((VAL/'polynomial-quad-nodes.json.gz').read_bytes()))
        pn=json.loads((VAL/'riccati-quad-nodes.json').read_text())
        assert verify(nodes,pn)==data
    parts=calculate(data)
    total=sum(parts,F(0))
    assert total < 56
    print(json.dumps({'scope':'Normalized benchmark only, conditional on the Riccati and polynomial residual comparison proofs.',
                      'polynomial_certificate_sha256':hashlib.sha256((VAL/'polynomial-residual-bound.json').read_bytes()).hexdigest(),
                      'parts_rational':[str(v) for v in parts],
                      'parts_decimal':[float(v) for v in parts],
                      'B2_upper_rational':str(total),'B2_upper_decimal':float(total),
                      'B2_strict_integer_upper':56},indent=2))


if __name__=='__main__':main()
