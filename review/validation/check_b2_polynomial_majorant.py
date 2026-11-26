"""Corroborate polynomial-based l1 and ladder bounds at 101 sampled times."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'simulation'))
import numpy as np

from b2_global_majorant import ladder_majorant
from b2_polynomial_majorant import R, certificate, coefficient_bounds
from check_b2_global_majorant import norms_by_degree
from grid_theorem import DEG, ladder, pad, polys
from model import riccati_P


def main():
    bounds=coefficient_bounds(certificate())
    max_p_ratio=max_ladder_ratio=0.0
    for t in np.linspace(0,1,101):
        p=riccati_P(t)
        actual=[np.ones((1,1)),*polys(t)]
        for k in (0,1,2):
            caps={d:float(v[0]) for d,v in bounds[k].items()}
            for degree,value in norms_by_degree(actual[k]).items():
                cap=caps.get(degree,0.0)
                assert value <= cap+1e-10
                if k and cap>1e-12:
                    max_p_ratio=max(max_p_ratio,value/cap)
            steps=8-2*k
            transformed={d:float(v[0]) for d,v in ladder_majorant(bounds[k],steps,R).items()}
            for axis in (0,1):
                q=ladder(pad(actual[k],DEG+steps+1),p,axis,steps)
                for degree,value in norms_by_degree(q).items():
                    cap=transformed.get(degree,0.0)
                    assert value<=cap+1e-10
                    if cap>1e-12:
                        max_ladder_ratio=max(max_ladder_ratio,value/cap)
    print(json.dumps({'sample_times':101,'max_polynomial_coefficient_ratio':max_p_ratio,
                      'max_ladder_coefficient_ratio':max_ladder_ratio,
                      'scope':'Floating-point samples only; the all-time claim uses the rational residual certificate.'},indent=2))


if __name__=='__main__':main()
