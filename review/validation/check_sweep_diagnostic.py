"""Sample-check one sweep certificate against independent numerical ODEs."""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate',required=True,type=Path)
    args=parser.parse_args()
    cert=json.loads(args.certificate.read_text())
    os.environ['EXEC_RHO']=str(float(F(cert['rho'])))
    os.environ['EXEC_LAM']=str(float(F(cert['lambda'])))
    sys.path.insert(0,str(ROOT/'simulation'))
    import numpy as np
    from model import riccati_P,RHO,LAM,T
    from grid_theorem import DEG,ladder,pad,polys
    from b2_global_majorant import ladder_majorant
    from check_b2_global_majorant import norms_by_degree
    assert (RHO,LAM,T)==(float(F(cert['rho'])),float(F(cert['lambda'])),1.)
    lower=float(F(cert['P_min_eigen_lower']))
    row_quarter=F(cert['P_row_quarter_upper'])
    blocks=[{0:[F(1)]}]
    for a,e in zip(cert['interpolant_degree_l1_bounds'],cert['coefficient_error_l1_bounds_at_T']):
        blocks.append({int(d):[F(v)+F(e[d])] for d,v in a.items()})
    min_eig=np.inf;max_row=0.;max_coeff_ratio=0.;max_ladder_ratio=0.
    for t in np.linspace(0,1,101):
        p=riccati_P(t)
        min_eig=min(min_eig,float(np.linalg.eigvalsh(p)[0]))
        max_row=max(max_row,float(np.max(np.sum(np.abs(p),axis=1))))
        assert min_eig>lower-1e-12 and max_row<float(4*row_quarter)+1e-12
        polynomials=[np.ones((1,1)),*polys(t)]
        for k,poly in enumerate(polynomials):
            caps={d:float(v[0]) for d,v in blocks[k].items()}
            for d,v in norms_by_degree(poly).items():
                cap=caps.get(d,0.)
                assert v<=cap+1e-10
                if k and cap>1e-12:max_coeff_ratio=max(max_coeff_ratio,v/cap)
            steps=8-2*k
            ladder_caps={d:float(v[0]) for d,v in ladder_majorant(blocks[k],steps,row_quarter).items()}
            for axis in (0,1):
                q=ladder(pad(poly,DEG+steps+1),p,axis,steps)
                for d,v in norms_by_degree(q).items():
                    cap=ladder_caps.get(d,0.)
                    assert v<=cap+1e-10
                    if cap>1e-12:max_ladder_ratio=max(max_ladder_ratio,v/cap)
    print(json.dumps({'rho':cert['rho'],'lambda':cert['lambda'],'sample_times':101,
                      'observed_min_eigen':min_eig,'observed_max_row':max_row,
                      'max_coefficient_ratio':max_coeff_ratio,'max_ladder_ratio':max_ladder_ratio,
                      'scope':'Floating-point sampled corroboration only.'}))


if __name__=='__main__':main()
