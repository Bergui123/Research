"""Rational terminal expansion enclosures for one selected sweep stencil.

The saved sweep residual certificate must first be reproduced independently.
See sweep-grid-bound.md for assumptions and limitations.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path

from check_expansion_enclosures import exp_negative,lower,upper
from verify_riccati_nodes import quadratic


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--nodes',required=True,type=Path)
    parser.add_argument('--certificate',required=True,type=Path)
    parser.add_argument('--selection',required=True,type=Path)
    args=parser.parse_args()
    raw=args.nodes.read_bytes()
    data=json.loads(gzip.decompress(raw))
    cert=json.loads(args.certificate.read_text())
    selected=json.loads(args.selection.read_text())
    ric,pol=data['riccati'],data['polynomial']
    assert hashlib.sha256(raw).hexdigest()==cert['input_sha256']
    assert (F(ric['rho']),F(ric['lambda']))==(F(cert['rho']),F(cert['lambda']))
    if 'case' in selected:
        case=selected['case']
        assert (F(str(case['rho'])),F(str(case['lambda'])))==(F(cert['rho']),F(cert['lambda']))
        best=case['alternative']['0.01']['best']
        assert case['B2_applied']==cert['B2_upper_integer']
    else:
        assert (F(selected['rho']),F(selected['lambda']))==(F(cert['rho']),F(cert['lambda']))
        assert selected['B2_applied']==cert['B2_upper_integer']
        best=selected['selection']['best']
    assert best is not None and best['target']==0.01
    n,den,pden=pol['segments'],pol['denominator'],ric['denominator']
    assert n==ric['segments']==cert['segments']
    nu=F(str(best['nu']));dx=F(4,round(4/best['dx']));h=best['h'];grid=2**best['n']
    assert abs(float(dx)-best['dx'])<1e-12 and abs(float(grid*dx/2)-best['L'])<1e-12
    integrals=[F(0)]*3
    for i in range(n):
        for k in range(3):
            a=quadratic(*[F(ric['nodes'][2*i+j][k],pden) for j in range(3)])
            integrals[k]+=(a[0]+a[1]/2+a[2]/3)/n
    terminal_p=[F(v,pden) for v in ric['nodes'][-1]]
    indices=[[tuple(x) for x in block] for block in pol['indices']]
    polys=[{ij:F(v,den) for ij,v in zip(indices[k],pol['nodes'][-1][k])} for k in range(2)]
    errors=[{int(d):F(v) for d,v in block.items()} for block in cert['coefficient_error_l1_bounds_at_T']]
    eps=F(cert['riccati_error_upper'])
    span=h*dx; result=[]
    for x,y in [(3+span,F(1)),(3-span,F(1)),(F(3),1+span),(F(3),1-span)]:
        r2=(x*x+y*y)/nu;delta2=dx*dx/nu
        factor,factor_error=F(1),F(0)
        for k in range(2):
            power=delta2**(k+1)
            value=sum(v*x**i*y**j/nu**((i+j)//2) for (i,j),v in polys[k].items())
            factor+=power*value
            factor_error+=power*sum(v*r2**(d//2) for d,v in errors[k].items())
        assert factor>factor_error>=0
        exponent=((terminal_p[0]*x*x+2*terminal_p[1]*x*y+terminal_p[2]*y*y)/(8*nu)
                  +(integrals[0]+integrals[2])/4)
        error=eps*(r2/8+F(1,2))
        assert exponent>error
        lo=lower(exp_negative(exponent+error)[0]*(factor-factor_error))
        hi=upper(exp_negative(exponent-error)[1]*(factor+factor_error))
        assert lo>0 and hi>lo
        result.append({'x':str(x),'y':str(y),'lower':str(lo),'upper':str(hi),
                       'lower_decimal':float(lo),'upper_decimal':float(hi),'width_upper':float(hi-lo)})
    print(json.dumps({'scope':'One selected sweep stencil, conditional on the separately verified Riccati and polynomial certificates.',
                      'rho':str(F(cert['rho'])),'lambda':str(F(cert['lambda'])),'n':grid,
                      'nu':str(nu),'dx':str(dx),'h':h,'nodes_sha256':cert['input_sha256'],
                      'nodes':result},indent=2))


if __name__=='__main__':main()
