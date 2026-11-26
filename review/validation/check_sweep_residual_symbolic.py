"""Spot-check sweep residual formulas by independent SymPy differentiation."""
from fractions import Fraction as F
import gzip
import json
from pathlib import Path

import sympy as s

from verify_riccati_nodes import quadratic
from verify_polynomial_nodes import residuals

ROOT=Path(__file__).resolve().parents[2]
IDS=['rho0-lam1','rho02-lam1','rho06-lam1','rho08-lam1','rho04-lam025','rho04-lam4']


def sym(value):return s.Rational(value.numerator,value.denominator)
def at(poly,z):return sum((v*z**i for i,v in enumerate(poly)),F(0))


def main():
    x,y=s.symbols('x y'); xs=[x,y];z=F(1,3);out=[]
    for ident in IDS:
        data=json.loads(gzip.decompress((ROOT/f'review/validation/sweep-nodes-{ident}.json.gz').read_bytes()))
        cert=json.loads((ROOT/f'review/validation/sweep-cert-{ident}.json').read_text())
        ric,pol=data['riccati'],data['polynomial']
        n,q,pq=pol['segments'],pol['denominator'],ric['denominator']
        indices=[[tuple(ij) for ij in block] for block in pol['indices']]
        rho,lam=F(ric['rho']),F(ric['lambda'])
        checked=[0,n//2,n-1]
        for i in checked:
            p=[quadratic(*[F(ric['nodes'][2*i+j][k],pq) for j in range(3)]) for k in range(3)]
            mat=[[p[0],p[1]],[p[1],p[2]]]
            blocks=[{ij:quadratic(*[F(pol['nodes'][2*i+j][k][a],q) for j in range(3)])
                     for a,ij in enumerate(indices[k])} for k in range(2)]
            actual=residuals(mat,blocks,n)
            ps=[[sym(at(mat[a][b],z)) for b in range(2)] for a in range(2)]
            pols=[sum(sym(at(v,z))*x**ij[0]*y**ij[1] for ij,v in block.items()) for block in blocks]
            dts=[sum(sym(n*(v[1]+2*z*v[2]))*x**ij[0]*y**ij[1] for ij,v in block.items()) for block in blocks]
            def ladder(f,a,count):
                for _ in range(count):
                    f=s.expand(s.diff(f,xs[a])-sum(ps[a][b]*xs[b] for b in range(2))*f/4)
                return f
            for k in range(2):
                gen=sum(s.diff(pols[k],xs[a],2)-sum(ps[a][b]*xs[b] for b in range(2))*s.diff(pols[k],xs[a])/2 for a in range(2))
                force=sum(ladder(s.Integer(1) if k==0 else pols[0],a,4)/12+
                          (ladder(s.Integer(1),a,6)/360 if k==1 else 0) for a in range(2))
                independent=s.Poly(s.expand(dts[k]-gen-force),x,y)
                got={key:sym(at(v,z)) for key,v in actual[k].items()}
                for key in set(got)|set(independent.monoms()):
                    assert got.get(key,0)==independent.coeff_monomial(x**key[0]*y**key[1])
            c=[[F(1),rho],[rho,F(1)]]
            entries=[]
            for a,b in ((0,0),(0,1),(1,1)):
                dp=n*(mat[a][b][1]+2*z*mat[a][b][2])
                residual=dp-2*lam*c[a][b]+sum(at(mat[a][h],z)*at(mat[h][b],z)/2 for h in (0,1))
                entries.append(abs(residual))
            # The all-interval integral bounds this interval's point value / n.
            assert (entries[0]+2*entries[1]+entries[2])/n <= F(cert['riccati']['residual_integral_upper'])
        out.append({'case':ident,'checked_segments':checked,'local_coordinate':str(z),'exact_residual_match':True})
    print(json.dumps({'scope':'Selected-interval diagnostics, not a proof of every interval.', 'cases':out},indent=2))


if __name__=='__main__':main()
