"""Exact polynomial residual and B2 majorant for one sweep case.

This uses a parameterized Riccati certificate and the same residual method as
verify_polynomial_nodes.py. See sweep-riccati.md for proof assumptions.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path

from b2_global_majorant import add, evaluate, integrate, ladder_majorant, scale
from sweep_riccati import verify as verify_riccati
from verify_polynomial_nodes import accum, bernstein, quadratic, residuals, round_up


def ladder_error(approx,error,steps,row_quarter,eps):
    a,e=approx,error
    for _ in range(steps):
        new_a,new_e={},{}
        for d in set(a)|set(e):
            if d:
                accum(new_a,d-1,scale(a.get(d,[]),d))
                accum(new_e,d-1,scale(e.get(d,[]),d))
            accum(new_a,d+1,scale(a.get(d,[]),row_quarter))
            accum(new_e,d+1,add(scale(e.get(d,[]),row_quarter),scale(a.get(d,[]),eps/2)))
        a,e=new_a,new_e
    return e


def error_bounds(interpolant,residual,row_quarter,eps):
    blocks=[]
    for k in range(2):
        source={d:[residual[k][d]+d*eps*interpolant[k][d]] for d in residual[k]}
        approx={0:[F(1)]} if k==0 else {d:[v] for d,v in interpolant[0].items()}
        error={} if k==0 else blocks[0]
        for d,value in ladder_error(approx,error,4,row_quarter,eps).items():
            accum(source,d,scale(value,F(1,6)))
        if k==1:
            for d,value in ladder_error({0:[F(1)]},{},6,row_quarter,eps).items():
                accum(source,d,scale(value,F(1,180)))
        block={}
        for d in range(4*(k+1),-1,-2):
            forcing=add(source.get(d,[]),scale(block.get(d+2,[]),(d+2)*(d+1)))
            block[d]=integrate(forcing)
        blocks.append(block)
    return blocks


def verify(data,min_hat):
    ric,pol=data['riccati'],data['polynomial']
    rc=verify_riccati(ric,min_hat)
    n,q,pq=pol['segments'],pol['denominator'],ric['denominator']
    assert n==ric['segments'] and q>0 and len(pol['nodes'])==2*n+1
    indices=[[(i,d-i) for d in range(0,4*k+1,2) for i in range(d+1)] for k in (1,2)]
    assert pol['indices']==[[list(x) for x in block] for block in indices]
    assert all(len(row)==2 and all(len(row[k])==len(indices[k]) for k in range(2)) for row in pol['nodes'])
    assert all(v==0 for block in pol['nodes'][0] for v in block)
    eps=round_up(F(rc['residual_integral_upper']),10**9)
    true_min=min_hat-eps
    assert true_min>0
    row_quarter=round_up(F(rc['true_row_upper'])/4,10**6)
    max_a=[{d:F(0) for d in range(0,4*k+1,2)} for k in (1,2)]
    max_r=[{d:F(0) for d in range(0,4*k+1,2)} for k in (1,2)]
    for i in range(n):
        ps=[quadratic(*[F(ric['nodes'][2*i+j][k],pq) for j in range(3)]) for k in range(3)]
        matrix=[[ps[0],ps[1]],[ps[1],ps[2]]]
        blocks=[]
        for k in range(2):
            block={ij:quadratic(*[F(pol['nodes'][2*i+j][k][z],q) for j in range(3)])
                   for z,ij in enumerate(indices[k])}
            blocks.append(block)
            for d in max_a[k]:
                bound=sum((max(map(abs,bernstein(v))) for ij,v in block.items() if sum(ij)==d),F(0))
                max_a[k][d]=max(max_a[k][d],bound)
        rs=residuals(matrix,blocks,n)
        for k in range(2):
            assert all(sum(ij) in max_r[k] for ij in rs[k])
            for d in max_r[k]:
                bound=sum((max(map(abs,bernstein(v))) for ij,v in rs[k].items() if sum(ij)==d),F(0))
                max_r[k][d]=max(max_r[k][d],bound)
    max_a=[{d:round_up(v,10**9) for d,v in b.items()} for b in max_a]
    max_r=[{d:round_up(v,10**12) for d,v in b.items()} for b in max_r]
    err=error_bounds(max_a,max_r,row_quarter,eps)
    assert all(value>=0 for block in err for poly in block.values() for value in poly)
    terminal=[{d:evaluate(v,F(1)) for d,v in b.items()} for b in err]
    coeff=[{0:[F(1)]}]+[{d:[v+terminal[k][d]] for d,v in max_a[k].items()} for k in (0,1)]
    gaussian=true_min/8
    def weight(d):
        assert d>=0 and d%2==0
        return F(1) if d==0 else (F(360*d,1957)/gaussian)**(d//2)
    parts=[]
    for k,block in enumerate(coeff):
        steps=8-2*k
        ladder=ladder_majorant(block,steps,row_quarter)
        parts.append(F(4,math.factorial(steps))*sum(p[0]*weight(d) for d,p in ladder.items()))
    total=sum(parts,F(0))
    return {'rho':ric['rho'],'lambda':ric['lambda'],'segments':n,'riccati':rc,
            'riccati_error_upper':str(eps),'P_min_eigen_lower':str(true_min),
            'P_row_quarter_upper':str(row_quarter),
            'interpolant_degree_l1_bounds':[{str(d):str(v) for d,v in b.items()} for b in max_a],
            'residual_degree_l1_bounds':[{str(d):str(v) for d,v in b.items()} for b in max_r],
            'coefficient_error_l1_bounds_at_T':[{str(d):str(v) for d,v in b.items()} for b in terminal],
            'B2_upper_rational':str(total),'B2_upper_decimal':float(total),
            'B2_upper_integer':math.floor(total)+1,
            'scope':'All-time rational residual/majorant conditional on stated Riccati and coefficient comparison proofs.'}
