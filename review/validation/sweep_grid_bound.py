"""Rational projected-box and grid bound for one selected sweep stencil.

Conditional on the saved sweep residual and boundary comparison proofs. See
sweep-grid-bound.md. This does not validate the amplitude or circuit models.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path

from check_boundary_bound import bound_for_a
from check_projected_box_bound import (log_interval,max_abs_interval,min_abs_interval,
                                       response_upper,round_down,round_up,sqrt_down,sqrt_up)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--nodes',required=True,type=Path)
    parser.add_argument('--certificate',required=True,type=Path)
    parser.add_argument('--expansion',required=True,type=Path)
    args=parser.parse_args()
    raw=args.nodes.read_bytes();data=json.loads(gzip.decompress(raw))
    cert=json.loads(args.certificate.read_text()); expansion=json.loads(args.expansion.read_text())
    ric=data['riccati'];rho,lam=F(cert['rho']),F(cert['lambda'])
    assert hashlib.sha256(raw).hexdigest()==cert['input_sha256']==expansion['nodes_sha256']
    assert (F(ric['rho']),F(ric['lambda']))==(rho,lam)==(F(expansion['rho']),F(expansion['lambda']))
    bound=cert['B2_upper_integer'];assert F(cert['B2_upper_rational'])<bound
    n,nu,dx,h=expansion['n'],F(expansion['nu']),F(expansion['dx']),expansion['h']
    assert n>0 and nu>0 and dx>0 and h>0 and len(expansion['nodes'])==4
    span=h*dx
    assert [(F(v['x']),F(v['y'])) for v in expansion['nodes']]==[(3+span,F(1)),(3-span,F(1)),(F(3),1+span),(F(3),1-span)]
    potential_coeff=lam*(1-rho*rho)
    assert potential_coeff>0
    r1=response_upper(F(3),dx,n,nu,h,potential_coeff=potential_coeff)
    r2=response_upper(F(1),dx,n,nu,h,potential_coeff=potential_coeff)
    period=n*dx
    def circ(x):
        y=x%period
        return min(y,period-y)
    radius=min(circ(r['left']) for r in (r1,r2))
    radius=min(radius,*(circ(r['right']) for r in (r1,r2)),
               *(abs(r['left']-dx) for r in (r1,r2)),
               *(abs(r['right']+dx) for r in (r1,r2)))
    assert 0<=radius<=period/2
    choices=[(bound_for_a(nu,dx,period,radius,F(k,1000)),F(k,1000)) for k in range(100,601,25)]
    ghost,a=min(choices)
    unit=[r1['plus']+r2['center'],r1['minus']+r2['center'],
          r1['center']+r2['plus'],r1['center']+r2['minus']]
    box=[ghost*u for u in unit]
    lattice=bound*(dx*dx/nu)**3
    slo=[F(v['lower']) for v in expansion['nodes']]
    shi=[F(v['upper']) for v in expansion['nodes']]
    assert all(lattice+u < s for u,s in zip(box,slo))
    slacks=[log_interval(s/(s-(lattice+u)))[1] for s,u in zip(slo,box)]
    den=ric['denominator'];p11,p12,p22=[F(v,den) for v in ric['nodes'][-1]]
    eps=F(cert['riccati_error_upper']);referr=2*eps
    reference=[((3*p11+p12)/2-referr,(3*p11+p12)/2+referr),
               ((3*p12+p22)/2-referr,(3*p12+p22)/2+referr)]
    ref_lower=sqrt_down(sum(min_abs_interval(v)**2 for v in reference))
    assert ref_lower>0
    factor=nu/span;diff=[]
    for k in range(2):
        plus,minus=2*k,2*k+1
        low=log_interval(slo[plus]/shi[minus])[0]
        high=log_interval(shi[plus]/slo[minus])[1]
        vs=(-factor*high,-factor*low);ref=reference[k]
        diff.append(max_abs_interval((vs[0]-ref[1],vs[1]-ref[0])))
    discrepancy=sqrt_up(sum(v*v for v in diff))
    remainder=factor*sqrt_up((slacks[0]+slacks[1])**2+(slacks[2]+slacks[3])**2)
    charge=round_up((discrepancy+remainder)/ref_lower,10**15)
    print(json.dumps({'scope':'One selected sweep stencil on an ideal rational grid; conditional on the stated boundary, B2, expansion and projected-box proofs.',
                      'rho':str(rho),'lambda':str(lam),'n':n,'nu':str(nu),'dx':str(dx),'h':h,
                      'B2_upper_applied':bound,'boundary_a':str(a),'radius':str(radius),'ghost_upper':str(ghost),
                      'box_upper':[str(round_up(v,10**15)) for v in box],
                      'lattice_upper':str(lattice),'S_lower':[str(v) for v in slo],
                      'E_below_S':True,'reference_norm_lower':str(round_down(ref_lower,10**15)),
                      'relative_grid_bound_upper':str(charge),'relative_grid_bound_upper_decimal':float(charge),
                      'half_target':'1/200','within_half_target':charge<F(1,200),
                      'uniformization_rates':[r1['uniformization_rate'],r2['uniformization_rate']]},indent=2))


if __name__=='__main__':main()
