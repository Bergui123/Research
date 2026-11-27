"""Floating-point diagnostics for one rational sweep grid comparison."""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--bounds',required=True,type=Path)
    parser.add_argument('--expansion',required=True,type=Path)
    parser.add_argument('--selection',required=True,type=Path)
    args=parser.parse_args()
    bound=json.loads(args.bounds.read_text());expansion=json.loads(args.expansion.read_text())
    selection=json.loads(args.selection.read_text())
    rho,lam=F(bound['rho']),F(bound['lambda'])
    assert (rho,lam)==(F(expansion['rho']),F(expansion['lambda']))
    if 'case' in selection:
        case=selection['case'];best=case['alternative']['0.01']['best']
        assert (rho,lam)==(F(str(case['rho'])),F(str(case['lambda'])))
    else:
        best=selection['selection']['best']
        assert (rho,lam)==(F(selection['rho']),F(selection['lambda']))
    assert best is not None
    n,nu,dx,h=bound['n'],float(F(bound['nu'])),float(F(bound['dx'])),bound['h']
    assert (n,nu,h)==(2**best['n'],best['nu'],best['h'])
    assert abs(dx-best['dx'])<1e-12 and abs(best['L']-n*dx/2)<1e-12
    os.environ['EXEC_RHO']=str(float(rho));os.environ['EXEC_LAM']=str(float(lam))
    sys.path.insert(0,str(ROOT/'simulation'))
    import numpy as np
    from scipy import sparse
    from scipy.sparse.linalg import expm_multiply
    from model import LAM,RHO,T,axis_grid,grid_operator
    from grid_theorem import expansion as numeric_expansion
    from step5_nu_sweep import neighbours
    assert (RHO,LAM,T)==(float(rho),float(lam),1.)
    for node in expansion['nodes']:
        coord=np.array([float(F(node['x'])),float(F(node['y']))])/np.sqrt(nu)
        value=numeric_expansion(coord,1,dx/np.sqrt(nu),2)
        assert float(F(node['lower']))<value<float(F(node['upper']))
    x1,_=axis_grid(n,n*dx/2,3.,dx);x2,_=axis_grid(n,n*dx/2,1.,dx)
    box,_=grid_operator(x1,x2,dx,nu,periodic=False)
    edge=np.zeros((n,n));edge[0]+=1;edge[-1]+=1;edge[:,0]+=1;edge[:,-1]+=1
    src=sparse.csr_matrix((nu/dx**2*float(F(bound['ghost_upper']))*edge).reshape(-1,1))
    aug=sparse.bmat([[-box,src],[None,sparse.csr_matrix((1,1))]],format='csr')
    start=np.zeros(n*n+1);start[-1]=1
    u=expm_multiply(aug,start)[:-1].reshape(n,n)
    actual=np.array([u[i,j] for i,j in neighbours(x1,x2,h)])
    caps=np.array([float(F(v)) for v in bound['box_upper']])
    assert np.all(actual<caps)
    measured=best['eps_grid_measured']
    assert measured<float(F(bound['relative_grid_bound_upper']))
    print(json.dumps({'rho':str(rho),'lambda':str(lam),'n':n,'nu':str(F(bound['nu'])),'h':h,
                      'floating_box_values':actual.tolist(),'box_upper':caps.tolist(),
                      'measured_rate_error':measured,
                      'relative_grid_bound_upper':bound['relative_grid_bound_upper_decimal'],
                      'within_half_target':bound['within_half_target'],
                      'expansion_values_inside_intervals':4,
                      'scope':'Floating-point corroboration only, not an enclosure.'},indent=2))


if __name__=='__main__':main()
