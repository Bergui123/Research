"""Re-cost one fixed sweep stencil using its rational grid-error bound.

The input must be a passing bound for the same grid. The search over saved
auxiliary profiles and the circuit model remain numerical or uncompiled.
"""
import argparse
from fractions import Fraction as F
import json
import math
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--bounds',required=True,type=Path)
    parser.add_argument('--selection',required=True,type=Path)
    args=parser.parse_args()
    bound=json.loads(args.bounds.read_text());selected=json.loads(args.selection.read_text())
    rho,lam=F(bound['rho']),F(bound['lambda'])
    if 'case' in selected:
        case=selected['case'];previous=case['alternative']['0.01']['best']
        assert (rho,lam)==(F(str(case['rho'])),F(str(case['lambda'])))
        assert case['B2_applied']==bound['B2_upper_applied']
        extension=False
    else:
        previous=selected['selection']['best']
        assert (rho,lam)==(F(selected['rho']),F(selected['lambda']))
        assert selected['B2_applied']==bound['B2_upper_applied']
        extension=True
    assert previous is not None and previous['target']==0.01
    assert (bound['n'],float(F(bound['nu'])),bound['h'])==(2**previous['n'],previous['nu'],previous['h'])
    assert abs(float(F(bound['dx']))-previous['dx'])<1e-12
    assert abs(previous['L']-bound['n']*float(F(bound['dx']))/2)<1e-12
    assert bound['E_below_S'] and bound['within_half_target']
    os.environ['EXEC_RHO']=str(float(rho));os.environ['EXEC_LAM']=str(float(lam))
    sys.path.insert(0,str(ROOT/'simulation'))
    import step8_resources as resources
    from model import LAM,RESULTS,RHO,T
    assert (RHO,LAM,T)==(float(rho),float(lam),1.)
    if extension:
        row=selected['row']
    else:
        data=json.loads((RESULTS/'step5.json').read_text())
        row=next(r for r in data['grid'] if (r['nu'],r['n'])==(previous['nu'],previous['n']))
    assert abs(row['dx']-previous['dx'])<1e-12 and abs(row['L']-previous['L'])<1e-12
    sp=next(s for s in row['spans'] if s['h']==previous['h'])
    profiles=[r for r in json.loads((RESULTS/'step7.json').read_text())['rows'] if r['safe_Lp']]
    for r in profiles:
        r['amp_err_env']=max(q['amp_err'] for q in profiles if (q['nu'],q['n'],q['kind'],q['p_a'])
                             ==(r['nu'],r['n'],r['kind'],r['p_a']) and q['eta_max']>=r['eta_max'])
    eps=float(F(bound['relative_grid_bound_upper']))
    if F.from_float(eps)<F(bound['relative_grid_bound_upper']):eps=math.nextafter(eps,math.inf)
    candidates=[]
    for pc in profiles:
        if (pc['nu'],pc['n'])!=(previous['nu'],min(previous['n'],6)):continue
        cfg={'nu':previous['nu'],'n':previous['n'],'L':row['L'],'dx':row['dx'],'h':previous['h'],
             'eps_grid':eps,'eps_grid_measured':sp['rel_err_v'],'p':sp['p'],
             'norm_ratio':row['norm_ratio'],'pgrid':pc}
        cost=resources.run_cost(cfg,0.01)
        if cost is not None:candidates.append(cost)
    best=min(candidates,key=lambda c:c['T_total']) if candidates else None
    print(json.dumps({'scope':'One fixed stencil; rational grid comparison is conditional on its proofs, while profile selection, amplitude inputs and circuit costs are numerical or uncompiled.',
                      'rho':str(rho),'lambda':str(lam),'grid_bound':str(F(bound['relative_grid_bound_upper'])),
                      'previous_grid_charge':previous['eps_grid'],'previous_T_total':previous['T_total'],
                      'profiles_costed':len(candidates),'best':best},indent=2))


if __name__=='__main__':main()
