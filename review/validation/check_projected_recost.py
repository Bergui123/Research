"""Re-cost only three exploratory stencil grids using the projected bounds."""
from fractions import Fraction
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'simulation'))
from model import LAM, RHO, T  # noqa: E402
import step8_resources as resources  # noqa: E402


def main():
    assert (RHO, LAM, T) == (0.4, 1.0, 1.0)
    bounds=json.loads((ROOT/'review/validation/projected-box-bound.json').read_text())['cases']
    old=json.loads((ROOT/'review/validation/b2-riccati-reselection.json').read_text())['case']
    selected=[old['alternative']['0.05']['best'],old['alternative']['0.01']['best'],old['larger_box']['selection']['best']]
    rows=json.loads((ROOT/'simulation/results/step5.json').read_text())['grid']
    profiles=[r for r in json.loads((ROOT/'simulation/results/step7.json').read_text())['rows'] if r['safe_Lp']]
    for r in profiles:
        r['amp_err_env']=max(q['amp_err'] for q in profiles if (q['nu'],q['n'],q['kind'],q['p_a'])
                             ==(r['nu'],r['n'],r['kind'],r['p_a']) and q['eta_max']>=r['eta_max'])
    result=[]
    for bound, previous in zip(bounds,selected):
        n=previous['n'];nu=previous['nu'];h=previous['h'];target=previous['target']
        if bound['target']=='0.001':
            row=old['larger_box']['row']
        else:
            row=next(r for r in rows if r['n']==n and r['nu']==nu)
        sp=next(s for s in row['spans'] if s['h']==h)
        eps=float(Fraction(bound['relative_grid_bound_upper']))
        # Converting the rational upper bound must not round downward.
        if Fraction.from_float(eps)<Fraction(bound['relative_grid_bound_upper']):
            import math
            eps=math.nextafter(eps,math.inf)
        candidates=[]
        for pc in profiles:
            if pc['nu']!=nu or pc['n']!=min(n,6):
                continue
            cfg={'nu':nu,'n':n,'L':row['L'],'dx':row['dx'],'h':h,
                 'eps_grid':eps,'eps_grid_measured':sp['rel_err_v'],
                 'p':sp['p'],'norm_ratio':row['norm_ratio'],'pgrid':pc}
            cost=resources.run_cost(cfg,target)
            if cost is not None:
                candidates.append(cost)
        best=min(candidates,key=lambda c:c['T_total']) if candidates else None
        result.append({'target':bound['target'],'profiles_costed':len(candidates),
                       'previous_grid_charge':previous['eps_grid'],
                       'new_grid_charge':eps,'previous_T_total':previous['T_total'],
                       'previous_qubits':previous['logical_qubits'],'best':best})
    print(json.dumps({'scope':'Re-costing of three specified benchmark stencils using the rational grid bounds, floating-point amplitudes and existing uncompiled cost model. No global grid search.', 'cases':result},indent=2))


if __name__=='__main__':main()
