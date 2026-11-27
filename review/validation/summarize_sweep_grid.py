"""Cross-check six sweep grid bounds, re-costings and two emulation seeds."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[2]
VAL=ROOT/'review/validation'
CASES=[('rho0-lam1',False),('rho02-lam1',False),('rho06-lam1',True),
       ('rho08-lam1',True),('rho04-lam025',False),('rho04-lam4',False)]


def read(name):return json.loads((VAL/name).read_text())


def main():
    rows=[]
    for ident,extended in CASES:
        nodes=VAL/f'sweep-nodes-{ident}.json.gz'
        cert=read(f'sweep-cert-{ident}.json')
        prefix='extended-' if extended else ''
        expansion=read(f'sweep-expansion-{prefix}{ident}.json')
        bound=read(f'sweep-grid-bound-{prefix}{ident}.json')
        recost=read(f'sweep-recost-{ident}.json');best=recost['best']
        assert best is not None and best['target']==0.01
        assert hashlib.sha256(nodes.read_bytes()).hexdigest()==cert['input_sha256']==expansion['nodes_sha256']
        assert (F(cert['rho']),F(cert['lambda']))==(F(bound['rho']),F(bound['lambda']))==(F(recost['rho']),F(recost['lambda']))
        assert (bound['n'],F(bound['nu']),F(bound['dx']),bound['h'])==(expansion['n'],F(expansion['nu']),F(expansion['dx']),expansion['h'])
        assert (bound['n'],float(F(bound['nu'])),bound['h'])==(2**best['n'],best['nu'],best['h'])
        assert abs(float(F(bound['dx']))-best['dx'])<1e-12 and abs(best['L']-bound['n']*float(F(bound['dx']))/2)<1e-12
        assert bound['B2_upper_applied']==cert['B2_upper_integer']
        assert F(cert['B2_upper_rational'])<cert['B2_upper_integer']
        assert bound['E_below_S'] and bound['within_half_target'] and F(bound['relative_grid_bound_upper'])<F(1,200)
        assert F(recost['grid_bound'])==F(bound['relative_grid_bound_upper'])
        assert F.from_float(best['eps_grid'])>=F(bound['relative_grid_bound_upper'])
        emus=[];sample_sets=[]
        for seed,name in [(20261122,f'sweep-emulation-{ident}.json'),
                          (20261123,f'sweep-emulation-independent-{ident}.json')]:
            emu=read(name);v=emu['validation']
            assert (F(emu['rho']),F(emu['lambda']),emu['grid'],emu['nu'],emu['h'])==(F(bound['rho']),F(bound['lambda']),bound['n'],best['nu'],best['h'])
            assert emu['profile']==best['kind'] and emu['eta_max']==best['eta_max']
            assert v['target']==0.01 and v['trials']==2000 and len(v['err_samples'])==2000
            successes=sum(error<=0.01 for error in v['err_samples'])
            assert v['frac_within_target']==successes/2000
            assert abs(v['lower_bound']-float(beta.ppf(0.001,successes,2001-successes)))<1e-14
            assert v['lower_bound']>=0.95 and v['amp_rel_err']<=v['amp_budget']
            sample_sets.append(v['err_samples'])
            if seed==20261123:assert emu['seed']==seed
            emus.append({'seed':seed,'trials':v['trials'],'fraction_within_target':v['frac_within_target'],
                         'one_sided_999_lower':v['lower_bound'],'amp_rel_err':v['amp_rel_err'],
                         'amp_budget':v['amp_budget']})
        assert sample_sets[0]!=sample_sets[1]
        rows.append({'case':ident,'rho':bound['rho'],'lambda':bound['lambda'],
                     'B2_upper':cert['B2_upper_integer'],'nu':bound['nu'],'n':bound['n'],
                     'dx':bound['dx'],'h':bound['h'],'grid_bound':bound['relative_grid_bound_upper'],
                     'grid_bound_decimal':bound['relative_grid_bound_upper_decimal'],
                     'T_total':best['T_total'],'logical_qubits':best['logical_qubits'],
                     'profile':best['kind'],'eta_max':best['eta_max'],'emulations':emus})
    print(json.dumps({'scope':'Six fixed 1% sweep candidates. Grid bounds depend on the stated rational comparison proofs; resource costs and amplitude inputs remain numerical or uncompiled. No global optimum is asserted.',
                      'cases':rows},indent=2))


if __name__=='__main__':main()
