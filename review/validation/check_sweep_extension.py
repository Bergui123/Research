"""Numerically test one specified larger box for one parameter sweep.

The saved B2 majorant and boundary comparison are analytic, subject to proof
review; box evolution, expansion, amplitude inputs and circuit costs are
numerical or uncompiled. Only this grid is costed, in a temporary directory.
"""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate',required=True,type=Path)
    parser.add_argument('--nu',required=True)
    parser.add_argument('--n',required=True,type=int)
    parser.add_argument('--half-width',required=True)
    args=parser.parse_args()
    cert=json.loads(args.certificate.read_text());rho,lam=F(cert['rho']),F(cert['lambda'])
    nu,half=F(args.nu),F(args.half_width)
    assert nu>0 and half>0 and args.n>0
    os.environ['EXEC_RHO']=str(float(rho));os.environ['EXEC_LAM']=str(float(lam))
    sys.path[:0]=[str(ROOT/'simulation'),str(ROOT/'review/validation')]
    import certify_grid
    import check_alternative_boundary_search as alternative
    from model import LAM,RESULTS,RHO,T
    from step5_nu_sweep import evaluate
    assert (RHO,LAM,T)==(float(rho),float(lam),1.)
    bound=cert['B2_upper_integer'];assert F(cert['B2_upper_rational'])<bound
    row=evaluate((float(nu),args.n,float(half)))
    assert row is not None
    certify_grid.truncation=alternative.analytic_truncation
    estimates=certify_grid.estimate((row,bound))
    available=[v for v in estimates if v['numerical_error_estimate'] is not None]
    assert all(v['measured']<=v['numerical_error_estimate']*(1+1e-6) for v in available)
    with tempfile.TemporaryDirectory(prefix='research-sweep-extension-') as directory:
        temp=Path(directory);shutil.copyfile(RESULTS/'step7.json',temp/'step7.json')
        step5=json.loads((RESULTS/'step5.json').read_text());step5['grid']=[row]
        (temp/'step5.json').write_text(json.dumps(step5))
        (temp/'grid_estimates.json').write_text(json.dumps({'B2':bound,'rows':estimates}))
        alternative.resources.RESULTS=temp
        selection=alternative.select(0.01)
    print(json.dumps({'scope':'One specified larger box, not a global grid search. B2 and boundary bounds depend on the stated proofs; box evolution, expansion, amplitudes and circuit costs remain numerical or uncompiled.',
                      'rho':str(rho),'lambda':str(lam),'B2_applied':bound,
                      'requested_nu':str(nu),'requested_n':args.n,'requested_half_width':str(half),
                      'row':row,'estimates':estimates,'selection':selection},indent=2))


if __name__=='__main__':main()
