"""Numerical sweep re-selection using a certified B2 majorant and boundary bound.

All modified inputs are temporary; box evolution, expansion and circuit costs
remain numerical. An absent candidate reflects this search and its bound only.
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
    args=parser.parse_args()
    cert=json.loads(args.certificate.read_text())
    rho,lam=F(cert['rho']),F(cert['lambda'])
    os.environ['EXEC_RHO']=str(float(rho));os.environ['EXEC_LAM']=str(float(lam))
    sys.path[:0]=[str(ROOT/'simulation'),str(ROOT/'review/validation')]
    import check_alternative_boundary_search as alternative
    from model import RHO,LAM,T
    assert (RHO,LAM,T)==(float(rho),float(lam),1.)
    bound=cert['B2_upper_integer']
    assert F(cert['B2_upper_rational'])<bound
    source=alternative.RESULTS
    sampled=json.loads((source/'grid_estimates.json').read_text())['B2']
    with tempfile.TemporaryDirectory(prefix='research-sweep-b2-') as directory:
        temp=Path(directory)
        for name in ('step5.json','step7.json','step8.json'):
            shutil.copyfile(source/name,temp/name)
        data=json.loads((source/'grid_estimates.json').read_text())
        data['B2']=bound
        (temp/'grid_estimates.json').write_text(json.dumps(data))
        alternative.RESULTS=temp
        case=alternative.case(['0.01'],extension=False)
    case['B2_sampled_charge']=sampled
    case['B2_applied']=bound
    print(json.dumps({'scope':'One sweep point; rational B2 and analytic boundary bounds depend on the documented proofs. Box evolution, expansion, amplitude inputs and circuit costs remain numerical.',
                      'case':case},indent=2))


if __name__=='__main__':main()
