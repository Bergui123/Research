"""Test one larger high-lambda box with the sweep B2 majorant.

Only the specified grid is costed, using temporary inputs. The box evolution,
expansion, amplitudes and circuit costs are numerical or uncompiled.
"""
from fractions import Fraction as F
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
os.environ['EXEC_RHO']='0.4';os.environ['EXEC_LAM']='4'
sys.path[:0]=[str(ROOT/'simulation'),str(ROOT/'review/validation')]
import certify_grid  # noqa: E402
import check_alternative_boundary_search as alternative  # noqa: E402
from model import LAM,RESULTS,RHO,T  # noqa: E402
from step5_nu_sweep import evaluate  # noqa: E402


def main():
    assert (RHO,LAM,T)==(0.4,4.,1.)
    cert=json.loads((ROOT/'review/validation/sweep-cert-rho04-lam4.json').read_text())
    bound=cert['B2_upper_integer']
    assert F(cert['rho'])==F(2,5) and F(cert['lambda'])==4 and F(cert['B2_upper_rational'])<bound
    previous=json.loads((ROOT/'review/validation/sweep-reselection-rho04-lam4.json').read_text())['case']
    assert previous['alternative']['0.01']['best'] is None and previous['B2_applied']==bound
    row=evaluate((16.,8,128/7))
    assert row is not None and abs(row['dx']-1/7)<1e-14
    certify_grid.truncation=alternative.analytic_truncation
    estimates=certify_grid.estimate((row,bound))
    available=[x for x in estimates if x['numerical_error_estimate'] is not None]
    assert all(x['measured']<=x['numerical_error_estimate']*(1+1e-6) for x in available)
    with tempfile.TemporaryDirectory(prefix='research-high-lambda-') as directory:
        temp=Path(directory)
        for name in ('step7.json',):shutil.copyfile(RESULTS/name,temp/name)
        step5=json.loads((RESULTS/'step5.json').read_text())
        step5['grid']=[row]
        (temp/'step5.json').write_text(json.dumps(step5))
        (temp/'grid_estimates.json').write_text(json.dumps({'B2':bound,'rows':estimates}))
        alternative.resources.RESULTS=temp
        selection=alternative.select(0.01)
    assert selection['best'] is not None
    print(json.dumps({'scope':'One specified larger high-lambda box, not a global grid search. B2 and boundary bounds depend on stated proofs; box evolution, expansion, amplitudes and circuit costs remain numerical or uncompiled.',
                      'rho':'2/5','lambda':'4','B2_applied':bound,'row':row,
                      'estimates':estimates,'selection':selection},indent=2))


if __name__=='__main__':main()
