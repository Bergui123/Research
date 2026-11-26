"""Re-run the benchmark exploratory selection with B2=56, in temporary files."""
from fractions import Fraction
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import check_alternative_boundary_search as alternative
from model import LAM, RHO, T

ROOT=Path(__file__).resolve().parents[2]


def main():
    assert (RHO,LAM,T)==(0.4,1.0,1.0)
    calculated=json.loads(subprocess.check_output(
        [sys.executable,str(ROOT/'review/validation/b2_polynomial_majorant.py')],text=True))
    saved=json.loads((ROOT/'review/validation/b2-polynomial-majorant.json').read_text())
    assert calculated==saved and Fraction(saved['B2_upper_rational'])<56
    source=alternative.RESULTS
    sampled=json.loads((source/'grid_estimates.json').read_text())['B2']
    with tempfile.TemporaryDirectory(prefix='research-b2-polynomial-') as directory:
        temporary=Path(directory)
        for name in ('step5.json','step7.json','step8.json'):
            shutil.copyfile(source/name,temporary/name)
        data=json.loads((source/'grid_estimates.json').read_text())
        data['B2']=56
        (temporary/'grid_estimates.json').write_text(json.dumps(data))
        alternative.RESULTS=temporary
        result=alternative.case(['0.05','0.01','0.001'],extension=True)
    result['B2_sampled_charge']=sampled
    result['B2_applied']=56
    print(json.dumps({'scope':'Benchmark only; B2=56 and the boundary bound depend on stated proofs; box evolution, expansion and circuit costs remain numerical. This tests the existing grid family plus one specified larger box.',
                      'case':result},indent=2))


if __name__=='__main__':main()
