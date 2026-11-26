"""Compare repeated selection and emulation with saved benchmark results."""
import argparse
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def compare(a,b,path,changes):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for key in a:
            compare(a[key],b[key],path+'/'+key,changes)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        # With a fixed random seed, all trial-error samples must match exactly.
        if path.endswith('/err_samples'):
            assert a==b,path
        for i,(left,right) in enumerate(zip(a,b)):
            compare(left,right,path+'/'+str(i),changes)
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        assert isinstance(b,(int,float)) and not isinstance(b,bool),path
        assert math.isfinite(a) and math.isfinite(b),path
        assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(path,a,b)
        if a!=b:
            changes.append((abs(a-b),path))
    else:
        assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--selection',type=Path,required=True)
    parser.add_argument('--emulation',type=Path,required=True)
    args=parser.parse_args()
    changes=[]
    for actual,saved in [(args.selection,ROOT/'review/validation/b2-polynomial-reselection.json'),
                         (args.emulation,ROOT/'review/validation/projected-emulation-56.json')]:
        compare(json.loads(actual.read_text()),json.loads(saved.read_text()),actual.name,changes)
    print(json.dumps({'selected_configs_and_trial_samples_match':True,
                      'numeric_fields_changed':len(changes),
                      'max_absolute_change':max((x[0] for x in changes),default=0.0),
                      'tolerances':{'relative':1e-10,'absolute':1e-12}},indent=2))


if __name__=='__main__':main()
