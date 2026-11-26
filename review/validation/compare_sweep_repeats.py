"""Compare sweep re-selections without mistaking roundoff for a changed result."""
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
IDS=['rho0-lam1','rho02-lam1','rho06-lam1','rho08-lam1','rho04-lam025','rho04-lam4']


def compare(a,b,path,changes):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k,changes)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i),changes)
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        assert isinstance(b,(int,float)) and not isinstance(b,bool),path
        assert math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(path,a,b)
        if a!=b:changes.append((abs(a-b),path))
    else:assert a==b,(path,a,b)


def main():
    changes=[]
    for ident in IDS:
        saved=ROOT/f'review/validation/sweep-reselection-{ident}.json'
        repeated=Path(f'/tmp/sweep-reselection-repeat-{ident}.json')
        compare(json.loads(saved.read_text()),json.loads(repeated.read_text()),ident,changes)
    print(json.dumps({'cases':len(IDS),'identifiers_and_decisions_match':True,
                      'numeric_fields_changed':len(changes),
                      'max_absolute_change':max((v for v,_ in changes),default=0.),
                      'tolerances':{'relative':1e-10,'absolute':1e-12}},indent=2))


if __name__=='__main__':main()
