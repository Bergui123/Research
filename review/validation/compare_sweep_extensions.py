"""Compare repeated larger-box selections with saved numerical evidence."""
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
PAIRS=[('sweep-extension-rho06-lam1.json','/tmp/sweep-extension-repeat-rho06-lam1.json'),
       ('sweep-extension-rho08-lam1.json','/tmp/sweep-extension-repeat-rho08-lam1.json'),
       ('high-lambda-extension.json','/tmp/high-lambda-extension-repeat.json')]


def compare(a,b,path,changes):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for key in a:compare(a[key],b[key],path+'/'+key,changes)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i),changes)
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        assert isinstance(b,(int,float)) and not isinstance(b,bool),path
        assert math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(path,a,b)
        if a!=b:changes.append({'path':path,'absolute':abs(a-b),'relative':abs(a-b)/max(abs(a),abs(b),1e-300)})
    else:assert a==b,(path,a,b)


def main():
    changes=[]
    for saved,repeat in PAIRS:
        compare(json.loads((ROOT/'review/validation'/saved).read_text()),json.loads(Path(repeat).read_text()),saved,changes)
    max_abs=max(changes,key=lambda x:x['absolute'],default=None)
    max_rel=max(changes,key=lambda x:x['relative'],default=None)
    print(json.dumps({'cases':len(PAIRS),'identifiers_and_decisions_match':True,
                      'numeric_fields_changed':len(changes),'max_absolute_change':max_abs,
                      'max_relative_change':max_rel,
                      'tolerances':{'relative':1e-10,'absolute':1e-12}},indent=2))


if __name__=='__main__':main()
