"""Count controlled SELECT at all previously selected benchmark/sweep sizes.

This runs the gate builder on one basis input, which suffices for static counts;
it is not an exhaustive correctness test at those sizes.
"""
import json
from pathlib import Path
import sys

import numpy as np
from qualtran.bloqs.mcmt import MultiControlZ

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'simulation'))
from select_circuit import Circuit, count, select, set_qubit  # noqa: E402
from step8_resources import tcount  # noqa: E402


def controlled_count(n,frac,m,w10):
    c=Circuit(1)
    reg=lambda k:[set_qubit(c,np.zeros(1,bool)) for _ in range(k)]
    regs=(reg(n),reg(n),*reg(3),reg(m))
    control=reg(1)[0]
    base=c.peak=c.live
    b=n+frac+1
    kbits=max(m-2*frac,2*b+3-2*frac)
    select(c,n,frac,(0,0),kbits,m,regs,w10,control=control)
    assert c.live==base and not c.sign.any()
    return c.ands,c.peak-base


def main():
    results=ROOT/'simulation/results'
    cases=[('benchmark',4,results)]
    for rho,lam in [('0','1'),('0.2','1'),('0.6','1'),('0.8','1'),('0.4','0.25'),('0.4','4')]:
        cases.append((f'rho={rho},lambda={lam}',round(10*float(rho)),results/f'rho{rho}_lam{lam}'))
    rows=[]
    for label,w10,directory in cases:
        for row in json.loads((directory/'compile_select.json').read_text()):
            n,m=row['n'],row['m']
            frac=2 if row['readout']=='A' else 0
            normal=count(n,frac,m,w10)
            controlled,work=controlled_count(n,frac,m,w10)
            assert 4*normal[0]==row['T_select'] and normal[2]==row['work']
            assert controlled==normal[0]+2
            reflected=tcount(MultiControlZ((1,)*(3+m)))
            model_refl=2*tcount(MultiControlZ((1,)*(4+m)))
            component=4*controlled+reflected
            allowance=4*normal[0]+model_refl
            assert component<=allowance
            rows.append({'case':label,'readout':row['readout'],'target':row['target'],
                         'n':n,'m':m,'uncontrolled_ANDs':normal[0],
                         'controlled_ANDs':controlled,'uncontrolled_work':normal[2],
                         'controlled_work':work,'T_component':component,
                         'T_model_allowance':allowance})
    assert len(rows)==18
    print(json.dumps({'scope':'Static gate-builder counts on one basis input; not a correctness or synthesis proof at these sizes.',
                      'max_extra_work':max(r['controlled_work']-r['uncontrolled_work'] for r in rows),
                      'min_T_slack':min(r['T_model_allowance']-r['T_component'] for r in rows),
                      'rows':rows},indent=2))


if __name__=='__main__':main()
