"""Actual-grid stochastic emulation of one re-costed sweep configuration.

This excludes GQSP truncation and gate synthesis and does not validate their
error budgets. It uses a fixed seed for reproducibility.
"""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--recost',required=True,type=Path)
    parser.add_argument('--seed',type=int,default=20261122)
    args=parser.parse_args()
    data=json.loads(args.recost.read_text());rho,lam=F(data['rho']),F(data['lambda'])
    cfg=data['best'];assert cfg is not None and cfg['target']==0.01
    os.environ['EXEC_RHO']=str(float(rho));os.environ['EXEC_LAM']=str(float(lam))
    sys.path.insert(0,str(ROOT/'simulation'))
    import numpy as np
    from model import LAM,RHO,T
    import step10_validation
    assert (RHO,LAM,T)==(float(rho),float(lam),1.)
    step10_validation.RNG=np.random.default_rng(args.seed)
    result=step10_validation.emulate(cfg)
    assert result['amp_rel_err']<=result['amp_budget']
    assert result['lower_bound']>=0.95
    output={'scope':'Actual-grid emulation with fixed seed; GQSP truncation and gate synthesis are excluded.',
            'rho':str(rho),'lambda':str(lam),'grid':2**cfg['n'],'nu':cfg['nu'],
            'h':cfg['h'],'profile':cfg['kind'],'eta_max':cfg['eta_max'],
            'validation':result}
    if args.seed!=20261122:output['seed']=args.seed
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
