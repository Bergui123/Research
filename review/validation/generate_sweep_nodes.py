"""Choose rational Riccati and expansion nodes for a selected sweep case.

Floating-point solutions choose candidate interpolants only; correctness comes
from the exact residual checks. Output is deterministic gzip JSON.
"""
import argparse
from fractions import Fraction as F
import gzip
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rho',required=True)
    parser.add_argument('--lambda',dest='lam',required=True)
    parser.add_argument('--segments',type=int,default=256)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    rho,lam=F(args.rho),F(args.lam)
    assert 0<=rho<1 and lam>0 and args.segments>0
    os.environ['EXEC_RHO']=str(float(rho));os.environ['EXEC_LAM']=str(float(lam))
    sys.path.insert(0,str(ROOT/'simulation'))
    import numpy as np
    from grid_theorem import solution
    from model import LAM,RHO,T,riccati_P
    assert (RHO,LAM,T)==(float(rho),float(lam),1.0)
    n,q=args.segments,10**15
    times=np.arange(2*n+1)/(2*n)
    pnodes=[]
    for t in times:
        matrix=riccati_P(t)
        pnodes.append([round(float(matrix[a,b])*q) for a,b in ((0,0),(0,1),(1,1))])
    pnodes[0]=[2*q,3*q//10,q]
    indices=[[(i,d-i) for d in range(0,4*k+1,2) for i in range(d+1)] for k in (1,2)]
    values=solution()(times).reshape(2,9,9,-1)
    polys=[[[round(float(values[k,i,j,t])*q) for i,j in indices[k]]
            for k in range(2)] for t in range(2*n+1)]
    data={'riccati':{'segments':n,'denominator':q,'rho':str(rho),'lambda':str(lam),
                     'components':['P11','P12','P22'],'nodes':pnodes},
          'polynomial':{'segments':n,'denominator':q,'indices':indices,'nodes':polys}}
    raw=json.dumps(data,separators=(',',':')).encode()
    args.output.write_bytes(gzip.compress(raw,mtime=0))
    print(f'{args.output}: {n} intervals, {len(raw)} uncompressed bytes')


if __name__=='__main__':main()
