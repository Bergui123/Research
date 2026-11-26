"""Compile SELECT and check a controlled small walk and inverse as matrices.

This counts the walk component only, not QROM, QFT, GQSP phases or controlled
amplitude-estimation state preparation. See controlled-walk.md.
"""
import json
import sys
from pathlib import Path

import numpy as np
from qualtran.bloqs.mcmt import MultiControlZ

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'simulation'))
from model import X0, axis_grid  # noqa: E402
from select_circuit import Circuit, count, select, set_qubit  # noqa: E402
from step8_resources import W10, potential_scale, tcount  # noqa: E402


def main():
    n, frac, m, dx, nu = 2, 2, 2, 0.5, 2.0
    npts=2**n
    axes=[axis_grid(npts,npts*dx/2,X0[i],dx)[0] for i in range(2)]
    vscale,kbits,q5=potential_scale(nu,*axes,dx)
    offsets=[round(2**frac*x[0]/dx) for x in axes]
    dim=2**(3+m)*npts**2
    states=np.arange(dim)
    anc=states//npts**2
    j1=(states//npts)%npts
    j2=states%npts
    t=(anc>>(m+2))&1
    k1=(anc>>(m+1))&1
    k0=(anc>>m)&1
    r=anc%(2**m)
    a_kin=4*nu/dx**2
    alpha=a_kin+vscale/2
    ry=np.array([[np.sqrt(a_kin/alpha),-np.sqrt(vscale/(2*alpha))],
                 [np.sqrt(vscale/(2*alpha)),np.sqrt(a_kin/alpha)]])
    had=np.array([[1.,1.],[1.,-1.]])/np.sqrt(2)
    prep=ry
    for _ in range(2+m): prep=np.kron(prep,had)
    prep=np.kron(prep,np.eye(npts**2))
    reflection=np.eye(dim)
    reflection[npts**2:,npts**2:]*=-1
    reflected=prep@reflection@prep.T
    q_ref=np.floor(q5*2**(m-kbits)).astype(int)[j1,j2]
    kinetic=t==0
    step=np.where(k0,-1,1)*kinetic
    out_j1=(j1+step*(1-k1))%npts
    out_j2=(j2+step*k1)%npts
    out_anc=anc^(kinetic.astype(int)<<m)
    output=out_anc*npts**2+out_j1*npts+out_j2
    sign=np.where(kinetic | ((t==1)&~(r<q_ref)),-1.,1.)
    ideal=np.zeros((dim,dim))
    ideal[output,states]=sign
    expected_walk=reflected@ideal
    measurements=[]
    for enabled in (0,1):
        c=Circuit(dim)
        load=lambda v,bits:[set_qubit(c,((v>>i)&1).astype(bool)) for i in range(bits)]
        regs=(load(j1,n),load(j2,n),*(set_qubit(c,v.astype(bool)) for v in (t,k1,k0)),load(r,m))
        ctrl=set_qubit(c,np.full(dim,bool(enabled)))
        base=c.live
        select(c,n,frac,offsets,kbits,m,regs,W10,control=ctrl)
        read=lambda reg:sum(c.bits[q].astype(int)<<i for i,q in enumerate(reg))
        newanc=(c.bits[regs[2]].astype(int)<<(m+2))+(c.bits[regs[3]].astype(int)<<(m+1))+(c.bits[regs[4]].astype(int)<<m)+read(regs[5])
        out=newanc*npts**2+read(regs[0])*npts+read(regs[1])
        matrix=np.zeros((dim,dim))
        matrix[out,states]=np.where(c.sign,-1.,1.)
        target=ideal if enabled else np.eye(dim)
        select_error=float(np.max(np.abs(matrix-target)))
        assert select_error == 0 and c.live==base
        controlled_reflection=reflected if enabled else np.eye(dim)
        walk=controlled_reflection@matrix
        walk_target=expected_walk if enabled else np.eye(dim)
        error=float(np.max(np.abs(walk-walk_target)))
        inverse=matrix.T@controlled_reflection.T
        inverse_error=float(np.max(np.abs(inverse-walk_target.T)))
        identity_error=float(np.max(np.abs(inverse@walk-np.eye(dim))))
        assert max(error,inverse_error,identity_error)<1e-12
        measurements.append({'control':enabled,'select_max_error':select_error,
                             'walk_max_error':error,'inverse_max_error':inverse_error,
                             'inverse_times_walk_max_error':identity_error,
                             'select_AND_count':c.ands,'peak_select_work':c.peak-base})
    normal=count(n,frac,m,W10)
    controlled_ands=measurements[1]['select_AND_count']
    assert controlled_ands==normal[0]+2
    actual_reflection=tcount(MultiControlZ((1,)*(3+m)))
    model_reflections=2*tcount(MultiControlZ((1,)*(4+m)))
    assert 4*controlled_ands+actual_reflection <= 4*normal[0]+model_reflections
    print(json.dumps({'n':n,'m':m,'block_dimension':dim,'control_values':measurements,
                      'uncontrolled_select_ANDs':normal[0],
                      'controlled_select_T':4*controlled_ands,
                      'controlled_reflection_T':actual_reflection,
                      'model_select_plus_two_reflections_T':4*normal[0]+model_reflections,
                      'actual_PREP_rotations_per_walk':2,'model_rotations_per_walk':4,
                      'scope':'Small controlled walk component, with ideal rotations and standard AND uncomputation; not the full GQSP/QAE circuit.'},indent=2))


if __name__=='__main__': main()
