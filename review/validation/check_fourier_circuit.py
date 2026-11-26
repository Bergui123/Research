"""Check the half-integer transform and evolution from elementary ideal gates."""
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'simulation'))
from step7_schrodinger import g_hat_fft, in_window, p_offset, profile  # noqa: E402


def circuit_matrix(n, ops):
    """Apply the gate list to all computational basis columns, little endian."""
    size=2**n
    rows=np.arange(size)
    matrix=np.eye(size,dtype=complex)
    for op in ops:
        kind=op[0]
        if kind=='H':
            q=op[1]; a=rows[(rows & (1<<q))==0]; b=a|(1<<q)
            x,y=matrix[a].copy(),matrix[b].copy()
            matrix[a],matrix[b]=(x+y)/np.sqrt(2),(x-y)/np.sqrt(2)
        elif kind=='P':
            q,angle=op[1:]
            matrix[((rows>>q)&1)==1]*=np.exp(1j*angle)
        elif kind=='CP':
            q,r,angle=op[1:]
            matrix[(((rows>>q)&1)&((rows>>r)&1))==1]*=np.exp(1j*angle)
        elif kind=='SWAP':
            q,r=op[1:]
            pair=rows[(((rows>>q)^(rows>>r))&1)==1]
            perm=rows.copy();perm[pair]^=(1<<q)|(1<<r)
            matrix=matrix[perm]
        elif kind=='GLOBAL':
            matrix*=np.exp(1j*op[1])
        else:
            raise ValueError(op)
    return matrix


def inverse_ops(ops):
    return [(op[0],*op[1:-1],-op[-1]) if op[0] in ('P','CP','GLOBAL') else op
            for op in reversed(ops)]


def transform_ops(n):
    size=2**n
    c=(size-1)/2
    ops=[('P',j,2*np.pi*c*2**j/size) for j in range(n)]
    for j in range(n-1,-1,-1):
        ops.append(('H',j))
        for k in range(j-1,-1,-1):
            ops.append(('CP',j,k,-np.pi/2**(j-k)))
    ops.extend(('SWAP',j,n-1-j) for j in range(n//2))
    ops.extend([('P',0,np.pi),('GLOBAL',-np.pi*c)])
    return ops


def evolution_ops(n, energy, duration, l_p):
    size=2**n
    c=(size-1)/2
    d_eta=np.pi/l_p
    return [('GLOBAL',-c*d_eta*energy*duration)]+[
        ('P',j,2**j*d_eta*energy*duration) for j in range(n)]


def main():
    results=[]
    l_p=4.0
    p_a=-0.43
    for n in (2,3,4,5):
        size=2**n
        dp=2*l_p/size
        p0=-l_p+dp*np.arange(size)
        p=p0+p_offset(l_p,dp,p_a)
        eta=(np.arange(size)-(size-1)/2)*np.pi/l_p
        ops=transform_ops(n)
        F=circuit_matrix(n,ops)
        inv=circuit_matrix(n,inverse_ops(ops))
        direct=np.exp(-1j*np.outer(eta,p0))/np.sqrt(size)
        transform_error=float(np.max(np.abs(F-direct)))
        inverse_error=float(np.max(np.abs(inv-direct.conj().T)))
        assert max(transform_error,inverse_error)<1e-12
        samples=profile(p,'smooth_w1')
        fft_error=float(np.max(np.abs(F@samples/np.sqrt(size)-g_hat_fft(samples,size))))
        assert fft_error<1e-12
        assert np.min(np.abs(p-p_a))<1e-12 and np.any(in_window(p,p_a))
        phase_errors=[]
        for energy in (0., 0.37*dp, dp, 2*dp, (size+1)*dp):
            phase=circuit_matrix(n,evolution_ops(n,energy,1.,l_p))
            actual=inv@phase@F
            reference=direct.conj().T@np.diag(np.exp(1j*eta*energy))@direct
            err=float(np.max(np.abs(actual-reference)))
            assert err<1e-12
            phase_errors.append(err)
        # A two-node translation wraps antiperiodically across the right boundary.
        shift=inv@circuit_matrix(n,evolution_ops(n,2*dp,1.,l_p))@F
        expected=np.zeros((size,size),complex)
        for j in range(size):
            unwrapped=j+2
            expected[j,unwrapped%size]=(-1)**(unwrapped//size)
        wrap_error=float(np.max(np.abs(shift-expected)))
        assert wrap_error<1e-12
        window=np.flatnonzero(in_window(p,p_a))
        valid=window[(window+2<size) & (p[window]+2*dp>=0)]
        assert len(valid)>0
        recovered=shift@samples
        recovery_error=float(np.max(np.abs(recovered[valid]-np.exp(-(p[valid]+2*dp)))))
        assert recovery_error<1e-12
        results.append({'n_p':n,'grid_size':size,'transform_max_error':transform_error,
                        'inverse_max_error':inverse_error,'fft_max_error':fft_error,
                        'evolution_max_error':max(phase_errors),'two_node_wrap_max_error':wrap_error,'shifted_window_recovery_max_error':recovery_error,
                        'qft_H':n,'qft_controlled_phases':n*(n-1)//2,
                        'qft_swaps':n//2,'input_phase_gates':n,'output_Z':1})
    radius=2.0
    y=np.linspace(-2*l_p+radius,2*l_p-radius,121)
    wrapped=np.zeros_like(y)
    for r in range(-50,51):
        if r:
            wrapped+=(-1.0)**r*profile(y+2*l_p*r,'smooth_w1')
    wrap_bound=2*np.exp(-radius)/(1-np.exp(-2*l_p))
    assert np.max(np.abs(wrapped))<=wrap_bound
    print(json.dumps({'scope':'Ideal elementary-gate simulation for scalar eigenvalues and shifted profile samples; no phase synthesis or finite-grid alias enclosure.',
                      'wrap_lemma_diagnostic':{'R':radius,'L':l_p,'sampled_y':len(y),
                          'finite_image_terms':100,'max_abs_finite_sum':float(np.max(np.abs(wrapped))),
                          'analytic_bound':float(wrap_bound)},'cases':results},indent=2))


if __name__=='__main__':main()
