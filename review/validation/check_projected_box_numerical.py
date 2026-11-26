"""Floating-point diagnostic of the rational box and rate-error bounds."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'simulation'))
from model import LAM, RHO, T, axis_grid, grid_operator  # noqa: E402
from step5_nu_sweep import neighbours  # noqa: E402


def main():
    assert (RHO, LAM, T) == (0.4, 1.0, 1.0)
    data = json.loads((ROOT/'review/validation/projected-box-bound.json').read_text())
    saved = json.loads((ROOT/'review/validation/b2-riccati-reselection.json').read_text())['case']
    selections = [saved['alternative']['0.05']['best'], saved['alternative']['0.01']['best'],
                  saved['larger_box']['selection']['best']]
    out = []
    for case, selection in zip(data['cases'], selections):
        n, nu, dx, h = case['n'], float(F(case['nu'])), float(F(case['dx'])), case['h']
        x1, _ = axis_grid(n,n*dx/2,3.,dx)
        x2, _ = axis_grid(n,n*dx/2,1.,dx)
        box, _ = grid_operator(x1,x2,dx,nu,periodic=False)
        edge = np.zeros((n,n))
        edge[0] += 1; edge[-1] += 1; edge[:,0] += 1; edge[:,-1] += 1
        source = sparse.csr_matrix((nu/dx**2*float(F(case['ghost_upper']))*edge).reshape(-1,1))
        aug = sparse.bmat([[-box,source],[None,sparse.csr_matrix((1,1))]],format='csr')
        start=np.zeros(n*n+1);start[-1]=1
        u=expm_multiply(aug,start)[:-1].reshape(n,n)
        actual=np.array([u[i,j] for i,j in neighbours(x1,x2,h)])
        upper=np.array([float(F(v)) for v in case['box_upper']])
        assert np.all(actual < upper)
        measured=selection['eps_grid_measured']
        assert measured < float(F(case['relative_grid_bound_upper']))
        out.append({'target':case['target'],'floating_box_values':actual.tolist(),
                    'box_upper':upper.tolist(),'measured_rate_error':measured,
                    'relative_grid_bound_upper':case['relative_grid_bound_upper_decimal']})
    print(json.dumps({'scope':'Floating-point diagnostic, not an independent proof of the exact bound.', 'cases':out},indent=2))


if __name__ == '__main__':
    main()
