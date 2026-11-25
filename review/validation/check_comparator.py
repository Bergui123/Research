"""Check saved comparator polynomials against an independent integer expression.

Run once per parameter pair, for example:
EXEC_RHO=0.6 EXEC_LAM=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_comparator.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'simulation'))

import numpy as np
from model import LAM, RHO, RESULTS, X0, axis_grid
from step8_resources import potential_scale

rows_a = json.loads((RESULTS / 'step8.json').read_text())['best'].values()
rows_b = json.loads((RESULTS / 'step9_sampling.json').read_text())['best'].values()
rows = []
for mode, configurations in (('A', rows_a), ('B', rows_b)):
    for cfg in configurations:
        n, dx, L = 2**cfg['n'], cfg['dx'], cfg['L']
        if mode == 'A':
            x, _ = axis_grid(n, L, X0[0], dx)
            y, _ = axis_grid(n, L, X0[1], dx)
            nu = cfg['nu']
        else:
            x = (np.arange(n) - n // 2) * dx
            y, nu = x, 1.0
        _, k, q = potential_scale(nu, x, y, dx)
        a, b = np.meshgrid(np.rint(4*x/dx).astype(np.int64),
                           np.rint(4*y/dx).astype(np.int64), indexing='ij')
        q_exact = (5*a*a + round(10*RHO)*a*b + 5*b*b) / 16.0
        m = cfg['m']
        actual_top = np.floor(q * 2.0**(m-k))
        expected_top = np.floor(q_exact * 2.0**(m-k))
        rows.append({'mode': mode, 'target': cfg['target'], 'grid': n, 'm': m, 'k': k,
                     'mismatched_values': int(np.count_nonzero(q != q_exact)),
                     'mismatched_comparator_integers': int(np.count_nonzero(actual_top != expected_top)),
                     'strict_scale': bool(q.max() < 2**k),
                     'representable': bool(expected_top.max() < 2**m)})
print(json.dumps({'rho': RHO, 'lambda': LAM, 'configurations': rows}, indent=2))
if any(r['mismatched_values'] or r['mismatched_comparator_integers'] or not r['strict_scale'] or not r['representable'] for r in rows):
    raise SystemExit(1)
