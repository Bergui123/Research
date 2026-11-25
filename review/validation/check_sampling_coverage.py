"""Independent-seed empirical coverage check for the saved 1% sampling configuration.

Run from the repository root:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_sampling_coverage.py

This is Monte Carlo evidence for the emulated distribution, not a mathematical
95% guarantee; GQSP truncation and gate synthesis are not emulated.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'simulation'))

import numpy as np
from scipy.stats import beta
import step9_sampling as sampling
from model import RESULTS

saved = json.loads((RESULTS / 'step9_sampling.json').read_text())
row = saved['best']['0.01']
cfg = sampling.configuration(0.01, row['n'], row['L'], row['kind'], row['eta_max'], saved['kappa95'])
if cfg is None or cfg['runs'] != row['runs']:
    raise RuntimeError('recreated configuration differs from saved run count')
sampling.RNG = np.random.default_rng(793007)
trials = 10000
r = sampling.validate(cfg, trials=trials)
def interval(fraction):
    successes = round(fraction * trials)
    return [float(beta.ppf(0.025, successes, trials - successes + 1)) if successes else 0.0,
            float(beta.ppf(0.975, successes + 1, trials - successes)) if successes < trials else 1.0]
print(json.dumps({
    'target': row['target'], 'seed': 793007, 'trials': trials,
    'runs_per_trial': row['runs'], 'acceptance_probability': row['P_succ'],
    'uniform_fraction': r['frac_uniform_within'],
    'uniform_binomial_95_percent_interval': interval(r['frac_uniform_within']),
    'at_X_fraction': r['frac_at_X_within'],
    'at_X_binomial_95_percent_interval': interval(r['frac_at_X_within']),
    'uniform_error_empirical_q95': r['q95_uniform'],
}, indent=2))
