"""Validate the saved readout-A sweep configuration against its actual rounded p-grid.

Run from the repository root, e.g.:
EXEC_RHO=0.6 EXEC_LAM=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python review/validation/check_sweep_readout.py

This validates the combined Schrödingerisation and rounding amplitude error
numerically. It does not emulate GQSP truncation or gate synthesis, and does
not make the finite-grid error a rigorous bound.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'simulation'))

import numpy as np
from model import LAM, RHO, RESULTS, reference_feedback
from step10_validation import feedback, readout_probabilities

cfg = json.loads((RESULTS / 'step8.json').read_text())['best']['0.01']
a, a_ideal = readout_probabilities(cfg)
amp_error = float(np.max(np.abs(np.sqrt(a / a_ideal) - 1)))
budget = 0.2 * cfg['delta_tot']
v = reference_feedback()[0]
rate_error = float(np.linalg.norm(feedback(a, cfg) - v) / np.linalg.norm(v))
print(json.dumps({
    'rho': RHO, 'lambda': LAM, 'target': cfg['target'], 'grid': 2**cfg['n'],
    'nu': cfg['nu'], 'span': cfg['h'], 'profile': cfg['kind'],
    'eta_max': cfg['eta_max'], 'n_p': cfg['n_p'],
    'combined_amplitude_relative_error': amp_error,
    'combined_amplitude_budget': budget,
    'actual_rate_relative_error': rate_error,
    'estimated_grid_error': cfg['eps_grid'],
    'passes_amplitude_budget': amp_error <= budget,
    'passes_rate_target_with_exact_amplitudes': rate_error <= cfg['target'],
}, indent=2))
