"""Diagnostic only: compare the rational intervals with the numerical expansion."""
from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'simulation'))
from grid_theorem import expansion  # noqa: E402


def main():
    data = json.loads((ROOT / 'review/validation/expansion-enclosures.json').read_text())
    checked = 0
    for case in data['cases']:
        nu, dx = float(Fraction(case['nu'])), float(Fraction(case['dx']))
        for node in case['nodes']:
            x = np.array([float(Fraction(node['x'])), float(Fraction(node['y']))]) / np.sqrt(nu)
            value = expansion(x, 1, dx / np.sqrt(nu), 2)
            assert float(Fraction(node['lower'])) < value < float(Fraction(node['upper']))
            checked += 1
    print(f'{checked} floating-point values inside rational intervals (diagnostic only)')


if __name__ == '__main__':
    main()
