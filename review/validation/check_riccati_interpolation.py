"""Compare the rational interpolant with a separate high-precision matrix exponential.

This is a sampled diagnostic, not the all-time residual proof.
"""

from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp

from verify_riccati_nodes import quadratic

ROOT = Path(__file__).resolve().parents[2]
mp.mp.dps = 60


def as_mp(value):
    q = Fraction(value)
    return mp.mpf(q.numerator) / q.denominator


def main():
    data = json.loads((ROOT / "review/validation/riccati-quad-nodes.json").read_text())
    certificate = json.loads((ROOT / "review/validation/riccati-residual-bound.json").read_text())
    n, den = data["segments"], data["denominator"]
    phi = mp.matrix([[2, mp.mpf(3) / 10], [mp.mpf(3) / 10, 1]])
    generator = mp.matrix([[0, 0, mp.mpf(1) / 2, 0], [0, 0, 0, mp.mpf(1) / 2],
                           [2, mp.mpf(4) / 5, 0, 0], [mp.mpf(4) / 5, 2, 0, 0]])
    max_norm = mp.mpf(0)
    samples = 201
    for i in range(samples):
        time = Fraction(2 * i + 1, 2 * samples)
        segment = min(int(time * n), n - 1)
        z = time * n - segment
        rows = data["nodes"][2 * segment:2 * segment + 3]
        entries = []
        for component in range(3):
            coeffs = quadratic(*(Fraction(row[component], den) for row in rows))
            entries.append(as_mp(sum(coeff * z**k for k, coeff in enumerate(coeffs))))
        interpolant = mp.matrix([[entries[0], entries[1]], [entries[1], entries[2]]])
        propagator = mp.expm(as_mp(time) * generator)
        initial = mp.matrix([[1, 0], [0, 1], [phi[0, 0], phi[0, 1]], [phi[1, 0], phi[1, 1]]])
        yz = propagator * initial
        y = mp.matrix([[yz[0, 0], yz[0, 1]], [yz[1, 0], yz[1, 1]]])
        zz = mp.matrix([[yz[2, 0], yz[2, 1]], [yz[3, 0], yz[3, 1]]])
        difference = zz * y**-1 - interpolant
        norm = mp.sqrt(sum(difference[a, b] ** 2 for a in range(2) for b in range(2)))
        max_norm = max(max_norm, norm)
        assert norm <= as_mp(certificate["residual_integral_upper_rational"])
    print(json.dumps({"samples": samples, "decimal_digits": mp.mp.dps,
                      "maximum_frobenius_difference": mp.nstr(max_norm, 15),
                      "residual_bound": certificate["residual_integral_upper"],
                      "scope": "Sampled high-precision cross-check, not a replacement for the residual proof."}, indent=2))


if __name__ == "__main__":
    main()
