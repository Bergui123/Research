"""Cross-check the saved rational boundary values at 80 decimal digits.

This floating-point check corroborates the exact rational evaluation; it is
not itself the proof of the boundary inequality.
"""

from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]
mp.mp.dps = 80

from check_boundary_bound import cosh_minus_one_upper, exp_negative_upper


def as_mp(value):
    q = Fraction(value)
    return mp.mpf(q.numerator) / q.denominator


def main():
    # Check both rational approximations against separate high-precision
    # evaluations, including zero, and check the periodic inequality on
    # shifted grids, at the origin and at the antipodal cusp.
    for x in ("0", "1/1000", "1/4", "1", "4", "12"):
        assert as_mp(exp_negative_upper(Fraction(x))) >= mp.exp(-as_mp(x))
    for x in ("0", "1/1000", "1/4", "1", "3"):
        assert as_mp(cosh_minus_one_upper(Fraction(x))) >= mp.cosh(as_mp(x)) - 1
    largest_residual = mp.mpf("-inf")
    for n in (4, 8, 16):
        dx, a = mp.mpf(1) / 3, mp.mpf(2) / 5
        period = n * dx
        def periodic_f(x):
            d = x % period
            return (mp.exp(-a * d) + mp.exp(-a * (period - d))) / (1 - mp.exp(-a * period))
        for offset in (mp.mpf(0), mp.mpf(1) / 7, dx / 2):
            for i in range(n):
                x = offset + i * dx
                residual = periodic_f(x + dx) + periodic_f(x - dx) - 2 * mp.cosh(a * dx) * periodic_f(x)
                largest_residual = max(largest_residual, residual)
                assert residual < mp.mpf("1e-70")
    saved = json.loads((ROOT / "review/validation/boundary-bound.json").read_text())
    comparisons = []
    for row in saved["cases"]:
        rho, lam = row["rho"], row["lambda"]
        results = ROOT / "simulation/results"
        if (rho, lam) != (0.4, 1.0):
            results /= f"rho{rho:g}_lam{lam:g}"
        cfg = json.loads((results / "step8.json").read_text())["best"][str(row["target"])]
        nu = mp.mpf(str(cfg["nu"]))
        dx = mp.mpf(4) / round(4 / cfg["dx"])
        period = 2**cfg["n"] * dx
        a, radius = as_mp(row["a"]), as_mp(row["radius"])
        exponent = 400 * nu * a**2 / 191 + 2 * nu / dx**2 * (mp.cosh(a * dx) - 1) - a * radius
        formula = mp.exp(exponent) * (1 + mp.exp(-a * (period - 2 * radius))) / (1 - mp.exp(-a * period))
        gap = as_mp(row["rational_boundary_upper"]) - min(mp.mpf(1), formula)
        assert gap >= 0
        comparisons.append({"rho": rho, "lambda": lam, "target": row["target"],
                            "upper_minus_formula": mp.nstr(gap, 12)})
    print(json.dumps({"precision_decimal_digits": mp.mp.dps,
                      "periodic_grid_sizes": [4, 8, 16],
                      "largest_periodic_residual": mp.nstr(largest_residual, 12),
                      "cases": comparisons}, indent=2))


if __name__ == "__main__":
    main()
