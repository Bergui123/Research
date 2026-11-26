"""Compare independent SymPy spatial differentiation with exact residual code."""
from fractions import Fraction as F
import gzip
import json
from pathlib import Path

import sympy as s

from verify_riccati_nodes import quadratic
from verify_polynomial_nodes import residuals

ROOT = Path(__file__).resolve().parents[2]


def as_sym(value):
    return s.Rational(value.numerator, value.denominator)


def at(poly, z):
    return sum((v * z**i for i, v in enumerate(poly)), F(0))


def main():
    data = json.loads(gzip.decompress((ROOT / "review/validation/polynomial-quad-nodes.json.gz").read_bytes()))
    pnodes = json.loads((ROOT / "review/validation/riccati-quad-nodes.json").read_text())
    n, q, pq = data['segments'], data['denominator'], pnodes['denominator']
    indices = [[tuple(ij) for ij in block] for block in data['indices']]
    x, y = s.symbols('x y')
    xs = [x, y]
    z = F(1, 3)
    for i in (0, 713, 1999, 3999):
        p = [quadratic(*[F(pnodes['nodes'][2*i+j][k], pq) for j in range(3)]) for k in range(3)]
        matrix = [[p[0], p[1]], [p[1], p[2]]]
        blocks = [{ij: quadratic(*[F(data['nodes'][2*i+j][k][a], q) for j in range(3)])
                   for a, ij in enumerate(indices[k])} for k in range(2)]
        numeric = residuals(matrix, blocks, n)
        ps = [[as_sym(at(matrix[a][b], z)) for b in range(2)] for a in range(2)]
        polys = [sum(as_sym(at(v, z)) * x**ij[0] * y**ij[1] for ij, v in block.items()) for block in blocks]
        dts = [sum(as_sym(n * (v[1] + 2 * z * v[2])) * x**ij[0] * y**ij[1]
                   for ij, v in block.items()) for block in blocks]

        def ladder(f, a, count):
            for _ in range(count):
                f = s.diff(f, xs[a]) - sum(ps[a][b] * xs[b] for b in range(2)) * f / 4
                f = s.expand(f)
            return f

        for k in range(2):
            gen = sum(s.diff(polys[k], xs[a], 2) - sum(ps[a][b] * xs[b] for b in range(2))
                      * s.diff(polys[k], xs[a]) / 2 for a in range(2))
            forcing = sum(ladder(s.Integer(1) if k == 0 else polys[0], a, 4) / 12
                          + (ladder(s.Integer(1), a, 6) / 360 if k == 1 else 0)
                          for a in range(2))
            independent = s.Poly(s.expand(dts[k] - gen - forcing), x, y)
            obtained = {key: as_sym(at(value, z)) for key, value in numeric[k].items()}
            for key in set(obtained) | set(independent.monoms()):
                assert obtained.get(key, 0) == independent.coeff_monomial(x**key[0] * y**key[1])
    print(json.dumps({'checked_segments': [0, 713, 1999, 3999], 'local_coordinate': str(z),
                      'result': 'exact agreement with SymPy spatial derivatives'}))


if __name__ == '__main__':
    main()
