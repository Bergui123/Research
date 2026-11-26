"""Exact residual enclosure of benchmark p1,p2; see polynomial-residual-bound.md."""
from fractions import Fraction as F
import gzip
import json
import math
from pathlib import Path

from b2_global_majorant import add, evaluate, integrate, scale
from verify_riccati_nodes import bernstein, quadratic, verify as verify_riccati

ROOT = Path(__file__).resolve().parents[2]
R = F(3, 5)
EPS = F(13, 10**9)


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def accum(target, key, value):
    target[key] = add(target.get(key, []), value)


def derivative(q, axis):
    out = {}
    for ij, value in q.items():
        if ij[axis]:
            key = list(ij)
            key[axis] -= 1
            accum(out, tuple(key), scale(value, ij[axis]))
    return out


def multiply_linear(q, matrix, axis, factor):
    out = {}
    for ij, value in q.items():
        for j in range(2):
            key = list(ij)
            key[j] += 1
            accum(out, tuple(key), scale(mul(value, matrix[axis][j]), factor))
    return out


def plus(a, b):
    out = dict(a)
    for key, value in b.items():
        accum(out, key, value)
    return out


def scaled(q, factor):
    return {key: scale(value, factor) for key, value in q.items()}


def ladder(q, matrix, axis, count):
    for _ in range(count):
        q = plus(derivative(q, axis), multiply_linear(q, matrix, axis, -F(1, 4)))
    return q


def generator(q, matrix):
    out = {}
    for axis in range(2):
        d = derivative(q, axis)
        out = plus(out, plus(derivative(d, axis),
                             multiply_linear(d, matrix, axis, -F(1, 2))))
    return out


def round_up(value, denominator):
    return F((value.numerator * denominator + value.denominator - 1) // value.denominator,
             denominator)


def residuals(matrix, blocks, segments):
    one = {(0, 0): [F(1)]}
    out = []
    for k in range(2):
        q = blocks[k]
        residual = {key: [segments * v[1], 2 * segments * v[2]]
                    for key, v in q.items()}
        residual = plus(residual, scaled(generator(q, matrix), -1))
        for axis in range(2):
            residual = plus(residual, scaled(ladder(one if k == 0 else blocks[0],
                                                    matrix, axis, 4), -F(1, 12)))
            if k == 1:
                residual = plus(residual, scaled(ladder(one, matrix, axis, 6), -F(1, 360)))
        out.append(residual)
    return out


def ladder_error(approx, error, steps):
    """Degreewise coefficient-l1 bound for D(P)^steps p - D(Phat)^steps phat."""
    a, e = approx, error
    for _ in range(steps):
        new_a, new_e = {}, {}
        for d in set(a) | set(e):
            if d:
                accum(new_a, d - 1, scale(a.get(d, []), d))
                accum(new_e, d - 1, scale(e.get(d, []), d))
            accum(new_a, d + 1, scale(a.get(d, []), R))
            accum(new_e, d + 1, add(scale(e.get(d, []), R),
                                     scale(a.get(d, []), EPS / 2)))
        a, e = new_a, new_e
    return e


def error_bounds(interpolant, residual):
    """Solve triangular coefficient comparison ODEs as polynomials in time."""
    blocks = []
    for k in range(2):
        source = {d: [residual[k][d] + d * EPS * interpolant[k][d]]
                  for d in residual[k]}
        approx = {0: [F(1)]} if k == 0 else {d: [v] for d, v in interpolant[0].items()}
        error = {} if k == 0 else blocks[0]
        for d, value in ladder_error(approx, error, 4).items():
            accum(source, d, scale(value, F(1, 6)))  # two axes, each 1/12
        if k == 1:
            for d, value in ladder_error({0: [F(1)]}, {}, 6).items():
                accum(source, d, scale(value, F(1, 180)))
        block = {}
        for d in range(4 * (k + 1), -1, -2):
            forcing = add(source.get(d, []),
                          scale(block.get(d + 2, []), (d + 2) * (d + 1)))
            block[d] = integrate(forcing)
        blocks.append(block)
    return blocks


def verify(nodes, riccati):
    assert nodes["segments"] == riccati["segments"] == 4000
    n, denom, pdenom = nodes["segments"], nodes["denominator"], riccati["denominator"]
    assert denom > 0 and pdenom > 0 and len(nodes["nodes"]) == len(riccati["nodes"]) == 2 * n + 1
    indices = [[(i, d - i) for d in range(0, 4 * k + 1, 2)
                for i in range(d + 1)] for k in (1, 2)]
    assert nodes["indices"] == [[list(x) for x in block] for block in indices]
    assert all(len(row) == 2 and all(len(row[k]) == len(indices[k]) for k in range(2))
               for row in nodes["nodes"])
    assert all(v == 0 for block in nodes["nodes"][0] for v in block)
    cert = verify_riccati(riccati)
    assert F(cert["residual_integral_upper_rational"]) < EPS
    max_a = [{d: F(0) for d in range(0, 4 * k + 1, 2)} for k in (1, 2)]
    max_r = [{d: F(0) for d in range(0, 4 * k + 1, 2)} for k in (1, 2)]
    for i in range(n):
        ps = [quadratic(*[F(riccati["nodes"][2*i+j][k], pdenom)
                          for j in range(3)]) for k in range(3)]
        matrix = [[ps[0], ps[1]], [ps[1], ps[2]]]
        blocks = []
        for k in range(2):
            block = {ij: quadratic(*[F(nodes["nodes"][2*i+j][k][z], denom)
                                     for j in range(3)])
                     for z, ij in enumerate(indices[k])}
            blocks.append(block)
            for d in max_a[k]:
                bound = sum((max(map(abs, bernstein(v))) for ij, v in block.items()
                             if sum(ij) == d), F(0))
                max_a[k][d] = max(max_a[k][d], bound)
        rs = residuals(matrix, blocks, n)
        for k in range(2):
            assert all(sum(ij) in max_r[k] for ij in rs[k])
            for d in max_r[k]:
                bound = sum((max(map(abs, bernstein(v))) for ij, v in rs[k].items()
                             if sum(ij) == d), F(0))
                max_r[k][d] = max(max_r[k][d], bound)
    max_a = [{d: round_up(v, 10**9) for d, v in block.items()} for block in max_a]
    max_r = [{d: round_up(v, 10**12) for d, v in block.items()} for block in max_r]
    errors = error_bounds(max_a, max_r)
    terminal = [{d: evaluate(v, F(1)) for d, v in block.items()} for block in errors]
    return {"segments": n, "riccati_uniform_frobenius_error_upper": str(EPS),
            "interpolant_degree_l1_bounds": [{str(d): str(v) for d, v in b.items()} for b in max_a],
            "residual_degree_l1_bounds": [{str(d): str(v) for d, v in b.items()} for b in max_r],
            "coefficient_error_l1_bounds_at_T": [{str(d): str(v) for d, v in b.items()} for b in terminal],
            "coefficient_error_l1_at_T_decimal": [{str(d): float(v) for d, v in b.items()} for b in terminal]}


def main():
    nodes = json.loads(gzip.decompress((ROOT / "review/validation/polynomial-quad-nodes.json.gz").read_bytes()))
    riccati = json.loads((ROOT / "review/validation/riccati-quad-nodes.json").read_text())
    print(json.dumps(verify(nodes, riccati), indent=2))


if __name__ == "__main__":
    main()
