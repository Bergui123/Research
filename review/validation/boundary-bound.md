# Boundary supersolution for the selected grids

This calculation bounds only the boundary source. The later estimate of the
grid error still uses floating-point values for the lattice remainder constant,
the expansion and the box evolution.

For the benchmark data, `Phi = [[2, 3/10], [3/10, 1]]`, so
`max_i (Phi^-1)[i,i] = 200/191`. By Cauchy--Schwarz in the `Phi` norm,
the initial Gaussian satisfies `u0(x) <= exp(-m x_i^2)` in either coordinate,
where `m = 191/(1600 nu)`. For any `a > 0`, completing the square gives

```
u0(x) <= exp(a^2/(4m)) exp(-a |x_i|).
```

Let `P = N dx` be the period and define
`f_P(x) = sum_{k in Z} exp(-a |x + k P|)`. Both `exp(-a |x|)` on the
infinite lattice and `f_P(x)` on the periodic lattice obey

```
nu D2 f <= lambda_a f,    lambda_a = 2 nu (cosh(a dx) - 1) / dx^2.
```

For the first function this follows by considering whether the two neighbours
cross zero; summing its translates proves the periodic case. Since the
potential is nonnegative, the discrete comparison principle bounds the
solutions by `exp(a^2/(4m) + lambda_a t)` times the respective function.

For each box face, use its perpendicular coordinate. Let `r` be the least
distance to zero among the wrapped face nodes measured on the periodic circle
and the first outside neighbours measured on the line. Then `0 <= r <= P/2`,
and the following bounds their supremum `g_*` for every `0 <= t <= T`:

```
g_* <= min(1, exp(a^2/(4m) + lambda_a T)
            * (exp(-a r) + exp(-a (P-r))) / (1-exp(-a P))).
```

The cap of one uses the maximum principle and `0 <= u0 <= 1`. The result applies
to the real-valued problem with the stated rational parameters, nonnegative
potential, positive `nu`, a square periodic grid and the same initial data on
the infinite lattice. It also applies when the value from the displayed
formula exceeds one.

`check_boundary_bound.py` chooses a positive rational `a`, evaluates an upper
bound on this expression using rational arithmetic, and rounds upward. It uses
`exp(-x) <= 1 / sum_{j=0}^{32} x^j/j!` for `x >= 0` and a geometric bound on
the positive tail of the cosh series. The floating-point optimizer only chooses
`a`; the inequality holds for every positive rational choice. An 80-decimal
evaluation in `check_boundary_bound_precision.py` independently checks the
formula's arithmetic but is not the proof of the inequality.

For the existing 0.1% benchmark, the rational upper value is 0.176100724357,
where the current pipeline charges a sampled value of 0.057447301062. Using
the former with the existing numerical `B2` and floating-point box evolution
gives an alternative grid-error estimate of 0.000909983. This exceeds that
configuration's 0.0005 grid allocation. It does not establish that the actual
error exceeds the allocation, nor that another configuration cannot meet it.
The production resource pipeline still uses the sampled value. A separate
temporary-file search with the conservative bound is recorded in
`alternative-boundary-selection.json`; it also tests one larger box. Its
remaining grid quantities are numerical estimates, and the production counts
have not been replaced.

Reproduce from the repository root:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 .venv/bin/python review/validation/check_boundary_bound.py --all \
 > review/validation/boundary-bound.json
.venv/bin/python review/validation/check_boundary_bound_precision.py \
 > review/validation/boundary-bound-precision.json
```
