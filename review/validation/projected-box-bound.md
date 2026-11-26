# Rational box-response comparison for three benchmark stencils

This calculation covers the three exploratory normalized-benchmark grids in
`expansion-enclosures.json`, with `nu = 8`, `rho = 2/5`, `lambda = T = 1`,
`H = I` and the fixed initial matrix in the paper. It is conditional on the
boundary, lattice-remainder and expansion derivations already recorded in
`boundary-bound.md`, `riccati-residual-bound.md` and
`polynomial-residual-bound.md`. These arguments, including the comparison
below, await independent review and are not Lean-formalized. The calculation
does not establish synthesized circuit costs or an end-to-end success rate.
The original production tables and parameter sweeps are unchanged.

## One-dimensional comparison

Write `a = nu/dx^2` and let `X_1, X_2` be independent continuous-time
nearest-neighbour random walks, each jumping in either direction at rate `a`.
Their first exit times from the respective box axes are `tau_1`, `tau_2`, and
`tau = min(tau_1, tau_2)`. The two-dimensional Dirichlet generator has
killing potential

```
V(x) = (x_1^2 + 2 rho x_1 x_2 + x_2^2)/(4 nu).
```

Completing the square shows, for either coordinate,
`V(x) >= v_i(x_i) = (1-rho^2) x_i^2/(4 nu)`. The forced box comparison
solution with constant boundary source `g` is, by the finite-state
Feynman--Kac formula,

```
U(x,T) = g E_x[ exp(-integral_0^tau V(X_s) ds) 1_{tau <= T} ].
```

Partition paths by their first exiting coordinate; enlarge each event to
`tau_i <= T` and replace `V` by `v_i`. The coordinate walks are independent,
so this yields

```
U(x,T) <= g (w_1(x_1,T) + w_2(x_2,T)),
w_i(x_i,T) = E_{x_i}[ exp(-integral_0^{tau_i} v_i(X_i(s)) ds)
                        1_{tau_i <= T} ].
```

The comparison is an upper bound, including corner contributions. It need not
be sharp. The existing all-time rational `g` bounds the periodic edge and the
first outside lattice values, whose difference is bounded by `g` because both
are nonnegative and at most `g`.

## Integer upper enclosure of each `w_i`

The one-dimensional equation has diagonal `-2a-v_i`, off-diagonal `a`,
zero initial data, and source `a` at each endpoint for a unit boundary. Add a
constant coordinate equal to one. If `G` is the augmented generator and
`Lambda >= max_x (2a+v_i(x))`, then `B = I + G/Lambda` is nonnegative and
substochastic. The last row of `B` is one. Thus

```
exp(h G) = exp(-z) sum_{j>=0} z^j B^j/j!,   z=Lambda h,
```

and every exact coordinate stays between zero and one. The script uses 128
equal substeps, truncates after `j=60`, and bounds the omitted tail by the
first omitted term times `1/(1-z/62)`, with a rational upper bound for
`exp(-z)`. It rounds every nonnegative matrix coefficient, Poisson weight,
matrix-vector product, partial sum and tail upward to multiples of `10^-18`.
Clamping an upper bound to one is sound because the exact coordinates are at
most one. Boundary coordinates are kept exactly one. The initial coordinates,
parameters, grid nodes, potential and time step are rational. The resulting
upper bounds are in `projected-box-bound.json`. The boundary parameter
`a=3/10` is chosen in advance; `bound_for_a` encloses its all-time boundary
expression with rational arithmetic, independently of any optimizer.

The independent diagnostic `check_projected_box_numerical.py` computes the full
two-dimensional forced system using SciPy for all three stencils. All 12
floating-point node values lie below their rational comparison bounds. A small
matrix-exponential test also checks the one-dimensional implementation. These
are diagnostics and do not prove the comparison or floating-point roundoff.

## Rate error and re-costing

Use `S` intervals from `expansion-enclosures.json`, `B_2 < 875` from the
Riccati certificate and the new upper bound on `U`. For each of the four nodes,
set `E = 875 (dx^2/nu)^3 + U_upper`; the script checks `E < S_lower` exactly.
It encloses the reference feedback from the terminal rational Riccati node and
its uniform error, and encloses the expansion's log differences using the
atanh series for `log`. Positive series tails and all final endpoints are
rounded outward. Integer square-root bounds give a lower bound on the reference
norm and an upper bound on the discrepancy and the paired log remainder.
Consequently the resulting relative grid bounds, subject to the stated proofs,
are

| target | rational grid upper bound | half target |
| --- | ---: | ---: |
| 5% | 0.007259339468832 | 0.025 |
| 1% | 0.001690797330240 | 0.005 |
| 0.1% | 0.000426924669709 | 0.0005 |

All three bounds cover the saved floating-point measured grid errors. This
establishes neither global grid optimality nor a validated circuit budget.
The earlier exploratory costings used smaller numerical grid charges. They
cannot simply inherit the larger bounds above. `check_projected_recost.py`
re-costs only these three fixed grids, scanning the saved safe auxiliary
profiles using the existing floating-point amplitude inputs and circuit-cost
model. The resulting best estimates are approximately `7.0437e13`,
`7.9327e15`, and `1.2582e17` T gates, with 196, 223 and 231 logical qubits.
No global grid-size search is claimed. The output and previous charges are in
`projected-recost.json`; all numerical cost-model and simulation assumptions
still apply. Each re-costed configuration met its target in 2,000 of 2,000
emulated trials; the one-sided 99.9% lower bound is 0.99655. GQSP truncation
and gate synthesis were not emulated. See `projected-emulation.json`.

## Reproduction

From the repository root, first reproduce the Riccati and polynomial
certificates as documented in the earlier derivations. Then run:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_projected_box_bound.py > /tmp/projected-box-bound.json
cmp /tmp/projected-box-bound.json review/validation/projected-box-bound.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s review/validation -p 'test_projected_box_bound.py' -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_projected_box_numerical.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_projected_recost.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_projected_emulation.py --all
```
