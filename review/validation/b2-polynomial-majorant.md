# Sharper benchmark lattice-remainder bound from polynomial certificates

This calculation applies only to the normalized benchmark in
`riccati-residual-bound.md` and `polynomial-residual-bound.md`: the fixed
terminal matrix, `rho=2/5`, `lambda=T=1`, and identity impact matrix. It
inherits their Riccati, coefficient-comparison and grid-expansion assumptions.
Those proofs and this deduction await independent review; they are not Lean
formalizations. The result does not by itself validate the amplitude or
circuit-cost models, the six parameter sweeps, or the production selections.

## Bound

The rational polynomial certificate bounds the coefficient-l1 norm of the
degree-`d` part of each interpolant by `A[k,d]`, and its true interpolation
error by a nonnegative comparison polynomial `W[k,d](t)`. All coefficients of
that comparison polynomial are nonnegative, so for every `0<=t<=1`,

```
||p[k,d](t)||_1 <= A[k,d] + W[k,d](1) = C[k,d].
```

The certified Riccati estimates give `||P||_infinity < 12/5` and
`P >= (899/1000) I`; the Gaussian determinant prefactor is at most one.
Writing `R=3/5` and `a=899/8000`, a single directional ladder operator
`D_a = partial_a - (P_a . xi)/4` transforms a coefficient-l1 degree bound
according to

```
(T C)[d-1] += d C[d],      (T C)[d+1] += R C[d].
```

For even `d`, the Gaussian radial weight is bounded using the rational
inequality `e > 1957/720`:

```
sup_{r>=0} r^d exp(-a r^2)
    <= ((360 d/1957)/a)^(d/2) = G[d],     G[0]=1.
```

The grid theorem has `m_k=8-2k` derivatives at order two. Bounding each
coefficient by `C[k,d]`, applying the ladder recurrence `m_k` times, summing
both coordinate directions, and integrating the time-independent bound over
`[0,1]` gives the exact rational expression

```
B2 <= sum_{k=0}^2 (4/m_k!) sum_d (T^m_k C[k])[d] G[d]
   = 14438613395725666200588957701465681808568920171762911100526401 /
     259484342549658374041982943826574158928750000000000000000000
   < 56.
```

The three terms are approximately `0.176319`, `5.777394`, and `49.689774`,
for a total `55.6434861`. This tightens the previous coefficient-majorant
bound `<875` using the already verified interpolant bounds, without assuming
that sampled maxima enclose unsampled times. The numerical mesh estimate used
by the production pipeline is much smaller and remains an estimate.

`b2_polynomial_majorant.py --verify-certificate` reruns the exact 4,000-piece
polynomial and Riccati residual checks and then evaluates the rational bound.
The result reproduced byte for byte. An independent SymPy generating-polynomial
calculation applies `D + R x` repeatedly and obtains the same exact rational
terms. A 101-time floating-point diagnostic checks the polynomial and ladder
coefficient inequalities; it is corroboration, not the all-time proof.

## Selection, enclosed grid bounds and re-costing

`check_b2_polynomial_reselection.py` repeats the benchmark resource search in
temporary files with `B2=56` and the analytic boundary charge. The existing
grid family has no 0.1% candidate under that numerical search; it additionally
tests the previously considered `256^2` larger box. It uses floating-point
box evolution, expansion values and resource inputs, and does not prove a
global optimum over box sizes or grids outside the searched family.

For the resulting three fixed stencils, `check_expansion_enclosures.py`
with `--case-set b2-56` encloses the 12 `S` values with the earlier rational
coefficient certificate. `check_projected_box_bound.py --b2 56 --case-set
b2-56` uses the projected one-dimensional killed-walk comparison in
`projected-box-bound.md` to bound the box response. It verifies `E<S` at every
node and obtains, subject to those proof assumptions:

| target | grid, spacing and span | relative grid upper bound | half target |
| --- | --- | ---: | ---: |
| 5% | `64^2`, `dx=1/2`, `h=8` | `0.022381564004643` | `0.025` |
| 1% | `128^2`, `dx=1/4`, `h=16` | `0.003253219793792` | `0.005` |
| 0.1% | `256^2`, `dx=2/15`, `h=16` | `0.000320128205863` | `0.0005` |

A separate full two-dimensional floating-point solve lies below each of the
12 projected box bounds, and the saved measured grid errors lie below these
relative bounds. This diagnostic does not prove the comparison.

The numerical search used smaller grid charges than the bounds in the table.
`check_projected_recost.py` therefore re-costs only these three fixed stencils,
scanning saved safe auxiliary profiles and using the rational bounds rounded
upward. It rejects bounds for mismatched grids. The result is approximately
`1.7588e13`, `1.0440e15`, and `1.0779e17` T gates, with 179, 203, and 231
logical qubits, respectively. These are estimates under the existing
floating-point amplitude inputs and uncompiled circuit model. They do not
establish global optimality, a complete error guarantee or the production
counts. Each re-costed configuration met its target in 2,000 of 2,000
actual-grid emulated trials; the one-sided 99.9% lower confidence bound is
`0.99655`. GQSP truncation and gate synthesis were not emulated. The saved
results are in `projected-emulation-56.json`. A second search reproduced the
same configurations, and a second emulation reproduced all trial-error samples
exactly. Two floating-point amplitude-error fields changed by at most
`2.91e-14`; the comparison uses `rtol=1e-10, atol=1e-12` for numeric fields.
Previous `<875` results remain saved separately.

## Reproduction

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/b2_polynomial_majorant.py --verify-certificate > /tmp/b2-polynomial-majorant.json
cmp /tmp/b2-polynomial-majorant.json review/validation/b2-polynomial-majorant.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_b2_polynomial_symbolic.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_b2_polynomial_majorant.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_expansion_enclosures.py --case-set b2-56 > /tmp/expansion-enclosures-56.json
cmp /tmp/expansion-enclosures-56.json review/validation/expansion-enclosures-56.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_projected_box_bound.py --b2 56 --case-set b2-56 > /tmp/projected-b2-56.json
cmp /tmp/projected-b2-56.json review/validation/projected-b2-56.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_projected_recost.py --bounds review/validation/projected-b2-56.json --selection review/validation/b2-polynomial-reselection.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_projected_box_numerical.py --bounds review/validation/projected-b2-56.json --selection review/validation/b2-polynomial-reselection.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s review/validation -p 'test_b2_polynomial_extension.py' -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_projected_emulation.py --all --selection review/validation/projected-recost-56.json > /tmp/projected-emulation-56-repeat.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_b2_polynomial_reselection.py > /tmp/b2-polynomial-reselection-repeat.json
.venv/bin/python review/validation/compare_b2_polynomial_repeats.py --selection /tmp/b2-polynomial-reselection-repeat.json --emulation /tmp/projected-emulation-56-repeat.json
```
