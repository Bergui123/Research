# Exact residual bounds for the benchmark grid-expansion polynomials

This concerns the same normalized benchmark and time interval as
`riccati-residual-bound.md`. The mathematics and implementation await
independent review; neither is formalized in Lean. The rational enclosures below
are conditional on the Riccati comparison and on the coefficient logarithmic
norm argument in `b2-global-majorant.md`. They do not enclose finite-box
propagation, establish the resource counts, or cover the six parameter sweeps.

## Interpolation and residual

The grid theorem defines `p_0 = 1`, `p_1`, `p_2` through

```
L_P q = sum_a [partial_a^2 q - (P_a . xi) partial_a q/2],
D_{P,a} q = partial_a q - (P_a . xi) q/4,
p_1' = L_P p_1 + sum_a D_{P,a}^4 1/12,
p_2' = L_P p_2 + sum_a (D_{P,a}^4 p_1/12 + D_{P,a}^6 1/360),
p_1(0) = p_2(0) = 0.
```

Their degrees are at most four and eight, with only even total degrees.
`generate_polynomial_nodes.py` computes numerical candidates at the same
8,001 endpoints and midpoints used for `P_hat`, and rounds their coefficients
to integers over `10^15`. Those numbers only choose the interpolants: the
bounds do not use an ODE-solver error estimate. The saved candidates are in
`polynomial-quad-nodes.json.gz`. Piecewise quadratic interpolation gives
continuous `p_hat_k`; the exact initial values are checked.

For each interval, `verify_polynomial_nodes.py` forms the two polynomial
residuals `r_k = p_hat_k' - F_k(P_hat,p_hat_1,p_hat_2)` with rational
arithmetic. It expands every spatial coefficient as a polynomial in the local
time variable. Bernstein convex-hull bounds enclose each residual coefficient
and each interpolant coefficient over the whole interval; summing absolute
values gives a bound by homogeneous degree. It takes the maximum over all
4,000 intervals and rounds each upper bound **outward** to a rational with
denominator `10^12` (residuals) or `10^9` (interpolants). All intermediate
operations use exact fractions. `check_polynomial_residual_symbolic.py`
independently differentiates the spatial polynomials with SymPy and finds
exact agreement for both residuals on four intervals at local time `1/3`.
This is a diagnostic, not an independent proof for every interval.

## Propagating coefficient errors

Let `A[k,d]` be the checked degree-`d` coefficient-l1 bound for `p_hat_k`,
`R[k,d]` the residual bound, and `W[k,d](t)` the l1 bound to construct for
the degree-`d` part of `p_k - p_hat_k`. The earlier rational Riccati
certificate gives `||P-P_hat||_F < e = 13/10^9` and both row-l1 norms below
`12/5`. For the true benchmark Riccati solution,
`P_11 >= P_22 >= P_12 >= 0`. Thus the coefficient-l1 logarithmic norm of
the drift on each homogeneous degree is nonpositive. The difference between
true and interpolated drift contributes at most `d e A[k,d]`, while the
Laplacian sends degree `d+2` to degree `d` with norm at most `(d+2)(d+1)`.

For each ladder application, if `A_d` bounds the interpolant input and `E_d`
the input difference, the next degree bounds are accumulated according to

```
A'_{d-1} += d A_d;       A'_{d+1} += (3/5) A_d,
E'_{d-1} += d E_d;       E'_{d+1} += (3/5) E_d + (e/2) A_d.
```

The last term follows from the row-l1 error `<= 2e`, divided by four. Apply
this recurrence four or six times to the appropriate interpolants, multiply
by the factors in the differential equations, and sum both axes. With those
forcing bounds `F[k,d](t)`, the nonnegative comparison polynomials satisfy

```
W[k,d](0) = 0,
W[k,d]' = R[k,d] + d e A[k,d] + F[k,d]
          + (d+2)(d+1) W[k,d+2].
```

Solving from the highest degree downward, and for `k=1` before `k=2`, gives
exact rational bounds at `t=1` in `polynomial-residual-bound.json`. The
largest degree-l1 bound for `p_1` is `2.511e-8`; that for `p_2` is `2.708e-7`.
These conservative bounds hold at every earlier time too: all comparison
polynomials have nonnegative coefficients. They bound coefficient errors, not
the truncation of the grid theorem's expansion.

## Expansion at the three exploratory stencils

For each stencil node `x`, `check_expansion_enclosures.py` bounds

```
S(x) = psi(x) [1 + delta^2 p_1(x/sqrt(nu))
                  + delta^4 p_2(x/sqrt(nu))],   delta = dx/sqrt(nu).
```

Since only even total degrees occur, each interpolant monomial at a rational
`x` and rational `nu` is rational. For each degree, its coefficient error is
bounded by `W[k,d](1) (|x|^2/nu)^(d/2)`. To bound the Gaussian prefactor, use
`(log det Y)' = tr(P)/2`; the interpolant integral is rational, and its
contribution to the exponent differs by at most `e/2`. The quadratic exponent
differs by at most `e |x|^2/(8 nu)`. The exponential is enclosed using a
positive Taylor sum for `exp(z)` and a geometric bound on its positive tail,
then reciprocated. Final endpoints are rounded outward to multiples of
`10^-15`. The denominator is checked positive before multiplication.

All twelve intervals for the three benchmark stencils are saved in
`expansion-enclosures.json`. Their largest width is `6.083e-9`. A separate
floating-point evaluation in `check_expansion_against_float.py` lies inside
each interval; this diagnostic is not itself an enclosure. The finite-box
boundary-response vector used in the grid budget is still produced by a
floating-point matrix exponential, and the production results still use
sampled error estimates. Therefore these expansion enclosures do not complete
the end-to-end validation or establish the stated resource counts.

## Reproduction

From the repository root, the saved nodes can be verified without regenerating
them. The first command takes several minutes on one core. Its output must
match the saved JSON, which is used by the later expansion check.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/verify_polynomial_nodes.py > /tmp/polynomial-check.json
cmp /tmp/polynomial-check.json review/validation/polynomial-residual-bound.json
.venv/bin/python -m unittest discover -s review/validation -p 'test_polynomial_certificate.py' -v
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_polynomial_residual_symbolic.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_expansion_enclosures.py > /tmp/expansion-check.json
cmp /tmp/expansion-check.json review/validation/expansion-enclosures.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_expansion_against_float.py
```
