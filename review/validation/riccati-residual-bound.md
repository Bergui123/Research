# Rational residual enclosure for the benchmark Riccati solution

This calculation concerns only the normalized benchmark with `rho = 2/5`,
`lambda = T = nu = 1`, the fixed `Phi` in the paper and `H = I`. It builds on
the Riccati and coefficient-majorant arguments in `b2-global-majorant.md`.
The proof and code have not received independent review or Lean formalization.

## Interpolant and residual

`generate_riccati_nodes.py` evaluates a candidate Riccati solution at the
endpoints and midpoint of each of 4,000 intervals. It rounds each symmetric
matrix entry to an integer over `10^15` and writes those integers to
`riccati-quad-nodes.json`. The initial node is the exact rational `Phi`.
The numerical solver is used only to choose this rational data. The result
below does not depend on its error estimate or on the nodes being accurate.

On each interval, interpolate the three rational matrices with a quadratic
`P_hat`, using a local coordinate `z` in `[0,1]`. The adjacent polynomials
share endpoint values, so `P_hat` is continuous. Its residual is

```
R = P_hat' - (2 C - P_hat^2/2).
```

Each residual entry is a rational polynomial of degree at most four in `z`.
For a power polynomial `q(z) = sum_j a_j z^j` of degree at most `n`, its
Bernstein coefficients are

```
b_k = sum_{j<=k} a_j binom(k,j)/binom(n,j).
```

Because the Bernstein basis is nonnegative and sums to one on `[0,1]`,
`|q(z)| <= max_k |b_k|`. The verifier uses this bound for each residual entry,
then bounds the Frobenius norm of the symmetric residual by
`|R_11| + 2|R_12| + |R_22|`. It sums the resulting rational upper bounds
times the interval width. It also checks, with exact rational arithmetic, that
every quadratic Bernstein coefficient matrix of `P_hat - 0.9 I` is positive
definite. Convexity then gives `P_hat > 0.9 I` throughout every interval.

Let `P` be the positive definite Riccati solution and `E = P - P_hat`. Since
`E' = -(P E + E P_hat)/2 - R`, positivity of `P` and `P_hat` implies

```
d ||E||_F / dt <= ||R||_F,      E(0) = 0.
```

This also holds across interval boundaries by integration on each interval.
The exact verifier obtains the uniform upper bound

```
||E||_F <= 603153375526238754063306787 /
           48000000000000000000000000000000000 < 13/10^9.
```

The Bernstein matrices also bound the interpolant's maximum row sum by
2.387691. A row's l1 perturbation is at most `sqrt(2) ||E||_F`; the verifier
uses the rational upper factor two. It checks the strict consequences

```
P > (899/1000) I,      ||P||_infinity < 12/5.
```

## Consequence for the coefficient majorant

The benchmark ordering `P_11 >= P_22 >= P_12 >= 0` and the hierarchy
recurrence are as in `b2-global-majorant.md`. The tighter row bound permits
`R = (12/5)/4 = 3/5` in its ladder recurrence. The Gaussian factor obeys
`psi <= exp(-(899/8000)|x|^2)` because the determinant prefactor is at most
one. For even degree `d`, the same rational lower bound `e > 1957/720` gives

```
sup_{r>=0} r^d exp(-(899/8000)r^2)
  <= ( (360 d / 1957) / (899/8000) )^(d/2).
```

Substitution in the coefficient recurrence gives the exact rational expression

```
283799853651087882644322188868100987850828541399789 /
324355428187072967552478679783217698660937500000 < 875.
```

It is an analytic upper expression for the normalized benchmark's `B_2`,
conditional on the stated comparison and hierarchy arguments. A separate SymPy
implementation reproduces it. A 101-time floating-point check corroborates
the coefficient inequalities; it does not prove them between sampled times.
For the same `C`, `lambda`, `Phi` and `T`, the change of variables
`x = sqrt(nu) xi` uses this normalized constant with `delta = dx/sqrt(nu)`.

## Checks and remaining scope

`test_riccati_certificate.py` checks interpolation, Bernstein conversion and
that corrupted nodes are rejected. `check_riccati_residual_symbolic.py`
compares seven residual polynomials with a separate SymPy construction.
`check_riccati_interpolation.py` compares the interpolant with a 60-digit
matrix-exponential calculation at 201 times; the largest observed Frobenius
difference is about `2.15e-13`. This sampled diagnostic is not the proof of
the all-time bound.

`check_b2_riccati_reselection.py` reruns the benchmark search using `B_2=875`
and the analytic boundary charge in temporary files, including the previously
tested larger box. Its box evolution, expansion, selection inputs and circuit
costs remain numerical or uncompiled. It does not prove the resulting resource
counts. It does not apply the benchmark-only `B_2` bound to the six parameter
sweeps. `check_b2_riccati_emulation.py` emulates the three changed benchmark
configurations: each meets its target in 2,000 of 2,000 trials, with a
one-sided 99.9% lower confidence bound of 0.99655. The emulation excludes GQSP
truncation and gate synthesis. The production pipeline and its tables remain
unchanged.

From the repository root, the saved certificate can be checked without
regenerating its nodes:

```
.venv/bin/python review/validation/verify_riccati_nodes.py
.venv/bin/python -m unittest discover -s review/validation -p 'test_riccati_certificate.py'
.venv/bin/python review/validation/check_riccati_residual_symbolic.py
.venv/bin/python review/validation/check_riccati_interpolation.py
.venv/bin/python review/validation/b2_riccati_majorant.py
.venv/bin/python review/validation/check_b2_riccati_symbolic.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 .venv/bin/python review/validation/check_b2_riccati_majorant.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 .venv/bin/python review/validation/check_b2_riccati_reselection.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 .venv/bin/python review/validation/check_b2_riccati_emulation.py --all
```
