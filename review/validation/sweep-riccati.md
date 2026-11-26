# Rational lattice-remainder majorants for six sweep points

This extends the benchmark residual certificates in
`riccati-residual-bound.md` and `b2-polynomial-majorant.md` to the six saved
parameter sweeps. It uses the same finite terminal matrix
`Phi=[[2,3/10],[3/10,1]]`, `T=1`, identity impact matrix, and normalized
`nu=1`; changing variables to `x=sqrt(nu) xi` gives the corresponding bound
with `delta=dx/sqrt(nu)` for other positive `nu`. The proof and implementation
have not been independently reviewed or formalized in Lean. The numerical
re-selections below do not validate the production resource tables.

## Riccati and polynomial checks

For `0<=rho<1` and `lambda>0`, the symmetric positive-definite solution
satisfies

```
P' = 2 lambda [[1,rho],[rho,1]] - P^2/2,    P(0)=Phi.
```

Writing `P=[[a,q],[q,b]]`, the benchmark's ordering extends to these sweep
parameters: `a-b` solves a scalar homogeneous equation with positive initial
value; `q' = 2 lambda rho - (a+b)q/2` preserves `q>=0`; and

```
(b-q)' = 2 lambda (1-rho) - (a+b)(b-q)/2 + det(P)/2
```

preserves `b>=q`. Positive definiteness is supplied by the Riccati comparison
in the manuscript. The ordering makes the coefficient-l1 logarithmic norm of
`-(P xi).grad/2` nonpositive on each homogeneous degree.

For each sweep, `generate_sweep_nodes.py` uses floating-point solvers only to
choose rational endpoint and midpoint values for 512 quadratic time pieces.
The saved gzip JSON files contain those integers, with denominator `10^15`.
`sweep_riccati.py` bounds each quartic residual entry on each entire interval
using rational Bernstein coefficients. It checks every Bernstein coefficient
matrix of `P_hat - m I` is positive definite, using `m=43/50` for
`rho=2/5, lambda=1/4` and `m=9/10` for the other five sweeps. The residual
comparison gives a uniform Frobenius error no larger than the exact sum of
interval residual bounds; this gives a positive lower eigenvalue and a row-l1
upper bound for the true `P`. The solver's tolerance is not used as evidence.

`sweep_polynomial.py` bounds both polynomial ODE residuals with the same
Bernstein construction. It propagates errors with the degreewise contraction,
Laplacian and ladder inequalities from `polynomial-residual-bound.md`, but
uses each sweep's certified Riccati error and row bound. Its comparison
polynomials have nonnegative coefficients, so their terminal values bound all
earlier times. The resulting coefficient-l1 bounds feed the same Gaussian
and ladder majorant as in `b2-polynomial-majorant.md`. Each output file saves
an exact rational expression and a strictly larger integer.

| rho | lambda | rational-majorant upper expression, rounded | integer charged |
| ---: | ---: | ---: | ---: |
| 0 | 1 | 24.126322 | 25 |
| 0.2 | 1 | 34.537804 | 35 |
| 0.6 | 1 | 91.366116 | 92 |
| 0.8 | 1 | 147.171763 | 148 |
| 0.4 | 0.25 | 9.057795 | 10 |
| 0.4 | 4 | 27652.335663 | 27653 |

The displayed values are rounded upward; exact values are in the six
`sweep-cert-*.json` files. The upper expressions can be far larger than the
actual constants. In particular, the global row and eigenvalue bounds make
the `lambda=4` estimate very conservative.

An independent SymPy generating-polynomial calculation reproduces each exact
B2 expression. A separate symbolic diagnostic checks both polynomial
residuals at local time `1/3` on three selected intervals per sweep and
checks the Riccati residual against the saved integral bound. At 101
floating-point times per sweep, independent ODE evaluations satisfy the
Riccati, coefficient and ladder inequalities. These samples do not prove the
bounds between times. Three focused tests cover a known constant Riccati
solution, rejected invalid data, and error propagation.

## Numerical re-selections, not validated resource counts

`check_sweep_reselection.py` repeats the saved 1% readout-A search for each
sweep using the corresponding integer `B2` and the analytic all-time boundary
charge. It copies inputs to a temporary directory; it does not modify
`simulation/results`. The boundary formula is valid for these nonnegative
quadratic potentials. The box evolution, expansion evaluation, amplitude
inputs and circuit costs in this exploratory search remain floating point or
uncompiled. No claim of an enclosed grid error, global optimum or validated
confidence follows from these numerical selections.

| rho | lambda | candidate in original tested grid family | numerical estimated T gates | qubits |
| ---: | ---: | --- | ---: | ---: |
| 0 | 1 | `nu=8, N=128, h=16` | `5.27e14` | 185 |
| 0.2 | 1 | `nu=16, N=128, h=16` | `7.29e14` | 204 |
| 0.6 | 1 | `nu=2, N=128, h=8` | `8.00e14` | 203 |
| 0.8 | 1 | `nu=4, N=128, h=8` | `8.82e14` | 205 |
| 0.4 | 0.25 | `nu=16, N=128, h=16` | `6.48e14` | 202 |
| 0.4 | 4 | none under this bound and family | — | — |

The absence of a candidate in the last row does not show that the actual
algorithm misses its target. A sharper time-dependent majorant or other grids
could give a different conclusion. The old production selections and all
previous output files remain unchanged. The re-selection data are saved as
`sweep-reselection-*.json`; their baseline entries retain the original counts.
A second run reproduced every saved numeric field and selection exactly;
`compare_sweep_repeats.py` also permits `rtol=1e-10, atol=1e-12` for future
floating-point differences.

## Reproduction

From the repository root, for each entry below, use the specified `m` value.
Each exact certificate check uses the saved candidate nodes and takes time;
`OPENBLAS_NUM_THREADS=1` is used for consistent diagnostics.

```sh
for spec in rho0-lam1:9/10 rho02-lam1:9/10 rho06-lam1:9/10 rho08-lam1:9/10 rho04-lam025:43/50 rho04-lam4:9/10; do
  id=${spec%:*}; m=${spec#*:}
  OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/verify_sweep_case.py --input "review/validation/sweep-nodes-$id.json.gz" --min-hat "$m" > "/tmp/sweep-cert-repeat-$id.json"
  cmp "/tmp/sweep-cert-repeat-$id.json" "review/validation/sweep-cert-$id.json"
  OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_sweep_symbolic.py --certificate "review/validation/sweep-cert-$id.json"
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_sweep_diagnostic.py --certificate "review/validation/sweep-cert-$id.json"
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_sweep_reselection.py --certificate "review/validation/sweep-cert-$id.json" > "/tmp/sweep-reselection-repeat-$id.json"
done
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_sweep_residual_symbolic.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s review/validation -p 'test_sweep_*.py' -v
.venv/bin/python review/validation/compare_sweep_repeats.py
```
