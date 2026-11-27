# Rational grid comparisons for six 1% sweep candidates

This extends the benchmark calculation in `projected-box-bound.md` to one
selected readout-A stencil for each of the six parameter sweeps. It relies on
the rational Riccati and polynomial certificates in `sweep-riccati.md`, the
all-time boundary supersolution in `boundary-bound.md`, the grid-expansion
remainder theorem and the projected killed-walk comparison. These arguments
and their implementations await independent review and are not formalized in
Lean. They do not validate the amplitude inputs, synthesized circuit cost,
whole-algorithm error or a global optimum.

## Stencil enclosures and box response

`sweep_expansion.py` checks that each saved node file matches its certificate
by SHA-256, and that the selected grid has the same `(rho,lambda,nu,N,dx,h)`
as the calculation. It encloses `S=psi(1+delta^2 p1+delta^4 p2)` at the four
stencil nodes using the rational polynomial and Riccati errors. The term
`(log det Y)'=tr(P)/2` is integrated from the piecewise-quadratic Riccati
interpolant. Exponentials are bounded by rational positive Taylor sums and a
geometric tail, then all endpoints are rounded outward. The grid is the
ideal rational grid `dx=4/k`; the script checks its numerical representation
and box width against the selected configuration.

For `C=[[1,rho],[rho,1]]` and `lambda>0`, completing the square gives, for
either coordinate,

```
V(x) = lambda x^T C x/(4 nu)
     >= lambda (1-rho^2) x_i^2/(4 nu).
```

As in `projected-box-bound.md`, a finite-state Feynman--Kac argument bounds
the two-dimensional forced box response by the sum of two independent
one-dimensional killed exit responses. `sweep_grid_bound.py` uses the
parameterized lower potential above, outward-rounded fixed-point
uniformization and the rational boundary supersolution. It chooses the
smallest upper bound among 21 rational boundary parameters from `0.100` to
`0.600` in steps of `0.025`; each is a valid upper bound independently, so
their minimum is also an upper bound. It combines the box response with the
case's integer `B2`, the `S` intervals and a rational enclosure of the
reference rate. Rational logarithm and integer square-root bounds give an
upper expression for relative grid error, and `E<S` is checked at all four
nodes.

The `rho=0.6` and `rho=0.8` candidates from the original numerical search
had projected upper expressions about `0.00977` and `0.00943`, above the
`0.005` allocation. These upper expressions do not show that the true errors
are too large. Two specified larger boxes were tested numerically and then
checked with the rational comparison. For `rho=0.4, lambda=4`, the original
family had no numerical candidate under the conservative majorant; one
specified larger box was tested and checked. None of these three extra boxes
is claimed to be optimal over possible box sizes.

## Results and re-costing

All six fixed candidates satisfy `E<S` at all 24 stencil nodes and have the
following relative grid upper expressions, conditional on the stated proofs.
Each is below the 1% target's `0.005` grid allocation. The displayed bounds
are rounded upward; exact rationals are in `sweep-grid-bound-*.json`.

| rho | lambda | nu | grid | dx | h | grid upper | estimated T gates | qubits |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | 8 | `128^2` | `1/4` | 16 | 0.002309 | `5.27e14` | 185 |
| 0.2 | 1 | 16 | `128^2` | `4/11` | 16 | 0.001324 | `7.95e14` | 204 |
| 0.6 | 1 | 2 | `256^2` | `1/12` | 16 | 0.000209 | `2.71e15` | 219 |
| 0.8 | 1 | 4 | `256^2` | `1/9` | 16 | 0.000254 | `3.26e15` | 223 |
| 0.4 | 0.25 | 16 | `128^2` | `2/5` | 16 | 0.001460 | `6.48e14` | 202 |
| 0.4 | 4 | 16 | `256^2` | `1/7` | 16 | 0.003466 | `1.15e16` | 227 |

`sweep_recost.py` rejects mismatched parameter cases or grids and re-costs
each fixed candidate with its rational grid bound rounded upward. It scans
the saved safe auxiliary profiles under the existing resource model. The
amplitude inputs and coherent gate costs remain numerical or uncompiled; the
values are not a validated total cost or a global minimum. They are not a
replacement for the production tables, whose sampled-bound selections remain
unchanged. The previous numerical searches and the two failing projected
upper expressions are also retained.

For each final stencil, an independent floating-point two-dimensional box
solve lies below all four projected box bounds, the four numerical expansion
values lie inside the rational `S` intervals, and the saved measured grid
error lies below the stated relative upper expression. These 24-node checks
are diagnostics, not proofs of the all-time comparison or floating-point
roundoff. The expansion, grid and re-cost outputs reproduced byte for byte on
a repeat in this environment. The three larger-box numerical selections
reproduced within `rtol=1e-10, atol=1e-12`; the largest absolute change was
2 T gates on an estimate above `3e15`, while the largest relative change was
in a near-zero measured error and had absolute size `2.80e-13`.

For each re-costed candidate, `sweep_emulate.py` ran 2,000 actual-grid trials
with seed `20261122` and 2,000 with seed `20261123`. Every trial in both runs
met the 1% target. Each run's one-sided 99.9% Clopper--Pearson lower bound is
`0.99655`. The sample vectors differ between seeds. The emulations include
amplitude-estimation sampling and the implemented Schrödingerisation and
potential rounding, but exclude GQSP truncation and gate synthesis. They are
empirical coverage checks, not end-to-end guarantees.

`summarize_sweep_grid.py` checks the six case identities and grid matches,
the `B2` and `E<S` inequalities, upward conversion of the grid charge,
emulation profile matches, the trial counts and the reported confidence
intervals. Its compact output is `sweep-grid-summary.json`.

## Reproduction

From the repository root, first reproduce the six rational residual
certificates using `sweep-riccati.md`. The following mapping identifies the
fixed candidates; `extended-` indicates the larger box used for that sweep.
The high-`lambda` selection is `high-lambda-extension.json`.

```sh
while read -r id prefix selection; do
  if [ "$prefix" = none ]; then prefix=; fi
  OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/sweep_expansion.py --nodes "review/validation/sweep-nodes-$id.json.gz" --certificate "review/validation/sweep-cert-$id.json" --selection "review/validation/$selection" > "/tmp/sweep-expansion-$id.json"
  cmp "/tmp/sweep-expansion-$id.json" "review/validation/sweep-expansion-$prefix$id.json"
  OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/sweep_grid_bound.py --nodes "review/validation/sweep-nodes-$id.json.gz" --certificate "review/validation/sweep-cert-$id.json" --expansion "review/validation/sweep-expansion-$prefix$id.json" > "/tmp/sweep-grid-$id.json"
  cmp "/tmp/sweep-grid-$id.json" "review/validation/sweep-grid-bound-$prefix$id.json"
  OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/sweep_recost.py --bounds "review/validation/sweep-grid-bound-$prefix$id.json" --selection "review/validation/$selection" > "/tmp/sweep-recost-$id.json"
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_sweep_grid_numerical.py --bounds "review/validation/sweep-grid-bound-$prefix$id.json" --expansion "review/validation/sweep-expansion-$prefix$id.json" --selection "review/validation/$selection"
done <<'CASES'
rho0-lam1 none sweep-reselection-rho0-lam1.json
rho02-lam1 none sweep-reselection-rho02-lam1.json
rho06-lam1 extended- sweep-extension-rho06-lam1.json
rho08-lam1 extended- sweep-extension-rho08-lam1.json
rho04-lam025 none sweep-reselection-rho04-lam025.json
rho04-lam4 none high-lambda-extension.json
CASES
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s review/validation -p 'test_sweep_enclosures.py' -v
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/summarize_sweep_grid.py
```

To repeat the two emulation runs,
call `sweep_emulate.py --recost review/validation/sweep-recost-ID.json` with
`--seed 20261122` and `--seed 20261123`; each run is computationally expensive.

To reproduce the three larger-box numerical selections and compare them with
the saved outputs:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_sweep_extension.py --certificate review/validation/sweep-cert-rho06-lam1.json --nu 2 --n 8 --half-width 32/3 > /tmp/sweep-extension-repeat-rho06-lam1.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_sweep_extension.py --certificate review/validation/sweep-cert-rho08-lam1.json --nu 4 --n 8 --half-width 128/9 > /tmp/sweep-extension-repeat-rho08-lam1.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python review/validation/check_high_lambda_extension.py > /tmp/high-lambda-extension-repeat.json
.venv/bin/python review/validation/compare_sweep_extensions.py
```
