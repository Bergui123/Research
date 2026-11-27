# Argument intervals for reconstructed Bessel tails

The previous [ball check](bessel-interval.md) evaluates two arguments per
segment. This check covers every argument in a narrow interval containing
both: the computed binary argument and the exact formula
`alpha_walk * eta_max * 2^j/(2^n_p-1)` with each saved binary `alpha_walk`
and `eta_max` independently varied by one adjacent floating-point value.
The interval does not represent a physical uncertainty model or enclose all
rounding in the calculation of those saved parameters.

For `T_n(x) = 2 sum_{k>=n} |J_k(x)|` and `0 < x <= b < n-1`, the
[Bessel derivative identity](https://dlmf.nist.gov/10.6.E1) and the
generating-function coefficient bounds in `bessel-tail.md` imply

```
|T_n(x) - T_n(x0)| <= |x-x0| [U_(n-1)(b) + U_(n+1)(b)] / 2,
U_m(b) >= 2 sum_{k>=m} |J_k(t)| for all 0 <= t <= b.
```

Each `U_m` is the exact rational bound returned by `omitted_upper(b,m)`.
Its generating-function bound increases with `t`, so using `b` bounds every
point in the interval. Uniform convergence follows from the same geometric
bound and justifies passing the derivative estimate through the sum.

For every segment, the checker first reproduces the saved binary argument and
the aggregate `sum(K+3)` walk count. It combines the rational perturbation
allowance with a fresh FLINT-ball enclosure at the binary argument, and checks
that their sum is below `0.999821` times the conservative reconstructed
threshold defined in `bessel-interval.md`. All 317 segments in 27 saved
configurations pass. The largest displayed ratio is approximately
`0.999820885`; exact rational endpoints are saved in
`bessel-perturbation.json`. Each enclosure requires at most 113 Bessel terms
before the generating-function remainder.

The statements depend on FLINT's ball arithmetic, the stated binary
rounding assumptions, and reconstruction from saved data. The original
individual orders were not saved, and the interval does not cover other
physical-parameter or modeling uncertainty. The result does not validate the
rest of the error budget or compile the circuit.

With the dependencies in `simulation/requirements.txt` and
`requirements-validation.txt` installed, run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s review/validation -p 'test_bessel_*.py' -v
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_bessel_perturbation.py > /tmp/bessel-perturbation.json
cmp /tmp/bessel-perturbation.json review/validation/bessel-perturbation.json
```
