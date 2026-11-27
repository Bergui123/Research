# Omitted Jacobi--Anger tail in the resource estimate

`simulation/step8_resources.py` computes `2 sum |J_k(x)|` for the 399 terms
`k=K+1,...,K+399` while selecting each order. The analytic cost argument
uses the infinite sum. The generating function in
[DLMF 10.12.1](https://dlmf.nist.gov/10.12.E1) gives, for any `r>1` and
`x>=0`, the Laurent coefficient bound

```
|J_k(x)| <= exp[x(r - 1/r)/2] r^(-k),
2 sum_{k>=n}|J_k(x)| <= 2 exp[-D] / (1 - 1/r),
D = n log(r) - x(r - 1/r)/2.
```

For the saved binary floating-point argument `x`, the checker chooses a
rational `r` and bounds `log(r)` below by 32 positive terms of its atanh
series. If this gives a positive rational lower bound `D0`, then
`exp(-D) <= 2^(-floor(D0))` because `e>2`. The resulting omitted-tail upper
bound is rational. The floating-point choice of `r` is merely a search for a
convenient positive rational; the inequality holds for any `r>1`.

The checker reconstructs the segment orders and thresholds for 27 saved
configurations: 18 production readout-A and readout-B benchmark/sweep
configurations, three benchmark alternatives and six sweep alternatives. The
sum of reconstructed `K+3` agrees with the saved walk-step count for each
configuration. Across 317 segments, each rational omitted-tail bound is less
than `3.41e-16` times the reconstructed threshold. The smallest observed
fractional slack of the 399-term floating-point sum is greater than `0.00023`.
The sum of the observed finite part and rational omitted bound is below the
reconstructed threshold in every case.

This does **not** enclose errors in SciPy's Bessel values, their summation,
the original physical parameters, or the recovered threshold. It uses the
binary value of `x` produced by the resource arithmetic as an exact input to
the omitted-tail bound. Individual selected orders were not saved separately;
the aggregate walk-step agreement does not prove they are identical if the
software or environment changes. Nor does this check synthesize GQSP phases or
validate the remaining circuit and amplitude assumptions. The full selected
infinite-tail inequality remains a numerical claim until those errors are
enclosed.

Reproduce from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s review/validation -p 'test_bessel_tail.py' -v
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_bessel_tail.py > /tmp/bessel-tail.json
cmp /tmp/bessel-tail.json review/validation/bessel-tail.json
```
