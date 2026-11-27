# Ball enclosures for reconstructed Bessel tails

The resource model needs the infinite sum
`2 sum_{k>K} |J_k(x)|` to fit its per-segment threshold. The companion
`bessel-tail.md` bounds the portion omitted after 399 computed terms, but the
remaining SciPy sum is not an enclosure. This check encloses the full infinite
sum for reconstructed arguments and orders, using python-flint 0.9.0 and its
[ball arithmetic](https://flintlib.org/doc/guide/types/ball-arithmetic.html).

For each segment, `check_bessel_interval.py` evaluates `J_(K+1)` and
`J_(K+2)` with FLINT's ball-valued Bessel function. At arguments at least 128,
it uses the [Bessel recurrence](https://dlmf.nist.gov/10.6.E1)
`J_(n+1)=2 n J_n/x-J_(n-1)` to enclose subsequent terms. For smaller
arguments it evaluates them directly, because upward recurrence can amplify
interval widths. It sums absolute values with ball arithmetic. Every eight
terms, it adds the rational generating-function remainder from
`check_bessel_tail.py` and checks the resulting upper endpoint against the
threshold. The algorithm stops only when the comparison is proved. None of
the 634 checks needed more than 113 explicit terms.

The two inputs checked per segment are (1) the reconstructed binary result
`x_float` from the resource calculation and (2) the exact rational expression
`alpha_walk * eta_max * 2^j/(2^n_p-1)`, with the saved binary `alpha_walk` and
`eta_max` treated as exact. The second does not enclose uncertainty in those
parameters or the interval between the two inputs. The segment order `K` is
reconstructed by the earlier SciPy-based checker; only its aggregate walk
count, not each historical order, can be compared with saved output.

The resource calculation stored `eps_rot = eps_state/(2 n_rot)`, rounded to
binary floating point. Under round-to-nearest arithmetic on positive normal
values, the script steps `eps_rot` down by one float, multiplies by `2 n_rot`,
divides by `4 n_p(n_p+1)`, and applies an additional factor `1-2^-52`. This
gives a conservative rational lower threshold even allowing rounding in the
original final division. It depends on the saved `eps_rot` being the output
of that calculation and on the saved walk count being the original count.
The checker also verifies that this threshold is no larger than the previous
binary reconstruction.

For all 317 segments in the 27 saved configurations, both argument choices
give a full infinite-tail upper bound below `0.999821` times this lower
threshold. The largest displayed ratio is about `0.999820881`; the rational
endpoints in `bessel-interval.json` determine the comparison. Results are
conditional on the FLINT implementation, the stated arithmetic assumptions
and the reconstructed data. They do not validate physical-input uncertainty,
GQSP phase synthesis, or the remainder of the circuit and error budget.

The later [perturbation check](bessel-perturbation.md) covers a narrow
continuous interval around the reconstructed arguments, with the further
limitations stated there.

From the repository root, install the pinned packages in
`simulation/requirements.txt` plus the additional validation dependency
`review/validation/requirements-validation.txt`, then run:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s review/validation -p 'test_bessel_interval.py' -v
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_bessel_interval.py > /tmp/bessel-interval.json
cmp /tmp/bessel-interval.json review/validation/bessel-interval.json
```
