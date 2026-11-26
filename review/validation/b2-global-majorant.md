# A coarse global bound for the benchmark lattice remainder

This argument applies to the normalized two-asset benchmark with `nu = 1`,
`T = lambda = 1`, `rho = 2/5`, and
`Phi = [[2, 3/10], [3/10, 1]]`. It bounds the `B_2` in the grid theorem.
It has not been independently reviewed or formalized. The bound is too large
to establish the paper's resource targets and is not used by the production
pipeline.

## Bounds on the Gaussian

Write the symmetric Riccati solution as `P = [[a, q], [q, b]]`. It satisfies

```
P' = 2 C - P^2/2,      P(0) = Phi,
C = [[1, 2/5], [2/5, 1]].
```

`Phi` has eigenvalues in `[7/10, 23/10]` by row bounds. Positivity of `P`
and the inverse equation `(P^-1)' = I/2 - 2 P^-1 C P^-1` give
`P >= (14/27) I` on `[0,1]`. At the upper barrier `(12/5) I`, the Riccati
vector field points inward because `2 lambda_max(C) - (12/5)^2/2 < 0`.
Thus `P <= (12/5) I`. In dimension two,
`||P||_infinity <= sqrt(2) ||P||_2 <= 18/5`.

The entries also satisfy `a >= b >= q >= 0`. Indeed, `a-b` solves a scalar
homogeneous equation with nonnegative initial value; `q' = 4/5 - (a+b)q/2`
preserves nonnegativity; and

```
(b-q)' = 6/5 - (a+b)(b-q)/2 + det(P)/2
```

preserves `b-q >= 0`. These statements use the positive definite Riccati
solution established in the manuscript. Since `Y' = P Y/2`, `det(Y) >= 1`.
Consequently the normalized Gaussian obeys
`psi(x,t) <= exp(-(7/108)|x|^2)` throughout the time interval.

## Coefficient recurrence

Let `C[k,d](t)` be the l1 norm of the coefficients of total degree `d` in
the polynomial `p_k`; set `U[0,0]=1` and other `U[0,d]=0`. On a homogeneous
degree, the l1 logarithmic norm of the drift
`-(P x).grad / 2` is nonpositive: the column for `x_1^i x_2^j` has sum
`-[i(a-q)+j(b-q)]/2`. The Laplacian maps degree `d+2` to `d`, with l1 norm
at most `(d+2)(d+1)` times the incoming norm.

For nonnegative degree bounds `U_d`, define the ladder majorant

```
(L_R U)_d = (d+1) U[d+1] + R U[d-1],    R = (18/5)/4 = 9/10.
```

This bounds the coefficients after one `D_a = partial_a - (P x)_a/4` for
either axis. With `c_j = 2/(2j+2)!`, recursively define, for `k = 1,2`,

```
U[k,d](t) = integral_0^t {
    (d+2)(d+1) U[k,d+2](s)
    + sum_{i<k} 2 c_{k-i} (L_R^(2(k-i)+2) U[i])_d(s)
} ds.
```

Start each block at its highest possible degree `4k` and work downwards. The
drift contraction, the Laplacian bound and the hierarchy equation imply
`C[k,d] <= U[k,d]`; all operations in the recurrence preserve nonnegative
coefficient bounds.

For each even `d`, `|x_1^i x_2^j| <= |x|^d` when `i+j=d`, and

```
sup_{r>=0} r^d exp(-(7/108) r^2)
    <= (38880 d / 13699)^(d/2),
```

with the value one for `d=0`. The inequality uses
`e > sum_{j=0}^6 1/j! = 1957/720`; it therefore uses rational arithmetic.
Combining these weights with `U[k]` and the three terms in the `B_2` integral
gives the exact rational upper expression

```
608469527952873377793049859602570313245647
------------------------------------------------  < 920674.
660896107550319725526780100000000000
```

The three term contributions are approximately 19.48, 6131.62 and 914522.39.
The proof relies on the stated Riccati and hierarchy equations and comparison
steps; exact arithmetic in the recurrence does not independently verify those
premises.

## Evidence and practical limit

`b2_global_majorant.py` evaluates the recurrence with exact rational arithmetic.
`check_b2_symbolic.py` reproduces the rational result with a separate SymPy
implementation. `check_b2_global_majorant.py` checks the Riccati, coefficient
and ladder inequalities at 101 floating-point time samples. The latter is a
diagnostic and does not prove behavior between samples.

The present pipeline charges a sampled value about 0.00553775. Replacing it by
920674 is therefore a vastly more conservative calculation. For the tested
larger-box 0.1% candidate, `delta^6 = 1/450^3` and the lattice remainder charge
alone becomes about 0.0101034 in amplitude. With the existing floating-point
evaluation of `S`, one positive term of the *resulting relative rate error
bound* is about 0.0655, above the 0.0005 grid allocation. This does not imply
that the actual rate error is that large. It shows that this coarse majorant
cannot validate the reported target. No resource counts were reselected from
it. A useful next step is a validated, sharper enclosure of the actual
time-dependent polynomial coefficients and their spatial suprema.

From the repository root:

```
.venv/bin/python review/validation/b2_global_majorant.py
.venv/bin/python review/validation/check_b2_symbolic.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 .venv/bin/python review/validation/check_b2_global_majorant.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 .venv/bin/python review/validation/check_b2_budget_impact.py
```
