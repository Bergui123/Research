# Elementary-gate check of the auxiliary Fourier convention

This is a floating-point simulation of ideal Hadamard, controlled-phase,
single-qubit phase and swap gates. It checks the convention for small powers
of two; it does not synthesize the phases into Clifford+T, validate their
precision or prove a finite-grid aliasing bound. The phase-estimation and
Hamiltonian simulation circuits are not included.

With `N=2^n`, `c=(N-1)/2`, `p_l=-L+2L l/N` and
`eta_k=(k-c) pi/L`, the forward unitary is

```
F[k,l] = exp(-i eta_k p_l)/sqrt(N)
       = exp(-i pi c) (-1)^k Qminus[k,l] exp(2 pi i c l/N),
Qminus[k,l] = exp(-2 pi i k l/N)/sqrt(N).
```

Use little-endian qubits for the integers `l` and `k`. The right-hand phase
is a product of single-qubit phase gates with angles `2 pi c 2^j/N`. The
left-hand `(-1)^k` is `Z` on qubit zero; `exp(-i pi c)` is a global phase.
The checker constructs `Qminus` from `n` Hadamards, `n(n-1)/2` controlled
phase gates and `floor(n/2)` swaps. It constructs the inverse by reversing
and adjointing the gate list. It compares both matrices with the direct
half-integer transform and compares the forward action with `g_hat_fft`
after accounting for that routine's `1/N` rather than `1/sqrt(N)`
normalization.

For a scalar Hamiltonian eigenvalue `E`, the intermediate diagonal is

```
exp(i eta_k E T) = exp(-i c pi E T/L)
                   product_j exp(i pi 2^j E T k_j/L).
```

The checker applies this diagonal as a global phase and `n` one-qubit phase
gates between the explicit forward and inverse circuits. At `N=4,8,16,32`,
with `L=4`, it compares five eigenvalue/time shifts per size against direct
matrices. The worst evolution entry error is `2.08e-14`; this is numerical
agreement, not a uniform error bound for all parameter values.

For an exact two-node translation, the circuit produces `g[l+2]` and changes
its sign when `l+2` crosses `N`, as required by antiperiodicity. The largest
matrix error over the four sizes is `7.08e-15`. The input is sampled at
`p_l+o`, where `o` puts one grid point at the specified negative window start;
the Fourier transform itself still uses `p_l`. For window points whose
translated index does not wrap and whose physical argument is nonnegative,
the same test compares the recovered sample with `exp(-(p+ET))`; the largest
error is `1.01e-15`. This exact-node test does not establish recovery for
arbitrary non-integer shifts or bound interpolation aliasing.

As a separate finite-sum diagnostic of the paper's antiperiodic wrap lemma,
the script evaluates 100 image terms at 121 values of `y` in `[-6,6]` for
`L=4`, `R=2` and the smooth-width-one profile. The largest observed absolute
sum is `0.135291`, below the lemma's analytic bound `0.270762`. The finite
floating-point diagnostic is not a proof of the infinite sum or all times;
the paper's geometric-series argument supplies the stated analytic bound.
It does not cover the aliasing error.

Reproduce from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_fourier_circuit.py > /tmp/fourier-circuit.json
cmp /tmp/fourier-circuit.json review/validation/fourier-circuit.json
PYTHONPATH=simulation OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s simulation/tests -p 'test_auxiliary_fourier.py' -v
```
