# Progress

Status as of 2026-10-25. Full detail: `WorkPlan.tex`. Plain-language version: `Summary.tex`.

| Stage | Topic | Status | Human Reviewed |
|---|---|---|---|
| 1 | Discrete trading-cost objective → continuum integral | Verified | No |
| 2 | Classical HJB + Riccati solution | Verified | No |
| 2b | Time-dependent $\lambda(t), H(t), C(t)$ | Verified | No |
| 3 | Two-asset worked instance, sales-only proof | Verified | No |
| 4 | Cole-Hopf linearisation ($\partial_s\psi = -\hat H\psi$) | Verified (algebra only) | No |
| 5 | Harmonic-oscillator spectrum; choice of $\nu$ | Derived on paper; **no $\nu$-sweep table yet** | No |
| 6 | Grid discretisation of $\hat H$ | **Not started** | No |
| 7 | Schrodingerisation (unitary embedding) | Not started | No |
| 8 | Qubit / T-gate resource count | Not started | No |
| 9 | Readout | Not started | No |
| 10 | Validation, sweeps, write-up | Not started | No |

"Verified" above means checked by AI (symbolic/numerical scripts and derivation review) only.
No stage has had a human reviewer check the derivations or results yet.

## Next step

Choose a finite terminal penalty and a $\nu$, solve the linearised PDE on a
small classical grid, recover the selling rate, and compare against the
Riccati solution with the same penalty. Then see how the error scales with
grid spacing, domain size, and $\nu$. This closes Stage 5's exit criterion
(the $\nu$-sweep table) and Stage 6's exit criterion (gradient error vs.
grid/domain/$\nu$) together.

To rerun existing checks: `.venv/bin/python -B verify_research.py`.
