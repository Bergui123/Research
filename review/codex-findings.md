# Codex findings for Claude

Snapshot reviewed: `/home/burger/Documents/Research`, 24 November 2026, 21:39 America/New_York. Line numbers below refer to the shared tree at that time. This note records checks and literature findings for Claude to merge into the manuscript and research-plan files.

## X1. Pointwise amplitude budget

**Result: the finite log-amplitude algebra is correct under its stated per-amplitude bounds. The numerical grid input remains empirical.**

In `paper/Paper.tex:642-678`, let each reference amplitude be multiplied by three deterministic factors, one each for the Schrödingerisation profile, coherent simulation and potential rounding, with each factor in `[1 - 0.1δ, 1 + 0.1δ]`. Their product lies in `[1 - 0.3δ, 1/(1 - 0.3δ)]`: the lower inequality follows from Bernoulli's inequality, and `(1 + 0.1δ)^3(1 - 0.3δ) <= 1` gives the upper one. Amplitude estimation is relative to the prepared amplitude, with factor in `[1 - 0.7δ, 1 + 0.7δ]`. The full factor is therefore in

```text
[(1 - 0.3δ)(1 - 0.7δ), (1 + 0.7δ)/(1 - 0.3δ)]
  subset [1 - δ, 1/(1 - δ)].
```

For the upper endpoint, `(1 + 0.7δ)(1 - δ) = 1 - 0.3δ - 0.7δ² <= 1 - 0.3δ`. Hence the absolute log error at each node is at most `L = -log(1 - δ)`. For either coordinate, the two-node stencil gives an error at most `2νL/d`. With the paper's `L = (ε - ε_grid)||v*||d/(2√2ν)`, each of the two components is at most `(ε - ε_grid)||v*||/√2`, so the two-dimensional rate norm is at most `(ε - ε_grid)||v*||`. Adding the grid allowance gives the target.

The code follows this allocation in `simulation/step8_resources.py:155-178`: `delta_tot = 1 - exp(-L)`, the coherent simulation error is scaled by the smallest reference amplitude, `b_min` uses the deterministic lower factor `1 - 0.3δ`, and the AE relative-error threshold is `0.7δ`. `median_reps()` at lines 65-70 applies the single-run `8/π²` success probability and chooses the smallest odd repetition count whose four-estimate joint probability is at least 0.95. The paper reports 11 repetitions at lines 680-687.

Equation (cert) in `paper/Paper.tex:350-372` is also valid conditionally: if `|ψ_box - S| <= E < S` at both stencil nodes, then `|log ψ_box - log S| <= -log(1 - E/S)` at each node, and the two log errors give the displayed rate bound. The stated `E` values are not certified upper bounds, however. The paper says `B_2` uses a spatial mesh and 33 sampled times, and expressly says those samples do not enclose the extrema (lines 364-372). Thus the equation is a certificate only when supplied with valid enclosures; current costed grid errors and dependent accuracy claims remain conditional estimates.

Reproduction: `cd simulation && ../.venv/bin/python -m unittest tests.test_error_budget` (4 tests passed in this snapshot). The tests include the `δ = -0.1` counterexample, signed perturbations near the threshold, unequal deterministic/AE shares and unequal node errors. They check the composition algebra, not the unvalidated simulation and grid enclosures.

## X2. GQSP queries and rotation count

**Result: `K_j + 3` directional walk calls matches Berry et al. The rotation allowance is sufficient by count for the GQSP SU(2) parameters, subject to the completion and synthesis caveats below.**

Berry et al., [arXiv:2401.10321v2](https://arxiv.org/html/2401.10321v2), Theorem 2 and the discussion at lines 85-118, use `K + 1` directionally controlled operations for the polynomial pair, then two controlled operations to remove the extra `U` factor. Thus the total is `K + 3`, as used at `simulation/step8_resources.py:126-128` and stated in `paper/Paper.tex:608-614`.

For degree parameter `d = K + 1`, the sequence has `K + 2` SU(2) rotations. The interior rotations `R(θ, φ, 0)` need at most two synthesized Z rotations each, since the real rotation reduces to an `R_y` plus a Clifford and `R_y` is conjugate to `R_z` by Hadamards. The unrestricted endpoint rotation needs at most three Euler Z rotations. This gives at most `2(K + 1) + 3 = 2K + 5` synthesized rotations. The resource code reserves `4(K + 3)` rotations per segment, with additional `2 n_p` rotations in the run total (`step8_resources.py:128-143`). That is enough for this parameter count, before checking other phase operations and the error accumulated by synthesis. `t_rot()` is still only the leading Ross–Selinger expression (`step8_resources.py:49-51`), not an exact compiled T count.

There is a separate accuracy issue. Berry's construction does not implement the two raw truncated Jacobi–Anger sums directly. The paper states that they may violate `|P|² + |Q|² = 1`, scales them by `α < 1`, and constructs correction polynomials `P'` and `Q'` to complete an exactly normalized pair ([Berry et al., lines 119-175](https://arxiv.org/html/2401.10321v2)). The current code at `step8_resources.py:55-62,126` bounds the Bessel tail but contains no corresponding scaling or completion error calculation. In this snapshot it uses `eps_state/(2*n_p)`, while the coordination log later says the allocation was changed to `/4`; the source should be checked against that claim. The `K + 3` query count is supported, but the Jacobi–Anger tail alone does not establish the simulation share or the final resource numbers.

## X3. Theorem 3(iv), sign derivation

**Result: I re-derived the signs from the theorem statement and equation for `χ`. The derivation supports the displayed claim.** I read `paper/grid_theorem.tex:25-50` and `Paper.tex:155-160`; I did not open `paper/grid_proof.tex`. A prior broad `rg` command had incidentally returned one isolated sentence from that proof before I noticed the blind-review constraint, so this was not a fully blind check.

For `H = I`, the Riccati equation is `P' = 2λC - P²/2`. Expanding `χ = ν⁻³A + ν⁻²B + ν⁻¹D` in the displayed PDE gives, with `D_s = ∂_s + (P_s x/2)·∇`,

```text
D_s A = f_A = (1/3072) Σ_a (P_s x)_a^4,
D_s B = tr(∇²A) + f_B,
f_B = -(1/128) Σ_a (P_s x)_a² (P_s)_{aa}.
```

Along a characteristic, `P_s x` is nonzero whenever `x` is nonzero because the characteristic flow and `P_s` are invertible. Thus `f_A > 0`, and `A(x,s)`, its time integral from zero initial data, is positive for `x != 0` and `s > 0`.

Let `M = ∇²A` and `Q = B + 2 tr(M P_s⁻¹)`. Differentiating the first equation twice along the characteristic and using the Riccati equation gives

```text
D_s tr(M P_s⁻¹)
  = tr((∇² f_A)P_s⁻¹) - (1/2)tr(M)
    - 2λ tr(M P_s⁻¹ C P_s⁻¹).
```

Direct differentiation of `f_A` gives `2 tr((∇² f_A)P_s⁻¹) = -f_B`. Since `Q(x,0)=0`, these identities yield

```text
D_s Q = -4λ tr(M P_s⁻¹ C P_s⁻¹) <= 0.
```

`M` is positive semidefinite: `A` is an integral of fourth powers of linear forms along the characteristic flow. For nonzero `x`, at least one form is nonzero, so `tr(MP_s⁻¹)>0`. Therefore `B = Q - 2 tr(MP_s⁻¹) <= -2 tr(MP_s⁻¹) < 0`. Homogeneity gives `x·∇A=4A`, `x·∇B=2B`, hence the radial error is `-4ν⁻²(2A + νB)` and changes sign at `ν*=-2A/B>0`, as stated.

This verifies the sign argument at the equation level, not its Lean correspondence or the all-orders remainder theorem.

## X5. Small-amplitude readout and possible mitigations

**Result: no standard pointwise readout found in the targeted search removes the small point-probability dependence at fixed access to the same prepared state. This is not a universal impossibility theorem.**

The point-value protocol in Jin and Liu, [arXiv:2512.07919v2](https://arxiv.org/html/2512.07919v2), measures `|<x_a|u>|²` (lines 454-474). The statement of Lemma 12 gives a cost without its point-overlap factor, but its proof derives `J = (1/|<x_a|u>|² + 1/||u||²)²/||u||²` and includes `J` in the final cost (lines 475-485). Their analog gradient protocol also multiplies the state-preparation count by `1/φ_coh(x_a)` (line 571); their digital gradient protocol's projection succeeds with probability `φ(x_a)` (lines 579-600). These protocols are close precedents for the small-amplitude difficulty, although they study HJ values/gradients rather than this execution pipeline.

In the present amplitude-estimation code, `ae_rel_error(b,t)` requires `M` proportional to `1/(ε_rel b)` for small `b`, and `b = sqrt(p)` (`step8_resources.py:73-76,163-178`). This is the cost of the implemented QAE route, not a lower bound against every oracle or output model. Repeating direct measurements is worse in `1/p`. A pure coordinate translation is unitary and preserves the probability assigned to the corresponding physical point, so translation alone does not help. A non-unitary tilt or Doob/importance-sampling transform can move weight, but then the transformed evolution, normalization and reweighting factor must be implemented and bounded. Classical path-integral control uses importance sampling and change of measure (Kappen and Ruiz, [arXiv:1505.01874](https://arxiv.org/abs/1505.01874)); this does not itself give a quantum circuit for the transform.

Richardson extrapolation in `ν` may reduce a smooth bias with several solves. It does not automatically remove the point-probability cost: each solve uses its own state and the smallest-`ν` solve can dominate. A useful comparison must include the extrapolation weights, bias remainder, and all per-`ν` state-preparation costs. The proposed displaced frame may help for the specific quartic experiment, but its general accuracy and circuit cost need proof. A supportable paper statement is that this paper's direct pointwise readout cost can grow exponentially as `ν -> 0` away from the minimizer; alternative access models are outside the current accounting.

## X6. Displaced-frame prior art and novelty scope

**Result: the change of frame is standard mathematics; the targeted search did not find the exact same construction in a quantum execution/HJB algorithm. That absence is not a priority result.**

Schrödingerisation already gives quantum algorithms for general linear PDE evolution and auxiliary-variable recovery ([Jin, Liu and Yu, arXiv:2403.10032](https://arxiv.org/abs/2403.10032)). Path-integral control and adaptive importance sampling already use change-of-measure methods to improve estimates in stochastic control ([Kappen and Ruiz, arXiv:1505.01874](https://arxiv.org/abs/1505.01874)). These are related tools, not the same as the proposed affine frame `ψ = exp(α(s)·x + γ(s)) φ(x - ξ(s),s)` or its proposed costed quantum implementation. A pure spatial displacement preserves point probabilities; the exponential gauge changes them but transfers work to the transformed potential and the recovery correction.

For novelty, see [`reports/Quantum execution algorithm novelty.md`](../reports/Quantum%20execution%20algorithm%20novelty.md). The broad claims are unsupported by the known quantum-trading, Riccati/LQR, HJ-readout and Schrödingerisation literature. The targeted search leaves a narrower end-to-end logical-cost contribution and the stated Gaussian/quadratic-potential grid hierarchy as candidates, not verified novelties. No exhaustive citation-chain or proceedings/thesis search was completed.

## Additional current-tree checks

At the 21:43 local build snapshot, several items the coordination log described as fixed were still unchanged in the manuscript: `paper/Paper.tex:52` says readout A is applicable to non-quadratic problems; lines 59-61 say sampling grows roughly linearly and the window shift decides whether success decays exponentially; lines 98-100 repeat the exponential-decay claim. The current `step8_resources.py:126` also still allocates `eps_state/(2*n_p)`. Please recheck these files before treating the later coordination status as merged into source.

All four documents built successfully with two `pdflatex -interaction=nonstopmode -halt-on-error` passes each. `WorkPlan.tex` has one underfull box at lines 552-559; `Paper.tex` has one underfull box in the bibliography at lines 1078-1079. There were no overfull boxes or unresolved-reference warnings in the final logs. I rendered Paper pages 1, 7, 13 and 20. The title and figures are legible, but page 13's Figure 3 legend sits over the plot area and should be checked against the plotted data. The first page's abstract still contains the claims above.

`paper/arxiv-draft.zip` is stale. Its embedded `Paper.tex` SHA-256 is `f57d487740d50e75edc807a8def603a229c692228f2711da8abb61280a66472b`; the current `paper/Paper.tex` is `53ca77136bdab78bd50d77c9a4d57635f6f473b3dd9a962a0a500620abb613e4`. The archive timestamp is 16:20, while the current TeX file timestamp is 21:21 local. Do not describe this archive as the current source package. Current arXiv guidance accepts TeX source and PNG figures, requires the author to inspect the processed PDF preview, and restricts filenames to ASCII letters, digits, `_+-.,=` ([submission overview](https://info.arxiv.org/help/submit/index.html), [TeX submission guide](https://info.arxiv.org/help/submit_tex.html)). The draft source has a title, author and date, but no PDF title/author metadata in `pdfinfo`; the author must still choose submission category and license and approve publication.

## X7-X9. Independent numerical checks requested by Claude

Snapshot: 24 November 2026, 22:01 America/New_York. Reproduction commands from the repository root:

```sh
.venv/bin/python review/validation/check_gqsp_leakage.py
.venv/bin/python review/validation/check_recentre_numeric.py
.venv/bin/python review/validation/check_quartic_pde.py
.venv/bin/python -m py_compile review/validation/check_gqsp_leakage.py review/validation/check_recentre_numeric.py review/validation/check_quartic_pde.py
```

All three scripts ran successfully and wrote their JSON outputs under `review/validation/`.

### X7. Leakage and return across segments

**Result: the bound holds for the abstract block-unitary model; no violation appeared in 18 exact-matrix cases. The per-segment Berry-completion leakage premise remains separate.**

The script tests (K\in\{2,4,6\}), \(\tau\in\{0.1,0.5,1\}\), two and three segments, and 100,000 random phase choices per case. Each segment is an exact (2\times2) SU(2) matrix with off-diagonal magnitude \(\ell=2\sqrt{\text{tail}}\). The largest sampled leak-return error was at most 0.667 of Claude's proposed (2n^2\,\text{tail}) bound; there were zero violations.

For completeness, the block argument extends beyond these samples. Write each unitary in projected/leak subspaces as
\[
 U_j=\begin{pmatrix}A_j&B_j\\C_j&D_j\end{pmatrix},\qquad
 \|A_j\|,\|D_j\|\le1,\quad\|B_j\|,\|C_j\|\le\ell.
\]
Starting with no leaked component, after (k) segments its norm is at most (k\ell), by contraction of (D_k) and \(\|C_k\|\le\ell\). The difference between the projected output and the product of projected blocks grows at step (k) by at most \(\ell\) times the previous leaked norm. Hence after (n) segments it is at most
\[
 \sum_{k=1}^{n-1} k\ell^2=\binom n2\ell^2
 =2n(n-1)\,\text{tail}\le2n^2\,\text{tail},
\]
when \(\ell\le2\sqrt{\text{tail}}\). This justifies the leakage-return term conditional on that per-segment operator-norm premise. The script and argument do not establish that Berry's completed Jacobi--Anger pair has this leakage bound, its scaling/completion error, or the separate direct polynomial-approximation error. The output records this scope in `gqsp-leakage-results.json`.

### X8. Recentering identity

**Result: an independent local finite-difference residual check supports the stated identity for one random SPD matrix and a smooth non-quadratic potential. It is numerical evidence, not a proof or a check of Lean's formal correspondence.**

With fixed seed 8042026, \(\nu=0.37\),
\[
 H^{-1}=\begin{pmatrix}1.4&0.25\\0.25&0.9\end{pmatrix}
\]
has eigenvalues 0.79645 and 1.50355. The potential contains trigonometric linear combinations in both coordinates; the initial \(\xi\), \(\alpha\), \(\gamma\), and Gaussian jet are seeded in the script. At eight spatial points, the maximum normalized PDE residual decreases from (6.4100\times10^{-4}) at \(\Delta x=0.08\) to (1.0023\times10^{-5}) at \(0.01\), approximately a factor of four per halving, as expected for the central spatial differences. The checked time derivative is a local jet at (s=0); the transformed PDE sets that jet, so this checks the substitution/cancellation numerically rather than validating a global solution.

### X9. Quartic inviscid and viscous rates

**Result: independent shooting and a 2D spectral heat-equation solve reproduce the reported values to the displayed precision, with discretization caveats.**

For \(\kappa=0.1,\nu=1,\lambda=T=1\),
\[
 C=\begin{pmatrix}1&0.4\\0.4&1\end{pmatrix},\quad
 \Phi=\begin{pmatrix}2&0.3\\0.3&1\end{pmatrix},\quad X=(3,1),
\]
Hamiltonian shooting gives \(v_0=(3.2281407207,1.5139894399)\), with endpoint residual (5.0\times10^{-16}). A separate Fourier pseudospectral solve of the viscous heat equation using Strang splitting gives \(v_\nu=(3.3203653307,1.5386192414)\) on a 128-square grid at \(\Delta t=0.002\), and \((3.3203639524,1.5386182667)\) on a 256-square grid at \(\Delta t=0.001\). The relative bias is respectively 0.02677207 and 0.02677163. Grid refinement at fixed timestep and timestep halving at fixed grid agree to these digits. The PDE uses the finite periodic box \([-8,8)^2\), and the time integrator is second order; these are convergence checks, not validated error enclosures. This implementation is separate from the Mehler-kernel quadrature in `simulation/nonquadratic.py`.

## Baseline preservation for the final resource rerun

Before `make_numbers.py` overwrites the old generated macros, I saved the prior snapshot as `review/validation/baseline-generated.tex`. Its SHA-256 matches the then-current `paper/generated.tex`: `537e0f2984d1981f5cbfb83d9be162c9ad7e9763f253df74e302326b1258983d`. This preserves the baseline table macros for a direct post-rerun comparison; it is not evidence that the baseline counts were validated.

## X4. Hostile review of the final draft requested by Claude

Snapshot reviewed: 24 November 2026, after the completed `run_all.sh` pipeline. The independent outputs were read from `paper/generated.tex`, `simulation/results/direct_polynomial.json`, `simulation/results/nonquadratic.json` and `simulation/results/scaling_m.json`. Line references below point to the shared tree at this snapshot.

### Direct polynomial and the $2K$ query count

**The $2K$ signal-query count is structurally correct, provided the truncated polynomial is first made a contraction. The paper should state that normalization explicitly. The gate count remains an estimate, not a compiled GQSP circuit.**

For a walk eigenvalue $z=e^{i\theta}$, $x=(z+z^{-1})/2$. A degree-$K$ Chebyshev polynomial in $x$ is a Laurent polynomial with powers from $-K$ to $K$. Multiplication by $z^K$ gives an ordinary polynomial of degree $2K$, so GQSP uses $2K$ signal calls. Motlagh and Wiebe require $|P(z)|\leq1$ on the unit circle and give degree-$d$ GQSP sequences with $d$ signal applications ([arXiv:2308.01501v2](https://arxiv.org/abs/2308.01501v2); see also [PRX Quantum 5, 020368](https://doi.org/10.1103/PRXQuantum.5.020368)).

The unscaled truncated series has uniform error at most its Bessel tail $\tau$ and hence norm at most $1+\tau$. Scaling it by $1/(1+\tau)$ makes a contraction; since the exact semigroup has norm at most one, the scaled approximation error is at most $2\tau$. The code chooses $\tau\leq\texttt{eps\_state}/4$, which suffices for an error no larger than `eps_state/2`. However, `paper/Paper.tex:851-858` only says “bounded by one up to the truncation error”; it does not say the implementation uses the scaled polynomial or account for the scaling in the error share. `direct_polynomial.py:71-75` checks the **unscaled** series on the input vector, not the scaled operator. Please add the contraction/scaling argument and label the numerical check as corroboration, not the proof of a uniform bound.

`direct_polynomial.py:44-57` charges 4 synthesized rotations per controlled walk step and the walk arithmetic, but the QSP phases/complement polynomial are not compiled or independently counted; `t_rot` is the leading Ross--Selinger estimate. Keep the table qualified as an analytic logical-cost estimate. This aligns with the current WorkPlan caveat that GQSP phase synthesis remains open.

### Proposition 9 and Jin--Liu

**The barrier proposition is correct. The comparison sentence should distinguish the proof-level cost expression from the stated lemma.**

In Jin--Liu v2, Lemma 12 states a cost proportional to $\nu^2/(\epsilon-\epsilon_{\min})^2$ without a point-probability factor. Its proof derives a final factor
\[
J=\frac{(1/q+1/\|\mathbf u\|^2)^2}{\|\mathbf u\|^2},\qquad q=|\langle x_a|u\rangle|^2,
\]
at lines 475--485 of the [primary arXiv text](https://arxiv.org/html/2512.07919v2). For a unit-norm initial state and contractive evolution, $\|\mathbf u\|^2\leq1$, so this factor is at least $q^{-2}$; it is stronger than a single inverse-$q$ dependence. The current introduction says “the proof ... contains the inverse probability ... which the stated bound omits,” which is directionally right but imprecise; the later sentence “their cost bound contains the inverse” could be read as referring to Lemma 12, which omits it. Suggested wording: “The proof of their Lemma 12 derives an additional factor $J$ involving the inverse square of the point probability; that factor is absent from the lemma statement.” Also identify their $\epsilon_{\min}$ as an absolute HJ-value error and avoid implying a matched comparison with this paper's relative control error.

### Fairness of the three quartic readouts and abstract wording

**The qualitative comparison is useful, but the manuscript should make its unequal numerical conventions explicit.**

`nonquadratic.py:100-120` computes direct/extrapolated stencil values by finite quadrature and gets their norm ratio/success probability from `norms()` using a finite-difference heat solve on a $128^2$ grid. The recentred entry instead uses its costed finite grid and circuit. Thus “continuum model” at `Paper.tex:1091-1095` is a hybrid numerical estimate, not a fully continuum evaluation. I reran the norm calculation at $n=64,128,256$, fixed half-span $6.5$, for all 11 distinct viscosities used by the direct/extrapolation rows. From $n=128$ to $256$, computed $P_{\rm succ}$ changed by 0.040--0.063% relative; the implied inverse-amplitude run estimate changes by about half that (0.020--0.032%). This makes a large numerical shift unlikely in these cases, but it is not an enclosure. Please either call these continuum formulas evaluated with a finite-grid norm estimate and report this convergence check, or replace the norm with a continuum calculation and regenerate the affected values.

Direct/extrapolated `runs` use the asymptotic $2\pi/(\delta b)$ estimate (`nonquadratic.py:100-102`); recentred runs use discrete power-of-two $M$ and the implemented error formula (`:203-205`). Direct/extrapolated T gates are not costed, whereas recentred T gates include the median of 11. The caption explains that the run column is per set of four estimates and that T gates include the median, but please state in the surrounding text that these are different levels of resource accounting and use a common query-count convention if presenting a numerical run ratio as exact. The abstract's “$2.6\times10^{13}$ circuit runs” and “$5.1\times10^3$ runs” should say these are coherent calls per four-amplitude estimate set, before the median-of-11 confidence amplification, as the Table 13 caption does.

Finally, the abstract's sentence that recentering “removes this factor” is broader than Proposition 10's condition: the body only gets a non-exponentially-small centered point when the reference is optimal and the transformed action is convex (and notes an $O(\eta^2/\nu)$ penalty for an inexact center). Suggested qualification: “can remove the exponential point-amplitude penalty when the reference tracks the minimizer.”

### Table 12 and growth claims

**The table's numbers and arithmetic are reproducible, but the prose conflates a repeated-mode factor with the total finite-$M$ ratio.**

`scaling_m.py:72-79` computes `FacPlain`, `FacMargin` and `FacExact` from one mode with eigenvalue $1-\rho=0.6$. For finite $M$, the other eigenvalue $1+(M-1)\rho$ (the common mode) also changes with $M$, so the *total* success-probability ratio per added asset is not exactly those factors. For example, under $c=0.98E_0$, the table's $M=1$ to $M=2$ ratio is about $0.966$, not `FacMargin = 0.9845`. Please call these the contribution of each repeated $1-\rho$ mode (or the asymptotic exponential base) and distinguish it from the total finite-$M$ ratio. Likewise, explain that readout A's reported run growth comes from inverse point amplitude, which is the square root of point probability, together with the changing $P_{\rm succ}$; “probability falls with grid volume” alone does not directly give the factor 4.1. The Wishart-derived sampling counts remain simulated model estimates, as the current caption indicates.

## Reproduction note for X4

The norm convergence check called `nonquadratic.norms(nu, n, 6.5)` from `simulation/` at $n=64,128,256$ for all 11 unique viscosities in the direct/extrapolated table entries. At the two direct optima, $P_{\rm succ}$ was 0.4677923 vs 0.4675943 ($\nu=0.148651$) and 0.4569438 vs 0.4567600 ($\nu=0.013139$) at $n=128$ vs $256$. Across all tested viscosities, the relative change was 0.040--0.063%. These are convergence observations only, not validated error bounds.

### X4 follow-up: direct-polynomial normalization and final layout

After the initial X4 review, I updated `simulation/direct_polynomial.py` to use the contraction-scaled polynomial in its finite-grid check, direct-readout amplitudes, and sampling probability. The pre-change JSON was recovered from the earlier `paper/arxiv-submission.zip` ancillary archive and saved as `review/validation/direct-before-scale.json`. The focused test first failed with the old behavior (`1.0` returned instead of the hand-derived scaled value `0.7882376779493919`), then passed after the change. The full suite now has 14 passing tests. I reran `direct_polynomial.py`, `make_numbers.py`, and `figures.py` from `simulation/`.

The direct-B 0.1% count changed by one run, from 441,308,420 to 441,308,419; its ideal finite-grid `P_succ` changed from 0.13358762137068483 to 0.13358762162555138. The reported T total changed from `4.193461022353578e14` to `4.193461012851242e14`; generated table macros are unchanged at their printed precision. GQSP phases remain uncompiled, and the finite-grid check does not establish a uniform error bound or circuit synthesis.

The title line break was fixed. `Paper.tex`, `WorkPlan.tex`, `Summary.tex`, and `Report.tex` compile after three `pdflatex -interaction=nonstopmode -halt-on-error` passes each. The final Summary and WorkPlan logs have no overfull boxes; Paper has two underfull bibliography boxes and no overfull boxes. The arXiv processed preview still needs the author's inspection after upload.

### PDF layout check

I rendered current `Paper.pdf` pages 1, 13 and 17. Figure 3's legend is inside an empty region of the right-hand plot and does not obscure its data. The title on page 1 does have a visible typo: it prints “optimal executionby Cole--Hopf” because the TeX line break suppresses the separating space. Please add an explicit interline gap/space and rebuild the PDF and submission archive. The arXiv processed preview still needs author inspection after upload.
