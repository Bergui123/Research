# Research project instructions

Every time you add or change content (new stage, new equation, new result),
write TWO versions, both updated together, same commit:

- **Long version** → `WorkPlan.tex`. Full derivation, source, caveats, exit
  criterion, per the existing stage format.
- **Short version** → `Summary.tex`. Equation + one-line status/next-step only.
  No derivation, no caveats prose.

Never update one without the other.

# Manuscript review and repair backlog

**Goal:** Resolve the documented inconsistencies and produce a reproducible,
accurately qualified manuscript suitable for an arXiv submission review.

**Status (24 November 2026, evening):** two reviews are recorded (`review/codex-findings.md`,
`review/claude-response-current.md`), their accepted findings are implemented, and the pipeline
(`simulation/run_all.sh`) was rerun in full. Open items are the unchecked boxes below. The second
review is by the same model as the one that made the revisions; a human referee has not read the paper.

**Working rules for the next agent:**

- Work through the tasks in order, marking a checkbox only after its stated
  verification passes. Record evidence, commands, remaining caveats and changed
  numerical results in `WorkPlan.tex`, and the corresponding equation and
  one-line status in `Summary.tex`, in the same commit as the research change.
- Preserve existing uncommitted work. Review `git status` and relevant diffs
  before editing; do not commit unrelated changes. If a claim cannot be proved,
  keep it explicitly conditional or empirical rather than presenting it as a
  guarantee. Do not assume that the current estimates survive a correction.
- Review `paper/Paper.tex`, `report/Report.tex`, `README.md`, generated tables,
  figures and standalone stage copies for consistency whenever a task changes
  their claims. If a `stageN/derivation.tex` file changes, update its copied
  `stageN/StageN.tex` as described in `README.md`.
- The numerical pipeline and build commands are documented in `README.md`.
  Reproduce baseline outputs before replacing them; save enough evidence to
  compare old and new results. Never invent a theorem, citation or test result.

## 1. Correct amplitude-error budgets

Progress: the finite allocation and the per-node estimation and rounding
budgets have been implemented and tested. The selected configurations still
depend on unvalidated numerical error estimates and circuit costs.

- [ ] **Issue:** The former bound used `|log(1+δ)| ≤ |δ|`, which fails for
  negative relative error. The paper and `simulation/step8_resources.py` now
  use a finite bound and the counts have been recalculated. The numerical
  assumptions behind the individual error shares still need validation.
- [x] **Repair:** Starting from the exact log bound in `paper/Paper.tex` and
  `WorkPlan.tex`, derive a sufficient bound for the sum of both stencil-node
  errors and both rate components. Account for how amplitude estimation,
  Schrödingerisation, potential rounding and coherent-simulation errors compose
  with respect to the same reference amplitude. Update `run_cost` and any
  affected downstream calculations; do not simply add relative errors without
  justifying the composition. Document the confidence event and assumptions.
- [x] **Verify** (`tests/test_error_budget.py`; `step10_validation.py` now asserts the amplitude share and
  `sqrt(a) >= b (1 - 0.3 delta)` on the actual grid for the benchmark and every sweep, which exposed and fixed the use of
  the 64^2 success probability on the 256^2 grid): Add focused checks for positive and negative perturbations,
  including the counterexample `δ = -0.1`, both unequal node errors, and a
  value near the allowed error threshold. Recalculate candidate configurations
  and counts, and demonstrate that the selected configurations satisfy the
  finite bound. If only an approximation is available, retain that qualification.

## 2. Repair the comparator endpoint

Progress: the scale now strictly exceeds the evaluated maximum, including at
an exact power of two; the small-unitary check uses the same scale. Endpoint
tests pass. All 18 selected benchmark and sweep configurations satisfy the
strict scale and representability checks in floating-point arithmetic; affected
resource and sampling results were regenerated. Circuit compilation remains
part of Task 7.

- [x] **Issue:** `simulation/step8_resources.py::potential_scale` chose a
  power-of-two scale using a ceiling. If the maximum is exactly that power of
  two, the ideal truncated value is `2**m`, which does not fit in `m` bits.
  The paper and implementation now use a strict scale inequality, and the
  selected inputs have been checked in floating point.
- [x] **Repair:** Use a strictly larger scale when required, or specify and
  implement a saturating encoding with a correctly derived error bound. Update
  resource accounting, `paper/Paper.tex`, and the corresponding Lean statement
  if its assumptions or stated correspondence change.
- [x] **Verify:** Test maxima exactly equal to a power of two and values on
  either side; check encoded values, representability and the rounding bound.
  Audit all selected benchmark and sweep configurations, then recalculate any
  affected costs and emulations.

## 3. Validate the grid-error and box estimates

- [ ] **Issue:** `simulation/certify_grid.py` samples 33 times and uses
  floating-point spatial and temporal extrema. These samples are not proven
  upper enclosures of the time-continuous boundary source or the supremum in
  the lattice remainder. The cancellation-free leading-order estimate is also
  not, by itself, a general finite-grid bound.
- [x] **Repair** (second option: retained as numerical values, charged at twice the sample, labelled in the paper): Derive
  and implement validated upper bounds for both the
  lattice constant and the periodic-box boundary source, including spatial,
  temporal and numerical errors, or explicitly retain the values as numerical
  estimates and make all dependent accuracy claims conditional. Track the
  theoretical bound separately from its numerical evaluation.
- [ ] **Verify:** Check the bound assumptions at every costed stencil and sweep
  point. Compare validated enclosures against finer independent evaluations,
  demonstrate the required `E < S`, and verify that each charged bound covers
  the measured error. Re-run the resource selection after any change. Never
  call sampled maxima or floating-point estimates certified without an enclosure.

## 4. Resolve the auxiliary Fourier convention

Progress: the phase-adjusted FFT agrees with explicit matrices for 4, 8, 16,
and 32 points. Circuit phase costs and a finite-grid wrap bound remain open.

- [ ] **Issue:** For even `N_p`, the half-integer frequencies in
  `simulation/step7_schrodinger.py::g_hat_fft` and
  `simulation/chebyshev.py::p_grid` are antiperiodic over `2 L_p`. The paper
  identifies this distinction, but the transform phases, circuit and
  wrap-around bound have not been independently verified together.
- [x] **Repair** (paper Section 5.2, `lean/ExecutionProofs/Wrap.lean`; phase rotations charged): Derive the exact discrete
  forward and inverse transforms,
  including all phase corrections and their circuit costs. Prove a wrap bound
  for the convention actually used, or change the implementation and paper to a
  consistently justified convention.
- [ ] **Verify** (transform and emulation checked; no circuit-level QFT was built): Compare the circuit-level transform and evolution with an
  explicit small matrix for multiple grid sizes and eigenvalues, including a
  shift across the boundary; check the recovery window and wrap bound. Re-run
  affected simulations, gate counts and qubit counts.

## 5. Fix the amplitude-estimation emulation

Progress: the sampler now uses all distinct outcomes; the endpoint regression
test and large-grid tail test pass, and the benchmark emulations, tables and
figures have been regenerated. The raw floating-point probability sum for the
largest selected `M = 2**24` differs from one by about `3.7e-14`; individual
probability enclosures and independent review remain open.

- [x] **Issue:** The former `simulation/step10_validation.py::sample_qae`
  truncated and renormalized the Fejér distribution for large `M` and
  duplicated an endpoint for smaller `M`. The revised sampler includes all
  distinct outcomes, and the reported benchmark validation was rerun.
- [x] **Repair:** Sample the exact `M` distinct outcomes efficiently, or
  provide and enforce a mathematically justified bound on the effect of any
  approximation on the reported success probabilities. Keep reproducible seeds.
- [x] **Verify:** For small even `M`, compare the full probability mass and
  empirical frequencies with the canonical QAE law, including phases at and
  between grid points. For saved large-`M` cases, check total and omitted mass
  and confirm no endpoint duplication. Re-run the affected validation and
  update the tables and figures.

## 6. Establish or qualify sampling confidence

- [ ] **Issue:** `simulation/step9_sampling.py` uses a Monte Carlo estimate of
  an asymptotic Gaussian/Wishart quantile for a discrete, perturbed output.
  `ceil(n_s / P_succ)` gives roughly `n_s` accepted samples in expectation, not
  with high probability. This fact alone does not show the reported target is
  missed. The 2,000-trial validation is empirical and excludes synthesis errors.
- [x] **Repair** (acceptance fluctuation by an exact binomial, coherent error charged, Wishart quantile asymptotic and
  labelled): Derive an end-to-end finite-sample bound for the actual
  distribution and acceptance fluctuations, including systematic and coherent
  errors, or describe confidence and counts as empirical estimates throughout.
  If using a fixed run count, budget the probability of too few accepted
  samples; if using a stopping rule, report expected and high-probability costs.
- [x] **Verify** (one-sided 99.9% Clopper-Pearson bounds asserted at least 95% for the benchmark and all sweeps): Check
  acceptance and error coverage for all target accuracies
  and sweep points with reproducible independent simulations and uncertainty
  intervals. Recompute counts and compare with the stated confidence claim.

## 7. Audit circuit costs and formalized claims

- [ ] **Issue:** The explicit small-unitary block-encoding check and the Lean
  axiom check do not establish the full circuit cost or that each paper
  statement exactly matches its formalization.
- [ ] **Repair** (SELECT compiled and simulated, GQSP count and tail rederived, QROM/QFT/GQSP phases still formulas): Audit PREPARE, SELECT, comparator arithmetic, QROM state
  preparation, Fourier phases, controlled operations, reflections, GQSP,
  rotation precision, ancillas, and measurement-based uncomputation, including
  operations used inside amplitude estimation. Check every paper theorem and
  hypothesis against its Lean declaration, especially comparator endpoints and
  shifted-window assumptions. Correct claims or proofs where they differ.
- [ ] **Verify:** Compile or otherwise independently cost at least one complete
  small walk step and its inverse/controlled use; compare with the resource
  model. Run `lake build` and `lake env lean scripts/AxiomCheck.lean` from
  `lean/`, and document exactly what is and is not formalized.

## 8. Reproduce all results and audit the manuscript

- [x] Run the numerical scripts in the order in `README.md`, including sweeps,
  `make_numbers.py` and `figures.py`; the documented sweeps take roughly
  40 minutes. Verify assertions, regenerated result files, table macros and
  figures. Record changed counts and investigate unexpected differences.
- [ ] (done for the paper, WorkPlan, Summary, report and README by the two agents; a human reading is still open)
  Compare every equation, assumption, benchmark value and limitation across
  the manuscript, `WorkPlan.tex`, `Summary.tex`, report, scripts and relevant
  stage documents. Check all citations against primary sources. The journal
  title for arXiv:2212.14703 was already corrected; the rest still need audit.
- [ ] Distinguish theoretical bounds, numerical estimates, empirical validation
  and unverified extrapolations throughout. Do not claim all issues were found
  merely because the listed checks pass.

## 9. Obtain an independent review and prepare for submission

- [ ] Obtain Claude Code's response to the final X4 review through
  `coordination.md`. The saved Claude review in
  `review/claude-response-current.md` covers an earlier snapshot, so it does
  not independently review the latest paper edits. Keep the final response and
  discussion in Markdown, then resolve or record each finding and rerun the
  relevant checks.
- [x] (`paper/make_submission.sh` builds and compiles it from a clean directory; the author still has to check arXiv's
  preview) Prepare a clean source package containing `paper/Paper.tex`, its four
  included TeX files, and all five figures with working paths. Provide a stable,
  accessible location or appropriate ancillary files for the simulation and Lean
  sources. Follow the current official arXiv submission instructions and check
  the processed preview before any submission.
- [ ] Build `Paper.tex`, `WorkPlan.tex`, `Summary.tex` and `report/Report.tex`
  with `pdflatex -interaction=nonstopmode -halt-on-error` from the directories
  specified in `README.md`, repeating until references resolve. Inspect the
  output for layout errors; the summary previously had a small overfull box.
  Confirm the author metadata, category, license and approval for publication
  with the author when those decisions are needed.

**Exit criterion:** Every completed task has reproducible evidence; all
remaining limitations are stated accurately; the two research-plan versions
and derived documents agree; the final counts are regenerated from the corrected
pipeline; and an independent reviewer has examined the resulting manuscript.
A local build does not guarantee arXiv acceptance or replace author approval.
