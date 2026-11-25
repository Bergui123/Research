# Schrödingerisation, HJ algorithms, and grid analysis

Search date: 2026-11-24. Scope: primary papers and preprints related to the
paper's Cole--Hopf/Hamilton--Jacobi route, Schrödingerisation, recovery from a
finite auxiliary coordinate, and the claimed all-order grid expansion. The
search is targeted, not an exhaustive priority search.

## Does prior work already combine Cole--Hopf, quantum PDE methods, and HJ readout?

### Takeaway

The high-level method is not new: quantum algorithms for Hamilton--Jacobi
solutions already recover point values and gradients, and recent work combines
Cole--Hopf with Schrödingerisation for HJ-related equations. This does not by
itself preempt the paper's domain-specific resource estimates or its particular
grid theorem, but the paper should not imply novelty for the general
Cole--Hopf-to-quantum-PDE or pointwise-gradient route.

### Cited Findings

- Jin, Liu, and Yu introduced Schrödingerisation for linear PDEs and gave
  examples including the heat, convection, Fokker--Planck, linear Boltzmann,
  and Black--Scholes equations. Their technical paper derives the auxiliary
  coordinate transport and recovery of the original solution. [Primary
  preprint](https://arxiv.org/abs/2212.14703).
- Jin and Liu's HJ paper develops quantum algorithms for viscosity solutions
  with convex Hamiltonians. It treats pointwise value and gradient estimation
  explicitly, including digital and analog gradient protocols. [Primary
  preprint](https://arxiv.org/abs/2512.07919).
- Jin, Liu, and Yu's 2026 paper applies Cole--Hopf and Schrödingerisation to
  the coupled heat equations for regularized Wasserstein proximal operators;
  its model includes a viscous HJ equation with quadratic Hamiltonian.
  [Primary preprint](https://arxiv.org/abs/2606.28752).
- Jin and Liu's 2025 technical work derives a smooth auxiliary profile and
  analyzes cut-off, polynomial, and Fourier initializations. It provides
  bounds for recovering the non-unitary solution from an auxiliary-coordinate
  interval. [Primary preprint](https://arxiv.org/abs/2505.00370).

### Inferences

- The defensible distinction is likely the concrete logical resource estimate
  for this optimal-execution benchmark, plus any genuinely new theorem about
  its discrete quadratic-potential heat equation. The broad combination of
  Cole--Hopf, HJ readout, and Schrödingerisation is already present in adjacent
  research.
- The manuscript already cites the HJ and Wasserstein papers. It should state
  what differs technically: model, discretization, estimator, error theorem,
  and level of resource accounting. The phrase "closest to this paper" should
  explain that distinction rather than suggest the prior work only treats
  asymptotics.
- The manuscript's readout A may remain useful for a broader class of models,
  but a quantum pointwise HJ-gradient algorithm is prior art. Any novelty
  claim should attach to the specific implementation and cost analysis.

### Gaps

- This search did not enumerate every HJ, control, or PDE algorithm paper, nor
  every follow-up and citation to the papers above. A final priority claim
  needs a broader citation-chain review.
- I did not compare the full algorithms' assumptions or error definitions
  line by line. The cited papers establish close subject overlap, not
  equivalence to the present execution problem.

## Are the finite auxiliary grid, profile, and shifted recovery window new?

### Takeaway

Smooth auxiliary profiles, finite periodic truncation, characteristic-based
recovery intervals, and bounds tied to spectral extrema all predate this
manuscript. The particular use of a lower-spectrum shift to move the recovery
window and improve the postselection probability may be a narrower contribution,
but the search did not establish priority for that formula. Present the broad
wrap and profile ideas as existing machinery, then isolate and compare the
specific shift result.

### Cited Findings

- Jin, Liu, and Yu's original Schrödingerisation paper extends the auxiliary
  coordinate periodically when its tails are small, uses Fourier
  discretization, and recovers heat-equation solutions on a positive-
  coordinate region. [Primary preprint](https://arxiv.org/abs/2212.14703).
- Jin, Liu, and Ma analyze auxiliary-coordinate truncation and spatial
  discretization errors for Schrödingerised ill-posed PDEs. Their recovery
  theorem selects a region beyond a spectral cutoff proportional to the
  maximum Fourier frequency squared times time, and their error analysis
  treats periodic truncation. They also smooth the initial auxiliary profile
  to improve convergence. [Primary preprint](https://arxiv.org/abs/2403.19123).
- Jin, Liu, Ma, Peng, and Yu give a detailed cut-off construction that is
  smooth, agrees with the required exponential profile on the recovery side,
  and has a finite interval chosen using eigenvalue bounds and time. Their
  Theorem 5.1 recovers the solution only in an interior interval after the
  characteristic shifts at both ends are accounted for. [Primary
  preprint](https://arxiv.org/abs/2505.00370).
- Nguyen and Tura's November 2026 paper gives a domain-of-dependence result
  for a periodic auxiliary interval and states the exact recovery slice and
  interval-width condition from the characteristic speed. [Primary
  preprint](https://arxiv.org/abs/2611.02357).

### Inferences

- The paper's claim that the auxiliary interval must cover the evolution
  range of the largest grid eigenvalue is an application of established
  characteristic/domain-of-dependence reasoning. The manuscript should cite
  the general finite-interval results and make its own contribution the
  problem-specific bound, profile, and verified cost tradeoff, if those are
  supported.
- A shift of a generator by a scalar changes only a global factor in the
  non-unitary solution, while it changes the auxiliary-coordinate transport
  velocities. The manuscript's proposition specializes this freedom to a
  positive heat operator, using a lower eigenvalue bound to move the recovery
  window and deriving a success-probability factor. This appears distinct in
  the sources checked, but it should be described as a spectral-offset
  optimization and compared directly with prior recovery-window results.
- The paper's smooth blend is a profile construction in a literature where
  smooth cut-off and interpolation profiles are established. Its newness is
  not supported by merely showing a direct matrix/FFT match.

### Gaps

- I found no exact match for the paper's formula using a window beginning at
  `-cT`, the condition `c <= E_min`, and the gain `exp(2 c T)` in the sources
  checked. That is not evidence that no earlier source contains the same
  scalar spectral shift.
- I did not verify the paper's discrete half-integer Fourier phases, wrap
  estimate, or gate cost. Literature search cannot substitute for the
  circuit-level and numerical checks in `AGENTS.md`.
- The 2026 preprint was submitted on 2026-11-02, before this search date. The
  manuscript's bibliography and related-work discussion should account for
  it if the result is relevant to the paper's final claims.

## Is the all-order lattice expansion established as a novel result?

### Takeaway

Modified-equation expansions for finite-difference schemes are longstanding,
so "expansion to every order" alone is not a novelty argument. The narrower
claim in the manuscript appears to be a uniform-on-the-infinite-lattice
remainder theorem for a multidimensional discrete heat equation with a
quadratic potential, followed by a rate-error expansion whose quartic and
quadratic terms have opposite signs. I found no direct prior match in this
targeted search, but cannot establish priority from that absence.

### Cited Findings

- Warming and Hyett's original paper introduced modified-equation analysis
  for the stability and accuracy of finite-difference approximations. It
  derives modified PDEs by Taylor expansion and uses them to analyze
  discretization errors. [Publisher record](https://www.sciencedirect.com/science/article/pii/0021999174900114).
- Jin and Liu's HJ algorithm paper includes discretization and gradient error
  estimates for its viscous HJ scheme, but the reviewed results use a
  different entropy-penalization construction. [Primary
  preprint](https://arxiv.org/abs/2512.07919).
- Jin, Liu, and Ma analyze truncation and spectral-discretization error for a
  Schrödingerised PDE in the auxiliary variable. That concerns the extra
  coordinate and does not match the manuscript's central-difference spatial
  grid theorem. [Primary preprint](https://arxiv.org/abs/2403.19123).

### Inferences

- The manuscript should make the theorem's exact scope its novelty statement:
  operator, initial data, norm, domain, order-uniform constants, and how the
  rate expansion follows. Avoid presenting generic modified-equation
  expansion as the new result.
- The sign/cancellation claim may be the sharper distinction, provided the
  proof's hypotheses and the transition from the infinite lattice to the
  periodic box remain explicit. The numerical box estimates are currently
  sampled, so they do not certify the resource budgets.

### Gaps

- The targeted searches did not find a prior theorem for the same
  multidimensional quadratic-potential lattice problem. Numerical PDE and
  difference-operator literature is broad, and a search-engine absence is not
  proof of originality.
- I did not audit the theorem's derivation or Lean formalization. This note is
  a literature comparison, not a correctness review of the stated expansion.
