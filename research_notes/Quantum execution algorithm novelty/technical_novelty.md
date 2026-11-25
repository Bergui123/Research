# Technical novelty in the Schrödingerisation and grid analysis

Search date: 2026-11-24. I read the paper's grid and Schrödingerisation sections
and compared their claims with primary arXiv papers. Searches included
"Schrödingerisation Hamilton-Jacobi quantum algorithm", "Schrödingerisation
finite p domain error", "shift Hamiltonian imaginary-time success
probability", "quadratic potential lattice heat equation error expansion",
and "uniform lattice heat-kernel expansion". The searches identify close
precedents but do not establish a comprehensive priority search.

## Is the shifted recovery window a distinct contribution?

### Takeaway

The recovery-window shift is useful for the reported success probability, but
the underlying operation is a scalar energy shift of imaginary-time evolution.
The paper should present the exact window identity as a problem-specific
derivation or implementation improvement, not as a new general principle,
unless it establishes a sharper distinction from known energy shifting.

### Cited Findings

- Jin, Liu, and Yu's Schrödingerisation paper encodes heat evolution as a
  unitary evolution in an auxiliary coordinate and explains that the recovered
  state-preparation cost depends on the norm ratio
  `||u(0)|| / ||u(t)||`. It also discusses imaginary-time ground-state
  preparation. [Jin, Liu, and Yu, arXiv:2212.13969](https://arxiv.org/abs/2212.13969)
- In its ground-state section, the same paper describes exponential filtering
  under `exp(-Ht)` and its norm decay, with the ground energy in the exponent.
  [Full text, Appendix C](https://arxiv.org/html/2212.13969v1)
- A 2026 primary preprint states that replacing `H` by `H+bI` in normalized
  imaginary-time evolution changes only an overall scalar, leaving the
  normalized state unchanged. It also studies choosing a constant shift to
  move the ground energy to zero. [Tang and Zhang, arXiv:2608.17792](https://arxiv.org/abs/2608.17792)
- The current manuscript's Proposition `shifted window` uses the bound
  `c <= E_min` so that every spectral component has entered the exponential
  part of the auxiliary profile throughout the selected interval. It then
  integrates that common exponential profile to obtain its stated success
  probability. (Paper, Section "Shift picture and recovery window".)

### Inferences

- Subtracting `c I` changes `exp(-H T)` by the scalar `exp(cT)` and changes the
  squared norm ratio by `exp(2cT)`. This explains the manuscript's success
  factor and follows from standard energy shifting; the exact lower-bounded
  interval is a tailored recovery statement.
- The shift appears to improve postselection probability without changing the
  normalized heat solution. Its novelty is therefore likely an application to
  the selected readout window, not the energy-shift mechanism itself.

### Gaps

- I did not find the manuscript's exact interval formula in the searched
  sources. Search absence does not prove priority.
- The exact implementation must still account for discrete normalization,
  finite auxiliary-grid boundaries, and errors in `E_min` and the profile.
  This literature check does not validate those resource or error claims.

## How much of the auxiliary-grid and wrap analysis is already known?

### Takeaway

Finite auxiliary-domain selection, periodic approximation, recovery-region
constraints, and smooth initial profiles all have close prior work. The paper's
specific finite-grid formula and half-integer Fourier phase convention may be
an adaptation to its discrete potential and shifted window, but the broad claim
that a finite periodic auxiliary grid needs a width bound is not new.

### Cited Findings

- Jin, Liu, and Ma analyze Schrödingerisation on a finite auxiliary interval.
  They choose the interval width using the maximum spectral frequency, time,
  and target error, describe the approximately periodic boundary condition,
  and derive an exponential boundary-error estimate. Their recovery theorem
  also restricts the recovery coordinate to lie past a spectral cutoff.
  [Jin, Liu, and Ma, arXiv:2403.19123](https://arxiv.org/abs/2403.19123),
  especially Sections 3–4.
- That paper constructs smooth auxiliary initial data that remain equal to the
  required exponential on the recovery half-line while improving convergence
  in the auxiliary coordinate. Its numerical section uses a cubic blend
  between `e^{-|p|}` and the positive-side exponential. The present paper's
  smooth blend is related in purpose and already cites this method as a
  variant. [Full text, Sections 4.3 and 6](https://arxiv.org/html/2403.19123v4)
- The 2022 Schrödingerisation paper uses ordinary integer Fourier modes for
  the auxiliary transform. The 2024 finite-domain analysis likewise defines
  frequencies as `(i-N/2)/L`, with integer spacing for even `N`.
  [Jin, Liu, and Yu, arXiv:2212.13969](https://arxiv.org/html/2212.13969v1);
  [Jin, Liu, and Ma, arXiv:2403.19123](https://arxiv.org/html/2403.19123v4)
- The current paper instead centers an even-size frequency grid at
  `(k-(N_p-1)/2) pi/L_p`, which is half-integer shifted. The paper identifies
  the resulting antiperiodicity and reports an explicit-matrix check at four
  grid sizes. This is a convention-specific phase correction, not a property
  established by the cited periodic-grid papers.

### Inferences

- Reframe the contribution as a derivation and validation of the wrap and
  phase corrections for this particular half-shifted, even-size grid and
  shifted recovery window. Cite the finite-domain work as prior art for the
  general truncation and smooth-profile ideas.
- Compare the paper's bound
  `L_p >= (alpha_H T + p_a + omega + 20)/2` term by term with the existing
  spectral-cutoff boundary estimate. Explain how it bounds the actual wrapped
  contribution for the chosen profile. The two derivations use different
  interval conventions, so matching symbols alone is not sufficient.
- Do not describe the exact circuit costs of the diagonal Fourier phases or
  the finite-grid wrap error as validated by the matrix agreement alone.
  The cited work uses a standard periodic transform; the manuscript uses an
  antiperiodic convention and says its phase costs remain open.

### Gaps

- This search did not prove whether the exact half-integer transform, its
  circuit implementation, or the manuscript's particular shifted-window wrap
  formula appeared elsewhere.
- The relevant manuscripts and implementation need a direct derivation-to-code
  comparison before any finite-grid error guarantee is claimed.

## Is the all-orders lattice expansion unprecedented?

### Takeaway

All-orders expansions with uniform lattice remainders predate this manuscript,
and all-orders modified-equation expansions are a standard numerical-analysis
idea. The more specific result may be the explicit expansion for a Gaussian
solution under a time-varying quadratic potential, together with the sign
structure of its recovered-control error. The search found no exact match, but
that narrower novelty claim still needs a wider priority review.

### Cited Findings

- Gurevich derives arbitrary-order time-asymptotic expansions for Green
  functions of spatially discrete parabolic equations and bounds remainders
  uniformly over the entire lattice. The operators are finite-difference
  approximations to constant-coefficient strongly elliptic operators. This is
  a close precedent for uniform lattice asymptotics, though it is not the same
  small-mesh expansion for a time-dependent quadratic potential and Gaussian
  initial data. [Gurevich, arXiv:1504.02673](https://arxiv.org/abs/1504.02673)
- The current paper constructs polynomial corrections `p_k(x,s) psi(x,s)` for
  every order, bounds the residual by global derivative suprema and applies a
  lattice maximum principle. It further derives a leading correction with
  quartic and quadratic components and a radial sign-change statement. (Paper,
  Theorem `grid error` and Appendix proof.)
- Jin and Liu's 2025/2026 preprint already gives quantum algorithms for
  extracting pointwise values and gradients of nonlinear Hamilton–Jacobi
  viscosity solutions using an entropy-penalisation framework. This is close
  to the broad task of quantum HJ gradient readout, though its method and
  solution class differ from the paper's compatible-diffusion Cole–Hopf
  treatment of a quadratic execution model. [Jin and Liu, arXiv:2512.07919](https://arxiv.org/abs/2512.07919)

### Inferences

- Avoid claiming that arbitrary-order lattice expansions or pointwise
  Hamilton–Jacobi gradient extraction are new in general. State the special
  setting and exact theorem proved here.
- A defensible candidate contribution is the explicit coefficient hierarchy
  and the opposite-sign radial terms for this quadratic-potential Gaussian
  model, if the derivation and Lean correspondence survive independent review.
- The paper's cancellation-free finite-grid error estimate depends on this
  theorem, but its costed numerical constants are sampled floating-point
  estimates according to the manuscript. This literature comparison does not
  certify those constants.

### Gaps

- I did not find a source with the same hierarchy, exact coefficient formula,
  or radial sign theorem. The search was targeted rather than exhaustive and
  does not support an absolute priority claim.
- The literature scan does not verify the paper's theorem, constants, Lean
  formalization, or correspondence to the implementation.
