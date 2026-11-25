# Quantum finance and control novelty

## Does quantum optimal execution already exist?

### Takeaway

Yes. Quantum annealing has been applied to a discrete, multi-period optimal
trading trajectory problem that includes transaction costs and market impact.
The narrow possible contribution here is a fault-tolerant logical resource
estimate for a continuous-time, HJB-based feedback calculation. The search did
not establish priority for that narrower claim.

### Cited Findings

- Rosenberg et al. formulated a multi-period portfolio optimization problem as
  a discrete optimization task and solved instances on a D-Wave quantum
  annealer. Their objective includes transaction costs, permanent and temporary
  market impact. The paper reports using a 512-qubit annealer and discusses
  physical-qubit embedding; a search of its full arXiv text found no T-gate
  resource estimate. [Original preprint](https://arxiv.org/abs/1508.06182)
- A 2026 preprint studies dynamic portfolio rebalancing with realistic
  transaction and allocation constraints, encodes asset weights as a QUBO,
  and uses a hybrid quantum-classical constrained optimization workflow. This
  is portfolio rebalancing, not a continuous-time execution PDE.
  [Original preprint](https://arxiv.org/abs/2607.03218)
- A 2026 preprint formulates a portfolio-rebalancing schedule as a QUBO and
  applies QAOA in a walk-forward evaluation. It concerns discrete schedule
  optimization and does not report fault-tolerant logical-qubit or T-gate
  counts. [Original preprint](https://arxiv.org/abs/2603.16904)

### Inferences

- The sentence at `paper/Paper.tex:128`, "not aware of one for optimal
  execution," should be bounded by the type of resource estimate. Read without
  the preceding word "resource," it could imply that no quantum optimal
  execution work exists, which the 2016 annealing paper contradicts.
- A defensible formulation would identify the particular gap: to the authors'
  knowledge, the cited search found no end-to-end fault-tolerant logical-qubit
  and T-gate estimate for a continuous-time, multi-asset execution feedback
  algorithm based on an HJB/PDE route. This is a scoped literature observation,
  not a proven priority result.

### Gaps

- Searches covered arXiv and publisher-indexed results using exact and varied
  phrases for quantum optimal execution, trading trajectory, HJB finance,
  market impact and resource estimates through 2026-11-24. This is not an
  exhaustive review of every database, conference proceeding, thesis, patent,
  non-English paper or unpublished project.
- The 2016 paper reports hardware and embedding resources, but I did not audit
  its entire experimental appendix for every hardware metric. The specific
  observation is only that the indexed full text contains no "T gate" match.

## How close are the quantum finance and control results?

### Takeaway

Quantum algorithms for the Riccati equation/LQR and for financial PDEs already
exist. A particularly relevant 2026 LQR preprint postdates the manuscript's
current bibliography. These works use different problem models or resource
accounting from the present paper, so they do not by themselves establish that
the manuscript's detailed benchmark cost is unoriginal.

### Cited Findings

- Krovi gives quantum algorithms for vector and matrix Riccati equations and
  applies them to the linear-quadratic regulator, described as a special case
  of the HJB equation. The paper gives query and gate complexity bounds, rather
  than the manuscript's concrete benchmark table of logical qubits and T gates.
  [Original preprint](https://arxiv.org/abs/2404.07303)
- Huang, Álvarez Yuste, Guseynov, Liu and González-Conde's 2026 preprint solves
  finite-horizon continuous-time LQR using a state-costate boundary-value
  system and a quantum linear-systems subroutine. It returns classical
  estimates of optimal cost and control at a selected time and gives query and
  gate complexity bounds under stated assumptions. It is not an optimal
  execution or finance-specific resource-count study.
  [Original preprint](https://arxiv.org/abs/2610.26843)
- Chakrabarti et al. give a detailed fault-tolerant resource estimate for
  derivative pricing, including complete circuit-level estimates for
  autocallable and TARF benchmarks. Their reparameterization examples report
  thousands of logical qubits and T-depths in the tens of millions.
  [Original paper](https://arxiv.org/abs/2012.03819)
- Wang et al.'s 2026 preprint uses quantum and quantum-inspired
  physics-informed neural networks to approximate a financial HJB/PDE for the
  Merton portfolio-optimization problem. It is a variational/PINN approach and
  does not present the same execution objective or a fault-tolerant logical
  resource count. [Original preprint](https://arxiv.org/abs/2604.03346)
- Fontanela, Jacquier and Oumgari apply a hybrid quantum imaginary-time method
  to Black--Scholes linear PDEs for European and Asian option pricing. Their
  construction is distinct from the manuscript's compatible-diffusion
  Cole--Hopf transformation for an execution-control HJB.
  [Original preprint](https://arxiv.org/abs/1912.02753),
  [journal record](https://doi.org/10.1137/21M1397878)

### Inferences

- The related-work sentence at `paper/Paper.tex:123` that describes Krovi is
  accurate at a high level, but it should be revisited to include the October
  2026 finite-horizon LQR preprint. Since that paper gives asymptotic query and
  gate bounds, avoid implying that only the current manuscript treats
  finite-horizon LQR or computes control outputs.
- The sentence at `paper/Paper.tex:127` about derivative-pricing resource
  estimates is supported by Chakrabarti et al. The sentence at line 128 can
  remain only if "resource estimate" is defined as a complete fault-tolerant
  logical-gate estimate for the specified continuous-time HJB/PDE method.
- The manuscript's strongest distinct claim in this comparison is the
  concrete resource accounting for its chosen execution pipeline, not the
  first quantum treatment of trading, the first quantum treatment of LQR, or
  the first quantum treatment of a financial PDE.

### Gaps

- I did not check every cited quantum-control or financial-PDE paper for
  hidden concrete hardware or fault-tolerant resource tables. The comparisons
  above are based on the linked primary sources and the parts surfaced by
  targeted searches.
- The manuscript's separate mathematical claims about the grid expansion,
  Schrödingerisation, auxiliary Fourier convention, and recovery window need
  the dedicated technical novelty review; this note does not establish their
  priority.

## What novelty wording is supportable?

### Takeaway

The current broad wording risks being read as a first quantum optimal
execution claim, despite prior annealing work. A narrower statement about a
specific fault-tolerant logical resource study is more accurate, but it should
be framed as the result of a documented search and not as proof that no such
paper exists.

### Cited Findings

- `paper/Paper.tex:127-128` states that derivative-pricing resource estimates
  exist and says the authors are "not aware of one for optimal execution."
  The cited derivative-pricing precedent is the complete resource analysis by
  Chakrabarti et al. [Source](https://arxiv.org/abs/2012.03819)
- The Rosenberg paper is explicitly about an "optimal trading trajectory"
  and includes market impact, but solves a discrete optimization problem on a
  quantum annealer rather than estimating fault-tolerant resources for an
  HJB/PDE feedback algorithm. [Source](https://arxiv.org/abs/1508.06182)
- Krovi's Riccati/LQR method and Huang et al.'s later finite-horizon LQR method
  already treat quantum optimal control and return control-related outputs;
  both present complexity bounds that are distinct from concrete benchmark
  logical-qubit/T-gate accounting. [Krovi](https://arxiv.org/abs/2404.07303),
  [Huang et al.](https://arxiv.org/abs/2610.26843)

### Inferences

- A suggested replacement is: "Quantum annealing has been applied to
  discrete multi-period trading trajectories. Fault-tolerant resource
  estimates have also been given for derivative pricing, and quantum
  algorithms for Riccati/LQR and financial PDEs have been studied. In the
  focused literature search described in the accompanying review notes, we
  found no end-to-end logical-qubit and T-gate estimate for a continuous-time,
  multi-asset optimal-execution feedback calculation using this HJB,
  Cole--Hopf and Schrödingerisation route. We therefore report a cost model for
  this benchmark, without claiming quantum advantage or priority beyond that
  scope."
- This wording should be used only after the broader novelty search and
  manuscript review finish. Search results alone cannot support an absolute
  "first" claim or certify novelty.

### Gaps

- Online search cannot establish absence of prior work. A scholar-led review
  should search bibliographic databases, proceedings and citations, and should
  examine the cited papers' reference lists before publication.
- No independent human literature review or author confirmation of the
  contribution statement has been completed in this search task.
