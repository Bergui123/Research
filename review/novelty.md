# Novelty record

Date: 24 November 2026. Method: arXiv API queries (titles, abstracts, full metadata), web
searches, and reading of the closest papers in full text. Codex searched independently
(`reports/Quantum execution algorithm novelty.md`, `research_notes/`) and reached the same
verdict. A search cannot prove priority; the paper says "we found no ..." and gives the date.

## Verdict per claim

| Claim in the paper | Verdict | Closest prior work |
|---|---|---|
| Quantum treatment of optimal execution | **not new**, not claimed | Rosenberg et al. 1508.06182 (annealing, discrete trajectories) |
| Quantum algorithm for LQ control / Riccati | **not new**, not claimed | Krovi 2404.07303; Huang et al. 2610.26843; Wang et al. 2609.37207; Dehaghani et al. 2507.09841 |
| Cole-Hopf + quantum PDE solver, HJ gradient readout | **not new**, not claimed | Jin, Liu 2512.07919; Jin, Liu, Yu 2606.28752; Uchida et al. 2412.17206 |
| Shifted recovery window, wide p domain, smooth profile | **not new**, not claimed (dropped from the contributions) | Gu, Jin 2501.10176 (spectral shift, domain rule); Jin, Liu, Ma 2403.19123; Jin et al. 2505.00370; Nguyen, Tura 2611.02357 |
| Grid error expansion to every order with uniform remainder | **not new in general** | Gyongy, Krylov 0905.3299 (bounded coefficients, powers of the mesh); Gurevich 1504.02673 (lattice heat kernels, powers of time) |
| Sign structure of the leading control error (Theorem 3(iv)) and its consequence for resource estimates | no earlier source found | |
| Logical qubit / T count for execution, or for an HJB equation by Schrodingerisation | no exact match found in the targeted search | adjacent full gate-level PDE algorithm: Wang, Hu and Zhang 2608.09675; derivative pricing: Chakrabarti 2012.03819, Stamatopoulos-Zeng 2307.14310; ODE solvers: Jennings 2309.07881, Pocrnic 2506.20760; fluids: Penuel 2406.06323; QHD: 2607.16996 |
| Quadratic saving of a direct polynomial for Hermitian semigroups | **not new** | Sachdeva-Vishnoi 2014; Gilyen et al. 1806.01838; Jin, Ma, Zuazua 2601.03616. New here: the measured factor for a complete pipeline |
| Amplitude barrier exp(-(W_j - min W)/2nu) | elementary; no earlier statement found for Cole-Hopf readout | Jin-Liu 2512.07919, Lemma 12: its proof derives $J=(1/q+1/\|u\|^2)^2/\|u\|^2$, where $q=|\langle x|u\rangle|^2$; the statement omits this factor (checked by Codex, X4) |
| Recentring along a classical extremal for a Cole-Hopf quantum algorithm | identity is standard; no earlier use in a quantum algorithm found | displaced frames: Heller 1975, Hagedorn 1980; importance sampling: Kappen-Ruiz 1505.01874 |
| Lean proofs | the formalisation itself | |

## Queries used

arXiv API: `cat:quant-ph` with "optimal execution", "trade execution", "Hamilton-Jacobi-Bellman",
"Cole-Hopf", "linearly solvable", "viscosity solution", "Riccati", "stochastic control",
"mean-field games", "Doob", "displaced frame", "moving frame" + imaginary time, "warm start" +
differential equation, "Feynman-Kac" + amplitude estimation, "resource estimate" + differential
equations; "Schrodingerization" with "resource", "T gates", "shift", "finance"; `q-fin.TR` and
`q-fin.CP` with "quantum"; "path integral control" + "importance sampling"; "linear combination of
Hamiltonian simulation" + resource. Web: the same topics plus Gyongy-Krylov expansions, quantum
limitations of Cole-Hopf / Burgers, quantum LQR, Schrodingerisation circuits and T counts.

## References checked against arXiv metadata

All arXiv identifiers in the bibliography were resolved through the API and their titles and
authors compared (two batches, 24 November 2026). Journal data of Fleming-Souganidis 1986,
Heller 1975 and Hagedorn 1980 were checked through secondary reference lists, not the originals.

## What was not done

No citation-chain crawl of the closest papers, no theses, patents or non-English sources. The
uniqueness of the "no earlier source found" rows rests on keyword searches and on reading the
closest papers.

Final web pass on 24 November 2026 searched exact combinations of "quantum optimal execution"
with Cole-Hopf, Schrödingerisation, HJB and logical resource estimation. It also returned a
general fault-tolerant algorithm-to-execution and resource-planning framework
([arXiv:2609.10965](https://arxiv.org/abs/2609.10965)); that paper is relevant to the resource
accounting context but does not match this execution/HJB pipeline. This pass did not include a
citation-chain crawl and does not establish priority.

## Follow-up primary-source check, 24 November 2026

A fresh search returned [Wang, Hu and Zhang, arXiv:2608.09675](https://arxiv.org/abs/2608.09675),
submitted 10 August 2026. I opened its arXiv HTML full text. It lifts scalar conservation laws
to a linear Liouville transport equation by a level-set method, discretizes by finite
differences, applies Schrödingerisation, and gives gate-level circuits, observable recovery,
error bounds and gate-complexity estimates for the complete algorithm. This establishes close
prior art for explicit end-to-end gate-level Schrödingerisation algorithms for nonlinear PDEs.
The formulation is distinct from the continuous-time optimal-execution HJB, its Cole-Hopf
heat equation, and the feedback readout costed here. Accordingly, the manuscript now cites it
as adjacent work and limits its no-match statement to optimal execution and HJB by
Schrödingerisation. The targeted search found no exact match for that narrower combination;
this is not a priority proof. No citation-chain crawl was done.

## Focused grid-theorem search, 24 November 2026

A second search targeted the paper's grid expansion and leading rate-error sign using phrases
combining Gaussian initial data, an unbounded quadratic potential, finite-difference expansion,
and optimal-control rate error. It returned Gurevich 1504.02673, already listed above, and
Chen and Li 2602.16339 on large-time asymptotics for lattice and fractional heat equations.
The latter concerns a different limit and model. The surfaced sources did not match the
paper's small-mesh expansion for this Gaussian/quadratic problem or its opposite-sign leading
rate-error terms. A broad query using only “radial error” returned unrelated image-registration
material and was discarded. This pass adds no priority evidence; it is a narrow search record,
not an exhaustive review or a proof that no source exists.

## Current arXiv search, 26 November 2026

I repeated searches on the official arXiv index, sorted by latest submission:
`"optimal execution" quantum` (8 results), `"Cole-Hopf" quantum` (15),
`quantum Hamilton-Jacobi optimal control` (15), and
`quantum Schrödingerisation gate complexity` (3). I inspected the result
titles and abstracts ([official result pages](https://arxiv.org/search/?query=%22optimal+execution%22+quantum&searchtype=all&abstracts=show&order=-announced_date_first&size=50),
[Cole-Hopf query](https://arxiv.org/search/?query=%22Cole-Hopf%22+quantum&searchtype=all&abstracts=show&order=-announced_date_first&size=50),
[HJB/control query](https://arxiv.org/search/?query=quantum+Hamilton-Jacobi+optimal+control&searchtype=all&abstracts=show&order=-announced_date_first&size=50),
[Schrödingerisation/gate query](https://arxiv.org/search/?query=quantum+Schr%C3%B6dingerisation+gate+complexity&searchtype=all&abstracts=show&order=-announced_date_first&size=50)).
The closest hits remain distinct from this paper's
execution-HJB pipeline: [Jin and Liu, 2512.07919](https://arxiv.org/abs/2512.07919)
gives quantum algorithms for pointwise values, gradients, minima and values at
minimizers of viscosity HJ solutions; [Jin, Liu and Yu, 2606.28752](https://arxiv.org/abs/2606.28752)
uses Cole-Hopf and Schrödingerisation for regularized Wasserstein proximal
operators; [Wang, Hu and Zhang, 2608.09675](https://arxiv.org/abs/2608.09675)
provides a gate-level Schrödingerisation algorithm for scalar conservation
laws. The targeted searches did not reveal an execution-specific HJB resource
estimate or the manuscript's small-mesh rate-error sign result. This supports
only a narrowly worded “no exact match found” statement, not priority: arXiv
full-text indexing, citation chains, theses, patents, non-English sources and
non-arXiv publications were not exhaustively searched. The manuscript's
contribution should remain its conditional logical resource estimate and its
stated Gaussian/quadratic grid-error result, not a general claim to the first
quantum treatment of HJB, Cole-Hopf or Schrödingerisation.
