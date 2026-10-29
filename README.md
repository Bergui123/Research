# Research work plan

**Current: Step 3 verification. Next: complete Step 4.**

We estimate the logical qubits and T gates needed to solve a reformulation of a two-asset execution problem. The classical Riccati solution provides the reference answer.

- Step 1 is not needed. Step 2 starts from [Equation 16 of Bergault, Drissi and Guéant (2022)](https://arxiv.org/pdf/2103.13773v4#page=7), with R = 0.
- Step 3 checks the notation, terminal conditions, sales-only proof and numerical benchmark.
- Step 4 applies Cole–Hopf, solves the linear PDE and compares its recovered control with Riccati using the same finite terminal penalty.
- Steps 5–10 cover grid choices, Schrödingerisation, resource counts, readout and validation.

[WorkPlan.tex](WorkPlan.tex) contains the derivations, sources and paper notes. [Summary.tex](Summary.tex) contains equations and one-line statuses. Update both together and include both in the same commit. The optional [symbolic HJB helper](stage2/BellmanPrinciple.py) requires SymPy.

Build each document with `pdflatex -interaction=nonstopmode -halt-on-error WorkPlan.tex` and the equivalent command for `Summary.tex`; repeat until references resolve. PDFs and build files are ignored by Git.
