# Research work plan

**Step 4 derivation recorded. Step 3 verification and Step 4 numerical PDE validation remain open.**

We estimate the logical qubits and T gates needed to solve a reformulation of a two-asset execution problem. The classical Riccati solution provides the reference answer.

- Step 1 is not needed. Step 2 starts from [Equation 16 of Bergault, Drissi and Guéant (2022)](https://arxiv.org/pdf/2103.13773v4#page=7), with R = 0.
- Step 3 checks the notation, terminal conditions, sales-only proof and numerical benchmark.
- Step 4 tests Chen's auxiliary-variable method, explains its missing nonlinear constraint, and derives a compatible-diffusion Cole–Hopf equation that preserves the quadratic model's feedback. Finite-domain PDE validation remains open.
- Steps 5–10 cover grid choices, Schrödingerisation, resource counts, readout and validation.

[WorkPlan.tex](WorkPlan.tex) contains the derivations, sources and paper notes. [Summary.tex](Summary.tex) contains equations and one-line statuses. Update both together and include both in the same commit. The optional [symbolic HJB helper](stage2/BellmanPrinciple.py) requires SymPy.

Each stage also has a standalone document: [stage2/Stage2.tex](stage2/Stage2.tex), [stage3/Stage3.tex](stage3/Stage3.tex) and [stage4/Stage4.tex](stage4/Stage4.tex). Each one's `derivation.tex` is also `\input` by `WorkPlan.tex`, so the standalone document and the work plan share the same content. Stage 4 cites [Chen (2001)](https://arxiv.org/pdf/math/0110005) and distinguishes that method from the compatible-diffusion transformation.

Build each document with `pdflatex -interaction=nonstopmode -halt-on-error WorkPlan.tex` and the equivalent command for `Summary.tex`; repeat until references resolve. PDFs and build files are ignored by Git.

For a standalone stage document, run `pdflatex -interaction=nonstopmode -halt-on-error Stage2.tex` (or `Stage3.tex`/`Stage4.tex`) twice from its own folder.
