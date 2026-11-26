# Research work plan

**The numerical pipeline has run; manuscript review and technical validation remain open.** For a target of 1% error in the initial optimal selling rate, the current provisional counts with Schrödingerisation are in the table of [Summary.tex](Summary.tex) (they are generated, so they are not repeated here). Readout A is amplitude estimation; readout B is sampling, which works only for quadratic models. A direct polynomial for the Hermitian generator needs fewer T gates than Schrödingerisation, and off the quadratic case readout A is exponentially expensive in 1/ν unless the state is recentred on a classical extremal. The full write-up is [paper/Paper.tex](paper/Paper.tex); [report/Report.tex](report/Report.tex) is a short summary. [coordination.md](coordination.md) and [review/](review/) hold the review record, [review/novelty.md](review/novelty.md) the literature check. The first figures (172 qubits, 4.0×10¹³ T) were withdrawn; WorkPlan.tex says why.

We estimate the logical qubits and T gates needed to solve a reformulation of a two-asset execution problem. The classical Riccati solution provides the reference answer.

An [alternative boundary-bound search](review/validation/alternative-boundary-selection.json) records how a conservative all-time boundary bound changes the readout-A selections, including one larger-box candidate. The production counts remain conditional on the numerical grid estimates.
For the normalized benchmark, a [rational Riccati residual calculation](review/validation/riccati-residual-bound.md) tightens an analytic lattice-remainder majorant. The remaining grid calculations and circuit costs still need validation.

- Step 1 is not needed. Step 2 starts from [Equation 16 of Bergault, Drissi and Guéant (2022)](https://arxiv.org/pdf/2103.13773v4#page=7), with R = 0.
- Step 3 verifies the notation, terminal conditions, sales-only proof and numerical benchmark.
- Step 4 tests Chen's auxiliary-variable method, explains its missing nonlinear constraint, and derives a compatible-diffusion Cole–Hopf equation that preserves the quadratic model's feedback. A finite-domain PDE solve validates it (second-order convergence).
- Steps 5–10 cover the choice of ν and the grid, Schrödingerisation, the block-encoding circuit and resource count, amplitude-estimation readout, and an end-to-end emulation of the costed pipeline.

[WorkPlan.tex](WorkPlan.tex) contains the derivations, sources and paper notes. [Summary.tex](Summary.tex) contains equations and one-line statuses. Update both together and include both in the same commit. The optional [symbolic HJB helper](stage2/BellmanPrinciple.py) requires SymPy.

Each stage also has a standalone, single-file document for Overleaf: [stage2/Stage2.tex](stage2/Stage2.tex), [stage3/Stage3.tex](stage3/Stage3.tex) and [stage4/Stage4.tex](stage4/Stage4.tex). Each contains its own merged copy of that folder's `derivation.tex`, the same content `WorkPlan.tex` separately `\input`s. **If you edit `derivation.tex`, re-merge it into the matching `StageN.tex` by hand** — they are copies, not shared live. Stage 4 cites [Chen (2001)](https://arxiv.org/pdf/math/0110005) and distinguishes that method from the compatible-diffusion transformation.

The numerics are in [simulation/](simulation/). Install the pinned packages from `simulation/requirements.txt` into `.venv`. Then, from `simulation/`, run `./run_all.sh` (the whole pipeline in order, about two hours on 16 cores), `../.venv/bin/python make_numbers.py` and `../.venv/bin/python figures.py`. Unit tests: `PYTHONPATH=. ../.venv/bin/python -m unittest discover -s tests`. `make_numbers.py` writes `paper/generated.tex`, which the paper, WorkPlan, Summary and report all read. Each script asserts its own checks and writes `results/<step>.json`. Figures go to `figures/`.

Every proposition, lemma and theorem of the paper has a formal counterpart proved in Lean 4 in [lean/](lean/), with Mathlib v4.36.0; Appendix B of the paper lists where the formal statement is weaker than the text. From `lean/`, run `lake exe cache get`, `lake build`, and `lake env lean scripts/AxiomCheck.lean`. The last command fails if a declaration uses `sorry` or any non-standard axiom. It also maps each statement to its Lean theorem. `paper/grid_theorem.tex`, `paper/grid_proof.tex` and `paper/grid_box.tex` are shared by the paper and WorkPlan.

Build each document with `pdflatex -interaction=nonstopmode -halt-on-error WorkPlan.tex` and the equivalent command for `Summary.tex`; repeat until references resolve. PDFs and build files are ignored by Git. Build the report from `report/` the same way; it reads its figures from `simulation/figures/`. Build the paper from `paper/`.

For a standalone stage document, upload just that one `StageN.tex` file to Overleaf (or run `pdflatex -interaction=nonstopmode -halt-on-error Stage2.tex`, etc., twice, locally from its own folder) — no other files needed.
