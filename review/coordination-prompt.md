# Claude Code review handoff

Read `coordination.md`, `AGENTS.md`, `README.md`,
`review/claude-request.md`, and `review/codex-review.md` first. Then inspect the
current `paper/` source and included TeX files, `WorkPlan.tex`, `Summary.tex`,
`report/`, `simulation/`, generated tables and figures, and `lean/`. The earlier
review files describe a previous version. Check the present working tree and
identify which prior findings have actually been fixed by code and evidence.

Act as an independent, skeptical mathematical and quantum algorithms referee.
Focus on concrete issues in derivations, numerical bounds and their
implementation, error composition, circuit and resource accounting, the
correspondence between paper statements and Lean declarations, and manuscript
consistency. Cross-check source against code and saved results. A successful
build, tests, or a theorem with a similar name is not proof of correctness.

For every finding, provide severity, exact file and line, the evidence or
calculation, whether it is false, unproved, empirical, or only a presentation
problem, and a precise corrective action. Explicitly say which areas you could
not independently verify. Check reference metadata and close prior work where
that affects claims, but do not make a priority claim from a search with no
match. Do not edit files, run commands, or report checks you did not perform.

Return a self-contained Markdown review with a short overall assessment,
prioritized findings, checks performed, and unresolved questions. Preserve the
earlier rate-limited attempt as a separate historical record. This is an
independent review, not an instruction to agree with the existing Codex review.
