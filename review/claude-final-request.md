# Final independent review request for Claude Code

Please review the latest committed state on branch codex/arxiv-review-coordination, including paper/Paper.tex, WorkPlan.tex, Summary.tex, review/novelty.md, and this request.

Focus on:

1. The X4 issues and their fixes in review/codex-findings.md: direct-polynomial degree/cost, the Jin–Liu attribution and probability factor, fairness of the readout comparison, Table 12 interpretation, and abstract consistency.
2. The new related-work discussion of Wang, Hu and Zhang, arXiv:2608.09675. Check the primary paper at https://arxiv.org/abs/2608.09675 and confirm whether the distinction from this execution HJB pipeline is accurate.
3. The grid-expansion/sign-structure novelty claim, including whether the comparison with Gurevich 1504.02673 is accurate and whether arXiv:2602.16339 changes the assessment.
4. Whether the novelty claim is appropriately narrow and supported by the cited sources. Do not treat a targeted search as proof of priority.
5. Whether the refreshed arXiv package and code-path references follow the current official source/ancillary-file guidance.
6. The 26 November grid-estimate terminology correction: confirm no code, generated table, manuscript or plan now presents sampled floating-point maxima as certified enclosures, and check the baseline-equivalence evidence in `review/validation/check_grid_estimate_rename.py`.
7. Whether the new abstract and limitations language correctly makes readout-A counts conditional on the sampled grid maxima covering the true spatial and temporal extrema, without qualifying unrelated results too broadly.
8. The separated amplitude-share check in `review/validation/check_amplitude_shares.py`: verify the common reference, decomposition and interpretation, especially the 5% benchmark's 0.964561 fraction of the Schrödingerisation share. The saved values are in `review/validation/amplitude-shares.json`.
9. The periodised-exponential boundary bound and its rational evaluation in `review/validation/boundary-bound.md` and `check_boundary_bound.py`: check the discrete comparison at the periodic cusp, the Gaussian majorant, exact arithmetic directions and applicability to first outside nodes. In particular, check the 0.1% benchmark comparison and its qualification as a numerical grid-error estimate rather than an enclosure.
10. The temporary-file alternative search and emulation in `check_alternative_boundary_search.py` and `check_alternative_boundary_emulation.py`: check that all original candidates are treated consistently, the larger-box substitution and its limited search scope are clear, and the new estimates are not presented as certified counts.
11. The coarse analytic B2 majorant in `review/validation/b2-global-majorant.md`: check the Riccati ordering, l1 logarithmic-norm argument, hierarchy recurrence, Gaussian weights, exact arithmetic and its normalized-variable scaling. Its large size does not validate the resource counts.
12. The rational Riccati residual certificate and tightened B2 majorant in `review/validation/riccati-residual-bound.md`: check Bernstein conversion, interval positivity, the Frobenius error inequality, row and eigenvalue bounds, coefficient recurrence and normalized scaling. Check the benchmark-only reselection and qualifications; the original search used floating-point expansion and box evolution.


13. The polynomial residual calculation in `review/validation/polynomial-residual-bound.md`: independently check the coefficient-l1 logarithmic norm, all forcing factors and error recurrences, the Bernstein interval arithmetic, exponential tail, rational stencil coordinates and saved endpoints. Its 12 expansion enclosures do not by themselves validate the search's costs.


14. The projected box comparison and re-costing in `review/validation/projected-box-bound.md`: verify the Feynman–Kac/union comparison, boundary-value factor, substochastic uniformization, fixed-point rounding and tail, rational log/square-root rate calculation, exact grid coordinates and re-costed amplitude budgets. Check whether the evidence supports the three stated grid bounds while preserving the circuit and sweep limitations.


15. The optional controlled SELECT and small controlled-walk check in `review/validation/controlled-walk.md`: audit the control-zero behavior, signs, inverse, uncomputation, reflected-walk decomposition, control arity, Qualtran cost and what the matrix check does and does not establish. Check the 18 static controlled-SELECT counts in `controlled-profiles.json`, including whether workspace and allowance comparisons are like-for-like. Do not infer that GQSP or the complete QAE state-preparation circuit is compiled.


16. The ideal auxiliary Fourier-gate construction in `review/validation/fourier-circuit.md`: check bit ordering, both diagonal phase factors, global phases, inverse gate order, evolution sign, shifted window, antiperiodic wrap and the stated distinction between the analytic image bound and unvalidated interpolation aliasing.


17. Check the added classical execution-HJB citation to Forsyth (2011), its description, and the updated novelty search record; assess the breadth of the remaining no-match claim.


18. Review `review/validation/b2-polynomial-majorant.md`: check that the terminal nonnegative coefficient-error polynomials bound all earlier times, the Gaussian weights and both directional factors give the stated rational `B2<56`, and the new stencils match the rational box bounds and re-costed profiles. Check the search's numerical and global-optimality qualifications and whether the earlier `B2<875` evidence is preserved.


19. Check `review/validation/sweep-riccati.md` and its six rational certificates: the parameterized Riccati entry ordering, Bernstein positivity/residual checks, generalized coefficient error propagation, Gaussian weights and scaling with nu. Verify the numerical re-selection scope, especially the absence of a lambda=4 candidate, and check that the production sweep claims remain conditional.

Give concrete findings with severity, file and line, evidence, and a corrective action. If you find no issue, say so and state what you could not verify. Do not edit the source files or run commands that write to simulation/results/. Save the response as Markdown, preferably review/claude-final-response.md, and append a short status to coordination.md.
