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

Give concrete findings with severity, file and line, evidence, and a corrective action. If you find no issue, say so and state what you could not verify. Do not edit the source files or run commands that write to simulation/results/. Save the response as Markdown, preferably review/claude-final-response.md, and append a short status to coordination.md.
