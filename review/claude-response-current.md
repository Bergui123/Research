# Independent review of the optimal-execution manuscript

**Reviewer:** Claude Code (Opus 5.5), read-only, 24 November 2026.

**Independence.** `coordination.md` assigns `paper/`, `simulation/` and `lean/` to another Claude session. I have no memory of that session, but this is the same model, so it may not meet the `AGENTS.md` exit criterion. I read F1–F3 in `coordination.md` before checking them, so my confirmations are not blind.

**Snapshot.** The tree changed while I worked. Line numbers were rechecked after the 21:25 log entry. `nonquadratic.py` appeared last, and I saw no saved output for it.

**Record.** `review/claude-review.md` and `review/claude-response.md` should stay unchanged as the record of the rate-limited attempt. Save this review separately.

## Overall assessment

The mathematics is in good shape. I re-derived Proposition 1, Lemma 2, Theorem 3, the box and window results, and the finite amplitude budget, and found no error. The Lean statements match the paper within Appendix B's qualifications.

It is not ready for arXiv, for three reasons:
- The headline counts carry avoidable factors: about 1.7–2.0× in both readouts, and about 1.9× more in readout B at 1%.
- Several prose claims contradict the saved results.
- The non-quadratic claim for readout A is unsupported.

The counts are otherwise conservative. I found only small undercounts (findings 8 and 9).

## Status of prior findings

| Finding | Status in current tree |
|---|---|
| Codex 1, amplitude composition | **Fixed.** `step8_resources.py:79-81,159-160,175-177`; algebra re-derived. Still conditional on the shares (finding 9). |
| Codex 2, grid certification | **Relabelled, not fixed.** 33 times and mesh 0.04 remain. |
| Codex 3, sampling confidence | **Qualified, not fixed.** |
| Codex 4, QAE sampler | **Fixed.** `step10_validation.py:43-52`; law re-derived. The test reuses the same formula. |
| Codex 5, Fourier convention | **Partly fixed.** See checks below. |
| Codex 6, comparator endpoint | **Fixed in code.** `WorkPlan.tex:785` still says `2^k ≥ max`. |
| Codex 7, bibliography | **Fixed in paper.** `WorkPlan.tex:1026-1027` still puts the preprint title on the PRA reference. |
| Codex 8, packaging | Not verifiable; I could not open `arxiv-draft.zip`. |
| F1, F2 (`coordination.md`) | **Confirmed, still open in `Paper.tex`** (findings 3, 4). |

## Major findings

**1. The LCU carries identity terms that double α_H.** *Presentation/accounting; "cheapest" is unsupported.*
- **Where:** `Paper.tex:547-590`; `step8_resources.py:122,171`; `step9_sampling.py:77-78`.
- **Evidence:** H_d = (α_H/2)I − (ν/Δx²)Σ(S_a+S_a†) − ½V′B. The non-identity part has LCU norm 4ν/Δx² + V′/2 = α_H/2.
  - The scalar contributes e^{iηα_H T/2}, a translation in p. That is the same free relabelling used at `Paper.tex:464-465`.
  - The Chebyshev remark at 536-538 already needs only degree ½α_Hη_maxT.
  - Keeping today's overhead, A at 1% drops from 6057 to about 3240 steps and B from 3332 to about 1860. A saves 1.7× (5%) to 2.0× (0.1%). Qubits `k2` and `s` disappear.
- **Fix:** use α = 4ν/Δx² + V′/2 and shift the window by α_H T/2 mod 2L_p. Keep the L_p rule in terms of the physical spectrum, then regenerate.

**2. The "6.4×" window gain is not realised in any costed configuration.** *False.*
- **Where:** `Paper.tex:468-469`; `WorkPlan.tex:695-696`; `Summary.tex:134`; `Report.tex:58-59`.
- **Evidence:** the window starts at the first multiple of Δp above −0.930.
  - Saved ratios against p_a = 0 are 0.24690/0.05136 = **4.81** (`step7.json:16286,16305`) and 0.28311/0.08071 = **3.51** (`:16001,16020`).
  - By the same calculation the gain is 3.5 for A 5%, B 5% and B 1%; 4.8 for A 1% and 0.1%; 5.6 for B 0.1%.
- **Cost effect:** for B at 1%, −3Δp = −0.9416 is admissible because −E_min T = −0.9780 (`window_check.json:8`). Only the 5% margin excludes it.
  - Adding it raises the window sum from 7.50 to 14.08 and P_succ from 0.281 to 0.528.
  - Runs fall from 729,685 to about 389,000, and T from 6.4×10¹² to about **3.4×10¹²**.
- **Fix:** start at the smallest grid multiple ≥ −E_min T, which `check_window.py` already verifies. Report realised gains.

**3. Success probability still decays exponentially in M, and sampling is not roughly linear.** *False.*
- **Where:** `Paper.tex:60-61,98-100,886-887`; `Report.tex:73-74`; `WorkPlan.tex:983-984`.
- **Evidence:** each mode μ = 0.6 contributes e^{0.95κ}(1/y)√(φ/p) = 0.96189.
  - With the collective mode (0.78718), P(10) = ½·0.78718·0.96189⁹ = 0.27744, matching `scaling_m.json`.
  - Almost all the decay is the 5% margin; an exact shift gives 0.9999 per mode here.
  - Runs B grow 13.98× from M = 1 to 10. The local exponent from 8 to 10 is 1.54.
- **Fix:** state the per-asset factor and say the shift changes the base of the exponential.

**4. "Applicable to non-quadratic problems" is unsupported.** *Unproved; misleading beside the quoted count.*
- **Where:** `Paper.tex:51-53,714,908-909`.
- **Evidence:** p_j ≤ exp(−(W_j−W_min)/2ν), so runs scale at least as exp((W_j−W_min)/4ν).
  - With W(X) = 11.0005 this is 1.4 at ν = 8, 16 at ν = 1, 9×10⁵ at ν = 0.2 and 9×10¹¹ at ν = 0.1.
  - All nine selected readout-A configurations use ν = 4–16, which is allowed only because the LQ feedback does not depend on ν.
  - The budget also needs the exact solution: `grid_error.py:63-67`, `certify_grid.py:75-87`, `step8_resources.py:158`.
- **Fix:** say the circuit carries over, but the count and the error budget do not. Add the barrier proposition.

**5. The abstract says "our budgets" avoid known-solution selection, but readout B does not.** *False as written; B's systematic error is empirical.*
- **Where:** `Paper.tex:45-47` against `step9_sampling.py:72-75`. Line 919 concedes B's error is "measured".
- **Evidence:** B's grid and domain are chosen by `rel(p_pde, P)` at the exact P, with no cancellation-free analogue.
- **Fix:** write "readout-A budgets", and either derive an estimate for B or disclose the selection.

## Moderate findings

**6. The factor 42 compares unlike guarantees.** *Presentation.* `Paper.tex:770-771`, `Summary.tex:168`.
- A is a conditional worst-case bound times 11. B is an expected run count at an asymptotic quantile.
- A single A set already succeeds in 96.0, 97.5 and 97.9% of trials. The single-set ratio is 3.8 at 1% and 2.0 at 0.1%.
- **Fix:** report both ratios.

**7. "At least 98% of trials" is contradicted by Table S.** *False.* `Paper.tex:834-836`.
- Saved values are 0.9775 (`rho0.4_lam4/step9_sampling.json:46`) and 0.9795 (`rho0.6_lam1`, printed as 98.0).
- **Fix:** generate this number by macro.

**8. Readout B omits its coherent error from the budget.** *Unproved.* `step9_sampling.py:86-91`.
- `eps_state` allows a second-moment error of 0.01εσ_min, at most 0.0126ε in the stated metric (√cond P = 1.264). The shot count still uses (ε−ρ).
- The worst-case undercount is 3.1% at 5% and 1%, and 5.0% at 0.1%. Section 7.2 of the paper never mentions this budget.
- **Fix:** subtract and document it.

**9. The 10% simulation share is not a proven bound.** *Unproved.* `step8_resources.py:126-130`.
- GQSP needs |P| ≤ 1. Rescaling the truncated series can double the tail, so the total is up to 1.5 times the share. Grover–Rudolph angle errors are also outside it, though negligible.
- **Fix:** use a tail of `eps_state/(4 n_p)`.

**10. K_j+3 is cited, not justified.** *Unproved here; I could not open the source.* `Paper.tex:610-614`.
- My own construction supports K+O(1). K directional steps give G = D e^{iτcosθ X}D† with D = diag(1, W†), so K+2 controlled uses with the signal in |+⟩.
- The truncated cos/sin pair is not exactly normalised, and no argument covers the completion.
- If directional control failed, the cost would revert to the 2K column.
- **Fix:** state the polynomial pair and the exact lemma.

**11. Lean scope is overstated outside Appendix B.** *False.*
- `README.md:18` says every result is proved in Lean; Appendix B lists exceptions.
- `Summary.tex:101-112` labels Eq. (cert) machine-checked, but I found no such declaration.
- `norm_ratio` takes arbitrary P and D, not P_T and det Y(T). P_T ≻ 0 is formalised only for H = I.
- **Fix:** align the wording and add both gaps to Appendix B.

**12. The search space is misdescribed, and some optima sit on its boundary.** *False caption.*
- `Paper.tex:730` says h ∈ {1,…,16}; the code uses {1,2,4,8,16} (`step5_nu_sweep.py:27`).
- The 0.1% optimum is at the maximum N and h. Sweep points (0.2, 1) and (0.4, 4) are at the maximum ν. Only three p-grids are tried (`step8_resources.py:273`).
- **Fix:** widen the search or write "best found".

**13. Reproducibility gaps.** *Presentation.*
- `boundary_check.json` has no generating script; only `make_numbers.py:139` reads it.
- `research-sources.zip` (`Paper.tex:1036`) is not in the tree.
- "Tens of microseconds" (`:83`) has no timing script.
- `Paper.tex:558` says the arithmetic is uncompiled, but `compile_select.json` exists and is not in the run list.

## Minor findings

- **Stale numbers, false:** `README.md:3` says 6.5×10¹² (now 6.4×10¹²). `WorkPlan.tex:345-347` has old R, P and n_s; saved values are 729,685, 0.28144 and 205,366.
- **Domain rule, false:** `Paper.tex:412-414` says a box 1.6 times larger. The code compares with one fixed box, 1.6 times the largest candidate (`step5_nu_sweep.py:79-99`); at ν = 1, L = 6 was tested against 21.6.
- **Spectral weight, false if sweeps are included:** `Paper.tex:489` says up to 10⁻⁵; λ = 4 gives 7.7×10⁻⁵.
- **Table 1, presentation:** the "sup" is over |ξ| < 9 on a periodic box (`grid_theorem.py:151-158`).
- **Notation, presentation:** `c` is the window shift (465), then (N_p−1)/2 (478), then p_a = −cT (497). ρ is both the correlation and B's systematic error.
- **Novelty, presentation:** contribution 2 lists the shift; `WorkPlan.tex:758-760` claims no novelty for it.
- **Figures and abstract, presentation:** the `tcount.png` legend covers a data point, and the Figure 1 legend says "bound" where the text says "estimate". The abstract says the counts are "provisional pending" checks.

## Checks performed (reading and hand calculation only)

- **Derivations:** constants 1/3072, 1/128 and 1/64; Eq. (mono); the sign argument; box, window, comparator and rounding results; Eq. (verification).
- **Readout:** median of 11 gives 0.965 joint (9 gives 0.939). The QAE bound follows from |sin θ̃ − sin θ| ≤ π/M.
- **A at 1%, recomputed:** δ, all four t_j, eps_rot, b_rot, oracle 3176 T, T_run, total T and 190 qubits. **B at 1%:** n_s, eps_rot and total T. All match the saved values.
- **Circuit:** four rotations and two reflections per step are right, because PREP·PREP† cancels between steps.
- **Fourier phases:** the η-side phases cancel and the final p-phase commutes with the readout. The 2n_p rotations at `step8_resources.py:143` suffice.
- **Wrap error:** at most e⁻²⁰√(3/‖g‖²) ≤ 3.6×10⁻⁹ per amplitude, about 7% of the Schrödingerisation share at worst (A, 0.1%). `Wrap.lean` bounds images of g, not the interpolant the circuit uses.
- **Lean:** read all 22 modules; no `sorry`, `native_decide` or added axiom found.
- **`compile_select.py`:** read, no error found. It covers SELECT only, at 0.6–5.8% below the charge.

## Not independently verified

- No script, test, `lake build`, axiom check or LaTeX build was run.
- **No reference metadata or prior work was checked.** WebFetch was redirected by a hook, and the indexed fetch and WebSearch were not permitted. From memory only, arXiv:2403.10032 gives explicit Schrödingerisation circuits, so check it against `Paper.tex:125`. I make no priority claim.
- QROM, QFT, reflection and GQSP compilation; true grid suprema; Ross–Selinger constants; archive contents; a line-by-line match of the stage copies.

## Unresolved questions

1. What exactly does arXiv:2401.10321 prove, and for which polynomial pair?
2. Do the optima move after findings 1 and 2?
3. Recentring leaves the generator unchanged only for quadratic running cost, which is the case `nonquadratic.py` tests. Otherwise it is time dependent and may be unbounded below, so the cost model and the V ≥ 0 lemmas do not transfer. That script also selects by error at the known answer (lines 109-111, 136-168).
4. `barrier_of_bias` assumes the bias is at least βν. Is there a proof for any non-quadratic case?

Skipped: everything needing execution or network access. Main risk: the tree is still changing, so recheck line numbers before acting, and do not quote the T counts until findings 1 and 2 are settled.
