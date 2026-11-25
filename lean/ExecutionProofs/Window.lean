import Mathlib

/-!
# The shifted recovery window (paper, Proposition 5)
-/

open MeasureTheory intervalIntegral

namespace ExecutionProofs

/-- One mode: a profile equal to `exp (-p)` on `p ≥ 0` shifts exactly for `p ≥ -Emin s`. -/
lemma window_mode {g : ℝ → ℝ} (hg : ∀ p, 0 ≤ p → g p = Real.exp (-p)) {E Emin s p : ℝ}
    (hE : Emin ≤ E) (hs : 0 ≤ s) (hp : -(Emin * s) ≤ p) :
    g (p + E * s) = Real.exp (-p) * Real.exp (-(E * s)) := by
  have h : 0 ≤ p + E * s := by
    have : Emin * s ≤ E * s := mul_le_mul_of_nonneg_right hE hs
    linarith
  rw [hg _ h, ← Real.exp_add]
  congr 1
  ring

/-- **Shifted window.** If every mode shifts, `w_E(s, p) = g (p + E s) a_E`, then
`w(s, p) = exp (-p) ψ(s)` for all `p ≥ -Emin s`. -/
theorem shifted_window {M W : Type*} [Fintype M] [AddCommGroup W] [Module ℝ W] {g : ℝ → ℝ}
    (hg : ∀ p, 0 ≤ p → g p = Real.exp (-p)) (E : M → ℝ) (φ : M → W) {Emin s p : ℝ}
    (hE : ∀ m, Emin ≤ E m) (hs : 0 ≤ s) (hp : -(Emin * s) ≤ p) :
    ∑ m, g (p + E m * s) • φ m = Real.exp (-p) • ∑ m, Real.exp (-(E m * s)) • φ m := by
  rw [Finset.smul_sum]
  refine Finset.sum_congr rfl fun m _ => ?_
  rw [window_mode hg (hE m) hs hp, smul_smul]

/-- Weight of the window `[-cT, -cT + ω]`: the numerator of the success probability. -/
theorem window_weight (c T ω : ℝ) :
    ∫ p in (-(c * T))..(-(c * T) + ω), Real.exp (-p) ^ 2
      = Real.exp (2 * c * T) * (1 - Real.exp (-2 * ω)) / 2 := by
  have hderiv : ∀ p ∈ Set.uIcc (-(c * T)) (-(c * T) + ω),
      HasDerivAt (fun p => -(Real.exp (-2 * p)) / 2) (Real.exp (-p) ^ 2) p := by
    intro p _
    have h := (((hasDerivAt_id p).const_mul (-2 : ℝ)).exp.neg).div_const 2
    refine h.congr_deriv ?_
    rw [← Real.exp_nat_mul]
    simp only [id, mul_one]
    push_cast
    ring_nf
  rw [integral_eq_sub_of_hasDerivAt hderiv]
  · have h1 : Real.exp (-2 * (-(c * T) + ω)) = Real.exp (2 * c * T) * Real.exp (-2 * ω) := by
      rw [← Real.exp_add]
      congr 1
      ring
    have h2 : Real.exp (-2 * -(c * T)) = Real.exp (2 * c * T) := by
      congr 1
      ring
    rw [h1, h2]
    ring
  · exact (by fun_prop : Continuous fun p : ℝ => Real.exp (-p) ^ 2).intervalIntegrable _ _

end ExecutionProofs
