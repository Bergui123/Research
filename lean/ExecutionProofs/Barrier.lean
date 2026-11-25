import Mathlib

/-!
# The amplitude barrier (paper, Proposition 9)

With `ψ = exp (-W / 4ν)` the normalised probability of a grid point is exponentially small in
`1 / ν` unless the point is close to the minimum of `W`. If moreover the recovered control has a
relative bias of at least `β ν`, a bias below `ε` makes it exponentially small in `1 / ε`.
-/

namespace ExecutionProofs

/-- **Amplitude barrier.** -/
theorem amplitude_barrier {G : Type*} [Fintype G] {ν : ℝ} (hν : 0 < ν) (W : G → ℝ) (j k : G) :
    Real.exp (-W j / (4 * ν)) ^ 2 / ∑ i, Real.exp (-W i / (4 * ν)) ^ 2
      ≤ Real.exp (-(W j - W k) / (2 * ν)) := by
  have hk : Real.exp (-W k / (4 * ν)) ^ 2 ≤ ∑ i, Real.exp (-W i / (4 * ν)) ^ 2 :=
    Finset.single_le_sum (f := fun i => Real.exp (-W i / (4 * ν)) ^ 2)
      (fun i _ => sq_nonneg _) (Finset.mem_univ k)
  have hpos : 0 < Real.exp (-W k / (4 * ν)) ^ 2 := by positivity
  calc Real.exp (-W j / (4 * ν)) ^ 2 / ∑ i, Real.exp (-W i / (4 * ν)) ^ 2
      ≤ Real.exp (-W j / (4 * ν)) ^ 2 / Real.exp (-W k / (4 * ν)) ^ 2 :=
        div_le_div_of_nonneg_left (sq_nonneg _) hpos hk
    _ = Real.exp (-(W j - W k) / (2 * ν)) := by
        rw [← Real.exp_nat_mul, ← Real.exp_nat_mul, ← Real.exp_sub]
        congr 1
        field_simp
        ring

/-- A relative bias of at least `β ν` that must stay below `ε` turns the barrier into
`exp (-w β / 2ε)`. -/
theorem barrier_of_bias {β ε ν w p : ℝ} (hβ : 0 < β) (hν : 0 < ν) (hw : 0 ≤ w)
    (hbias : β * ν ≤ ε) (hp : p ≤ Real.exp (-w / (2 * ν))) :
    p ≤ Real.exp (-(w * β) / (2 * ε)) := by
  have hε : 0 < ε := lt_of_lt_of_le (mul_pos hβ hν) hbias
  refine hp.trans (Real.exp_le_exp.mpr ?_)
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith [mul_le_mul_of_nonneg_left hbias hw]

end ExecutionProofs
