import Mathlib

/-!
# Taylor bound for the central second difference
-/

open Set Finset

namespace ExecutionProofs

/-- The coefficients `2 / (2 j + 2)!` of the expansion of the second difference. -/
noncomputable def c (j : ℕ) : ℝ := 2 / ((2 * j + 2).factorial : ℝ)

lemma c_pos (j : ℕ) : 0 < c j := by
  unfold c
  positivity

@[simp] lemma c_zero : c 0 = 1 := by
  simp [c]

lemma c_one : c 1 = 1 / 12 := by
  norm_num [c, Nat.factorial]

lemma c_two : c 2 = 1 / 360 := by
  norm_num [c, Nat.factorial]

/-- Lagrange form of Taylor's theorem at `0`, for either sign of `h`. -/
lemma exists_taylor {L : ℝ → ℝ} {m : ℕ} (hL : ContDiff ℝ (m + 1 : ℕ) L) {h : ℝ} (hh : h ≠ 0) :
    ∃ ξ, L h = ∑ k ∈ range (m + 1), iteratedDeriv k L 0 * h ^ k / (k.factorial : ℝ)
      + iteratedDeriv (m + 1) L ξ * h ^ (m + 1) / ((m + 1).factorial : ℝ) := by
  have hu : UniqueDiffOn ℝ (uIcc 0 h) := uniqueDiffOn_uIcc hh.symm
  have heq : ∀ k ≤ m + 1, ∀ x ∈ uIcc 0 h,
      iteratedDerivWithin k L (uIcc 0 h) x = iteratedDeriv k L x := fun k hk x hx =>
    iteratedDerivWithin_eq_iteratedDeriv hu
      ((hL.of_le (by exact_mod_cast hk)).contDiffAt) hx
  have hdiff : DifferentiableOn ℝ (iteratedDerivWithin m L (uIcc 0 h)) (uIoo 0 h) := by
    have h1 : Differentiable ℝ (iteratedDeriv m L) :=
      hL.differentiable_iteratedDeriv m (by exact_mod_cast Nat.lt_succ_self m)
    have h2 : DifferentiableOn ℝ (iteratedDerivWithin m L (uIcc 0 h)) (uIcc 0 h) :=
      h1.differentiableOn.congr fun x hx => heq m (Nat.le_succ m) x hx
    exact h2.mono uIoo_subset_uIcc_self
  obtain ⟨ξ, hξ, hT⟩ := taylor_mean_remainder_lagrange (f := L) (n := m) hh.symm
    ((hL.of_le (by exact_mod_cast Nat.le_succ m)).contDiffOn) hdiff
  refine ⟨ξ, ?_⟩
  rw [taylor_within_apply, heq (m + 1) le_rfl ξ (uIoo_subset_uIcc_self hξ), sub_zero] at hT
  have hsum : ∑ k ∈ range (m + 1), ((k.factorial : ℝ)⁻¹ * h ^ k) •
      iteratedDerivWithin k L (uIcc 0 h) 0
      = ∑ k ∈ range (m + 1), iteratedDeriv k L 0 * h ^ k / (k.factorial : ℝ) := by
    refine sum_congr rfl fun k hk => ?_
    rw [heq k (mem_range.mp hk).le 0 left_mem_uIcc, smul_eq_mul]
    ring
  rw [hsum] at hT
  linarith

lemma sum_range_two_mul (F : ℕ → ℝ) (m : ℕ) :
    ∑ k ∈ range (2 * m), F k = ∑ j ∈ range m, (F (2 * j) + F (2 * j + 1)) := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [Nat.mul_succ, sum_range_succ, sum_range_succ, ih, sum_range_succ]
    ring

/-- **Second difference.** The central second difference equals its Taylor sum up to the next
term, with the next coefficient. -/
theorem second_difference_bound {L : ℝ → ℝ} {J : ℕ} (hL : ContDiff ℝ (2 * J + 4 : ℕ) L) {M : ℝ}
    (hM : ∀ t, |iteratedDeriv (2 * J + 4) L t| ≤ M) {Δ : ℝ} (hΔ : Δ ≠ 0) :
    |(L Δ + L (-Δ) - 2 * L 0) / Δ ^ 2
      - ∑ j ∈ range (J + 1), c j * Δ ^ (2 * j) * iteratedDeriv (2 * j + 2) L 0|
      ≤ c (J + 1) * Δ ^ (2 * (J + 1)) * M := by
  obtain ⟨ξ, hξ⟩ := exists_taylor (m := 2 * J + 3) hL hΔ
  obtain ⟨η, hη⟩ := exists_taylor (m := 2 * J + 3) hL (neg_ne_zero.mpr hΔ)
  have hrange : 2 * J + 3 + 1 = 2 * (J + 2) := by ring
  rw [hrange, sum_range_two_mul, sum_range_succ'] at hξ hη
  set A := iteratedDeriv (2 * (J + 2)) L ξ with hA
  set B := iteratedDeriv (2 * (J + 2)) L η with hB
  have hΔ2 : Δ ^ 2 ≠ 0 := pow_ne_zero 2 hΔ
  have heven : ∀ k : ℕ, (-Δ) ^ (2 * k) = Δ ^ (2 * k) := fun k => by
    rw [pow_mul, neg_sq, ← pow_mul]
  have hodd : ∀ k : ℕ, (-Δ) ^ (2 * k + 1) = -Δ ^ (2 * k + 1) := fun k => by
    rw [pow_succ, heven, pow_succ]
    ring
  have hkey : (L Δ + L (-Δ) - 2 * L 0) / Δ ^ 2
      - ∑ j ∈ range (J + 1), c j * Δ ^ (2 * j) * iteratedDeriv (2 * j + 2) L 0
      = (A + B) * Δ ^ (2 * (J + 1)) / ((2 * (J + 2)).factorial : ℝ) := by
    have hadd : L Δ + L (-Δ) - 2 * L 0
        = Δ ^ 2 * ∑ j ∈ range (J + 1), c j * Δ ^ (2 * j) * iteratedDeriv (2 * j + 2) L 0
          + Δ ^ 2 * ((A + B) * Δ ^ (2 * (J + 1)) / ((2 * (J + 2)).factorial : ℝ)) := by
      rw [hξ, hη, mul_sum]
      simp only [heven, hodd]
      have hs : ∑ j ∈ range (J + 1), (iteratedDeriv (2 * (j + 1)) L 0 * Δ ^ (2 * (j + 1))
            / ((2 * (j + 1)).factorial : ℝ)
            + iteratedDeriv (2 * (j + 1) + 1) L 0 * Δ ^ (2 * (j + 1) + 1)
            / ((2 * (j + 1) + 1).factorial : ℝ))
          + ∑ j ∈ range (J + 1), (iteratedDeriv (2 * (j + 1)) L 0 * Δ ^ (2 * (j + 1))
            / ((2 * (j + 1)).factorial : ℝ)
            + iteratedDeriv (2 * (j + 1) + 1) L 0 * -Δ ^ (2 * (j + 1) + 1)
            / ((2 * (j + 1) + 1).factorial : ℝ))
          = ∑ j ∈ range (J + 1),
              Δ ^ 2 * (c j * Δ ^ (2 * j) * iteratedDeriv (2 * j + 2) L 0) := by
        rw [← sum_add_distrib]
        refine sum_congr rfl fun j _ => ?_
        have h2 : 2 * (j + 1) = 2 * j + 2 := by ring
        rw [h2, c]
        ring
      have hpow : Δ ^ (2 * (J + 2)) = Δ ^ 2 * Δ ^ (2 * (J + 1)) := by ring
      rw [hpow]
      simp only [mul_zero, pow_zero, zero_add, pow_one, Nat.factorial_zero, Nat.factorial_one,
        Nat.cast_one, div_one, mul_one, iteratedDeriv_zero]
      linear_combination hs
    rw [hadd, add_div, mul_div_cancel_left₀ _ hΔ2, mul_div_cancel_left₀ _ hΔ2]
    ring
  have hc : c (J + 1) = 2 / ((2 * (J + 2)).factorial : ℝ) := by
    unfold c
    rfl
  rw [hkey, hc, abs_div, abs_mul, Nat.abs_cast, heven_abs]
  · have hAM : |A + B| ≤ 2 * M := by
      have h1 := hM ξ
      have h2 := hM η
      have h3 := abs_add_le A B
      have e : 2 * (J + 2) = 2 * J + 4 := by ring
      rw [hA, hB, e] at h3 ⊢
      linarith
    have hp : 0 ≤ Δ ^ (2 * (J + 1)) := by
      rw [pow_mul]
      positivity
    calc |A + B| * Δ ^ (2 * (J + 1)) / ((2 * (J + 2)).factorial : ℝ)
        ≤ 2 * M * Δ ^ (2 * (J + 1)) / ((2 * (J + 2)).factorial : ℝ) := by gcongr
      _ = 2 / ((2 * (J + 2)).factorial : ℝ) * Δ ^ (2 * (J + 1)) * M := by ring
where
  heven_abs {Δ : ℝ} {k : ℕ} : |Δ ^ (2 * k)| = Δ ^ (2 * k) := by
    rw [pow_mul]
    exact abs_of_nonneg (by positivity)

end ExecutionProofs
