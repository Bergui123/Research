import Mathlib

/-!
# Wrap-around on the antiperiodic auxiliary grid (paper, Lemma 6)

With half-integer frequencies the interpolated profile is antiperiodic over `2 L`. Its images are
`(-1)^r g (y + 2 L r)`, `r ≠ 0`; a profile that decays like `exp (-|p|)` bounds them.
-/

namespace ExecutionProofs

/-- One pair of images. -/
lemma image_pair_le {g : ℝ → ℝ} {L R y : ℝ} (hL : 0 < L)
    (hg : ∀ p, R ≤ |p| → |g p| ≤ Real.exp (-|p|)) (hy : |y| ≤ 2 * L - R) (r : ℕ) :
    |g (y + 2 * L * (r + 1)) + g (y - 2 * L * (r + 1))|
      ≤ 2 * Real.exp (-R) * Real.exp (-(2 * L)) ^ r := by
  have hr : (0 : ℝ) ≤ r := Nat.cast_nonneg r
  have hLr : 0 ≤ 2 * L * r := by positivity
  obtain ⟨hy1, hy2⟩ := abs_le.mp hy
  have hexp : Real.exp (-R) * Real.exp (-(2 * L)) ^ r = Real.exp (-(R + 2 * L * r)) := by
    rw [← Real.exp_nat_mul, ← Real.exp_add]
    congr 1
    ring
  have hplus : R + 2 * L * r ≤ |y + 2 * L * (r + 1)| := by
    refine le_trans ?_ (le_abs_self _)
    nlinarith
  have hminus : R + 2 * L * r ≤ |y - 2 * L * (r + 1)| := by
    refine le_trans ?_ (neg_le_abs _)
    nlinarith
  have h1 : |g (y + 2 * L * (r + 1))| ≤ Real.exp (-(R + 2 * L * r)) :=
    (hg _ (by linarith)).trans (Real.exp_le_exp.mpr (by linarith))
  have h2 : |g (y - 2 * L * (r + 1))| ≤ Real.exp (-(R + 2 * L * r)) :=
    (hg _ (by linarith)).trans (Real.exp_le_exp.mpr (by linarith))
  calc |g (y + 2 * L * (r + 1)) + g (y - 2 * L * (r + 1))|
      ≤ |g (y + 2 * L * (r + 1))| + |g (y - 2 * L * (r + 1))| := abs_add_le _ _
    _ ≤ 2 * Real.exp (-R) * Real.exp (-(2 * L)) ^ r := by
        rw [mul_assoc 2 (Real.exp (-R)), hexp]
        linarith

/-- **Wrap bound.** The images of an antiperiodic profile contribute at most
`2 exp (-R) / (1 - exp (-2L))` at every `y` with `|y| ≤ 2L - R`. -/
theorem wrap_bound {g : ℝ → ℝ} {L R y : ℝ} (hL : 0 < L)
    (hg : ∀ p, R ≤ |p| → |g p| ≤ Real.exp (-|p|)) (hy : |y| ≤ 2 * L - R) :
    Summable (fun r : ℕ =>
        (-1 : ℝ) ^ (r + 1) * (g (y + 2 * L * (r + 1)) + g (y - 2 * L * (r + 1)))) ∧
      |∑' r : ℕ, (-1 : ℝ) ^ (r + 1) * (g (y + 2 * L * (r + 1)) + g (y - 2 * L * (r + 1)))|
        ≤ 2 * Real.exp (-R) / (1 - Real.exp (-(2 * L))) := by
  have hq1 : Real.exp (-(2 * L)) < 1 := Real.exp_lt_one_iff.mpr (by linarith)
  have hgeo : HasSum (fun r : ℕ => 2 * Real.exp (-R) * Real.exp (-(2 * L)) ^ r)
      (2 * Real.exp (-R) / (1 - Real.exp (-(2 * L)))) := by
    rw [div_eq_mul_inv]
    exact (hasSum_geometric_of_lt_one (Real.exp_nonneg _) hq1).mul_left _
  have hterm : ∀ r : ℕ,
      ‖(-1 : ℝ) ^ (r + 1) * (g (y + 2 * L * (r + 1)) + g (y - 2 * L * (r + 1)))‖
        ≤ 2 * Real.exp (-R) * Real.exp (-(2 * L)) ^ r := by
    intro r
    rw [norm_mul, norm_pow, norm_neg, norm_one, one_pow, one_mul, Real.norm_eq_abs]
    exact image_pair_le hL hg hy r
  refine ⟨Summable.of_norm_bounded hgeo.summable hterm, ?_⟩
  rw [← Real.norm_eq_abs]
  exact tsum_of_norm_bounded hgeo hterm

end ExecutionProofs
