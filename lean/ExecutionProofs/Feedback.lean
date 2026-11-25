import Mathlib

/-!
# Exact feedback and scale invariance (paper, Proposition 1 and Lemma 2)

General impact: `A = H⁻¹` symmetric. `ψ = exp (-(x · P x) / (8 ν) - θ)` with the Riccati equation
`P' = 2 λ C - ½ P A P` and `θ' = tr (A P) / 4` solves the linear equation, and its logarithmic
gradient returns the feedback.
-/

open Set Matrix

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The Gaussian of Proposition 1. -/
noncomputable def psiG (ν : ℝ) (P : Matrix ι ι ℝ) (θ : ℝ) (x : ι → ℝ) : ℝ :=
  Real.exp (-(x ⬝ᵥ (P *ᵥ x)) / (8 * ν) - θ)

/-- First partial derivative along axis `b`, as a genuine derivative. -/
noncomputable def dd (b : ι) (f : (ι → ℝ) → ℝ) (x : ι → ℝ) : ℝ :=
  deriv (fun t : ℝ => f (x + Pi.single b t)) 0

lemma mulVec_add_single (P : Matrix ι ι ℝ) (x : ι → ℝ) (a b : ι) (t : ℝ) :
    (P *ᵥ (x + Pi.single b t)) a = (P *ᵥ x) a + P a b * t := by
  rw [Matrix.mulVec_add, Pi.add_apply]
  congr 1
  simp [Matrix.mulVec, dotProduct, Pi.single_apply]

lemma hasDerivAt_psiG {ν θ : ℝ} {P : Matrix ι ι ℝ} (hP : P.IsSymm) (x : ι → ℝ) (b : ι) :
    HasDerivAt (fun t : ℝ => psiG ν P θ (x + Pi.single b t))
      (-((P *ᵥ x) b) / (4 * ν) * psiG ν P θ x) 0 := by
  have hquad : ∀ t : ℝ, (x + Pi.single b t) ⬝ᵥ (P *ᵥ (x + Pi.single b t))
      = x ⬝ᵥ (P *ᵥ x) + 2 * t * (P *ᵥ x) b + t ^ 2 * P b b := by
    intro t
    have h1 : P *ᵥ Pi.single b t = fun i => P i b * t := by
      ext i
      simp [Matrix.mulVec, dotProduct, Pi.single_apply]
    have h2 : x ⬝ᵥ (fun i => P i b * t) = t * (P *ᵥ x) b := by
      simp only [dotProduct, Matrix.mulVec, Finset.mul_sum]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [← hP.apply b i]
      ring
    have h3 : ∀ y : ι → ℝ, (Pi.single b t : ι → ℝ) ⬝ᵥ y = t * y b := fun y => by
      simp [dotProduct, Pi.single_apply]
    rw [Matrix.mulVec_add, add_dotProduct, dotProduct_add, dotProduct_add, h1, h2, h3, h3]
    ring
  unfold psiG
  simp only [hquad]
  have hq : HasDerivAt (fun t : ℝ => x ⬝ᵥ (P *ᵥ x) + 2 * t * (P *ᵥ x) b + t ^ 2 * P b b)
      (2 * (P *ᵥ x) b) 0 := by
    have h1 : HasDerivAt (fun t : ℝ => 2 * t * (P *ᵥ x) b) (2 * (P *ᵥ x) b) 0 := by
      simpa using ((hasDerivAt_id (0 : ℝ)).const_mul 2).mul_const ((P *ᵥ x) b)
    have h2 : HasDerivAt (fun t : ℝ => t ^ 2 * P b b) 0 0 := by
      simpa using (hasDerivAt_pow 2 (0 : ℝ)).mul_const (P b b)
    simpa using (h1.const_add (x ⬝ᵥ (P *ᵥ x))).fun_add h2
  refine ((hq.fun_neg.div_const (8 * ν)).sub_const θ).exp.congr_deriv ?_
  simp only [mul_zero, zero_mul, add_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow]
  ring

lemma dd_psiG {ν θ : ℝ} {P : Matrix ι ι ℝ} (hP : P.IsSymm) (b : ι) :
    dd b (psiG ν P θ) = fun x => -((P *ᵥ x) b) / (4 * ν) * psiG ν P θ x := by
  funext x
  exact (hasDerivAt_psiG hP x b).deriv

lemma dd_dd_psiG {ν θ : ℝ} {P : Matrix ι ι ℝ} (hP : P.IsSymm) (a b : ι) (x : ι → ℝ) :
    dd a (dd b (psiG ν P θ)) x
      = ((P *ᵥ x) a * (P *ᵥ x) b / (4 * ν) ^ 2 - P b a / (4 * ν)) * psiG ν P θ x := by
  rw [dd_psiG hP]
  unfold dd
  have h1 : HasDerivAt (fun t : ℝ => -((P *ᵥ (x + Pi.single a t)) b) / (4 * ν))
      (-(P b a) / (4 * ν)) 0 := by
    simp only [mulVec_add_single]
    have := ((((hasDerivAt_id (0 : ℝ)).const_mul (P b a)).const_add ((P *ᵥ x) b)).fun_neg).div_const
      (4 * ν)
    simpa using this
  have h2 := (h1.fun_mul (hasDerivAt_psiG (ν := ν) (θ := θ) hP x a)).deriv
  rw [h2]
  simp only [Pi.single_zero, add_zero]
  ring

/-- **Proposition 1, the equation.** `ψ` solves `∂ₛ ψ = ν ∑ A_ab ∂ₐ ∂_b ψ - (λ / 4ν) (x · C x) ψ`. -/
theorem psiG_solves {ν lam s : ℝ} {I : Set ℝ} (hν : ν ≠ 0) {A C : Matrix ι ι ℝ} {P : ℝ → Matrix ι ι ℝ}
    {θ : ℝ → ℝ} (hPs : (P s).IsSymm)
    (hP : ∀ i j, HasDerivWithinAt (fun s => P s i j)
      (2 * lam * C i j - 1 / 2 * (P s * A * P s) i j) I s)
    (hθ : HasDerivWithinAt θ (1 / 4 * (A * P s).trace) I s) (x : ι → ℝ) :
    HasDerivWithinAt (fun s => psiG ν (P s) (θ s) x)
      (ν * ∑ a, ∑ b, A a b * dd a (dd b (psiG ν (P s) (θ s))) x
        - lam / (4 * ν) * (x ⬝ᵥ (C *ᵥ x)) * psiG ν (P s) (θ s) x) I s := by
  have hform : HasDerivWithinAt (fun s => x ⬝ᵥ (P s *ᵥ x))
      (2 * lam * (x ⬝ᵥ (C *ᵥ x)) - 1 / 2 * (x ⬝ᵥ ((P s * A * P s) *ᵥ x))) I s := by
    simp only [dotProduct, Matrix.mulVec]
    refine (HasDerivWithinAt.fun_sum fun i _ => (HasDerivWithinAt.fun_sum fun j _ =>
      (hP i j).mul_const (x j)).const_mul (x i)).congr_deriv ?_
    simp only [Finset.mul_sum, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring
  have hPAP : x ⬝ᵥ ((P s * A * P s) *ᵥ x) = ∑ a, ∑ b, A a b * ((P s *ᵥ x) a * (P s *ᵥ x) b) := by
    have h1 : x ⬝ᵥ ((P s * A * P s) *ᵥ x) = (P s *ᵥ x) ⬝ᵥ (A *ᵥ (P s *ᵥ x)) := by
      rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
        ← Matrix.mulVec_transpose, hPs]
    rw [h1]
    simp only [dotProduct, Matrix.mulVec, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    ring
  have htr : (A * P s).trace = ∑ a, ∑ b, A a b * P s b a := by
    simp [Matrix.trace, Matrix.mul_apply]
  have hd : HasDerivWithinAt (fun s => psiG ν (P s) (θ s) x)
      ((-(2 * lam * (x ⬝ᵥ (C *ᵥ x)) - 1 / 2 * (x ⬝ᵥ ((P s * A * P s) *ᵥ x))) / (8 * ν)
        - 1 / 4 * (A * P s).trace) * psiG ν (P s) (θ s) x) I s := by
    unfold psiG
    refine (((hform.fun_neg.div_const (8 * ν)).fun_sub hθ).exp).congr_deriv ?_
    ring
  refine hd.congr_deriv ?_
  have hsum : ν * ∑ a, ∑ b, A a b * dd a (dd b (psiG ν (P s) (θ s))) x
      = (1 / (16 * ν) * (∑ a, ∑ b, A a b * ((P s *ᵥ x) a * (P s *ᵥ x) b))
        - 1 / 4 * ∑ a, ∑ b, A a b * P s b a) * psiG ν (P s) (θ s) x := by
    simp only [dd_dd_psiG hPs, Finset.mul_sum, ← Finset.sum_sub_distrib, Finset.sum_mul]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    field_simp
    ring
  rw [hsum, hPAP, htr]
  field_simp
  ring

omit [DecidableEq ι] in
lemma log_psiG (ν θ : ℝ) (P : Matrix ι ι ℝ) (x : ι → ℝ) :
    Real.log (psiG ν P θ x) = -(x ⬝ᵥ (P *ᵥ x)) / (8 * ν) - θ := Real.log_exp _

/-- **Proposition 1, the feedback.** `-2 ν A ∇ log ψ = ½ A P x` for every `ν > 0`. -/
theorem feedback {ν θ : ℝ} (hν : ν ≠ 0) {A P : Matrix ι ι ℝ} (hP : P.IsSymm) (x : ι → ℝ) (a : ι) :
    -2 * ν * ∑ b, A a b * dd b (fun x => Real.log (psiG ν P θ x)) x
      = 1 / 2 * (A *ᵥ (P *ᵥ x)) a := by
  have h : ∀ b, dd b (fun x => Real.log (psiG ν P θ x)) x = -((P *ᵥ x) b) / (4 * ν) := by
    intro b
    unfold dd
    have h1 := (hasDerivAt_psiG (ν := ν) (θ := θ) hP x b).log
      (by simp [psiG, Real.exp_ne_zero])
    rw [h1.deriv]
    simp only [Pi.single_zero, add_zero]
    have : psiG ν P θ x ≠ 0 := Real.exp_ne_zero _
    field_simp
  simp only [h, Matrix.mulVec, dotProduct, Finset.mul_sum]
  refine Finset.sum_congr rfl fun b _ => ?_
  field_simp
  simp only [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  ring

omit [DecidableEq ι] in
/-- The quadratic value solves the original Hamilton–Jacobi equation. -/
theorem hjb {lam s : ℝ} {I : Set ℝ} {A C : Matrix ι ι ℝ} {P : ℝ → Matrix ι ι ℝ} (hPs : (P s).IsSymm)
    (hP : ∀ i j, HasDerivWithinAt (fun s => P s i j)
      (2 * lam * C i j - 1 / 2 * (P s * A * P s) i j) I s) (x : ι → ℝ) :
    HasDerivWithinAt (fun s => 1 / 2 * (x ⬝ᵥ (P s *ᵥ x)))
      (lam * (x ⬝ᵥ (C *ᵥ x)) - 1 / 4 * ((P s *ᵥ x) ⬝ᵥ (A *ᵥ (P s *ᵥ x)))) I s := by
  have hform : HasDerivWithinAt (fun s => x ⬝ᵥ (P s *ᵥ x))
      (2 * lam * (x ⬝ᵥ (C *ᵥ x)) - 1 / 2 * (x ⬝ᵥ ((P s * A * P s) *ᵥ x))) I s := by
    simp only [dotProduct, Matrix.mulVec]
    refine (HasDerivWithinAt.fun_sum fun i _ => (HasDerivWithinAt.fun_sum fun j _ =>
      (hP i j).mul_const (x j)).const_mul (x i)).congr_deriv ?_
    simp only [Finset.mul_sum, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring
  refine (hform.const_mul (1 / 2)).congr_deriv ?_
  have h1 : x ⬝ᵥ ((P s * A * P s) *ᵥ x) = (P s *ᵥ x) ⬝ᵥ (A *ᵥ (P s *ᵥ x)) := by
    rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
      ← Matrix.mulVec_transpose, hPs]
  rw [h1]
  ring

omit [DecidableEq ι] in
/-- **Lemma 2 (scale invariance).** `ψ_ν(x, s) = ψ_1(x / √ν, s)`. -/
theorem scale_invariance {ν : ℝ} (hν : 0 < ν) (P : Matrix ι ι ℝ) (θ : ℝ) (x : ι → ℝ) :
    psiG ν P θ x = psiG 1 P θ ((Real.sqrt ν)⁻¹ • x) := by
  unfold psiG
  congr 1
  rw [Matrix.mulVec_smul, dotProduct_smul, smul_dotProduct, smul_eq_mul, smul_eq_mul]
  have h : (Real.sqrt ν)⁻¹ * (Real.sqrt ν)⁻¹ = ν⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt hν.le]
  rw [← mul_assoc, h]
  field_simp

end ExecutionProofs
