import ExecutionProofs.Riccati
import ExecutionProofs.Feedback

/-!
# Proposition 1 for a general impact matrix, and the norm ratio of Lemma 2

`A = H⁻¹` positive definite. The substitution `Y = S Ỹ S⁻¹`, `Z = S⁻¹ Z̃ S⁻¹` with `S = A^{1/2}` reduces
the flow `Y' = ½ A Z`, `Z' = 2 λ C Y` to the case `H = 1`. Jacobi's formula identifies the prefactor of
`ψ` as `(det Y)^{-1/2}`.
-/

open Set Matrix MeasureTheory

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-! ### Jacobi's formula -/

/-- Derivative of the determinant of a matrix path. -/
lemma hasDeriv_det {M : ℝ → Matrix ι ι ℝ} {M' : Matrix ι ι ℝ} {I : Set ℝ} {s : ℝ}
    (hM : ∀ i j, HasDerivWithinAt (fun r => M r i j) (M' i j) I s) :
    HasDerivWithinAt (fun r => (M r).det) (((M s).adjugate * M').trace) I s := by
  simp only [Matrix.det_apply']
  refine (HasDerivWithinAt.fun_sum fun σ _ =>
    (HasDerivWithinAt.fun_finsetProd fun k _ => hM (σ k) k).const_mul _).congr_deriv ?_
  have h1 : ((M s).adjugate * M').trace = ∑ k, ((M s).updateCol k fun i => M' i k).det := by
    simp only [Matrix.trace, Matrix.diag_apply]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [← Matrix.cramer_apply, Matrix.cramer_eq_adjugate_mulVec]
    rfl
  rw [h1]
  simp only [Matrix.det_apply', Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun σ _ => ?_
  rw [← Finset.mul_prod_erase Finset.univ _ (Finset.mem_univ k)]
  simp only [Matrix.updateCol_self, smul_eq_mul]
  have : ∏ j ∈ Finset.univ.erase k, (M s).updateCol k (fun i => M' i k) (σ j) j
      = ∏ j ∈ Finset.univ.erase k, M s (σ j) j :=
    Finset.prod_congr rfl fun j hj => by rw [Matrix.updateCol_ne (Finset.ne_of_mem_erase hj)]
  rw [this]
  ring

/-- If `Y' = ½ A Z` and `Y` stays invertible from `Y 0 = 1`, then `det Y > 0` and
`(½ log det Y)' = ¼ tr (A Z Y⁻¹)`. -/
lemma logDet_deriv {A : Matrix ι ι ℝ} {Y Z : ℝ → Matrix ι ι ℝ} {T : ℝ} (hY0 : Y 0 = 1)
    (hY : ∀ s ∈ Icc 0 T, ∀ i j, HasDerivWithinAt (fun s => Y s i j) (1 / 2 * (A * Z s) i j) (Icc 0 T) s)
    (hdet : ∀ s ∈ Icc 0 T, (Y s).det ≠ 0) {s : ℝ} (hs : s ∈ Icc 0 T) :
    0 < (Y s).det ∧ HasDerivWithinAt (fun s => 1 / 2 * Real.log (Y s).det)
      (1 / 4 * (A * (Z s * (Y s)⁻¹)).trace) (Icc 0 T) s := by
  have hD : ∀ r ∈ Icc 0 T, HasDerivWithinAt (fun r => (Y r).det)
      (((Y r).adjugate * ((1 / 2 : ℝ) • (A * Z r))).trace) (Icc 0 T) r := fun r hr =>
    hasDeriv_det fun i j => (hY r hr i j).congr_deriv (by simp [Matrix.smul_apply])
  have hpos : 0 < (Y s).det := by
    by_contra hneg
    have hcont : ContinuousOn (fun r => (Y r).det) (Icc 0 s) := fun r hr =>
      ((hD r ⟨hr.1, hr.2.trans hs.2⟩).continuousWithinAt).mono (Icc_subset_Icc_right hs.2)
    have hmem : (0 : ℝ) ∈ Icc (Y s).det (Y 0).det := ⟨not_lt.mp hneg, by simp [hY0]⟩
    obtain ⟨r, hr, h0⟩ := intermediate_value_Icc' hs.1 hcont hmem
    exact hdet r ⟨hr.1, hr.2.trans hs.2⟩ h0
  refine ⟨hpos, (((hD s hs).log hpos.ne').const_mul (1 / 2)).congr_deriv ?_⟩
  have hadj : (Y s).adjugate = (Y s).det • (Y s)⁻¹ := by
    rw [Matrix.inv_def, Ring.inverse_eq_inv', smul_smul, mul_inv_cancel₀ hpos.ne', one_smul]
  have htr : ((Y s)⁻¹ * (A * Z s)).trace = (A * (Z s * (Y s)⁻¹)).trace := by
    rw [Matrix.trace_mul_comm, Matrix.mul_assoc]
  rw [hadj, Matrix.mul_smul, Matrix.smul_mul, Matrix.trace_smul, Matrix.trace_smul, htr,
    smul_eq_mul, smul_eq_mul]
  field_simp
  ring

/-- For `H = 1` the phase `θ` of the grid theorem is `½ log det Y`. -/
theorem theta_eq {C Φ : Matrix ι ι ℝ} {lam T : ℝ} {Y Z : ℝ → Matrix ι ι ℝ}
    (h : RiccatiFlow C Φ lam T Y Z) (hdet : ∀ s ∈ Icc 0 T, (Y s).det ≠ 0) {θ : ℝ → ℝ} (h0 : θ 0 = 0)
    (hθ : ∀ s ∈ Icc 0 T, HasDerivWithinAt θ (1 / 4 * (RiccatiFlow.P Y Z s).trace) (Icc 0 T) s)
    {s : ℝ} (hs : s ∈ Icc 0 T) : θ s = 1 / 2 * Real.log (Y s).det := by
  have hd : ∀ r ∈ Icc 0 T, HasDerivWithinAt (fun r => θ r - 1 / 2 * Real.log (Y r).det) 0
      (Icc 0 T) r := by
    intro r hr
    have h1 := (logDet_deriv (A := 1) h.Y_zero (fun r hr i j => by
      simpa using h.Y_deriv r hr i j) hdet hr).2
    refine ((hθ r hr).fun_sub h1).congr_deriv ?_
    rw [Matrix.one_mul, RiccatiFlow.P, sub_self]
  have := constant_of_derivWithin_zero (f := fun r => θ r - 1 / 2 * Real.log (Y r).det)
    (a := 0) (b := T) (fun r hr => (hd r hr).differentiableWithinAt)
    (fun r hr => (hd r ⟨hr.1, hr.2.le⟩).derivWithin
      (uniqueDiffOn_Icc (hr.1.trans_lt hr.2) r ⟨hr.1, hr.2.le⟩)) s hs
  simp only [h0, h.Y_zero, Matrix.det_one, Real.log_one] at this
  linarith

/-! ### General impact -/

lemma hasDeriv_sandwich {M : ℝ → Matrix ι ι ℝ} {M' : Matrix ι ι ℝ} {I : Set ℝ} {s : ℝ}
    (hM : ∀ i j, HasDerivWithinAt (fun r => M r i j) (M' i j) I s) (L R : Matrix ι ι ℝ) (i j : ι) :
    HasDerivWithinAt (fun r => (L * M r * R) i j) ((L * M' * R) i j) I s := by
  simp only [Matrix.mul_apply]
  exact HasDerivWithinAt.fun_sum fun k _ =>
    (HasDerivWithinAt.fun_sum fun l _ => (hM l k).const_mul _).mul_const _

/-- **Riccati flow for general impact.** -/
theorem exists_general_flow {A C Φ : Matrix ι ι ℝ} {lam T : ℝ} (hT : 0 < T) (hlam : 0 ≤ lam)
    (hA : A.PosDef) (hCs : C.IsSymm) (hΦs : Φ.IsSymm) (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v))
    (hΦ : ∀ v, 0 ≤ v ⬝ᵥ (Φ *ᵥ v)) :
    ∃ Y Z : ℝ → Matrix ι ι ℝ, Y 0 = 1 ∧ Z 0 = Φ ∧
      (∀ s ∈ Icc 0 T, ∀ i j,
        HasDerivWithinAt (fun s => Y s i j) (1 / 2 * (A * Z s) i j) (Icc 0 T) s) ∧
      (∀ s ∈ Icc 0 T, ∀ i j,
        HasDerivWithinAt (fun s => Z s i j) (2 * lam * (C * Y s) i j) (Icc 0 T) s) ∧
      ∀ s ∈ Icc 0 T, (Y s).det ≠ 0 ∧ (Z s * (Y s)⁻¹).IsSymm ∧ ∀ i j,
        HasDerivWithinAt (fun s => (Z s * (Y s)⁻¹) i j)
          (2 * lam * C i j - 1 / 2 * (Z s * (Y s)⁻¹ * A * (Z s * (Y s)⁻¹)) i j) (Icc 0 T) s := by
  obtain ⟨S, hSs, hSS⟩ : ∃ S : Matrix ι ι ℝ, S.IsSymm ∧ S * S = A := by
    open scoped MatrixOrder in
    refine ⟨CFC.sqrt A, ?_, CFC.sqrt_mul_sqrt_self A (Matrix.nonneg_iff_posSemidef.mpr hA.posSemidef)⟩
    have := (Matrix.nonneg_iff_posSemidef.mp (CFC.sqrt_nonneg A)).isHermitian
    simpa [Matrix.IsSymm, Matrix.IsHermitian] using this
  have hSdet : IsUnit S.det := by
    refine isUnit_iff_ne_zero.mpr fun h0 => hA.det_pos.ne' ?_
    rw [← hSS, Matrix.det_mul, h0, zero_mul]
  have hSi : S * S⁻¹ = 1 := Matrix.mul_nonsing_inv _ hSdet
  have hiS : S⁻¹ * S = 1 := Matrix.nonsing_inv_mul _ hSdet
  have hform : ∀ M : Matrix ι ι ℝ, (∀ v, 0 ≤ v ⬝ᵥ (M *ᵥ v)) → ∀ v, 0 ≤ v ⬝ᵥ ((S * M * S) *ᵥ v) := by
    intro M hM v
    rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
      ← Matrix.mulVec_transpose, hSs]
    exact hM _
  have hsymm : ∀ M : Matrix ι ι ℝ, M.IsSymm → (S * M * S).IsSymm := by
    intro M hM
    rw [Matrix.IsSymm, Matrix.transpose_mul, Matrix.transpose_mul, hSs, hM, Matrix.mul_assoc]
  have hSis : S⁻¹.IsSymm := by rw [Matrix.IsSymm, Matrix.transpose_nonsing_inv, hSs]
  obtain ⟨Y, Z, h⟩ := exists_riccatiFlow (S * C * S) (S * Φ * S) lam hT.le
  have hdet : ∀ s ∈ Icc 0 T, (Y s).det ≠ 0 := fun s hs =>
    h.det_ne_zero_of_nonneg hlam (hform C hC) (hform Φ hΦ) hs
  have hP : ∀ s, S⁻¹ * Z s * S⁻¹ * (S * Y s * S⁻¹)⁻¹ = S⁻¹ * RiccatiFlow.P Y Z s * S⁻¹ := by
    intro s
    rw [Matrix.mul_inv_rev, Matrix.mul_inv_rev, Matrix.nonsing_inv_nonsing_inv _ hSdet, RiccatiFlow.P]
    calc S⁻¹ * Z s * S⁻¹ * (S * ((Y s)⁻¹ * S⁻¹))
        = S⁻¹ * Z s * (S⁻¹ * S) * ((Y s)⁻¹ * S⁻¹) := by simp only [Matrix.mul_assoc]
      _ = S⁻¹ * (Z s * (Y s)⁻¹) * S⁻¹ := by rw [hiS, Matrix.mul_one]; simp only [Matrix.mul_assoc]
  refine ⟨fun s => S * Y s * S⁻¹, fun s => S⁻¹ * Z s * S⁻¹, ?_, ?_, fun s hs i j => ?_, fun s hs i j => ?_,
    fun s hs => ⟨?_, ?_, fun i j => ?_⟩⟩
  · simp only [h.Y_zero, Matrix.mul_one, hSi]
  · simp only [h.Z_zero]
    calc S⁻¹ * (S * Φ * S) * S⁻¹ = (S⁻¹ * S) * Φ * (S * S⁻¹) := by simp only [Matrix.mul_assoc]
      _ = Φ := by rw [hiS, hSi, Matrix.one_mul, Matrix.mul_one]
  · refine (hasDeriv_sandwich (M' := (1 / 2 : ℝ) • Z s) (fun i j =>
      (h.Y_deriv s hs i j).congr_deriv (by simp [Matrix.smul_apply])) S S⁻¹ i j).congr_deriv ?_
    have : A * (S⁻¹ * Z s * S⁻¹) = S * Z s * S⁻¹ := by
      calc A * (S⁻¹ * Z s * S⁻¹) = S * (S * S⁻¹) * Z s * S⁻¹ := by rw [← hSS]; simp only [Matrix.mul_assoc]
        _ = S * Z s * S⁻¹ := by rw [hSi, Matrix.mul_one]
    rw [this, Matrix.mul_smul, Matrix.smul_mul, Matrix.smul_apply, smul_eq_mul]
  · refine (hasDeriv_sandwich (M' := (2 * lam) • (S * C * S * Y s)) (fun i j =>
      (h.Z_deriv s hs i j).congr_deriv (by simp [Matrix.smul_apply])) S⁻¹ S⁻¹ i j).congr_deriv ?_
    have : S⁻¹ * (S * C * S * Y s) * S⁻¹ = C * (S * Y s * S⁻¹) := by
      calc S⁻¹ * (S * C * S * Y s) * S⁻¹ = (S⁻¹ * S) * C * (S * Y s * S⁻¹) := by
            simp only [Matrix.mul_assoc]
        _ = C * (S * Y s * S⁻¹) := by rw [hiS, Matrix.one_mul]
    rw [Matrix.mul_smul, Matrix.smul_mul, Matrix.smul_apply, smul_eq_mul, this]
  · rw [Matrix.det_mul, Matrix.det_mul]
    exact mul_ne_zero (mul_ne_zero hSdet.ne_zero (hdet s hs))
      (by rw [Matrix.det_nonsing_inv, Ring.inverse_eq_inv']; exact inv_ne_zero hSdet.ne_zero)
  · rw [hP, Matrix.IsSymm, Matrix.transpose_mul, Matrix.transpose_mul, hSis,
      h.P_symm_of_det (hsymm C hCs) (hsymm Φ hΦs) hs (hdet s hs), Matrix.mul_assoc]
  · simp only [hP]
    refine (hasDeriv_sandwich (M' := (2 * lam) • (S * C * S)
      - (1 / 2 : ℝ) • (RiccatiFlow.P Y Z s * RiccatiFlow.P Y Z s)) (fun i j =>
      (h.P_deriv_of_det hT hdet hs i j).congr_deriv
        (by simp [Matrix.smul_apply, Matrix.sub_apply])) S⁻¹ S⁻¹ i j).congr_deriv ?_
    generalize RiccatiFlow.P Y Z s = P
    have e1 : S⁻¹ * (S * C * S) * S⁻¹ = C := by
      calc S⁻¹ * (S * C * S) * S⁻¹ = (S⁻¹ * S) * C * (S * S⁻¹) := by simp only [Matrix.mul_assoc]
        _ = C := by rw [hiS, hSi, Matrix.one_mul, Matrix.mul_one]
    have e2 : S⁻¹ * P * S⁻¹ * A * (S⁻¹ * P * S⁻¹) = S⁻¹ * (P * P) * S⁻¹ := by
      calc S⁻¹ * P * S⁻¹ * A * (S⁻¹ * P * S⁻¹)
          = S⁻¹ * P * (S⁻¹ * S) * (S * S⁻¹) * P * S⁻¹ := by rw [← hSS]; simp only [Matrix.mul_assoc]
        _ = S⁻¹ * (P * P) * S⁻¹ := by
            rw [hiS, hSi, Matrix.mul_one, Matrix.mul_one]; simp only [Matrix.mul_assoc]
    rw [e2, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_smul,
      Matrix.smul_mul, e1]
    simp [Matrix.smul_apply, Matrix.sub_apply]

/-- **Proposition 1.** For every positive definite `A = H⁻¹`, positive semidefinite `C`, `Φ`, and every
`ν > 0`, the function `ψ = (det Y)^{-1/2} exp (-(x · P x) / 8ν)` solves the linear equation with data
`exp (-(x · Φ x) / 8ν)`, returns the feedback `½ A P x`, and `½ x · P x` solves the Hamilton–Jacobi
equation. -/
theorem proposition_one {A C Φ : Matrix ι ι ℝ} {lam T : ℝ} (hT : 0 < T) (hlam : 0 ≤ lam)
    (hA : A.PosDef) (hCs : C.IsSymm) (hΦs : Φ.IsSymm) (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v))
    (hΦ : ∀ v, 0 ≤ v ⬝ᵥ (Φ *ᵥ v)) :
    ∃ Y P : ℝ → Matrix ι ι ℝ, P 0 = Φ ∧
      (∀ s ∈ Icc 0 T, 0 < (Y s).det ∧ (P s).IsSymm ∧
        ∀ x, HasDerivWithinAt (fun s => 1 / 2 * (x ⬝ᵥ (P s *ᵥ x)))
          (lam * (x ⬝ᵥ (C *ᵥ x)) - 1 / 4 * ((P s *ᵥ x) ⬝ᵥ (A *ᵥ (P s *ᵥ x)))) (Icc 0 T) s) ∧
      ∀ ν, 0 < ν → ∃ ψ : (ι → ℝ) → ℝ → ℝ,
        (∀ x, ∀ s ∈ Icc 0 T, ψ x s = (Y s).det ^ (-(1 / 2) : ℝ)
          * Real.exp (-(x ⬝ᵥ (P s *ᵥ x)) / (8 * ν))) ∧
        (∀ x, ψ x 0 = Real.exp (-(x ⬝ᵥ (Φ *ᵥ x)) / (8 * ν))) ∧
        (∀ x, ∀ s ∈ Icc 0 T, HasDerivWithinAt (ψ x)
          (ν * ∑ a, ∑ b, A a b * dd a (dd b fun x => ψ x s) x
            - lam / (4 * ν) * (x ⬝ᵥ (C *ᵥ x)) * ψ x s) (Icc 0 T) s) ∧
        ∀ x, ∀ s ∈ Icc 0 T, ∀ a, -2 * ν * ∑ b, A a b * dd b (fun x => Real.log (ψ x s)) x
          = 1 / 2 * (A *ᵥ (P s *ᵥ x)) a := by
  obtain ⟨Y, Z, hY0, hZ0, hY, -, hmain⟩ := exists_general_flow hT hlam hA hCs hΦs hC hΦ
  have hlog : ∀ s ∈ Icc 0 T, _ := fun s hs =>
    logDet_deriv (A := A) hY0 hY (fun s hs => (hmain s hs).1) hs
  have hP0 : Z 0 * (Y 0)⁻¹ = Φ := by rw [hY0, hZ0, inv_one, Matrix.mul_one]
  refine ⟨Y, fun s => Z s * (Y s)⁻¹, hP0,
    fun s hs => ⟨(hlog s hs).1, (hmain s hs).2.1, fun x => hjb (hmain s hs).2.1 (hmain s hs).2.2 x⟩,
    fun ν hν => ⟨fun x s => psiG ν (Z s * (Y s)⁻¹) (1 / 2 * Real.log (Y s).det) x, ?_, ?_, ?_, ?_⟩⟩
  · intro x s hs
    change psiG ν (Z s * (Y s)⁻¹) (1 / 2 * Real.log (Y s).det) x = _
    rw [psiG, sub_eq_add_neg, Real.exp_add, mul_comm, Real.rpow_def_of_pos (hlog s hs).1]
    congr 2
    ring
  · intro x
    change psiG ν (Z 0 * (Y 0)⁻¹) (1 / 2 * Real.log (Y 0).det) x = _
    rw [psiG, hP0, hY0, Matrix.det_one, Real.log_one, mul_zero, sub_zero]
  · intro x s hs
    exact psiG_solves hν.ne' (hmain s hs).2.1 (hmain s hs).2.2 (hlog s hs).2 x
  · intro x s hs a
    exact feedback hν.ne' (hmain s hs).2.1 x a

/-! ### The norm ratio -/

/-- Gaussian integral of `ψ²`. -/
theorem integral_psiG_sq {ν : ℝ} (hν : 0 < ν) (θ : ℝ) {P : Matrix ι ι ℝ} (hP : P.PosDef) :
    ∫ x, psiG ν P θ x ^ 2
      = Real.exp (-(2 * θ)) * √(4 * Real.pi * ν) ^ Fintype.card ι / √P.det := by
  obtain ⟨S, hSs, hSS, hS0⟩ : ∃ S : Matrix ι ι ℝ, S.IsSymm ∧ S * S = P ∧ 0 ≤ S.det := by
    open scoped MatrixOrder in
    have hnn := Matrix.nonneg_iff_posSemidef.mp (CFC.sqrt_nonneg P)
    refine ⟨CFC.sqrt P, ?_, CFC.sqrt_mul_sqrt_self P (Matrix.nonneg_iff_posSemidef.mpr hP.posSemidef),
      hnn.det_nonneg⟩
    simpa [Matrix.IsSymm, Matrix.IsHermitian] using hnn.isHermitian
  have hdet : S.det * S.det = P.det := by rw [← Matrix.det_mul, hSS]
  have hSne : S.det ≠ 0 := fun h0 => hP.det_pos.ne' (by rw [← hdet, h0, zero_mul])
  have hSsqrt : √P.det = S.det := by rw [← hdet, Real.sqrt_mul_self hS0]
  set f : (ι → ℝ) → ℝ := fun y => ∏ i, Real.exp (-(1 / (4 * ν)) * y i ^ 2) with hf
  have hfm : AEStronglyMeasurable f (Measure.map (toLin' S) volume) := by
    refine Continuous.aestronglyMeasurable ?_
    fun_prop
  have hval : ∀ x, psiG ν P θ x ^ 2 = Real.exp (-(2 * θ)) * f (toLin' S x) := by
    intro x
    have hq : x ⬝ᵥ (P *ᵥ x) = ∑ i, (S *ᵥ x) i ^ 2 := by
      rw [← hSS, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hSs]
      simp only [dotProduct, sq]
    rw [psiG, ← Real.exp_nat_mul, hf]
    simp only [toLin'_apply]
    rw [← Real.exp_sum, ← Real.exp_add, ← Finset.mul_sum, ← hq]
    congr 1
    push_cast
    field_simp
    ring
  have hmeas : Measurable (toLin' S) := (LinearMap.continuous_of_finiteDimensional _).measurable
  simp_rw [hval]
  rw [integral_const_mul, ← integral_map hmeas.aemeasurable hfm,
    Real.map_matrix_volume_pi_eq_smul_volume_pi hSne, integral_smul_measure, hf,
    integral_fintype_prod_volume_eq_pow (fun t : ℝ => Real.exp (-(1 / (4 * ν)) * t ^ 2)),
    integral_gaussian, hSsqrt, ENNReal.toReal_ofReal (abs_nonneg _), abs_inv, abs_of_nonneg hS0,
    smul_eq_mul]
  have hb : Real.pi / (1 / (4 * ν)) = 4 * Real.pi * ν := by field_simp
  rw [hb]
  ring

/-- **Lemma 2, norm ratio.** It does not depend on `ν`. -/
theorem norm_ratio {ν : ℝ} (hν : 0 < ν) {P Φ : Matrix ι ι ℝ} (hP : P.PosDef) (hΦ : Φ.PosDef) {D : ℝ}
    (hD : 0 < D) :
    (∫ x, psiG ν P (1 / 2 * Real.log D) x ^ 2) / ∫ x, psiG ν Φ 0 x ^ 2
      = 1 / D * √(Φ.det / P.det) := by
  rw [integral_psiG_sq hν _ hP, integral_psiG_sq hν _ hΦ, Real.sqrt_div' _ hP.det_pos.le]
  have h1 : Real.exp (-(2 * (1 / 2 * Real.log D))) = 1 / D := by
    rw [show -(2 * (1 / 2 * Real.log D)) = -Real.log D by ring, Real.exp_neg, Real.exp_log hD, one_div]
  have h2 : (0 : ℝ) < √(4 * Real.pi * ν) ^ Fintype.card ι := by
    have := Real.pi_pos
    positivity
  have h3 := Real.sqrt_pos.mpr hP.det_pos
  have h4 := Real.sqrt_pos.mpr hΦ.det_pos
  rw [h1, mul_zero, neg_zero, Real.exp_zero]
  field_simp

end ExecutionProofs
