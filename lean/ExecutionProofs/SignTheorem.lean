import ExecutionProofs.Sign

/-!
# Strict signs, homogeneity, uniqueness and the radial sign change (paper, Theorem 3 (iv))
-/

open Set Matrix MvPolynomial

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Derivative of the inverse of a differentiable matrix family. -/
lemma hasDeriv_inv {M : ℝ → Matrix ι ι ℝ} {M' : Matrix ι ι ℝ} {T s : ℝ} (hT : 0 < T)
    (hs : s ∈ Icc 0 T) (hdet : ∀ r ∈ Icc 0 T, (M r).det ≠ 0)
    (hM : ∀ i j, HasDerivWithinAt (fun r => M r i j) (M' i j) (Icc 0 T) s) (i j : ι) :
    HasDerivWithinAt (fun r => (M r)⁻¹ i j) ((-((M s)⁻¹ * M' * (M s)⁻¹)) i j) (Icc 0 T) s := by
  have huniq : UniqueDiffWithinAt ℝ (Icc 0 T) s := uniqueDiffOn_Icc hT s hs
  have hNd : ∀ a b, DifferentiableWithinAt ℝ (fun r => (M r)⁻¹ a b) (Icc 0 T) s :=
    differentiableWithinAt_inv (fun a b => (hM a b).differentiableWithinAt) (hdet s hs)
  set N' : Matrix ι ι ℝ := Matrix.of fun a b =>
    derivWithin (fun r => (M r)⁻¹ a b) (Icc 0 T) s with hN'
  have hN : ∀ a b, HasDerivWithinAt (fun r => (M r)⁻¹ a b) (N' a b) (Icc 0 T) s := fun a b =>
    (hNd a b).hasDerivWithinAt
  have hid : N' * M s + (M s)⁻¹ * M' = 0 := by
    ext a b
    have hd : HasDerivWithinAt (fun r => ((M r)⁻¹ * M r) a b)
        ((N' * M s + (M s)⁻¹ * M') a b) (Icc 0 T) s := by
      simp only [Matrix.mul_apply, Matrix.add_apply]
      refine (HasDerivWithinAt.fun_sum fun k _ => (hN a k).mul (hM k b)).congr_deriv ?_
      rw [← Finset.sum_add_distrib]
    have hc : HasDerivWithinAt (fun r => ((M r)⁻¹ * M r) a b) 0 (Icc 0 T) s := by
      refine (hasDerivWithinAt_const s (Icc 0 T) ((1 : Matrix ι ι ℝ) a b)).congr ?_ ?_
      · intro r hr
        rw [Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.mpr (hdet r hr))]
      · rw [Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.mpr (hdet s hs))]
    exact huniq.eq_deriv _ hd hc
  have hunit := isUnit_iff_ne_zero.mpr (hdet s hs)
  have h1 : N' * M s * (M s)⁻¹ = N' := by
    rw [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunit, Matrix.mul_one]
  have h2 : N' * M s = -((M s)⁻¹ * M') := eq_neg_of_add_eq_zero_left hid
  have h3 : N' = -((M s)⁻¹ * M' * (M s)⁻¹) := by
    rw [← h1, h2, Matrix.neg_mul]
  rw [← h3]
  exact hN i j

omit [DecidableEq ι] in
/-- Derivative of a quadratic form with moving matrix and vector. -/
lemma hasDeriv_quad {T s : ℝ} {w : ℝ → ι → ℝ} {w' : ι → ℝ} {M : ℝ → Matrix ι ι ℝ}
    {M' : Matrix ι ι ℝ} (hw : ∀ i, HasDerivWithinAt (fun r => w r i) (w' i) (Icc 0 T) s)
    (hM : ∀ i j, HasDerivWithinAt (fun r => M r i j) (M' i j) (Icc 0 T) s) :
    HasDerivWithinAt (fun r => w r ⬝ᵥ (M r *ᵥ w r))
      (w' ⬝ᵥ (M s *ᵥ w s) + w s ⬝ᵥ (M' *ᵥ w s) + w s ⬝ᵥ (M s *ᵥ w')) (Icc 0 T) s := by
  simp only [dotProduct, Matrix.mulVec]
  refine (HasDerivWithinAt.fun_sum fun i _ => (hw i).fun_mul (HasDerivWithinAt.fun_sum fun j _ =>
    (hM i j).fun_mul (hw j))).congr_deriv ?_
  simp only [Finset.mul_sum, ← Finset.sum_add_distrib, mul_add]
  exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring

/-- Hessian of a polynomial at a point. -/
noncomputable def hess (A : MvPolynomial ι ℝ) (y : ι → ℝ) : Matrix ι ι ℝ :=
  Matrix.of fun i j => eval y (pderiv i (pderiv j A))

omit [Fintype ι] [DecidableEq ι] in
lemma hess_symm (A : MvPolynomial ι ℝ) (y : ι → ℝ) : (hess A y).IsSymm := by
  ext i j
  simp only [hess, Matrix.transpose_apply, Matrix.of_apply]
  rw [pderiv_comm]

/-- Second derivatives of the transport operator. -/
lemma eval_pderiv2_polyGen (P : Matrix ι ι ℝ) (hP : P.IsSymm) (A : MvPolynomial ι ℝ) (y : ι → ℝ)
    (i j : ι) :
    eval y (pderiv i (pderiv j (polyGen 0 P A)))
      = -(∑ a, (P *ᵥ y) a / 2 * eval y (pderiv a (pderiv i (pderiv j A))))
        - 1 / 2 * (hess A y * P) i j - 1 / 2 * (P * hess A y) i j := by
  have hterm : ∀ a b, eval y (pderiv i (pderiv j ((P a b / 2) • (X b * pderiv a A))))
      = P a b / 2 * ((if b = i then 1 else 0) * eval y (pderiv j (pderiv a A)))
        + P a b / 2 * ((if b = j then 1 else 0) * eval y (pderiv i (pderiv a A)))
        + P a b / 2 * (y b * eval y (pderiv a (pderiv i (pderiv j A)))) := by
    intro a b
    have hX : ∀ k : ι, pderiv k (X b : MvPolynomial ι ℝ) = MvPolynomial.C (if b = k then 1 else 0) := by
      intro k
      rw [pderiv_X]
      by_cases h : b = k
      · subst h
        simp
      · simp [h]
    simp only [Derivation.map_smul, pderiv_mul, map_add, hX, pderiv_C, zero_mul, zero_add,
      smul_eval, map_mul, eval_C, eval_X]
    have hcomm : pderiv i (pderiv j (pderiv a A)) = pderiv a (pderiv i (pderiv j A)) := by
      rw [pderiv_comm j a, pderiv_comm i a]
    rw [hcomm]
    ring
  have hsum : eval y (pderiv i (pderiv j (polyGen 0 P A)))
      = -∑ a, ∑ b, eval y (pderiv i (pderiv j ((P a b / 2) • (X b * pderiv a A)))) := by
    simp only [polyGen, zero_smul, zero_sub, map_neg, map_sum]
  rw [hsum]
  simp only [hterm, Finset.sum_add_distrib]
  have e1 : ∑ a, ∑ b, P a b / 2 * ((if b = i then 1 else 0) * eval y (pderiv j (pderiv a A)))
      = 1 / 2 * (P * hess A y) i j := by
    rw [Matrix.mul_apply, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Finset.sum_eq_single i]
    · rw [hess, Matrix.of_apply, ← hP.apply i a, pderiv_comm j a]
      simp
      ring
    · intro b _ hb
      simp [hb]
    · simp
  have e2 : ∑ a, ∑ b, P a b / 2 * ((if b = j then 1 else 0) * eval y (pderiv i (pderiv a A)))
      = 1 / 2 * (hess A y * P) i j := by
    rw [Matrix.mul_apply, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Finset.sum_eq_single j]
    · rw [hess, Matrix.of_apply]
      simp
      ring
    · intro b _ hb
      simp [hb]
    · simp
  have e3 : ∑ a, ∑ b, P a b / 2 * (y b * eval y (pderiv a (pderiv i (pderiv j A))))
      = ∑ a, (P *ᵥ y) a / 2 * eval y (pderiv a (pderiv i (pderiv j A))) := by
    refine Finset.sum_congr rfl fun a _ => ?_
    simp only [Matrix.mulVec, dotProduct, Finset.sum_div, Finset.sum_mul]
    exact Finset.sum_congr rfl fun b _ => by ring
  rw [e1, e2, e3]
  ring

lemma mulVec_line' (P : Matrix ι ι ℝ) (x : ι → ℝ) (a j : ι) (t : ℝ) :
    (P *ᵥ (x + Pi.single j t)) a = (P *ᵥ x) a + P a j * t := by
  rw [Matrix.mulVec_add, Pi.add_apply]
  congr 1
  simp [Matrix.mulVec, dotProduct, Pi.single_apply]

lemma hasDerivAt_lin_line (P : Matrix ι ι ℝ) (x : ι → ℝ) (a j : ι) :
    HasDerivAt (fun t : ℝ => (P *ᵥ (x + Pi.single j t)) a) (P a j) 0 := by
  simp only [mulVec_line']
  simpa using ((hasDerivAt_id (0 : ℝ)).const_mul (P a j)).const_add ((P *ᵥ x) a)

/-- Hessian of the quartic source. -/
lemma eval_pderiv2_srcA (P : Matrix ι ι ℝ) (y : ι → ℝ) (i j : ι) :
    eval y (pderiv i (pderiv j (srcA P)))
      = 1 / 256 * ∑ a, ((P *ᵥ y) a) ^ 2 * (P a i * P a j) := by
  have h1 : ∀ x, eval x (pderiv j (srcA P)) = 1 / 768 * ∑ a, ((P *ᵥ x) a) ^ 3 * P a j := by
    intro x
    refine eval_pderiv_of_eval (eval_srcA P) ?_
    refine ((HasDerivAt.fun_sum fun a _ =>
      (hasDerivAt_lin_line P x a j).fun_pow 4).const_mul (1 / 3072)).congr_deriv ?_
    simp only [Pi.single_zero, add_zero, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    push_cast
    ring
  refine eval_pderiv_of_eval h1 ?_
  refine ((HasDerivAt.fun_sum fun a _ =>
    ((hasDerivAt_lin_line P y a i).fun_pow 3).mul_const (P a j)).const_mul (1 / 768)).congr_deriv ?_
  simp only [Pi.single_zero, add_zero, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  push_cast
  ring

/-- Source matrix of the Hessian along a characteristic. -/
noncomputable def srcHess (P : Matrix ι ι ℝ) (z : ι → ℝ) : Matrix ι ι ℝ :=
  Matrix.of fun i j => 1 / 256 * ∑ a, (z a) ^ 2 * (P a i * P a j)

/-! ### The setting -/

/-- Hypotheses on the model data, bundled. -/
structure Model (C Φ : Matrix ι ι ℝ) (lam T φ : ℝ) (Y Z : ℝ → Matrix ι ι ℝ) : Prop where
  hT : 0 < T
  hlam : 0 ≤ lam
  hCs : C.IsSymm
  hΦs : Φ.IsSymm
  hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v)
  hφ : 0 < φ
  hΦ : ∀ v, φ * (v ⬝ᵥ v) ≤ v ⬝ᵥ (Φ *ᵥ v)
  flow : RiccatiFlow C Φ lam T Y Z

namespace Model

variable {C Φ : Matrix ι ι ℝ} {lam T φ : ℝ} {Y Z : ℝ → Matrix ι ι ℝ}

lemma detY (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) : (Y s).det ≠ 0 :=
  h.flow.det_ne_zero h.hlam h.hC h.hφ h.hΦ hs

lemma symm (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) :
    (RiccatiFlow.P Y Z s).IsSymm := h.flow.P_symm h.hlam h.hCs h.hΦs h.hC h.hφ h.hΦ hs

lemma deriv (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) (i j : ι) :
    HasDerivWithinAt (fun s => RiccatiFlow.P Y Z s i j)
      (2 * lam * C i j - 1 / 2 * (RiccatiFlow.P Y Z s * RiccatiFlow.P Y Z s) i j) (Icc 0 T) s :=
  h.flow.P_deriv h.hT h.hlam h.hC h.hφ h.hΦ hs i j

lemma pos_form (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) {x : ι → ℝ} (hx : x ≠ 0) :
    0 < x ⬝ᵥ (RiccatiFlow.P Y Z s *ᵥ x) := by
  obtain ⟨m, hm, hc⟩ := h.flow.exists_coercive h.hlam h.hC h.hφ h.hΦ
  have h1 : 0 < x ⬝ᵥ x := by
    rcases (Finset.sum_nonneg fun i _ => mul_self_nonneg (x i)).lt_or_eq with h' | h'
    · exact h'
    · exact absurd (dotProduct_self_eq_zero.mp h'.symm) hx
  exact (mul_pos hm h1).trans_le (hc s hs x)

lemma nonneg_form (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) (x : ι → ℝ) :
    0 ≤ x ⬝ᵥ (RiccatiFlow.P Y Z s *ᵥ x) := by
  by_cases hx : x = 0
  · simp [hx]
  · exact (h.pos_form hs hx).le

lemma detP (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) :
    (RiccatiFlow.P Y Z s).det ≠ 0 := by
  intro hdet
  obtain ⟨v, hv, hPv⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet
  have := h.pos_form hs hv
  rw [hPv, dotProduct_zero] at this
  exact lt_irrefl _ this

lemma mul_inv (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) :
    RiccatiFlow.P Y Z s * (RiccatiFlow.P Y Z s)⁻¹ = 1 :=
  Matrix.mul_nonsing_inv _ (isUnit_iff_ne_zero.mpr (h.detP hs))

lemma inv_mul (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) :
    (RiccatiFlow.P Y Z s)⁻¹ * RiccatiFlow.P Y Z s = 1 :=
  Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.mpr (h.detP hs))

lemma inv_symm (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) :
    ((RiccatiFlow.P Y Z s)⁻¹).IsSymm := by
  rw [Matrix.IsSymm, Matrix.transpose_nonsing_inv, h.symm hs]

lemma P_mulVec (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) (ξ : ι → ℝ) :
    RiccatiFlow.P Y Z s *ᵥ (Y s *ᵥ ξ) = Z s *ᵥ ξ := by
  rw [Matrix.mulVec_mulVec, RiccatiFlow.P_mul_Y (h.detY hs)]

lemma Z_mulVec_ne (h : Model C Φ lam T φ Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) {ξ : ι → ℝ}
    (hξ : ξ ≠ 0) : Z s *ᵥ ξ ≠ 0 := by
  intro hZ
  have h1 := h.flow.le_form h.hlam h.hC hs ξ
  rw [hZ, dotProduct_zero] at h1
  have h2 := h.hΦ ξ
  have h3 : 0 < ξ ⬝ᵥ ξ := by
    rcases (Finset.sum_nonneg fun i _ => mul_self_nonneg (ξ i)).lt_or_eq with h' | h'
    · exact h'
    · exact absurd (dotProduct_self_eq_zero.mp h'.symm) hξ
  nlinarith [h.hφ]

variable {A A' B B' D D' : ℝ → MvPolynomial ι ℝ}

/-- The Hessian of `A` along a characteristic. -/
lemma hess_deriv (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D') {s : ℝ}
    (hs : s ∈ Icc 0 T) (ξ : ι → ℝ) (i j : ι) :
    HasDerivWithinAt (fun s => hess (A s) (Y s *ᵥ ξ) i j)
      ((-(1 / 2 : ℝ) • (hess (A s) (Y s *ᵥ ξ) * RiccatiFlow.P Y Z s)
        - (1 / 2 : ℝ) • (RiccatiFlow.P Y Z s * hess (A s) (Y s *ᵥ ξ))
        + srcHess (RiccatiFlow.P Y Z s) (Z s *ᵥ ξ)) i j) (Icc 0 T) s := by
  have hd := (g.A_diff.map ((pderiv i).toLinearMap ∘ₗ (pderiv j).toLinearMap)).hasDeriv_along hs
    (fun a => h.flow.hasDeriv_Y_mulVec hs ξ a)
  refine hd.congr_deriv ?_
  simp only [LinearMap.comp_apply, Derivation.coeFn_coe]
  rw [g.A_eq s hs, map_add, map_add, map_add, eval_pderiv2_polyGen _ (h.symm hs),
    eval_pderiv2_srcA, h.P_mulVec hs]
  simp only [Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul, srcHess,
    Matrix.of_apply]
  have e : ∑ a, 1 / 2 * (Z s *ᵥ ξ) a * eval (Y s *ᵥ ξ) (pderiv a (pderiv i (pderiv j (A s))))
      = ∑ a, (Z s *ᵥ ξ) a / 2 * eval (Y s *ᵥ ξ) (pderiv a (pderiv i (pderiv j (A s)))) :=
    Finset.sum_congr rfl fun a _ => by ring
  rw [e]
  ring

omit [DecidableEq ι] in
lemma form_srcHess (P : Matrix ι ι ℝ) (z w : ι → ℝ) :
    w ⬝ᵥ (srcHess P z *ᵥ w) = 1 / 256 * ∑ a, (z a) ^ 2 * ((P *ᵥ w) a) ^ 2 := by
  have h1 : ∀ a, ((P *ᵥ w) a) ^ 2 = ∑ i, ∑ j, w i * (P a i * P a j) * w j := by
    intro a
    rw [sq, Matrix.mulVec, dotProduct, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring
  have hR : 1 / 256 * ∑ a, (z a) ^ 2 * ((P *ᵥ w) a) ^ 2
      = ∑ a, ∑ i, ∑ j, w i * (1 / 256 * (z a ^ 2 * (P a i * P a j)) * w j) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [h1 a]
    simp only [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring
  rw [hR]
  simp only [dotProduct, Matrix.mulVec, srcHess, Matrix.of_apply, Finset.mul_sum, Finset.sum_mul]
  calc ∑ i, ∑ j, ∑ a, w i * (1 / 256 * (z a ^ 2 * (P a i * P a j)) * w j)
      = ∑ i, ∑ a, ∑ j, w i * (1 / 256 * (z a ^ 2 * (P a i * P a j)) * w j) :=
        Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ a, ∑ i, ∑ j, w i * (1 / 256 * (z a ^ 2 * (P a i * P a j)) * w j) := Finset.sum_comm

/-- Second derivative of `A` in a transported direction, along a characteristic. -/
lemma hess_form_deriv (h : Model C Φ lam T φ Y Z)
    (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D') {s : ℝ} (hs : s ∈ Icc 0 T) (ξ η : ι → ℝ) :
    HasDerivWithinAt (fun s => (Y s *ᵥ η) ⬝ᵥ (hess (A s) (Y s *ᵥ ξ) *ᵥ (Y s *ᵥ η)))
      (1 / 256 * ∑ a, ((Z s *ᵥ ξ) a) ^ 2 * ((Z s *ᵥ η) a) ^ 2) (Icc 0 T) s := by
  have hw : ∀ i, HasDerivWithinAt (fun r => (Y r *ᵥ η) i)
      (((1 / 2 : ℝ) • (RiccatiFlow.P Y Z s *ᵥ (Y s *ᵥ η))) i) (Icc 0 T) s := by
    intro i
    rw [h.P_mulVec hs]
    exact h.flow.hasDeriv_Y_mulVec hs η i
  refine (hasDeriv_quad hw (h.hess_deriv g hs ξ)).congr_deriv ?_
  set P := RiccatiFlow.P Y Z s with hP
  set M := hess (A s) (Y s *ᵥ ξ) with hM
  set w := Y s *ᵥ η with hwdef
  have hPs : P.IsSymm := h.symm hs
  have e1 : (P *ᵥ w) ⬝ᵥ (M *ᵥ w) = w ⬝ᵥ ((P * M) *ᵥ w) := by
    rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec w, ← Matrix.mulVec_transpose, hPs]
  have e2 : w ⬝ᵥ (M *ᵥ (P *ᵥ w)) = w ⬝ᵥ ((M * P) *ᵥ w) := by
    rw [Matrix.mulVec_mulVec]
  rw [Matrix.add_mulVec, Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.smul_mulVec,
    dotProduct_add, dotProduct_sub, dotProduct_smul, dotProduct_smul, smul_dotProduct,
    Matrix.mulVec_smul, dotProduct_smul, e1, e2, form_srcHess, h.P_mulVec hs]
  simp only [smul_eq_mul]
  ring

/-- The Hessian of `A` is positive semidefinite. -/
lemma hess_nonneg (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D')
    {s : ℝ} (hs : s ∈ Icc 0 T) (x v : ι → ℝ) : 0 ≤ v ⬝ᵥ (hess (A s) x *ᵥ v) := by
  have hunit := isUnit_iff_ne_zero.mpr (h.detY hs)
  have hY : ∀ x : ι → ℝ, Y s *ᵥ ((Y s)⁻¹ *ᵥ x) = x := fun x => by
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hunit, Matrix.one_mulVec]
  have hmono := monotoneOn_of_hasDerivWithinAt_nonneg (convex_Icc 0 T)
    (f := fun r => (Y r *ᵥ ((Y s)⁻¹ *ᵥ v)) ⬝ᵥ (hess (A r) (Y r *ᵥ ((Y s)⁻¹ *ᵥ x))
      *ᵥ (Y r *ᵥ ((Y s)⁻¹ *ᵥ v))))
    (fun r hr => (h.hess_form_deriv g hr _ _).continuousWithinAt)
    (fun r hr => (h.hess_form_deriv g (interior_subset hr) _ _).mono interior_subset)
    (fun r _ => mul_nonneg (by norm_num)
      (Finset.sum_nonneg fun a _ => mul_nonneg (sq_nonneg _) (sq_nonneg _)))
  have h0 := hmono ⟨le_rfl, h.hT.le⟩ hs hs.1
  simp only [hY] at h0
  refine le_trans (le_of_eq ?_) h0
  simp [g.A_zero, hess, Matrix.mulVec, dotProduct]

/-- Strict positivity in the radial direction. -/
lemma hess_radial_pos (h : Model C Φ lam T φ Y Z)
    (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D') {s : ℝ} (hs : s ∈ Icc 0 T) (hs0 : 0 < s)
    {x : ι → ℝ} (hx : x ≠ 0) : 0 < x ⬝ᵥ (hess (A s) x *ᵥ x) := by
  have hunit := isUnit_iff_ne_zero.mpr (h.detY hs)
  have hY : Y s *ᵥ ((Y s)⁻¹ *ᵥ x) = x := by
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hunit, Matrix.one_mulVec]
  have hξ : (Y s)⁻¹ *ᵥ x ≠ 0 := fun h0 => hx (by rw [← hY, h0, Matrix.mulVec_zero])
  have hmono := strictMonoOn_of_hasDerivWithinAt_pos (convex_Icc 0 T)
    (f := fun r => (Y r *ᵥ ((Y s)⁻¹ *ᵥ x)) ⬝ᵥ (hess (A r) (Y r *ᵥ ((Y s)⁻¹ *ᵥ x))
      *ᵥ (Y r *ᵥ ((Y s)⁻¹ *ᵥ x))))
    (fun r hr => (h.hess_form_deriv g hr _ _).continuousWithinAt)
    (fun r hr => (h.hess_form_deriv g (interior_subset hr) _ _).mono interior_subset)
    (fun r hr => by
      have hne := h.Z_mulVec_ne (interior_subset hr) hξ
      obtain ⟨a, ha⟩ := Function.ne_iff.mp hne
      refine mul_pos (by norm_num) (Finset.sum_pos' (fun a _ =>
        mul_nonneg (sq_nonneg _) (sq_nonneg _)) ⟨a, Finset.mem_univ a, ?_⟩)
      have : 0 < ((Z r *ᵥ ((Y s)⁻¹ *ᵥ x)) a) ^ 2 := by
        have : (Z r *ᵥ ((Y s)⁻¹ *ᵥ x)) a ≠ 0 := ha
        positivity
      exact mul_pos this this)
  have h0 := hmono ⟨le_rfl, h.hT.le⟩ hs hs0
  simp only [hY] at h0
  refine lt_of_eq_of_lt ?_ h0
  simp [g.A_zero, hess, Matrix.mulVec, dotProduct]

/-- **Quartic part is positive.** -/
theorem A_pos (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D')
    {s : ℝ} (hs : s ∈ Icc 0 T) (hs0 : 0 < s) {x : ι → ℝ} (hx : x ≠ 0) : 0 < eval x (A s) := by
  have hunit := isUnit_iff_ne_zero.mpr (h.detY hs)
  have hY : Y s *ᵥ ((Y s)⁻¹ *ᵥ x) = x := by
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hunit, Matrix.one_mulVec]
  have hξ : (Y s)⁻¹ *ᵥ x ≠ 0 := fun h0 => hx (by rw [← hY, h0, Matrix.mulVec_zero])
  have hd : ∀ r ∈ Icc 0 T, HasDerivWithinAt (fun r => eval (Y r *ᵥ ((Y s)⁻¹ *ᵥ x)) (A r))
      (1 / 3072 * ∑ a, ((Z r *ᵥ ((Y s)⁻¹ *ᵥ x)) a) ^ 4) (Icc 0 T) r := by
    intro r hr
    refine (transport_along h.flow hr (h.detY hr) (f := fun r => srcA (RiccatiFlow.P Y Z r))
      g.A_diff (g.A_eq r hr) _).congr_deriv ?_
    rw [eval_srcA, h.P_mulVec hr]
  have hmono := strictMonoOn_of_hasDerivWithinAt_pos (convex_Icc 0 T)
    (fun r hr => (hd r hr).continuousWithinAt)
    (fun r hr => (hd r (interior_subset hr)).mono interior_subset)
    (fun r hr => by
      have hne := h.Z_mulVec_ne (interior_subset hr) hξ
      obtain ⟨a, ha⟩ := Function.ne_iff.mp hne
      refine mul_pos (by norm_num) (Finset.sum_pos' (fun a _ => by positivity)
        ⟨a, Finset.mem_univ a, ?_⟩)
      have : (Z r *ᵥ ((Y s)⁻¹ *ᵥ x)) a ≠ 0 := ha
      positivity)
  have h0 := hmono ⟨le_rfl, h.hT.le⟩ hs hs0
  simp only [hY, g.A_zero, map_zero] at h0
  exact h0

/-- Derivative of `tr (∇²A · P⁻¹)` along a characteristic. -/
lemma trace_deriv (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D')
    {s : ℝ} (hs : s ∈ Icc 0 T) (ξ : ι → ℝ) :
    HasDerivWithinAt
      (fun s => (hess (A s) (Y s *ᵥ ξ) * (RiccatiFlow.P Y Z s)⁻¹).trace)
      (-(1 / 2) * (hess (A s) (Y s *ᵥ ξ)).trace
        + 1 / 256 * ∑ a, ((Z s *ᵥ ξ) a) ^ 2 * RiccatiFlow.P Y Z s a a
        - 2 * lam * (hess (A s) (Y s *ᵥ ξ)
          * ((RiccatiFlow.P Y Z s)⁻¹ * C * (RiccatiFlow.P Y Z s)⁻¹)).trace) (Icc 0 T) s := by
  have hRd := hasDeriv_inv h.hT hs (fun r hr => h.detP hr)
    (M := fun r => RiccatiFlow.P Y Z r)
    (M' := (2 * lam) • C - (1 / 2 : ℝ) • (RiccatiFlow.P Y Z s * RiccatiFlow.P Y Z s))
    (fun i j => by
      refine (h.deriv hs i j).congr_deriv ?_
      simp [Matrix.smul_apply, Matrix.sub_apply])
  have hd : HasDerivWithinAt
      (fun s => (hess (A s) (Y s *ᵥ ξ) * (RiccatiFlow.P Y Z s)⁻¹).trace)
      (((-(1 / 2 : ℝ) • (hess (A s) (Y s *ᵥ ξ) * RiccatiFlow.P Y Z s)
        - (1 / 2 : ℝ) • (RiccatiFlow.P Y Z s * hess (A s) (Y s *ᵥ ξ))
        + srcHess (RiccatiFlow.P Y Z s) (Z s *ᵥ ξ)) * (RiccatiFlow.P Y Z s)⁻¹
        + hess (A s) (Y s *ᵥ ξ) * (-((RiccatiFlow.P Y Z s)⁻¹ * ((2 * lam) • C
          - (1 / 2 : ℝ) • (RiccatiFlow.P Y Z s * RiccatiFlow.P Y Z s))
          * (RiccatiFlow.P Y Z s)⁻¹))).trace) (Icc 0 T) s := by
    simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply (M := hess _ _)]
    refine (HasDerivWithinAt.fun_sum fun i _ => HasDerivWithinAt.fun_sum fun k _ =>
      (h.hess_deriv g hs ξ i k).fun_mul (hRd k i)).congr_deriv ?_
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Matrix.add_apply, Matrix.mul_apply, Matrix.mul_apply, ← Finset.sum_add_distrib]
  refine hd.congr_deriv ?_
  set P := RiccatiFlow.P Y Z s with hP
  set M := hess (A s) (Y s *ᵥ ξ) with hM
  have hPR : P * P⁻¹ = 1 := h.mul_inv hs
  have hRP : P⁻¹ * P = 1 := h.inv_mul hs
  have hPs : P.IsSymm := h.symm hs
  have t1 : (M * P * P⁻¹).trace = M.trace := by
    rw [Matrix.mul_assoc, hPR, Matrix.mul_one]
  have t2 : (P * M * P⁻¹).trace = M.trace := by
    rw [Matrix.trace_mul_comm, ← Matrix.mul_assoc, hRP, Matrix.one_mul]
  have t3 : (srcHess P (Z s *ᵥ ξ) * P⁻¹).trace
      = 1 / 256 * ∑ a, ((Z s *ᵥ ξ) a) ^ 2 * P a a := by
    have hid : ∀ a, ∑ i, ∑ k, P a i * P a k * P⁻¹ k i = P a a := by
      intro a
      have h1 : ∀ i, ∑ k, P a i * P a k * P⁻¹ k i = P a i * (1 : Matrix ι ι ℝ) a i := by
        intro i
        rw [← hPR, Matrix.mul_apply, Finset.mul_sum]
        exact Finset.sum_congr rfl fun k _ => by ring
      simp only [h1, Matrix.one_apply]
      simp
    simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, srcHess, Matrix.of_apply,
      Finset.mul_sum, Finset.sum_mul]
    calc ∑ i, ∑ k, ∑ a, 1 / 256 * ((Z s *ᵥ ξ) a ^ 2 * (P a i * P a k)) * P⁻¹ k i
        = ∑ i, ∑ a, ∑ k, 1 / 256 * ((Z s *ᵥ ξ) a ^ 2 * (P a i * P a k)) * P⁻¹ k i :=
          Finset.sum_congr rfl fun i _ => Finset.sum_comm
      _ = ∑ a, ∑ i, ∑ k, 1 / 256 * ((Z s *ᵥ ξ) a ^ 2 * (P a i * P a k)) * P⁻¹ k i :=
          Finset.sum_comm
      _ = ∑ a, 1 / 256 * ((Z s *ᵥ ξ) a ^ 2 * P a a) := by
          refine Finset.sum_congr rfl fun a _ => ?_
          rw [← hid a, Finset.mul_sum, Finset.mul_sum]
          refine Finset.sum_congr rfl fun i _ => ?_
          rw [Finset.mul_sum, Finset.mul_sum]
          exact Finset.sum_congr rfl fun k _ => by ring
  have t4 : (P⁻¹ * ((2 * lam) • C - (1 / 2 : ℝ) • (P * P)) * P⁻¹)
      = (2 * lam) • (P⁻¹ * C * P⁻¹) - (1 / 2 : ℝ) • 1 := by
    rw [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_smul, Matrix.mul_smul, Matrix.smul_mul,
      Matrix.smul_mul, ← Matrix.mul_assoc P⁻¹ P P, hRP, Matrix.one_mul, hPR]
  rw [t4, Matrix.trace_add, Matrix.add_mul, Matrix.sub_mul, Matrix.smul_mul, Matrix.smul_mul,
    Matrix.trace_add, Matrix.trace_sub, Matrix.trace_smul, Matrix.trace_smul, t1, t2, t3,
    Matrix.mul_neg, Matrix.trace_neg, Matrix.mul_sub, Matrix.trace_sub, Matrix.mul_smul,
    Matrix.mul_smul, Matrix.trace_smul, Matrix.trace_smul, Matrix.mul_one]
  simp only [smul_eq_mul]
  ring

/-- Derivative of `B + 2 tr (∇²A · P⁻¹)` along a characteristic. -/
lemma B_trace_deriv (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D')
    {r : ℝ} (hr : r ∈ Icc 0 T) (ξ : ι → ℝ) :
    HasDerivWithinAt
      (fun r => eval (Y r *ᵥ ξ) (B r) + 2 * (hess (A r) (Y r *ᵥ ξ) * (RiccatiFlow.P Y Z r)⁻¹).trace)
      (-(4 * lam * (hess (A r) (Y r *ᵥ ξ)
        * ((RiccatiFlow.P Y Z r)⁻¹ * C * (RiccatiFlow.P Y Z r)⁻¹)).trace)) (Icc 0 T) r := by
  have hB := transport_along h.flow hr (h.detY hr)
    (f := fun r => srcB (RiccatiFlow.P Y Z r) (A r)) g.B_diff (g.B_eq r hr) ξ
  have hT := (h.trace_deriv g hr ξ).const_mul 2
  refine (hB.fun_add hT).congr_deriv ?_
  rw [eval_srcB, h.P_mulVec hr]
  have e1 : ∑ a, eval (Y r *ᵥ ξ) (pderiv a (pderiv a (A r))) = (hess (A r) (Y r *ᵥ ξ)).trace := by
    simp [Matrix.trace, hess]
  have e2 : ∑ a, RiccatiFlow.P Y Z r a a * ((Z r *ᵥ ξ) a) ^ 2
      = ∑ a, ((Z r *ᵥ ξ) a) ^ 2 * RiccatiFlow.P Y Z r a a :=
    Finset.sum_congr rfl fun a _ => mul_comm _ _
  rw [e1, e2]
  generalize (hess (A r) (Y r *ᵥ ξ)).trace = t
  generalize (hess (A r) (Y r *ᵥ ξ)
    * ((RiccatiFlow.P Y Z r)⁻¹ * C * (RiccatiFlow.P Y Z r)⁻¹)).trace = t'
  ring

/-- The derivative of `B_trace_deriv` is nonpositive. -/
lemma trace_RCR_nonneg (h : Model C Φ lam T φ Y Z)
    (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D') {r : ℝ} (hr : r ∈ Icc 0 T) (y : ι → ℝ) :
    0 ≤ (hess (A r) y * ((RiccatiFlow.P Y Z r)⁻¹ * C * (RiccatiFlow.P Y Z r)⁻¹)).trace := by
  have hRs : ((RiccatiFlow.P Y Z r)⁻¹).IsSymm := h.inv_symm hr
  have hN0' := h.hess_nonneg g hr y
  generalize (RiccatiFlow.P Y Z r)⁻¹ = R at *
  have hN : (R * C * R).IsSymm := by
    rw [Matrix.IsSymm, Matrix.transpose_mul, Matrix.transpose_mul, hRs, h.hCs, Matrix.mul_assoc]
  have hN0 : ∀ v, 0 ≤ v ⬝ᵥ ((R * C * R) *ᵥ v) := by
    intro v
    have e : v ⬝ᵥ ((R * C * R) *ᵥ v) = (R *ᵥ v) ⬝ᵥ (C *ᵥ (R *ᵥ v)) := by
      rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec v R,
        ← Matrix.mulVec_transpose, hRs]
    rw [e]
    exact h.hC _
  have e : (hess (A r) y * (R * C * R)).trace
      = ∑ i, ∑ j, (R * C * R) i j * hess (A r) y i j := by
    simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply (M := hess _ _)]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [hN.apply i j, mul_comm]
  rw [e]
  exact sum_mul_nonneg hN hN0 hN0'

/-- `B + 2 tr (∇²A · P⁻¹) ≤ 0`. -/
lemma B_le (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D') {s : ℝ}
    (hs : s ∈ Icc 0 T) (x : ι → ℝ) :
    eval x (B s) + 2 * (hess (A s) x * (RiccatiFlow.P Y Z s)⁻¹).trace ≤ 0 := by
  have hunit := isUnit_iff_ne_zero.mpr (h.detY hs)
  obtain ⟨ξ, hξ⟩ : ∃ ξ, Y s *ᵥ ξ = x := ⟨(Y s)⁻¹ *ᵥ x, by
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hunit, Matrix.one_mulVec]⟩
  have hanti := antitoneOn_of_hasDerivWithinAt_nonpos (convex_Icc 0 T)
    (fun r hr => (h.B_trace_deriv g hr ξ).continuousWithinAt)
    (fun r hr => (h.B_trace_deriv g (interior_subset hr) ξ).mono interior_subset)
    (fun r hr => by
      have ht := h.trace_RCR_nonneg g (interior_subset hr) (Y r *ᵥ ξ)
      have := mul_nonneg (mul_nonneg (by norm_num : (0 : ℝ) ≤ 4) h.hlam) ht
      linarith)
  have h0 := hanti ⟨le_rfl, h.hT.le⟩ hs hs.1
  simp only [hξ] at h0
  refine h0.trans (le_of_eq ?_)
  simp [g.A_zero, g.B_zero, hess, Matrix.trace, Matrix.mul_apply]

/-- **Quadratic part is negative.** -/
theorem B_neg (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D')
    {s : ℝ} (hs : s ∈ Icc 0 T) (hs0 : 0 < s) {x : ι → ℝ} (hx : x ≠ 0) : eval x (B s) < 0 := by
  have h1 := h.B_le g hs x
  have h2 := h.hess_radial_pos g hs hs0 hx
  have h3 := h.pos_form hs hx
  have h4 := form_le_trace (h.symm hs) (h.nonneg_form hs) (h.mul_inv hs) (hess_symm (A s) x)
    (h.hess_nonneg g hs x) x
  have e : ∑ i, ∑ j, hess (A s) x i j * (RiccatiFlow.P Y Z s)⁻¹ i j
      = (hess (A s) x * (RiccatiFlow.P Y Z s)⁻¹).trace := by
    simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply (M := hess _ _)]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [(h.inv_symm hs).apply j i]
  rw [e] at h4
  have h5 : 0 < (hess (A s) x * (RiccatiFlow.P Y Z s)⁻¹).trace := by
    by_contra hcon
    have := mul_nonpos_of_nonneg_of_nonpos h3.le (not_lt.mp hcon)
    linarith
  linarith

end Model

end ExecutionProofs
