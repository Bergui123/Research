import ExecutionProofs.LinearSystem

/-!
# The Riccati flow

`Y' = Z / 2`, `Z' = 2 λ C Y`, `Y(0) = 1`, `Z(0) = Φ`. Existence, invertibility of `Y`, and for
`P = Z Y⁻¹`: symmetry, coercivity and `P' = 2 λ C - P² / 2`. Everything follows from
`Yᵀ Z ⪰ Φ`.
-/

open Set Matrix

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The block system of the paper (Eq. 3) with `H = 1`, entrywise on `[0, T]`. -/
structure RiccatiFlow (C Φ : Matrix ι ι ℝ) (lam T : ℝ) (Y Z : ℝ → Matrix ι ι ℝ) : Prop where
  Y_zero : Y 0 = 1
  Z_zero : Z 0 = Φ
  Y_deriv : ∀ s ∈ Icc 0 T, ∀ i j,
    HasDerivWithinAt (fun s => Y s i j) (1 / 2 * Z s i j) (Icc 0 T) s
  Z_deriv : ∀ s ∈ Icc 0 T, ∀ i j,
    HasDerivWithinAt (fun s => Z s i j) (2 * lam * (C * Y s) i j) (Icc 0 T) s

/-- The flow exists. -/
theorem exists_riccatiFlow (C Φ : Matrix ι ι ℝ) (lam : ℝ) {T : ℝ} (hT : 0 ≤ T) :
    ∃ Y Z, RiccatiFlow C Φ lam T Y Z := by
  classical
  -- unknowns `(false, i, j) ↦ Y i j`, `(true, i, j) ↦ Z i j`
  let κ : Bool × ι → Bool × ι × ι → ℝ → ℝ := fun d n _ =>
    match d, n with
    | (true, k), (false, i, _) => if k = i then 1 / 2 else 0
    | (false, k), (true, i, _) => 2 * lam * C i k
    | _, _ => 0
  let u₀ : Bool × ι × ι → ℝ := fun n => if n.1 then Φ n.2.1 n.2.2 else (1 : Matrix ι ι ℝ) n.2.1 n.2.2
  obtain ⟨Kb, hKb⟩ : ∃ Kb : ℝ, ∀ d n, ∀ s ∈ Icc 0 T, |κ d n s| ≤ Kb := by
    refine ⟨∑ d, ∑ n, |κ d n 0|, fun d n s _ => ?_⟩
    have h1 : |κ d n s| = |κ d n 0| := rfl
    rw [h1]
    exact (Finset.single_le_sum (f := fun n => |κ d n 0|) (fun _ _ => abs_nonneg _)
      (Finset.mem_univ n)).trans (Finset.single_le_sum (f := fun d => ∑ n, |κ d n 0|)
      (fun _ _ => Finset.sum_nonneg fun _ _ => abs_nonneg _) (Finset.mem_univ d))
  obtain ⟨A, hA⟩ : ∃ A : ℝ, ∀ n, |u₀ n| ≤ A :=
    ⟨∑ n, |u₀ n|, fun n => Finset.single_le_sum (f := fun n => |u₀ n|)
      (fun _ _ => abs_nonneg _) (Finset.mem_univ n)⟩
  obtain ⟨u, h0, -, hu⟩ := exists_linear_solution
    (nb := fun (d : Bool × ι) (n : Bool × ι × ι) => (d.1, d.2, n.2.2)) hT (κ := κ)
    (fun d n => continuousOn_const) hKb (V := fun _ => 0) (fun _ => le_rfl)
    (b := fun _ _ => 0) (fun _ => continuousOn_const) (β := 0) (fun _ _ _ => by simp) hA
  refine ⟨fun s => Matrix.of fun i j => u (false, i, j) s,
    fun s => Matrix.of fun i j => u (true, i, j) s, ?_, ?_, ?_, ?_⟩
  · ext i j
    simpa [u₀] using h0 (false, i, j)
  · ext i j
    simpa [u₀] using h0 (true, i, j)
  · intro s hs i j
    refine (hu (false, i, j) s hs).congr_deriv ?_
    simp [κ, Fintype.sum_prod_type]
  · intro s hs i j
    refine (hu (true, i, j) s hs).congr_deriv ?_
    simp [κ, Fintype.sum_prod_type, Matrix.mul_apply, Finset.mul_sum, mul_assoc]

/-- Entries of the inverse of a differentiable family are differentiable. -/
lemma differentiableWithinAt_inv {M : ℝ → Matrix ι ι ℝ} {S : Set ℝ} {s : ℝ}
    (hM : ∀ i j, DifferentiableWithinAt ℝ (fun s => M s i j) S s) (hdet : (M s).det ≠ 0) (i j : ι) :
    DifferentiableWithinAt ℝ (fun s => (M s)⁻¹ i j) S s := by
  have hdetd : ∀ N : ℝ → Matrix ι ι ℝ, (∀ i j, DifferentiableWithinAt ℝ (fun s => N s i j) S s) →
      DifferentiableWithinAt ℝ (fun s => (N s).det) S s := by
    intro N hN
    simp only [Matrix.det_apply']
    refine DifferentiableWithinAt.fun_sum fun σ _ => ?_
    exact (DifferentiableWithinAt.fun_finsetProd fun k _ => hN _ _).const_mul _
  have hadj : DifferentiableWithinAt ℝ (fun s => (M s).adjugate i j) S s := by
    simp only [Matrix.adjugate_apply]
    refine hdetd _ fun a b => ?_
    by_cases h : a = j
    · subst h
      simp only [Matrix.updateRow_self]
      exact differentiableWithinAt_const _
    · simp only [Matrix.updateRow_ne h]
      exact hM a b
  have h : ∀ s, (M s)⁻¹ i j = ((M s).det)⁻¹ * (M s).adjugate i j := by
    intro s
    rw [Matrix.inv_def, Matrix.smul_apply, Ring.inverse_eq_inv', smul_eq_mul]
  simp only [h]
  exact ((hdetd M hM).inv hdet).mul hadj

namespace RiccatiFlow

variable {C Φ : Matrix ι ι ℝ} {lam T : ℝ} {Y Z : ℝ → Matrix ι ι ℝ}

/-- Derivative of `(Y v)_k`. -/
lemma hasDeriv_Y_mulVec (h : RiccatiFlow C Φ lam T Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) (v : ι → ℝ)
    (k : ι) : HasDerivWithinAt (fun s => (Y s *ᵥ v) k) (1 / 2 * (Z s *ᵥ v) k) (Icc 0 T) s := by
  simp only [Matrix.mulVec, dotProduct]
  refine (HasDerivWithinAt.fun_sum fun j _ => (h.Y_deriv s hs k j).mul_const (v j)).congr_deriv ?_
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- Derivative of `(Z v)_k`. -/
lemma hasDeriv_Z_mulVec (h : RiccatiFlow C Φ lam T Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) (v : ι → ℝ)
    (k : ι) :
    HasDerivWithinAt (fun s => (Z s *ᵥ v) k) (2 * lam * (C *ᵥ (Y s *ᵥ v)) k) (Icc 0 T) s := by
  rw [Matrix.mulVec_mulVec]
  simp only [Matrix.mulVec, dotProduct]
  refine (HasDerivWithinAt.fun_sum fun j _ => (h.Z_deriv s hs k j).mul_const (v j)).congr_deriv ?_
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- Derivative of the bilinear form `(Y v) · (Z w)`. -/
lemma hasDeriv_form (h : RiccatiFlow C Φ lam T Y Z) {s : ℝ} (hs : s ∈ Icc 0 T) (v w : ι → ℝ) :
    HasDerivWithinAt (fun s => (Y s *ᵥ v) ⬝ᵥ (Z s *ᵥ w))
      (1 / 2 * ((Z s *ᵥ v) ⬝ᵥ (Z s *ᵥ w)) + 2 * lam * ((Y s *ᵥ v) ⬝ᵥ (C *ᵥ (Y s *ᵥ w))))
      (Icc 0 T) s := by
  simp only [dotProduct]
  refine (HasDerivWithinAt.fun_sum fun k _ =>
    (h.hasDeriv_Y_mulVec hs v k).mul (h.hasDeriv_Z_mulVec hs w k)).congr_deriv ?_
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  ring

/-- The form `(Y v) · (Z w)` is symmetric. -/
lemma form_symm (h : RiccatiFlow C Φ lam T Y Z) (hC : C.IsSymm) (hΦ : Φ.IsSymm) {s : ℝ}
    (hs : s ∈ Icc 0 T) (v w : ι → ℝ) :
    (Y s *ᵥ v) ⬝ᵥ (Z s *ᵥ w) = (Y s *ᵥ w) ⬝ᵥ (Z s *ᵥ v) := by
  have hsymm : ∀ {M : Matrix ι ι ℝ}, M.IsSymm → ∀ a b : ι → ℝ, a ⬝ᵥ (M *ᵥ b) = b ⬝ᵥ (M *ᵥ a) := by
    intro M hM a b
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hM, dotProduct_comm]
  set f : ℝ → ℝ := fun s => (Y s *ᵥ v) ⬝ᵥ (Z s *ᵥ w) - (Y s *ᵥ w) ⬝ᵥ (Z s *ᵥ v) with hf
  have hderiv : ∀ r ∈ Icc 0 T, HasDerivWithinAt f 0 (Icc 0 T) r := by
    intro r hr
    refine ((h.hasDeriv_form hr v w).sub (h.hasDeriv_form hr w v)).congr_deriv ?_
    rw [hsymm hC (Y r *ᵥ v) (Y r *ᵥ w), dotProduct_comm (Z r *ᵥ v)]
    ring
  have hconst : f s = f 0 := by
    have := constant_of_derivWithin_zero (f := f) (a := 0) (b := T)
      (fun r hr => (hderiv r hr).differentiableWithinAt)
      (fun r hr => by
        rcases (Icc_subset_Icc_left le_rfl : Icc 0 T ⊆ Icc 0 T) with _
        have hr' : r ∈ Icc 0 T := ⟨hr.1, hr.2.le⟩
        exact (hderiv r hr').derivWithin
          (uniqueDiffOn_Icc (hr.1.trans_lt hr.2) r hr'))
    exact this s hs
  have h0 : f 0 = 0 := by
    simp only [hf, h.Y_zero, h.Z_zero, Matrix.one_mulVec]
    rw [hsymm hΦ v w, sub_self]
  have : f s = 0 := hconst.trans h0
  linarith [this, show f s = (Y s *ᵥ v) ⬝ᵥ (Z s *ᵥ w) - (Y s *ᵥ w) ⬝ᵥ (Z s *ᵥ v) from rfl]

/-- `Yᵀ Z ⪰ Φ`. -/
lemma le_form (h : RiccatiFlow C Φ lam T Y Z) (hlam : 0 ≤ lam)
    (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v)) {s : ℝ} (hs : s ∈ Icc 0 T) (v : ι → ℝ) :
    v ⬝ᵥ (Φ *ᵥ v) ≤ (Y s *ᵥ v) ⬝ᵥ (Z s *ᵥ v) := by
  have hmono : MonotoneOn (fun s => (Y s *ᵥ v) ⬝ᵥ (Z s *ᵥ v)) (Icc 0 T) := by
    refine monotoneOn_of_hasDerivWithinAt_nonneg (convex_Icc 0 T)
      (fun r hr => (h.hasDeriv_form hr v v).continuousWithinAt)
      (fun r hr => (h.hasDeriv_form (interior_subset hr) v v).mono interior_subset)
      (fun r hr => ?_)
    have h1 : 0 ≤ (Z r *ᵥ v) ⬝ᵥ (Z r *ᵥ v) := Finset.sum_nonneg fun i _ => mul_self_nonneg _
    have h2 := hC (Y r *ᵥ v)
    positivity
  have := hmono ⟨le_rfl, hs.1.trans hs.2⟩ hs hs.1
  simpa [h.Y_zero, h.Z_zero] using this

variable {φ : ℝ}

/-- `Y` stays invertible. -/
lemma det_ne_zero (h : RiccatiFlow C Φ lam T Y Z) (hlam : 0 ≤ lam)
    (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v)) (hφ : 0 < φ) (hΦ : ∀ v, φ * (v ⬝ᵥ v) ≤ v ⬝ᵥ (Φ *ᵥ v)) {s : ℝ}
    (hs : s ∈ Icc 0 T) : (Y s).det ≠ 0 := by
  intro hdet
  obtain ⟨v, hv, hYv⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet
  have h1 := h.le_form hlam hC hs v
  rw [hYv, zero_dotProduct] at h1
  have h2 := hΦ v
  have h3 : 0 < v ⬝ᵥ v := by
    rcases (Finset.sum_nonneg fun i _ => mul_self_nonneg (v i)).lt_or_eq with h' | h'
    · exact h'
    · exact absurd (dotProduct_self_eq_zero.mp h'.symm) hv
  nlinarith

/-- `Y` stays invertible already for `Φ ⪰ 0`. -/
lemma det_ne_zero_of_nonneg (h : RiccatiFlow C Φ lam T Y Z) (hlam : 0 ≤ lam)
    (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v)) (hΦ : ∀ v, 0 ≤ v ⬝ᵥ (Φ *ᵥ v)) {s : ℝ} (hs : s ∈ Icc 0 T) :
    (Y s).det ≠ 0 := by
  intro hdet
  obtain ⟨v, hv, hYv⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet
  have hsub : Icc 0 s ⊆ Icc 0 T := Icc_subset_Icc_right hs.2
  rcases hs.1.eq_or_lt with h0 | hpos
  · rw [← h0, h.Y_zero, Matrix.one_mulVec] at hYv
    exact hv hYv
  -- the form is nonnegative, nondecreasing and vanishes at `s`, so it vanishes on `[0, s]`
  have hmono : MonotoneOn (fun r => (Y r *ᵥ v) ⬝ᵥ (Z r *ᵥ v)) (Icc 0 T) := by
    refine monotoneOn_of_hasDerivWithinAt_nonneg (convex_Icc 0 T)
      (fun r hr => (h.hasDeriv_form hr v v).continuousWithinAt)
      (fun r hr => (h.hasDeriv_form (interior_subset hr) v v).mono interior_subset)
      (fun r hr => ?_)
    have h1 : 0 ≤ (Z r *ᵥ v) ⬝ᵥ (Z r *ᵥ v) := Finset.sum_nonneg fun i _ => mul_self_nonneg _
    have h2 := hC (Y r *ᵥ v)
    positivity
  have hform0 : ∀ r ∈ Icc 0 s, (Y r *ᵥ v) ⬝ᵥ (Z r *ᵥ v) = 0 := by
    intro r hr
    refine le_antisymm ?_ ((hΦ v).trans (h.le_form hlam hC (hsub hr) v))
    have := hmono (hsub hr) hs hr.2
    simpa [hYv] using this
  have hZ0 : ∀ r ∈ Icc 0 s, Z r *ᵥ v = 0 := by
    intro r hr
    have hd := (h.hasDeriv_form (hsub hr) v v).mono hsub
    have hc : HasDerivWithinAt (fun r => (Y r *ᵥ v) ⬝ᵥ (Z r *ᵥ v)) 0 (Icc 0 s) r :=
      (hasDerivWithinAt_const r (Icc 0 s) (0 : ℝ)).congr hform0 (hform0 r hr)
    have heq := (uniqueDiffOn_Icc hpos r hr).eq_deriv _ hd hc
    have h1 : 0 ≤ (Z r *ᵥ v) ⬝ᵥ (Z r *ᵥ v) := Finset.sum_nonneg fun i _ => mul_self_nonneg _
    have h2 : 0 ≤ lam * ((Y r *ᵥ v) ⬝ᵥ (C *ᵥ (Y r *ᵥ v))) := mul_nonneg hlam (hC _)
    exact dotProduct_self_eq_zero.mp (by linarith)
  refine hv (funext fun k => ?_)
  have hconst := constant_of_derivWithin_zero (f := fun r => (Y r *ᵥ v) k) (a := 0) (b := s)
    (fun r hr => ((h.hasDeriv_Y_mulVec (hsub hr) v k).mono hsub).differentiableWithinAt)
    (fun r hr => by
      have hr' : r ∈ Icc 0 s := ⟨hr.1, hr.2.le⟩
      have hd := (h.hasDeriv_Y_mulVec (hsub hr') v k).mono hsub
      rw [hZ0 r hr', Pi.zero_apply, mul_zero] at hd
      exact hd.derivWithin (uniqueDiffOn_Icc hpos r hr')) s ⟨hs.1, le_rfl⟩
  simpa [hYv, h.Y_zero] using hconst.symm

/-- The Riccati matrix. -/
noncomputable def P (Y Z : ℝ → Matrix ι ι ℝ) (s : ℝ) : Matrix ι ι ℝ := Z s * (Y s)⁻¹

lemma P_mul_Y {s : ℝ} (hdet : (Y s).det ≠ 0) : P Y Z s * Y s = Z s := by
  rw [P, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.mpr hdet), Matrix.mul_one]

/-- `x · P x = (Y v) · (Z v)` with `v = Y⁻¹ x`. -/
lemma P_form {s : ℝ} (hdet : (Y s).det ≠ 0) (x y : ι → ℝ) :
    x ⬝ᵥ (P Y Z s *ᵥ y) = (Y s *ᵥ ((Y s)⁻¹ *ᵥ x)) ⬝ᵥ (Z s *ᵥ ((Y s)⁻¹ *ᵥ y)) := by
  have hY : Y s *ᵥ ((Y s)⁻¹ *ᵥ x) = x := by
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ (isUnit_iff_ne_zero.mpr hdet),
      Matrix.one_mulVec]
  rw [hY, P, ← Matrix.mulVec_mulVec]

/-- `P` is symmetric. -/
lemma P_symm_of_det (h : RiccatiFlow C Φ lam T Y Z) (hCs : C.IsSymm) (hΦs : Φ.IsSymm) {s : ℝ}
    (hs : s ∈ Icc 0 T) (hdet : (Y s).det ≠ 0) : (P Y Z s).IsSymm := by
  ext i j
  have key := h.form_symm hCs hΦs hs ((Y s)⁻¹ *ᵥ Pi.single j 1) ((Y s)⁻¹ *ᵥ Pi.single i 1)
  rw [← P_form hdet, ← P_form hdet] at key
  simpa [Matrix.transpose_apply] using key

lemma P_symm (h : RiccatiFlow C Φ lam T Y Z) (hlam : 0 ≤ lam) (hCs : C.IsSymm) (hΦs : Φ.IsSymm)
    (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v)) (hφ : 0 < φ) (hΦ : ∀ v, φ * (v ⬝ᵥ v) ≤ v ⬝ᵥ (Φ *ᵥ v)) {s : ℝ}
    (hs : s ∈ Icc 0 T) : (P Y Z s).IsSymm :=
  h.P_symm_of_det hCs hΦs hs (h.det_ne_zero hlam hC hφ hΦ hs)

/-- Uniform coercivity of `P`. -/
lemma exists_coercive (h : RiccatiFlow C Φ lam T Y Z) (hlam : 0 ≤ lam)
    (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v)) (hφ : 0 < φ) (hΦ : ∀ v, φ * (v ⬝ᵥ v) ≤ v ⬝ᵥ (Φ *ᵥ v)) :
    ∃ m, 0 < m ∧ ∀ s ∈ Icc 0 T, ∀ x, m * (x ⬝ᵥ x) ≤ x ⬝ᵥ (P Y Z s *ᵥ x) := by
  -- bound the entries of `Y` on `[0, T]`
  have hbdd : ∀ i j, ∃ B, ∀ s ∈ Icc 0 T, |Y s i j| ≤ B := by
    intro i j
    obtain ⟨B, hB⟩ := isCompact_Icc.exists_bound_of_continuousOn
      (fun s hs => (h.Y_deriv s hs i j).continuousWithinAt)
    exact ⟨B, hB⟩
  choose B hB using hbdd
  set K : ℝ := 1 + ∑ i, ∑ j, |B i j| with hK
  have hK1 : 0 < K := by
    have : 0 ≤ ∑ i, ∑ j, |B i j| :=
      Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => abs_nonneg _
    linarith
  have hYK : ∀ s ∈ Icc 0 T, ∀ i j, |Y s i j| ≤ K := by
    intro s hs i j
    have h1 : |B i j| ≤ ∑ j, |B i j| :=
      Finset.single_le_sum (f := fun j => |B i j|) (fun _ _ => abs_nonneg _) (Finset.mem_univ j)
    have h2 : ∑ j, |B i j| ≤ ∑ i, ∑ j, |B i j| :=
      Finset.single_le_sum (f := fun i => ∑ j, |B i j|)
        (fun _ _ => Finset.sum_nonneg fun _ _ => abs_nonneg _) (Finset.mem_univ i)
    have h3 := (hB i j s hs).trans (le_abs_self _)
    linarith
  set n : ℝ := Fintype.card ι + 1 with hn
  have hn0 : 0 < n := by positivity
  refine ⟨φ / (K ^ 2 * n ^ 2), by positivity, fun s hs x => ?_⟩
  have hdet := h.det_ne_zero hlam hC hφ hΦ hs
  set v := (Y s)⁻¹ *ᵥ x with hv
  have hYv : Y s *ᵥ v = x := by
    rw [hv, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ (isUnit_iff_ne_zero.mpr hdet),
      Matrix.one_mulVec]
  have hform : φ * (v ⬝ᵥ v) ≤ x ⬝ᵥ (P Y Z s *ᵥ x) := by
    rw [P_form hdet]
    exact (hΦ v).trans (h.le_form hlam hC hs v)
  have hvv : 0 ≤ v ⬝ᵥ v := Finset.sum_nonneg fun i _ => mul_self_nonneg _
  have hcomp : ∀ i, x i * x i ≤ K ^ 2 * n * (v ⬝ᵥ v) := by
    intro i
    have h1 : x i = ∑ j, Y s i j * v j := by
      rw [← hYv]
      rfl
    have h2 : (∑ j, Y s i j * v j) ^ 2 ≤ (∑ j, (Y s i j) ^ 2) * ∑ j, (v j) ^ 2 :=
      Finset.sum_mul_sq_le_sq_mul_sq _ _ _
    have h3 : ∑ j, (Y s i j) ^ 2 ≤ n * K ^ 2 := by
      calc ∑ j, (Y s i j) ^ 2 ≤ ∑ _j : ι, K ^ 2 := by
            refine Finset.sum_le_sum fun j _ => ?_
            rw [← sq_abs]
            exact pow_le_pow_left₀ (abs_nonneg _) (hYK s hs i j) 2
        _ = Fintype.card ι * K ^ 2 := by
            rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
        _ ≤ n * K ^ 2 := by
            have : 0 ≤ K ^ 2 := by positivity
            rw [hn]
            nlinarith
    have h4 : ∑ j, (v j) ^ 2 = v ⬝ᵥ v := by
      simp only [dotProduct, sq]
    rw [← sq, h1]
    rw [h4] at h2
    nlinarith
  have hsum : x ⬝ᵥ x ≤ K ^ 2 * n ^ 2 * (v ⬝ᵥ v) := by
    calc x ⬝ᵥ x = ∑ i, x i * x i := rfl
      _ ≤ ∑ _i : ι, K ^ 2 * n * (v ⬝ᵥ v) := Finset.sum_le_sum fun i _ => hcomp i
      _ = Fintype.card ι * (K ^ 2 * n * (v ⬝ᵥ v)) := by
          rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
      _ ≤ n * (K ^ 2 * n * (v ⬝ᵥ v)) := by
          refine mul_le_mul_of_nonneg_right ?_ (by positivity)
          rw [hn]
          linarith
      _ = K ^ 2 * n ^ 2 * (v ⬝ᵥ v) := by ring
  calc φ / (K ^ 2 * n ^ 2) * (x ⬝ᵥ x) ≤ φ / (K ^ 2 * n ^ 2) * (K ^ 2 * n ^ 2 * (v ⬝ᵥ v)) :=
        mul_le_mul_of_nonneg_left hsum (by positivity)
    _ = φ * (v ⬝ᵥ v) := by field_simp
    _ ≤ x ⬝ᵥ (P Y Z s *ᵥ x) := hform

/-- **Riccati equation** `P' = 2 λ C - P² / 2`, entrywise. -/
theorem P_deriv_of_det (h : RiccatiFlow C Φ lam T Y Z) (hT : 0 < T)
    (hdet : ∀ r ∈ Icc 0 T, (Y r).det ≠ 0) {s : ℝ} (hs : s ∈ Icc 0 T) (i j : ι) :
    HasDerivWithinAt (fun s => P Y Z s i j) (2 * lam * C i j - 1 / 2 * (P Y Z s * P Y Z s) i j)
      (Icc 0 T) s := by
  have huniq : UniqueDiffWithinAt ℝ (Icc 0 T) s := uniqueDiffOn_Icc hT s hs
  have hNd : ∀ a b, DifferentiableWithinAt ℝ (fun s => (Y s)⁻¹ a b) (Icc 0 T) s :=
    differentiableWithinAt_inv (fun a b => (h.Y_deriv s hs a b).differentiableWithinAt) (hdet s hs)
  set N' : Matrix ι ι ℝ := Matrix.of fun a b =>
    derivWithin (fun s => (Y s)⁻¹ a b) (Icc 0 T) s with hN'
  have hN : ∀ a b, HasDerivWithinAt (fun s => (Y s)⁻¹ a b) (N' a b) (Icc 0 T) s := fun a b =>
    (hNd a b).hasDerivWithinAt
  -- differentiate `Y⁻¹ Y = 1`
  have hid : N' * Y s + (Y s)⁻¹ * ((1 / 2 : ℝ) • Z s) = 0 := by
    ext a b
    have hd : HasDerivWithinAt (fun r => ((Y r)⁻¹ * Y r) a b)
        ((N' * Y s + (Y s)⁻¹ * ((1 / 2 : ℝ) • Z s)) a b) (Icc 0 T) s := by
      simp only [Matrix.mul_apply, Matrix.add_apply]
      refine (HasDerivWithinAt.fun_sum fun k _ =>
        (hN a k).mul (h.Y_deriv s hs k b)).congr_deriv ?_
      rw [← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun k _ => ?_
      simp [Matrix.smul_apply]
    have hc : HasDerivWithinAt (fun r => ((Y r)⁻¹ * Y r) a b) 0 (Icc 0 T) s := by
      refine (hasDerivWithinAt_const s (Icc 0 T) ((1 : Matrix ι ι ℝ) a b)).congr ?_ ?_
      · intro r hr
        rw [Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.mpr (hdet r hr))]
      · rw [Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.mpr (hdet s hs))]
    exact huniq.eq_deriv _ hd hc
  have hunit := isUnit_iff_ne_zero.mpr (hdet s hs)
  have hN'eq : N' = -((1 / 2 : ℝ) • ((Y s)⁻¹ * Z s * (Y s)⁻¹)) := by
    have h1 : N' * Y s * (Y s)⁻¹ = N' := by
      rw [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunit, Matrix.mul_one]
    have h2 : N' * Y s = -((Y s)⁻¹ * ((1 / 2 : ℝ) • Z s)) := eq_neg_of_add_eq_zero_left hid
    rw [← h1, h2, Matrix.neg_mul, Matrix.mul_smul, Matrix.smul_mul]
  have hd : HasDerivWithinAt (fun r => P Y Z r i j)
      (((2 * lam) • (C * Y s) * (Y s)⁻¹ + Z s * N') i j) (Icc 0 T) s := by
    simp only [P, Matrix.mul_apply, Matrix.add_apply]
    refine (HasDerivWithinAt.fun_sum fun k _ =>
      (h.Z_deriv s hs i k).mul (hN k j)).congr_deriv ?_
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    simp [Matrix.smul_apply]
  refine hd.congr_deriv ?_
  have h3 : (2 * lam) • (C * Y s) * (Y s)⁻¹ = (2 * lam) • C := by
    rw [Matrix.smul_mul, Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunit, Matrix.mul_one]
  have h4 : Z s * N' = -((1 / 2 : ℝ) • (P Y Z s * P Y Z s)) := by
    rw [hN'eq, Matrix.mul_neg, Matrix.mul_smul]
    simp only [P, Matrix.mul_assoc]
  rw [h3, h4]
  simp only [Matrix.add_apply, Matrix.neg_apply, Matrix.smul_apply, smul_eq_mul]
  ring

theorem P_deriv (h : RiccatiFlow C Φ lam T Y Z) (hT : 0 < T) (hlam : 0 ≤ lam)
    (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v)) (hφ : 0 < φ) (hΦ : ∀ v, φ * (v ⬝ᵥ v) ≤ v ⬝ᵥ (Φ *ᵥ v)) {s : ℝ}
    (hs : s ∈ Icc 0 T) (i j : ι) :
    HasDerivWithinAt (fun s => P Y Z s i j) (2 * lam * C i j - 1 / 2 * (P Y Z s * P Y Z s) i j)
      (Icc 0 T) s :=
  h.P_deriv_of_det hT (fun _ hr => h.det_ne_zero hlam hC hφ hΦ hr) hs i j

end RiccatiFlow

end ExecutionProofs
