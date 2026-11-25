import ExecutionProofs.GridExpansion
import ExecutionProofs.Riccati

/-!
# Polynomials times a Gaussian

`Dop` is the action of `∂ₐ` on the polynomial factor of `p · gauss`. `PolyFam` is the class of
polynomial families with continuous coefficients on `[0, T]`; it is closed under all operations
used below and has uniform polynomial growth.
-/

open Set Matrix MvPolynomial

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The Gaussian factor `exp (-(x · P x) / (8 ν) - θ)`. -/
noncomputable def gauss (ν : ℝ) (P : Matrix ι ι ℝ) (θ : ℝ) (x : ι → ℝ) : ℝ :=
  Real.exp (-(x ⬝ᵥ (P *ᵥ x)) / (8 * ν) - θ)

omit [DecidableEq ι] in
lemma gauss_pos (ν : ℝ) (P : Matrix ι ι ℝ) (θ : ℝ) (x : ι → ℝ) : 0 < gauss ν P θ x :=
  Real.exp_pos _

/-- `∂ₐ (p · gauss) = (Dop a p) · gauss`. -/
noncomputable def Dop (ν : ℝ) (P : Matrix ι ι ℝ) (a : ι) (p : MvPolynomial ι ℝ) :
    MvPolynomial ι ℝ :=
  pderiv a p - ∑ b, (P a b / (4 * ν)) • (X b * p)

lemma single_dot (y : ι → ℝ) (a : ι) (t : ℝ) :
    (Pi.single a t : ι → ℝ) ⬝ᵥ y = t * y a := by
  simp [dotProduct, Pi.single_apply]

/-- The quadratic form along a coordinate line. -/
lemma form_line {P : Matrix ι ι ℝ} (hP : P.IsSymm) (x : ι → ℝ) (a : ι) (t : ℝ) :
    (x + Pi.single a t) ⬝ᵥ (P *ᵥ (x + Pi.single a t))
      = x ⬝ᵥ (P *ᵥ x) + 2 * t * (P *ᵥ x) a + t ^ 2 * P a a := by
  have h1 : P *ᵥ Pi.single a t = fun i => P i a * t := by
    ext i
    simp [Matrix.mulVec, dotProduct, Pi.single_apply]
  have h2 : x ⬝ᵥ (fun i => P i a * t) = t * (P *ᵥ x) a := by
    simp only [dotProduct, Matrix.mulVec, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← hP.apply a i]
    ring
  rw [Matrix.mulVec_add, add_dotProduct, dotProduct_add, dotProduct_add, h1, h2,
    single_dot, single_dot]
  ring

omit [Fintype ι] in
/-- Evaluation of a polynomial along a coordinate line. -/
lemma hasDerivAt_eval_line (p : MvPolynomial ι ℝ) (x : ι → ℝ) (a : ι) (t : ℝ) :
    HasDerivAt (fun t : ℝ => eval (x + Pi.single a t) p) (eval (x + Pi.single a t) (pderiv a p)) t := by
  induction p using MvPolynomial.induction_on with
  | C r => simpa using hasDerivAt_const t r
  | add p q hp hq => simpa using hp.fun_add hq
  | mul_X p i hp =>
    have hX : HasDerivAt (fun t : ℝ => eval (x + Pi.single a t) (X i : MvPolynomial ι ℝ))
        (eval (x + Pi.single a t) (pderiv a (X i : MvPolynomial ι ℝ))) t := by
      simp only [eval_X, Pi.add_apply, pderiv_X]
      by_cases h : i = a
      · subst h
        simpa using (hasDerivAt_id t).const_add (x i)
      · simpa [Pi.single_eq_of_ne h, Pi.single_eq_of_ne (Ne.symm h)] using
          hasDerivAt_const t (x i)
    have h := hp.fun_mul hX
    simp only [map_mul, pderiv_mul, map_add]
    simpa [mul_comm, add_comm] using h

/-- The Gaussian along a coordinate line. -/
lemma hasDerivAt_gauss_line {ν : ℝ} {P : Matrix ι ι ℝ} (hP : P.IsSymm) (θ : ℝ) (x : ι → ℝ) (a : ι)
    (t : ℝ) : HasDerivAt (fun t : ℝ => gauss ν P θ (x + Pi.single a t))
      (-((P *ᵥ (x + Pi.single a t)) a) / (4 * ν) * gauss ν P θ (x + Pi.single a t)) t := by
  unfold gauss
  simp only [form_line hP]
  have hlin : (P *ᵥ (x + Pi.single a t)) a = (P *ᵥ x) a + P a a * t := by
    rw [Matrix.mulVec_add, Pi.add_apply]
    congr 1
    simp [Matrix.mulVec, dotProduct, Pi.single_apply]
  rw [hlin]
  have hq : HasDerivAt (fun t : ℝ => x ⬝ᵥ (P *ᵥ x) + 2 * t * (P *ᵥ x) a + t ^ 2 * P a a)
      (2 * (P *ᵥ x) a + 2 * t * P a a) t := by
    have h1 : HasDerivAt (fun t : ℝ => 2 * t * (P *ᵥ x) a) (2 * (P *ᵥ x) a) t := by
      simpa using ((hasDerivAt_id t).const_mul 2).mul_const ((P *ᵥ x) a)
    have h2 : HasDerivAt (fun t : ℝ => t ^ 2 * P a a) (2 * t * P a a) t := by
      simpa using (hasDerivAt_pow 2 t).mul_const (P a a)
    exact (h1.const_add _).add h2
  refine ((hq.fun_neg.div_const (8 * ν)).sub_const θ).exp.congr_deriv ?_
  ring

/-- One derivative of `p · gauss` along a coordinate line. -/
lemma hasDerivAt_line {ν : ℝ} {P : Matrix ι ι ℝ} (hP : P.IsSymm) (θ : ℝ) (p : MvPolynomial ι ℝ)
    (x : ι → ℝ) (a : ι) (t : ℝ) :
    HasDerivAt (fun t : ℝ => eval (x + Pi.single a t) p * gauss ν P θ (x + Pi.single a t))
      (eval (x + Pi.single a t) (Dop ν P a p) * gauss ν P θ (x + Pi.single a t)) t := by
  refine ((hasDerivAt_eval_line p x a t).fun_mul
    (hasDerivAt_gauss_line hP θ x a t)).congr_deriv ?_
  have h : eval (x + Pi.single a t) (∑ b, (P a b / (4 * ν)) • (X b * p))
      = (P *ᵥ (x + Pi.single a t)) a / (4 * ν) * eval (x + Pi.single a t) p := by
    simp only [map_sum, smul_eval, map_mul, eval_X, Matrix.mulVec, dotProduct, Finset.sum_div,
      Finset.sum_mul]
    refine Finset.sum_congr rfl fun b _ => ?_
    ring
  rw [Dop, map_sub, h]
  ring

/-- All derivatives of `p · gauss` along a coordinate line. -/
lemma iteratedDeriv_polyGauss {ν : ℝ} {P : Matrix ι ι ℝ} (hP : P.IsSymm) (θ : ℝ)
    (p : MvPolynomial ι ℝ) (x : ι → ℝ) (a : ι) (m : ℕ) :
    iteratedDeriv m (fun t : ℝ => eval (x + Pi.single a t) p * gauss ν P θ (x + Pi.single a t))
      = fun t => eval (x + Pi.single a t) ((Dop ν P a)^[m] p) * gauss ν P θ (x + Pi.single a t) := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [iteratedDeriv_succ, ih]
    funext t
    rw [Function.iterate_succ_apply']
    exact (hasDerivAt_line hP θ _ x a t).deriv

lemma pd_polyGauss {ν : ℝ} {P : Matrix ι ι ℝ} (hP : P.IsSymm) (θ : ℝ) (p : MvPolynomial ι ℝ)
    (x : ι → ℝ) (a : ι) (m : ℕ) :
    pd a m (fun x => eval x p * gauss ν P θ x) x
      = eval x ((Dop ν P a)^[m] p) * gauss ν P θ x := by
  unfold pd
  rw [iteratedDeriv_polyGauss hP]
  simp

lemma contDiff_polyGauss {ν : ℝ} {P : Matrix ι ι ℝ} (hP : P.IsSymm) (θ : ℝ)
    (p : MvPolynomial ι ℝ) (x : ι → ℝ) (a : ι) (n : ℕ) :
    ContDiff ℝ n (fun t : ℝ => eval (x + Pi.single a t) p * gauss ν P θ (x + Pi.single a t)) := by
  refine contDiff_of_differentiable_iteratedDeriv fun m _ => ?_
  rw [iteratedDeriv_polyGauss hP]
  exact fun t => (hasDerivAt_line hP θ _ x a t).differentiableAt

/-! ### Polynomial families with continuous coefficients -/

/-- Finite sums of fixed polynomials with coefficients continuous on `[0, T]`. -/
inductive PolyFam (T : ℝ) : (ℝ → MvPolynomial ι ℝ) → Prop
  | smul (α : ℝ → ℝ) (hα : ContinuousOn α (Icc 0 T)) (m : MvPolynomial ι ℝ) :
      PolyFam T fun s => α s • m
  | add {q r : ℝ → MvPolynomial ι ℝ} (hq : PolyFam T q) (hr : PolyFam T r) :
      PolyFam T fun s => q s + r s

namespace PolyFam

variable {T : ℝ} {q r : ℝ → MvPolynomial ι ℝ}

omit [Fintype ι] [DecidableEq ι] in
lemma const (m : MvPolynomial ι ℝ) : PolyFam T fun _ => m := by
  simpa using PolyFam.smul (fun _ => (1 : ℝ)) continuousOn_const m

omit [Fintype ι] [DecidableEq ι] in
lemma map (ℓ : MvPolynomial ι ℝ →ₗ[ℝ] MvPolynomial ι ℝ) (hq : PolyFam T q) :
    PolyFam T fun s => ℓ (q s) := by
  induction hq with
  | smul α hα m => simpa using PolyFam.smul α hα (ℓ m)
  | add _ _ ih₁ ih₂ => simpa using ih₁.add ih₂

omit [Fintype ι] [DecidableEq ι] in
lemma const_smul {β : ℝ → ℝ} (hβ : ContinuousOn β (Icc 0 T)) (hq : PolyFam T q) :
    PolyFam T fun s => β s • q s := by
  induction hq with
  | smul α hα m => simpa [smul_smul] using PolyFam.smul (fun s => β s * α s) (hβ.mul hα) m
  | add _ _ ih₁ ih₂ => simpa [smul_add] using ih₁.add ih₂

omit [Fintype ι] [DecidableEq ι] in
lemma sum {J : Type*} (S : Finset J) {f : J → ℝ → MvPolynomial ι ℝ}
    (hf : ∀ j ∈ S, PolyFam T (f j)) : PolyFam T fun s => ∑ j ∈ S, f j s := by
  classical
  induction S using Finset.induction_on with
  | empty => simpa using const (T := T) (0 : MvPolynomial ι ℝ)
  | insert a S ha ih =>
    simp only [Finset.sum_insert ha]
    exact (hf a (Finset.mem_insert_self a S)).add
      (ih fun j hj => hf j (Finset.mem_insert_of_mem hj))

omit [DecidableEq ι] in
lemma dop {ν : ℝ} {P : ℝ → Matrix ι ι ℝ} (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T))
    (a : ι) (hq : PolyFam T q) : PolyFam T fun s => Dop ν (P s) a (q s) := by
  have h1 : PolyFam T fun s => (pderiv a).toLinearMap (q s) := hq.map (pderiv a).toLinearMap
  have h2 : PolyFam T fun s => ∑ b, (-(P s a b / (4 * ν))) • (LinearMap.mulLeft ℝ (X b) (q s)) :=
    sum _ fun b _ => (hq.map (LinearMap.mulLeft ℝ (X b))).const_smul (((hP a b).div_const _).neg)
  have h := h1.add h2
  refine (congrArg (PolyFam T) (funext fun s => ?_)).mp h
  simp [Dop, sub_eq_add_neg, Finset.sum_neg_distrib]

omit [DecidableEq ι] in
lemma dop_iterate {ν : ℝ} {P : ℝ → Matrix ι ι ℝ}
    (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T)) (a : ι) (hq : PolyFam T q) (m : ℕ) :
    PolyFam T fun s => (Dop ν (P s) a)^[m] (q s) := by
  induction m with
  | zero => simpa using hq
  | succ m ih =>
    simp only [Function.iterate_succ_apply']
    exact ih.dop hP a

omit [DecidableEq ι] in
/-- A fixed polynomial grows at most polynomially. -/
lemma growth_const (m : MvPolynomial ι ℝ) :
    ∃ B : ℝ, ∃ k : ℕ, 0 ≤ B ∧ ∀ x : ι → ℝ, |eval x m| ≤ B * (1 + x ⬝ᵥ x) ^ k := by
  have hbase : ∀ x : ι → ℝ, 1 ≤ 1 + x ⬝ᵥ x := fun x => by
    have : 0 ≤ x ⬝ᵥ x := Finset.sum_nonneg fun i _ => mul_self_nonneg _
    linarith
  induction m using MvPolynomial.induction_on with
  | C r => exact ⟨|r|, 0, abs_nonneg _, fun x => by simp⟩
  | add p q hp hq =>
    obtain ⟨B₁, k₁, h₁, hb₁⟩ := hp
    obtain ⟨B₂, k₂, h₂, hb₂⟩ := hq
    refine ⟨B₁ + B₂, max k₁ k₂, add_nonneg h₁ h₂, fun x => ?_⟩
    rw [map_add, add_mul]
    refine (abs_add_le _ _).trans (add_le_add ((hb₁ x).trans ?_) ((hb₂ x).trans ?_))
    · exact mul_le_mul_of_nonneg_left (pow_le_pow_right₀ (hbase x) (le_max_left _ _)) h₁
    · exact mul_le_mul_of_nonneg_left (pow_le_pow_right₀ (hbase x) (le_max_right _ _)) h₂
  | mul_X p i hp =>
    obtain ⟨B, k, hB, hb⟩ := hp
    refine ⟨B, k + 1, hB, fun x => ?_⟩
    have hi : |x i| ≤ 1 + x ⬝ᵥ x := by
      have h1 : x i * x i ≤ x ⬝ᵥ x :=
        Finset.single_le_sum (f := fun j => x j * x j) (fun j _ => mul_self_nonneg _)
          (Finset.mem_univ i)
      rw [abs_le]
      constructor <;> nlinarith [sq_nonneg (x i - 1), sq_nonneg (x i + 1)]
    rw [map_mul, eval_X, abs_mul, pow_succ, ← mul_assoc]
    exact mul_le_mul (hb x) hi (abs_nonneg _)
      (mul_nonneg hB (pow_nonneg (zero_le_one.trans (hbase x)) _))

omit [DecidableEq ι] in
/-- Uniform polynomial growth of a family on `[0, T]`. -/
lemma growth (hq : PolyFam T q) : ∃ B : ℝ, ∃ k : ℕ, 0 ≤ B ∧
    ∀ s ∈ Icc 0 T, ∀ x : ι → ℝ, |eval x (q s)| ≤ B * (1 + x ⬝ᵥ x) ^ k := by
  have hbase : ∀ x : ι → ℝ, 1 ≤ 1 + x ⬝ᵥ x := fun x => by
    have : 0 ≤ x ⬝ᵥ x := Finset.sum_nonneg fun i _ => mul_self_nonneg _
    linarith
  induction hq with
  | smul α hα m =>
    obtain ⟨A, hA⟩ := isCompact_Icc.exists_bound_of_continuousOn hα
    obtain ⟨B, k, hB, hb⟩ := growth_const m
    refine ⟨|A| * B, k, by positivity, fun s hs x => ?_⟩
    rw [smul_eval, abs_mul, mul_assoc]
    exact mul_le_mul ((hA s hs).trans (le_abs_self A)) (hb x) (abs_nonneg _) (abs_nonneg _)
  | add _ _ ih₁ ih₂ =>
    obtain ⟨B₁, k₁, h₁, hb₁⟩ := ih₁
    obtain ⟨B₂, k₂, h₂, hb₂⟩ := ih₂
    refine ⟨B₁ + B₂, max k₁ k₂, add_nonneg h₁ h₂, fun s hs x => ?_⟩
    rw [map_add, add_mul]
    refine (abs_add_le _ _).trans (add_le_add ((hb₁ s hs x).trans ?_) ((hb₂ s hs x).trans ?_))
    · exact mul_le_mul_of_nonneg_left (pow_le_pow_right₀ (hbase x) (le_max_left _ _)) h₁
    · exact mul_le_mul_of_nonneg_left (pow_le_pow_right₀ (hbase x) (le_max_right _ _)) h₂

end PolyFam

/-- `(1 + r)^k` is dominated by `exp (κ r)`. -/
lemma pow_mul_exp_le {κ : ℝ} (hκ : 0 < κ) (k : ℕ) :
    ∃ C, 0 ≤ C ∧ ∀ r, 0 ≤ r → (1 + r) ^ k * Real.exp (-(κ * r)) ≤ C := by
  set η : ℝ := κ / (k + 1) with hη
  have hη0 : 0 < η := by positivity
  refine ⟨(1 + 1 / η) ^ k, by positivity, fun r hr => ?_⟩
  have h1 : 1 + r ≤ (1 + 1 / η) * Real.exp (η * r) := by
    have he : 1 + η * r ≤ Real.exp (η * r) := by
      have := Real.add_one_le_exp (η * r)
      linarith
    have h2 : 1 + r ≤ (1 + 1 / η) * (1 + η * r) := by
      have : (1 + 1 / η) * (1 + η * r) = 1 + r + (1 / η + η * r) := by
        field_simp
        ring
      have h3 : 0 ≤ 1 / η + η * r := by positivity
      linarith
    exact h2.trans (mul_le_mul_of_nonneg_left he (by positivity))
  have h3 : (1 + r) ^ k ≤ (1 + 1 / η) ^ k * Real.exp (κ * r) := by
    calc (1 + r) ^ k ≤ ((1 + 1 / η) * Real.exp (η * r)) ^ k :=
          pow_le_pow_left₀ (by linarith) h1 k
      _ = (1 + 1 / η) ^ k * Real.exp (k * (η * r)) := by
          rw [mul_pow, ← Real.exp_nat_mul]
      _ ≤ (1 + 1 / η) ^ k * Real.exp (κ * r) := by
          refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (by positivity)
          have h4 : (k : ℝ) * η ≤ κ := by
            rw [hη, mul_div_assoc']
            rw [div_le_iff₀ (by positivity)]
            nlinarith
          nlinarith
  have hexp : Real.exp (κ * r) * Real.exp (-(κ * r)) = 1 := by
    rw [← Real.exp_add]
    simp
  calc (1 + r) ^ k * Real.exp (-(κ * r))
      ≤ (1 + 1 / η) ^ k * Real.exp (κ * r) * Real.exp (-(κ * r)) :=
        mul_le_mul_of_nonneg_right h3 (Real.exp_pos _).le
    _ = (1 + 1 / η) ^ k := by
        rw [mul_assoc, hexp, mul_one]

omit [DecidableEq ι] in
/-- A polynomial family times a uniformly coercive Gaussian is bounded. -/
lemma PolyFam.bounded {T ν : ℝ} (hν : 0 < ν) {q : ℝ → MvPolynomial ι ℝ} (hq : PolyFam T q)
    {P : ℝ → Matrix ι ι ℝ} {θ : ℝ → ℝ} (hθ : ContinuousOn θ (Icc 0 T)) {m : ℝ} (hm : 0 < m)
    (hP : ∀ s ∈ Icc 0 T, ∀ x, m * (x ⬝ᵥ x) ≤ x ⬝ᵥ (P s *ᵥ x)) :
    ∃ N, 0 ≤ N ∧ ∀ s ∈ Icc 0 T, ∀ x, |eval x (q s) * gauss ν (P s) (θ s) x| ≤ N := by
  obtain ⟨B, k, hB, hb⟩ := hq.growth
  obtain ⟨Θ, hΘ⟩ := isCompact_Icc.exists_bound_of_continuousOn hθ
  obtain ⟨C, hC, hc⟩ := pow_mul_exp_le (κ := m / (8 * ν)) (by positivity) k
  refine ⟨B * C * Real.exp Θ, by positivity, fun s hs x => ?_⟩
  have hr : 0 ≤ x ⬝ᵥ x := Finset.sum_nonneg fun i _ => mul_self_nonneg _
  have hg : gauss ν (P s) (θ s) x ≤ Real.exp (-(m / (8 * ν) * (x ⬝ᵥ x))) * Real.exp Θ := by
    unfold gauss
    rw [← Real.exp_add]
    refine Real.exp_le_exp.mpr ?_
    have h1 : -(θ s) ≤ Θ := (neg_le_abs _).trans (hΘ s hs)
    have h2 : m / (8 * ν) * (x ⬝ᵥ x) ≤ x ⬝ᵥ (P s *ᵥ x) / (8 * ν) := by
      rw [div_mul_eq_mul_div]
      exact div_le_div_of_nonneg_right (hP s hs x) (by positivity)
    rw [neg_div]
    linarith
  rw [abs_mul, abs_of_pos (gauss_pos _ _ _ _)]
  calc |eval x (q s)| * gauss ν (P s) (θ s) x
      ≤ B * (1 + x ⬝ᵥ x) ^ k * (Real.exp (-(m / (8 * ν) * (x ⬝ᵥ x))) * Real.exp Θ) :=
        mul_le_mul (hb s hs x) hg (gauss_pos _ _ _ _).le (by positivity)
    _ = B * ((1 + x ⬝ᵥ x) ^ k * Real.exp (-(m / (8 * ν) * (x ⬝ᵥ x)))) * Real.exp Θ := by ring
    _ ≤ B * C * Real.exp Θ := by
        have := hc (x ⬝ᵥ x) hr
        gcongr

end ExecutionProofs
