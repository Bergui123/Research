import ExecutionProofs.Leading
import ExecutionProofs.PSD

/-!
# Sign structure of the leading coefficient (paper, Theorem 3 (iv))

`p₁ = ν⁻³ A + ν⁻² B + ν⁻¹ D` with `A`, `B`, `D` independent of `ν`, `A` a quartic form that is
positive away from the origin, `B` a quadratic form that is negative away from the origin, and `D`
constant in space.
-/

open Set Matrix MvPolynomial

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-! ### Derivatives along characteristics -/

omit [DecidableEq ι] in
lemma hasDeriv_eval_curve {T s : ℝ} {y : ℝ → ι → ℝ} {y' : ι → ℝ}
    (hy : ∀ a, HasDerivWithinAt (fun s => y s a) (y' a) (Icc 0 T) s) (m : MvPolynomial ι ℝ) :
    HasDerivWithinAt (fun s => eval (y s) m) (∑ a, y' a * eval (y s) (pderiv a m)) (Icc 0 T) s := by
  classical
  induction m using MvPolynomial.induction_on with
  | C r => simpa using hasDerivWithinAt_const s (Icc 0 T) r
  | add p q hp hq =>
    simp only [map_add, mul_add, Finset.sum_add_distrib]
    exact hp.fun_add hq
  | mul_X p i hp =>
    simp only [map_mul, eval_X]
    refine (hp.fun_mul (hy i)).congr_deriv ?_
    have h1 : ∀ a, y' a * eval (y s) (pderiv a (p * X i))
        = y' a * eval (y s) (pderiv a p) * y s i + y' a * (eval (y s) p * (if a = i then 1 else 0)) := by
      intro a
      rw [pderiv_mul, pderiv_X, map_add, map_mul, map_mul, eval_X]
      by_cases h : a = i
      · subst h
        simp
        ring
      · simp [h]
        ring
    simp only [h1, Finset.sum_add_distrib, ← Finset.sum_mul]
    congr 1
    simp [Finset.sum_ite_eq', mul_comm]

omit [DecidableEq ι] in
lemma DiffFam.hasDeriv_along {T s : ℝ} (hs : s ∈ Icc 0 T) {p p' : ℝ → MvPolynomial ι ℝ}
    (h : DiffFam T p p') {y : ℝ → ι → ℝ} {y' : ι → ℝ}
    (hy : ∀ a, HasDerivWithinAt (fun s => y s a) (y' a) (Icc 0 T) s) :
    HasDerivWithinAt (fun s => eval (y s) (p s))
      (eval (y s) (p' s) + ∑ a, y' a * eval (y s) (pderiv a (p s))) (Icc 0 T) s := by
  obtain ⟨k, m, α, α', hp, hp', hα⟩ := h
  rw [hp' s hs, hp s]
  simp only [hp, map_sum, smul_eval]
  refine (HasDerivWithinAt.fun_sum fun i _ =>
    (hα i s hs).fun_mul (hasDeriv_eval_curve hy (m i))).congr_deriv ?_
  rw [Finset.sum_add_distrib]
  congr 1
  simp only [Derivation.map_smul, smul_eval, Finset.mul_sum]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun i _ => by ring

omit [Fintype ι] [DecidableEq ι] in
lemma DiffFam.map {T : ℝ} {p p' : ℝ → MvPolynomial ι ℝ} (h : DiffFam T p p')
    (ℓ : MvPolynomial ι ℝ →ₗ[ℝ] MvPolynomial ι ℝ) :
    DiffFam T (fun s => ℓ (p s)) fun s => ℓ (p' s) := by
  obtain ⟨k, m, α, α', hp, hp', hα⟩ := h
  refine ⟨k, fun i => ℓ (m i), α, α', fun s => ?_, fun s hs => ?_, hα⟩
  · simp [hp]
  · simp [hp' s hs]

omit [DecidableEq ι] in
lemma eval_polyGen (ν : ℝ) (P : Matrix ι ι ℝ) (x : ι → ℝ) (p : MvPolynomial ι ℝ) :
    eval x (polyGen ν P p) = ν * ∑ a, eval x (pderiv a (pderiv a p))
      - ∑ a, (P *ᵥ x) a / 2 * eval x (pderiv a p) := by
  rw [polyGen, map_sub, smul_eval, map_sum, map_sum]
  congr 1
  refine Finset.sum_congr rfl fun a _ => ?_
  simp only [map_sum, smul_eval, map_mul, eval_X, Matrix.mulVec, dotProduct, Finset.sum_div,
    Finset.sum_mul]
  exact Finset.sum_congr rfl fun b _ => by ring

/-- Along `y = Y ξ` the transport operator is cancelled. -/
lemma transport_along {C Φ : Matrix ι ι ℝ} {lam T s : ℝ} {Y Z : ℝ → Matrix ι ι ℝ}
    (hYZ : RiccatiFlow C Φ lam T Y Z) (hs : s ∈ Icc 0 T) (hdet : (Y s).det ≠ 0)
    {p p' f : ℝ → MvPolynomial ι ℝ} (hp : DiffFam T p p')
    (hp' : p' s = polyGen 0 (RiccatiFlow.P Y Z s) (p s) + f s) (ξ : ι → ℝ) :
    HasDerivWithinAt (fun s => eval (Y s *ᵥ ξ) (p s)) (eval (Y s *ᵥ ξ) (f s)) (Icc 0 T) s := by
  refine (hp.hasDeriv_along hs (fun a => hYZ.hasDeriv_Y_mulVec hs ξ a)).congr_deriv ?_
  rw [hp', map_add, eval_polyGen, zero_mul, zero_sub]
  have h : RiccatiFlow.P Y Z s *ᵥ (Y s *ᵥ ξ) = Z s *ᵥ ξ := by
    rw [Matrix.mulVec_mulVec, RiccatiFlow.P_mul_Y hdet]
  rw [h]
  have h2 : ∑ a, 1 / 2 * (Z s *ᵥ ξ) a * eval (Y s *ᵥ ξ) (pderiv a (p s))
      = ∑ a, (Z s *ᵥ ξ) a / 2 * eval (Y s *ᵥ ξ) (pderiv a (p s)) :=
    Finset.sum_congr rfl fun a _ => by ring
  rw [h2]
  ring

/-! ### The three sources -/

/-- Linear form `x ↦ (P x)_a`. -/
noncomputable def lin (P : Matrix ι ι ℝ) (a : ι) : MvPolynomial ι ℝ := ∑ b, P a b • X b

omit [DecidableEq ι] in
lemma eval_lin (P : Matrix ι ι ℝ) (a : ι) (x : ι → ℝ) : eval x (lin P a) = (P *ᵥ x) a := by
  simp [lin, Matrix.mulVec, dotProduct, smul_eval]

omit [DecidableEq ι] in
lemma pderiv_lin (P : Matrix ι ι ℝ) (a i : ι) : pderiv i (lin P a) = MvPolynomial.C (P a i) := by
  classical
  simp only [lin, map_sum, Derivation.map_smul, pderiv_X]
  rw [Finset.sum_eq_single i]
  · rw [Pi.single_eq_same, MvPolynomial.smul_eq_C_mul, mul_one]
  · intro b _ hb
    rw [Pi.single_eq_of_ne hb, smul_zero]
  · simp

/-- Source of the quartic part. -/
noncomputable def srcA (P : Matrix ι ι ℝ) : MvPolynomial ι ℝ :=
  ∑ a, (1 / 3072 : ℝ) • lin P a ^ 4

/-- Source of the quadratic part. -/
noncomputable def srcB (P : Matrix ι ι ℝ) (A : MvPolynomial ι ℝ) : MvPolynomial ι ℝ :=
  ∑ a, pderiv a (pderiv a A) - ∑ a, (P a a / 128) • lin P a ^ 2

/-- Source of the constant part. -/
noncomputable def srcD (P : Matrix ι ι ℝ) (B : MvPolynomial ι ℝ) : MvPolynomial ι ℝ :=
  ∑ a, pderiv a (pderiv a B) + ∑ a, (P a a ^ 2 / 64) • (1 : MvPolynomial ι ℝ)

omit [DecidableEq ι] in
lemma eval_srcA (P : Matrix ι ι ℝ) (x : ι → ℝ) :
    eval x (srcA P) = 1 / 3072 * ∑ a, ((P *ᵥ x) a) ^ 4 := by
  simp only [srcA, map_sum, smul_eval, map_pow, eval_lin, Finset.mul_sum]

omit [DecidableEq ι] in
lemma eval_srcB (P : Matrix ι ι ℝ) (A : MvPolynomial ι ℝ) (x : ι → ℝ) :
    eval x (srcB P A) = ∑ a, eval x (pderiv a (pderiv a A))
      - 1 / 128 * ∑ a, P a a * ((P *ᵥ x) a) ^ 2 := by
  simp only [srcB, map_sub, map_sum, smul_eval, map_pow, eval_lin, Finset.mul_sum]
  congr 1
  exact Finset.sum_congr rfl fun a _ => by ring

omit [DecidableEq ι] in
lemma eval_srcD (P : Matrix ι ι ℝ) (B : MvPolynomial ι ℝ) (x : ι → ℝ) :
    eval x (srcD P B) = ∑ a, eval x (pderiv a (pderiv a B)) + 1 / 64 * ∑ a, P a a ^ 2 := by
  simp only [srcD, map_add, map_sum, smul_eval, map_one, Finset.mul_sum]
  congr 1
  exact Finset.sum_congr rfl fun a _ => by ring

omit [Fintype ι] [DecidableEq ι] in
lemma PolyFam.mul {T : ℝ} {q r : ℝ → MvPolynomial ι ℝ} (hq : PolyFam T q) (hr : PolyFam T r) :
    PolyFam T fun s => q s * r s := by
  induction hq with
  | smul α hα m =>
    have h := (hr.map (LinearMap.mulLeft ℝ m)).const_smul hα
    simpa [smul_mul_assoc] using h
  | add _ _ ih₁ ih₂ => simpa [add_mul] using ih₁.add ih₂

omit [Fintype ι] [DecidableEq ι] in
lemma PolyFam.pow {T : ℝ} {q : ℝ → MvPolynomial ι ℝ} (hq : PolyFam T q) (n : ℕ) :
    PolyFam T fun s => q s ^ n := by
  induction n with
  | zero => simpa using PolyFam.const (T := T) (1 : MvPolynomial ι ℝ)
  | succ n ih =>
    simp only [pow_succ]
    exact ih.mul hq

omit [DecidableEq ι] in
lemma polyFam_lin {T : ℝ} {P : ℝ → Matrix ι ι ℝ}
    (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T)) (a : ι) :
    PolyFam T fun s => lin (P s) a :=
  PolyFam.sum _ fun b _ => PolyFam.smul _ (hP a b) _

/-! ### Degrees -/

omit [Fintype ι] [DecidableEq ι] in
lemma pderiv_mem_lower {n : ℕ} (a : ι) {p : MvPolynomial ι ℝ}
    (hp : p ∈ restrictTotalDegree ι ℝ (n + 1)) : pderiv a p ∈ restrictTotalDegree ι ℝ n := by
  rw [mem_restrictTotalDegree] at hp
  rw [p.as_sum, map_sum]
  refine Submodule.sum_mem _ fun d hd => ?_
  have hdn : (d.sum fun _ e => e) ≤ n + 1 := (le_totalDegree hd).trans hp
  rw [pderiv_monomial, mem_restrictTotalDegree]
  by_cases ha : d a = 0
  · simp [ha]
  · refine (totalDegree_monomial_le _ _).trans ?_
    change ((d - Finsupp.single a 1).sum fun _ e => e) ≤ n
    have hle : Finsupp.single a 1 ≤ d := by
      rw [Finsupp.single_le_iff]
      omega
    have hsplit : (d.sum fun _ e => e)
        = ((d - Finsupp.single a 1).sum fun _ e => e) + 1 := by
      conv_lhs => rw [← tsub_add_cancel_of_le hle]
      rw [Finsupp.sum_add_index' (fun _ => rfl) (fun _ _ _ => rfl),
        Finsupp.sum_single_index rfl]
    omega

omit [Fintype ι] [DecidableEq ι] in
lemma pderiv_eq_zero_of_mem_zero (a : ι) {p : MvPolynomial ι ℝ}
    (hp : p ∈ restrictTotalDegree ι ℝ 0) : pderiv a p = 0 := by
  rw [mem_restrictTotalDegree, Nat.le_zero] at hp
  rw [totalDegree_eq_zero_iff_eq_C.mp hp, pderiv_C]

omit [DecidableEq ι] in
lemma lin_pow_mem (P : Matrix ι ι ℝ) (a : ι) (n : ℕ) :
    lin P a ^ n ∈ restrictTotalDegree ι ℝ n := by
  have h1 : lin P a ∈ restrictTotalDegree ι ℝ 1 := by
    refine Submodule.sum_mem _ fun b _ => Submodule.smul_mem _ _ ?_
    rw [mem_restrictTotalDegree, totalDegree_X]
  rw [mem_restrictTotalDegree] at h1 ⊢
  exact (totalDegree_pow _ n).trans (by nlinarith)

/-- The three families. -/
structure ABD (T : ℝ) (P : ℝ → Matrix ι ι ℝ) (A A' B B' D D' : ℝ → MvPolynomial ι ℝ) : Prop where
  A_zero : A 0 = 0
  B_zero : B 0 = 0
  D_zero : D 0 = 0
  A_fam : PolyFam T A
  B_fam : PolyFam T B
  D_fam : PolyFam T D
  D_mem : ∀ s, D s ∈ restrictTotalDegree ι ℝ 0
  A_diff : DiffFam T A A'
  B_diff : DiffFam T B B'
  D_diff : DiffFam T D D'
  A_eq : ∀ s ∈ Icc 0 T, A' s = polyGen 0 (P s) (A s) + srcA (P s)
  B_eq : ∀ s ∈ Icc 0 T, B' s = polyGen 0 (P s) (B s) + srcB (P s) (A s)
  D_eq : ∀ s ∈ Icc 0 T, D' s = polyGen 0 (P s) (D s) + srcD (P s) (B s)

theorem exists_ABD {T : ℝ} (hT : 0 ≤ T) {P : ℝ → Matrix ι ι ℝ}
    (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T)) :
    ∃ A A' B B' D D', ABD T P A A' B B' D D' := by
  have hfA : PolyFam T fun s => srcA (P s) :=
    PolyFam.sum _ fun a _ => ((polyFam_lin hP a).pow 4).const_smul continuousOn_const
  have hdA : ∀ s, srcA (P s) ∈ restrictTotalDegree ι ℝ 4 := fun s =>
    Submodule.sum_mem _ fun a _ => Submodule.smul_mem _ _ (lin_pow_mem _ a 4)
  obtain ⟨A, A', hA0, hAf, hAm, hAd, hAe⟩ := exists_poly_solution_diff (ν := 0) hT hP hfA hdA
  have hlap : ∀ {q : ℝ → MvPolynomial ι ℝ}, PolyFam T q →
      PolyFam T fun s => ∑ a, pderiv a (pderiv a (q s)) := fun hq =>
    PolyFam.sum _ fun a _ => (hq.map (pderiv a).toLinearMap).map (pderiv a).toLinearMap
  have hfB : PolyFam T fun s => srcB (P s) (A s) := by
    have h2 : PolyFam T fun s => ∑ a, (-(P s a a / 128)) • lin (P s) a ^ 2 :=
      PolyFam.sum _ fun a _ => ((polyFam_lin hP a).pow 2).const_smul (((hP a a).div_const _).neg)
    refine (congrArg (PolyFam T) (funext fun s => ?_)).mp ((hlap hAf).add h2)
    simp [srcB, sub_eq_add_neg, Finset.sum_neg_distrib]
  have hdB : ∀ s, srcB (P s) (A s) ∈ restrictTotalDegree ι ℝ 2 := fun s =>
    Submodule.sub_mem _ (Submodule.sum_mem _ fun a _ =>
      pderiv_mem_lower a (pderiv_mem_lower a (hAm s)))
      (Submodule.sum_mem _ fun a _ => Submodule.smul_mem _ _ (lin_pow_mem _ a 2))
  obtain ⟨B, B', hB0, hBf, hBm, hBd, hBe⟩ := exists_poly_solution_diff (ν := 0) hT hP hfB hdB
  have hfD : PolyFam T fun s => srcD (P s) (B s) := by
    have h2 : PolyFam T fun s => ∑ a, (P s a a ^ 2 / 64) • (1 : MvPolynomial ι ℝ) :=
      PolyFam.sum _ fun a _ => PolyFam.smul _ (((hP a a).pow 2).div_const _) _
    exact (hlap hBf).add h2
  have hdD : ∀ s, srcD (P s) (B s) ∈ restrictTotalDegree ι ℝ 0 := fun s =>
    Submodule.add_mem _ (Submodule.sum_mem _ fun a _ =>
      pderiv_mem_lower a (pderiv_mem_lower a (hBm s)))
      (Submodule.sum_mem _ fun a _ => Submodule.smul_mem _ _ (by
        rw [mem_restrictTotalDegree, totalDegree_one]))
  obtain ⟨D, D', hD0, hDf, hDm, hDd, hDe⟩ := exists_poly_solution_diff (ν := 0) hT hP hfD hdD
  exact ⟨A, A', B, B', D, D', hA0, hB0, hD0, hAf, hBf, hDf, hDm, hAd, hBd, hDd, hAe, hBe, hDe⟩

/-- The combination solves the `p₁` equation. -/
theorem ABD.isHierarchy {T ν : ℝ} (hν : ν ≠ 0) {P : ℝ → Matrix ι ι ℝ}
    {A A' B B' D D' : ℝ → MvPolynomial ι ℝ} (h : ABD T P A A' B B' D D') :
    IsHierarchy T ν P 1 fun k s =>
      if k = 0 then 1 else (ν⁻¹ ^ 3) • A s + (ν⁻¹ ^ 2) • B s + ν⁻¹ • D s := by
  refine ⟨fun s => by simp, fun k h1 h2 => ?_, fun k hk => ?_, fun k h1 h2 x s hs => ?_⟩
  · obtain rfl : k = 1 := by omega
    simp [h.A_zero, h.B_zero, h.D_zero]
  · rcases Nat.eq_zero_or_pos k with rfl | hpos
    · simpa using PolyFam.const (T := T) (1 : MvPolynomial ι ℝ)
    · obtain rfl : k = 1 := by omega
      simp only [one_ne_zero, ↓reduceIte]
      exact ((h.A_fam.const_smul continuousOn_const).add
        (h.B_fam.const_smul continuousOn_const)).add (h.D_fam.const_smul continuousOn_const)
  · obtain rfl : k = 1 := by omega
    simp only [one_ne_zero, ↓reduceIte, map_add, smul_eval]
    have hA := h.A_diff.hasDeriv_eval x hs
    have hB := h.B_diff.hasDeriv_eval x hs
    have hD := h.D_diff.hasDeriv_eval x hs
    refine (((hA.const_mul _).fun_add (hB.const_mul _)).fun_add (hD.const_mul _)).congr_deriv ?_
    rw [h.A_eq s hs, h.B_eq s hs, h.D_eq s hs]
    rw [eval_source_one ν (P s) (p := fun i => if i = 0 then 1 else
      (ν⁻¹ ^ 3) • A s + (ν⁻¹ ^ 2) • B s + ν⁻¹ • D s) (by simp) x]
    have hDzero : ∀ a, pderiv a (D s) = 0 := fun a => pderiv_eq_zero_of_mem_zero a (h.D_mem s)
    have hsrc : ν / 12 * ∑ a, ((-(P s *ᵥ x) a / (4 * ν)) ^ 4
        + 6 * (-(P s *ᵥ x) a / (4 * ν)) ^ 2 * (-P s a a / (4 * ν)) + 3 * (-P s a a / (4 * ν)) ^ 2)
        = ν⁻¹ ^ 3 * (1 / 3072 * ∑ a, ((P s *ᵥ x) a) ^ 4)
          - ν⁻¹ ^ 2 * (1 / 128 * ∑ a, P s a a * ((P s *ᵥ x) a) ^ 2)
          + ν⁻¹ * (1 / 64 * ∑ a, P s a a ^ 2) := by
      simp only [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun a _ => ?_
      field_simp
      ring
    rw [hsrc]
    simp only [polyGen_add, polyGen_smul, map_add, smul_eval]
    have hD2 : ∀ μ : ℝ, eval x (polyGen μ (P s) (D s)) = 0 := by
      intro μ
      rw [eval_polyGen]
      simp [hDzero]
    rw [hD2, hD2, eval_polyGen 0, eval_polyGen 0, eval_polyGen ν, eval_polyGen ν, eval_srcA, eval_srcB,
      eval_srcD]
    have hinv : ν⁻¹ * ν = 1 := inv_mul_cancel₀ hν
    generalize ∑ a, eval x (pderiv a (pderiv a (A s))) = SA
    generalize ∑ a, eval x (pderiv a (pderiv a (B s))) = SB
    linear_combination (-(ν⁻¹ ^ 2 * SA) - ν⁻¹ * SB) * hinv

end ExecutionProofs
