import ExecutionProofs.PolyGauss

/-!
# The polynomial hierarchy

Existence of polynomials `p_k(s)` with continuous coefficients such that `ψ_k = p_k · gauss` solves
the hierarchy of `grid_expansion`.
-/

open Set Matrix MvPolynomial

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The transport-diffusion operator on the polynomial factor:
`ν ∑ₐ ∂ₐ² p - ½ ∑ₐᵦ P_ab x_b ∂ₐ p`. -/
noncomputable def polyGen (ν : ℝ) (P : Matrix ι ι ℝ) (p : MvPolynomial ι ℝ) : MvPolynomial ι ℝ :=
  ν • ∑ a, pderiv a (pderiv a p) - ∑ a, ∑ b, (P a b / 2) • (X b * pderiv a p)

/-! ### The Gaussian solves the heat equation -/

omit [DecidableEq ι] in
lemma eval_linear (P : Matrix ι ι ℝ) (ν : ℝ) (a : ι) (x : ι → ℝ) (p : MvPolynomial ι ℝ) :
    eval x (∑ b, (P a b / (4 * ν)) • (X b * p)) = (P *ᵥ x) a / (4 * ν) * eval x p := by
  simp only [map_sum, smul_eval, map_mul, eval_X, Matrix.mulVec, dotProduct, Finset.sum_div,
    Finset.sum_mul]
  exact Finset.sum_congr rfl fun b _ => by ring

omit [DecidableEq ι] in
lemma eval_dop (P : Matrix ι ι ℝ) (ν : ℝ) (a : ι) (x : ι → ℝ) (p : MvPolynomial ι ℝ) :
    eval x (Dop ν P a p) = eval x (pderiv a p) - (P *ᵥ x) a / (4 * ν) * eval x p := by
  rw [Dop, map_sub, eval_linear]

lemma eval_pderiv_linear (P : Matrix ι ι ℝ) (ν : ℝ) (a : ι) (x : ι → ℝ) (p : MvPolynomial ι ℝ) :
    eval x (pderiv a (∑ b, (P a b / (4 * ν)) • (X b * p)))
      = P a a / (4 * ν) * eval x p + (P *ᵥ x) a / (4 * ν) * eval x (pderiv a p) := by
  have h : ∀ b, pderiv a ((P a b / (4 * ν)) • (X b * p))
      = (P a b / (4 * ν)) • (X b * pderiv a p) + (P a b / (4 * ν)) • (pderiv a (X b) * p) := by
    intro b
    rw [(pderiv a).map_smul, pderiv_mul, smul_add, add_comm]
  simp only [map_sum, h, map_add, Finset.sum_add_distrib]
  have h1 : ∑ b, eval x ((P a b / (4 * ν)) • (X b * pderiv a p))
      = (P *ᵥ x) a / (4 * ν) * eval x (pderiv a p) := by
    rw [← map_sum, eval_linear]
  have h2 : ∑ b, eval x ((P a b / (4 * ν)) • (pderiv a (X b) * p))
      = P a a / (4 * ν) * eval x p := by
    rw [Finset.sum_eq_single a]
    · simp [smul_eval]
    · intro b _ hb
      simp [pderiv_X, Pi.single_eq_of_ne hb]
    · simp
  rw [h1, h2, add_comm]

/-- Two applications of `Dop`. -/
lemma eval_dop_dop (P : Matrix ι ι ℝ) (ν : ℝ) (a : ι) (x : ι → ℝ) (p : MvPolynomial ι ℝ) :
    eval x (Dop ν P a (Dop ν P a p)) = eval x (pderiv a (pderiv a p))
      - 2 * ((P *ᵥ x) a / (4 * ν)) * eval x (pderiv a p)
      + (((P *ᵥ x) a / (4 * ν)) ^ 2 - P a a / (4 * ν)) * eval x p := by
  rw [eval_dop, eval_dop]
  have h : eval x (pderiv a (Dop ν P a p)) = eval x (pderiv a (pderiv a p))
      - (P a a / (4 * ν) * eval x p + (P *ᵥ x) a / (4 * ν) * eval x (pderiv a p)) := by
    rw [Dop, map_sub, map_sub, eval_pderiv_linear]
  rw [h]
  ring

/-- `ν ∑ₐ ∂ₐ²` on `p · gauss`, split into transport-diffusion and multiplication. -/
lemma sum_dop_dop {ν : ℝ} (hν : ν ≠ 0) (P : Matrix ι ι ℝ) (x : ι → ℝ) (p : MvPolynomial ι ℝ) :
    ν * ∑ a, eval x (Dop ν P a (Dop ν P a p)) = eval x (polyGen ν P p)
      + ν * (∑ a, (((P *ᵥ x) a / (4 * ν)) ^ 2 - P a a / (4 * ν))) * eval x p := by
  have h1 : eval x (polyGen ν P p) = ν * ∑ a, eval x (pderiv a (pderiv a p))
      - ∑ a, (P *ᵥ x) a / 2 * eval x (pderiv a p) := by
    rw [polyGen, map_sub, smul_eval, map_sum, map_sum]
    congr 1
    refine Finset.sum_congr rfl fun a _ => ?_
    simp only [map_sum, smul_eval, map_mul, eval_X, Matrix.mulVec, dotProduct, Finset.sum_div,
      Finset.sum_mul]
    exact Finset.sum_congr rfl fun b _ => by ring
  rw [h1]
  simp only [eval_dop_dop, Finset.sum_add_distrib, Finset.sum_sub_distrib, mul_add, mul_sub,
    Finset.mul_sum]
  have h2 : ∑ a, ν * (2 * ((P *ᵥ x) a / (4 * ν)) * eval x (pderiv a p))
      = ∑ a, (P *ᵥ x) a / 2 * eval x (pderiv a p) := by
    refine Finset.sum_congr rfl fun a _ => ?_
    field_simp
    ring
  rw [h2]
  congr 1
  rw [sub_mul, Finset.sum_mul, Finset.sum_mul, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun a _ => ?_
  ring

omit [DecidableEq ι] in
/-- Time derivative of the Gaussian along the Riccati flow. -/
lemma hasDeriv_gauss {C : Matrix ι ι ℝ} {lam ν T : ℝ} (hν : ν ≠ 0) {P : ℝ → Matrix ι ι ℝ}
    {θ : ℝ → ℝ} {s : ℝ} (hPs : (P s).IsSymm)
    (hP : ∀ i j, HasDerivWithinAt (fun s => P s i j)
      (2 * lam * C i j - 1 / 2 * (P s * P s) i j) (Icc 0 T) s)
    (hθ : HasDerivWithinAt θ (1 / 4 * (P s).trace) (Icc 0 T) s) (x : ι → ℝ) :
    HasDerivWithinAt (fun s => gauss ν (P s) (θ s) x)
      ((ν * (∑ a, (((P s *ᵥ x) a / (4 * ν)) ^ 2 - P s a a / (4 * ν)))
        - lam / (4 * ν) * (x ⬝ᵥ (C *ᵥ x))) * gauss ν (P s) (θ s) x) (Icc 0 T) s := by
  have hform : HasDerivWithinAt (fun s => x ⬝ᵥ (P s *ᵥ x))
      (2 * lam * (x ⬝ᵥ (C *ᵥ x)) - 1 / 2 * (x ⬝ᵥ ((P s * P s) *ᵥ x))) (Icc 0 T) s := by
    simp only [dotProduct, Matrix.mulVec]
    refine (HasDerivWithinAt.fun_sum fun i _ => (HasDerivWithinAt.fun_sum fun j _ =>
      (hP i j).mul_const (x j)).const_mul (x i)).congr_deriv ?_
    simp only [Finset.mul_sum, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring
  have hsq : x ⬝ᵥ ((P s * P s) *ᵥ x) = ∑ a, ((P s *ᵥ x) a) ^ 2 := by
    rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hPs]
    simp only [dotProduct, sq]
  unfold gauss
  refine (((hform.fun_neg.div_const (8 * ν)).sub hθ).exp).congr_deriv ?_
  have hcoef : -(2 * lam * (x ⬝ᵥ (C *ᵥ x)) - 1 / 2 * (x ⬝ᵥ ((P s * P s) *ᵥ x))) / (8 * ν)
      - 1 / 4 * (P s).trace = ν * (∑ a, (((P s *ᵥ x) a / (4 * ν)) ^ 2 - P s a a / (4 * ν)))
        - lam / (4 * ν) * (x ⬝ᵥ (C *ᵥ x)) := by
    rw [hsq, Finset.sum_sub_distrib, Matrix.trace]
    simp only [Matrix.diag_apply, div_pow, ← Finset.sum_div]
    field_simp
    ring
  rw [← hcoef]
  simp only [Pi.sub_apply]
  ring

/-! ### The polynomial ODE -/

omit [Fintype ι] [DecidableEq ι] in
lemma mem_restrict_of_monomials {n : ℕ} (ℓ : MvPolynomial ι ℝ →ₗ[ℝ] MvPolynomial ι ℝ)
    (hℓ : ∀ d : ι →₀ ℕ, (d.sum fun _ e => e) ≤ n → ∀ r : ℝ,
      ℓ (monomial d r) ∈ restrictTotalDegree ι ℝ n)
    {p : MvPolynomial ι ℝ} (hp : p ∈ restrictTotalDegree ι ℝ n) :
    ℓ p ∈ restrictTotalDegree ι ℝ n := by
  rw [mem_restrictTotalDegree] at hp
  rw [p.as_sum, map_sum]
  refine Submodule.sum_mem _ fun d hd => hℓ d ?_ _
  exact (le_totalDegree hd).trans hp

omit [Fintype ι] [DecidableEq ι] in
lemma X_mul_pderiv_mem {n : ℕ} (a b : ι) {p : MvPolynomial ι ℝ}
    (hp : p ∈ restrictTotalDegree ι ℝ n) : X b * pderiv a p ∈ restrictTotalDegree ι ℝ n := by
  refine mem_restrict_of_monomials
    (LinearMap.mulLeft ℝ (X b) ∘ₗ (pderiv a).toLinearMap) (fun d hd r => ?_) hp
  simp only [LinearMap.comp_apply, LinearMap.mulLeft_apply, Derivation.coeFn_coe]
  rw [pderiv_monomial]
  by_cases ha : d a = 0
  · simp [ha]
  · have hX : (X b : MvPolynomial ι ℝ) = monomial (Finsupp.single b 1) 1 := rfl
    rw [hX, monomial_mul_monomial, mem_restrictTotalDegree]
    refine (totalDegree_monomial_le _ _).trans ?_
    change ((Finsupp.single b 1 + (d - Finsupp.single a 1)).sum fun _ e => e) ≤ n
    have hle : Finsupp.single a 1 ≤ d := by
      rw [Finsupp.single_le_iff]
      omega
    have hsplit : (d.sum fun _ e => e)
        = ((d - Finsupp.single a 1).sum fun _ e => e) + 1 := by
      conv_lhs => rw [← tsub_add_cancel_of_le hle]
      rw [Finsupp.sum_add_index' (fun _ => rfl) (fun _ _ _ => rfl),
        Finsupp.sum_single_index rfl]
    rw [Finsupp.sum_add_index' (fun _ => rfl) (fun _ _ _ => rfl), Finsupp.sum_single_index rfl]
    omega

omit [Fintype ι] in
lemma pderiv_mem {n : ℕ} (a : ι) {p : MvPolynomial ι ℝ}
    (hp : p ∈ restrictTotalDegree ι ℝ n) : pderiv a p ∈ restrictTotalDegree ι ℝ n := by
  refine mem_restrict_of_monomials (pderiv a).toLinearMap (fun d hd r => ?_) hp
  simp only [Derivation.coeFn_coe]
  rw [pderiv_monomial, mem_restrictTotalDegree]
  refine (totalDegree_monomial_le _ _).trans (le_trans ?_ hd)
  change ((d - Finsupp.single a 1).sum fun _ e => e) ≤ d.sum fun _ e => e
  by_cases ha : d a = 0
  · have h0 : d - Finsupp.single a 1 = d := by
      ext i
      by_cases hi : i = a
      · subst hi
        simp [ha]
      · simp [Ne.symm hi]
    rw [h0]
  · have hle : Finsupp.single a 1 ≤ d := by
      rw [Finsupp.single_le_iff]
      omega
    conv_rhs => rw [← tsub_add_cancel_of_le hle]
    rw [Finsupp.sum_add_index' (fun _ => rfl) (fun _ _ _ => rfl), Finsupp.sum_single_index rfl]
    omega

omit [DecidableEq ι] in
lemma polyGen_add (ν : ℝ) (P : Matrix ι ι ℝ) (p q : MvPolynomial ι ℝ) :
    polyGen ν P (p + q) = polyGen ν P p + polyGen ν P q := by
  simp only [polyGen, map_add, mul_add, smul_add, Finset.sum_add_distrib]
  abel

omit [DecidableEq ι] in
lemma polyGen_smul (ν : ℝ) (P : Matrix ι ι ℝ) (r : ℝ) (p : MvPolynomial ι ℝ) :
    polyGen ν P (r • p) = r • polyGen ν P p := by
  simp only [polyGen, Derivation.map_smul, mul_smul_comm, smul_sub, Finset.smul_sum]
  congr 1
  · exact Finset.sum_congr rfl fun a _ => smul_comm _ _ _
  · exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => smul_comm _ _ _

/-- `polyGen` as a linear map. -/
noncomputable def polyGenₗ (ν : ℝ) (P : Matrix ι ι ℝ) :
    MvPolynomial ι ℝ →ₗ[ℝ] MvPolynomial ι ℝ where
  toFun := polyGen ν P
  map_add' := polyGen_add ν P
  map_smul' := polyGen_smul ν P

lemma polyGen_mem {n : ℕ} (ν : ℝ) (P : Matrix ι ι ℝ) {p : MvPolynomial ι ℝ}
    (hp : p ∈ restrictTotalDegree ι ℝ n) : polyGen ν P p ∈ restrictTotalDegree ι ℝ n := by
  refine Submodule.sub_mem _ (Submodule.smul_mem _ _ (Submodule.sum_mem _ fun a _ =>
    pderiv_mem a (pderiv_mem a hp))) (Submodule.sum_mem _ fun a _ => Submodule.sum_mem _ fun b _ =>
    Submodule.smul_mem _ _ (X_mul_pderiv_mem a b hp))

omit [Fintype ι] [DecidableEq ι] in
/-- A polynomial family has bounded degree. -/
lemma PolyFam.exists_degree {T : ℝ} {q : ℝ → MvPolynomial ι ℝ} (hq : PolyFam T q) :
    ∃ n, ∀ s, q s ∈ restrictTotalDegree ι ℝ n := by
  induction hq with
  | smul α hα m =>
    refine ⟨m.totalDegree, fun s => Submodule.smul_mem _ _ ?_⟩
    rw [mem_restrictTotalDegree]
  | add _ _ ih₁ ih₂ =>
    obtain ⟨n₁, h₁⟩ := ih₁
    obtain ⟨n₂, h₂⟩ := ih₂
    refine ⟨max n₁ n₂, fun s => Submodule.add_mem _ ?_ ?_⟩
    · have := h₁ s
      rw [mem_restrictTotalDegree] at this ⊢
      exact this.trans (le_max_left _ _)
    · have := h₂ s
      rw [mem_restrictTotalDegree] at this ⊢
      exact this.trans (le_max_right _ _)

omit [Fintype ι] [DecidableEq ι] in
/-- Linear functionals of a polynomial family are continuous. -/
lemma PolyFam.continuousOn_functional {T : ℝ} {q : ℝ → MvPolynomial ι ℝ} (hq : PolyFam T q)
    (ℓ : MvPolynomial ι ℝ →ₗ[ℝ] ℝ) : ContinuousOn (fun s => ℓ (q s)) (Icc 0 T) := by
  induction hq with
  | smul α hα m =>
    simp only [map_smul, smul_eq_mul]
    exact hα.mul continuousOn_const
  | add _ _ ih₁ ih₂ =>
    simp only [map_add]
    exact ih₁.add ih₂

/-- `p' ` is the coefficientwise derivative of `p` on `[0, T]`. -/
def DiffFam (T : ℝ) (p p' : ℝ → MvPolynomial ι ℝ) : Prop :=
  ∃ (k : ℕ) (m : Fin k → MvPolynomial ι ℝ) (α α' : Fin k → ℝ → ℝ),
    (∀ s, p s = ∑ i, α i s • m i) ∧ (∀ s ∈ Icc 0 T, p' s = ∑ i, α' i s • m i) ∧
    ∀ i, ∀ s ∈ Icc 0 T, HasDerivWithinAt (α i) (α' i s) (Icc 0 T) s

omit [Fintype ι] [DecidableEq ι] in
lemma DiffFam.hasDeriv_eval {T : ℝ} {p p' : ℝ → MvPolynomial ι ℝ} (h : DiffFam T p p')
    (x : ι → ℝ) {s : ℝ} (hs : s ∈ Icc 0 T) :
    HasDerivWithinAt (fun s => eval x (p s)) (eval x (p' s)) (Icc 0 T) s := by
  obtain ⟨k, m, α, α', hp, hp', hα⟩ := h
  rw [hp' s hs]
  simp only [hp, map_sum, smul_eval]
  exact HasDerivWithinAt.fun_sum fun i _ => (hα i s hs).mul_const _

/-- **Polynomial ODE.** `p' = polyGen p + f`, `p(0) = 0` has a polynomial family solution, of degree
at most that of the source. -/
theorem exists_poly_solution_diff {T ν : ℝ} (hT : 0 ≤ T) {P : ℝ → Matrix ι ι ℝ}
    (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T)) {f : ℝ → MvPolynomial ι ℝ}
    (hf : PolyFam T f) {n : ℕ} (hn : ∀ s, f s ∈ restrictTotalDegree ι ℝ n) :
    ∃ p p' : ℝ → MvPolynomial ι ℝ, p 0 = 0 ∧ PolyFam T p ∧
      (∀ s, p s ∈ restrictTotalDegree ι ℝ n) ∧ DiffFam T p p' ∧
      ∀ s ∈ Icc 0 T, p' s = polyGen ν (P s) (p s) + f s := by
  classical
  set W : Submodule ℝ (MvPolynomial ι ℝ) := restrictTotalDegree ι ℝ n with hW
  let e := Module.finBasis ℝ W
  obtain ⟨g, hg⟩ := W.subtype.exists_leftInverse_of_injective (Submodule.ker_subtype W)
  have hgW : ∀ w (hw : w ∈ W), g w = ⟨w, hw⟩ := fun w hw => by
    have := LinearMap.congr_fun hg ⟨w, hw⟩
    simpa using this
  let ρ : Fin (Module.finrank ℝ W) → MvPolynomial ι ℝ →ₗ[ℝ] ℝ := fun i => e.coord i ∘ₗ g
  have hρ : ∀ w ∈ W, ∑ i, ρ i w • (e i : MvPolynomial ι ℝ) = w := by
    intro w hw
    have h1 := congrArg (W.subtype) (e.sum_repr ⟨w, hw⟩)
    rw [map_sum] at h1
    simp only [map_smul, Submodule.subtype_apply] at h1
    refine Eq.trans (Finset.sum_congr rfl fun i _ => ?_) h1
    simp only [ρ, LinearMap.comp_apply, hgW w hw, Module.Basis.coord_apply]
  -- the matrix of the generator in the basis
  let M : Fin (Module.finrank ℝ W) → Fin (Module.finrank ℝ W) → ℝ → ℝ := fun j i s =>
    ρ i (polyGen ν (P s) (e j : MvPolynomial ι ℝ))
  have hgenfam : ∀ q : MvPolynomial ι ℝ, PolyFam T fun s => polyGen ν (P s) q := by
    intro q
    have h1 : PolyFam T fun _ => ν • ∑ a, pderiv a (pderiv a q) := PolyFam.const _
    have h2 : PolyFam T fun s => ∑ a, ∑ b, (-(P s a b / 2)) • (X b * pderiv a q) :=
      PolyFam.sum _ fun a _ => PolyFam.sum _ fun b _ =>
        PolyFam.smul _ (((hP a b).div_const _).neg) _
    refine (congrArg (PolyFam T) (funext fun s => ?_)).mp (h1.add h2)
    simp [polyGen, sub_eq_add_neg, Finset.sum_neg_distrib]
  have hM : ∀ j i, ContinuousOn (M j i) (Icc 0 T) := fun j i =>
    (hgenfam _).continuousOn_functional (ρ i)
  have hβ : ∀ i, ContinuousOn (fun s => ρ i (f s)) (Icc 0 T) := fun i =>
    hf.continuousOn_functional (ρ i)
  obtain ⟨Kb, hKb⟩ : ∃ Kb, ∀ j i, ∀ s ∈ Icc 0 T, |M j i s| ≤ Kb := by
    have h := fun j i => isCompact_Icc.exists_bound_of_continuousOn (hM j i)
    choose B hB using h
    refine ⟨∑ j, ∑ i, |B j i|, fun j i s hs => ((hB j i s hs).trans (le_abs_self _)).trans ?_⟩
    exact (Finset.single_le_sum (f := fun i => |B j i|) (fun _ _ => abs_nonneg _)
      (Finset.mem_univ i)).trans (Finset.single_le_sum (f := fun j => ∑ i, |B j i|)
      (fun _ _ => Finset.sum_nonneg fun _ _ => abs_nonneg _) (Finset.mem_univ j))
  obtain ⟨β, hβb⟩ : ∃ β, ∀ i, ∀ s ∈ Icc 0 T, |ρ i (f s)| ≤ β := by
    have h := fun i => isCompact_Icc.exists_bound_of_continuousOn (hβ i)
    choose B hB using h
    refine ⟨∑ i, |B i|, fun i s hs => ((hB i s hs).trans (le_abs_self _)).trans ?_⟩
    exact Finset.single_le_sum (f := fun i => |B i|) (fun _ _ => abs_nonneg _) (Finset.mem_univ i)
  obtain ⟨α, hα0, -, hα⟩ := exists_linear_solution
    (nb := fun (j : Fin (Module.finrank ℝ W)) (_ : Fin (Module.finrank ℝ W)) => j) hT
    (κ := M) hM hKb (V := fun _ => 0) (fun _ => le_rfl) hβ hβb (u₀ := fun _ => 0) (A := 0)
    (fun _ => by simp)
  refine ⟨fun s => ∑ i, α i s • (e i : MvPolynomial ι ℝ),
    fun s => ∑ i, (∑ j, M j i s * α j s - 0 * α i s + ρ i (f s)) • (e i : MvPolynomial ι ℝ),
    ?_, ?_, fun s => ?_, ⟨_, fun i => (e i : MvPolynomial ι ℝ), α, _, fun _ => rfl,
      fun _ _ => rfl, fun i s hs => hα i s hs⟩, fun s _ => ?_⟩
  · simp [hα0]
  · exact PolyFam.sum _ fun i _ =>
      PolyFam.smul _ (fun s hs => (hα i s hs).continuousWithinAt) _
  · exact Submodule.sum_mem _ fun i _ => Submodule.smul_mem _ _ (e i).2
  · have hgen : polyGen ν (P s) (∑ j, α j s • (e j : MvPolynomial ι ℝ))
        = ∑ j, α j s • polyGen ν (P s) (e j : MvPolynomial ι ℝ) := by
      have h := map_sum (polyGenₗ ν (P s)) (fun j => α j s • (e j : MvPolynomial ι ℝ)) Finset.univ
      simp only [map_smul] at h
      exact h
    have hexp : ∀ j, polyGen ν (P s) (e j : MvPolynomial ι ℝ)
        = ∑ i, M j i s • (e i : MvPolynomial ι ℝ) := fun j =>
      (hρ _ (polyGen_mem ν (P s) (e j).2)).symm
    have hfexp : f s = ∑ i, ρ i (f s) • (e i : MvPolynomial ι ℝ) := (hρ _ (hn s)).symm
    rw [hgen]
    conv_rhs => rw [hfexp]
    simp only [hexp, Finset.smul_sum, smul_smul]
    rw [Finset.sum_comm, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [zero_mul, sub_zero, add_smul, Finset.sum_smul]
    congr 1
    exact Finset.sum_congr rfl fun j _ => by rw [mul_comm]

/-- Evaluated form, with the degree left implicit. -/
theorem exists_poly_solution {T ν : ℝ} (hT : 0 ≤ T) {P : ℝ → Matrix ι ι ℝ}
    (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T)) {f : ℝ → MvPolynomial ι ℝ}
    (hf : PolyFam T f) :
    ∃ p : ℝ → MvPolynomial ι ℝ, p 0 = 0 ∧ PolyFam T p ∧
      ∀ x, ∀ s ∈ Icc 0 T, HasDerivWithinAt (fun s => eval x (p s))
        (eval x (polyGen ν (P s) (p s)) + eval x (f s)) (Icc 0 T) s := by
  obtain ⟨n, hn⟩ := hf.exists_degree
  obtain ⟨p, p', h0, hfam, -, hd, hp'⟩ := exists_poly_solution_diff (ν := ν) hT hP hf hn
  refine ⟨p, h0, hfam, fun x s hs => ?_⟩
  rw [← map_add, ← hp' s hs]
  exact hd.hasDeriv_eval x hs

/-- Source of the `k`-th equation of the hierarchy. -/
noncomputable def source (ν : ℝ) (P : Matrix ι ι ℝ) (p : ℕ → MvPolynomial ι ℝ) (k : ℕ) :
    MvPolynomial ι ℝ :=
  ν • ∑ i ∈ Finset.range k, c (k - i) • ∑ a, (Dop ν P a)^[2 * (k - i) + 2] (p i)

omit [DecidableEq ι] in
lemma source_congr {ν : ℝ} {P : Matrix ι ι ℝ} {p q : ℕ → MvPolynomial ι ℝ} {k : ℕ}
    (h : ∀ i < k, p i = q i) : source ν P p k = source ν P q k := by
  have e : ∀ i ∈ Finset.range k,
      c (k - i) • ∑ a, (Dop ν P a)^[2 * (k - i) + 2] (p i)
        = c (k - i) • ∑ a, (Dop ν P a)^[2 * (k - i) + 2] (q i) := fun i hi => by
    rw [h i (Finset.mem_range.mp hi)]
  exact congrArg (fun z => ν • z) (Finset.sum_congr rfl e)

omit [DecidableEq ι] in
lemma PolyFam.source {T ν : ℝ} {P : ℝ → Matrix ι ι ℝ}
    (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T)) {p : ℕ → ℝ → MvPolynomial ι ℝ} {k : ℕ}
    (hp : ∀ i < k, PolyFam T (p i)) : PolyFam T fun s => source ν (P s) (fun i => p i s) k := by
  have h1 : ∀ i < k, PolyFam T fun s =>
      c (k - i) • ∑ a, (Dop ν (P s) a)^[2 * (k - i) + 2] (p i s) := by
    intro i hi
    have h2 : PolyFam T fun s => ∑ a, (Dop ν (P s) a)^[2 * (k - i) + 2] (p i s) :=
      PolyFam.sum _ fun a _ => (hp i hi).dop_iterate hP a _
    exact h2.const_smul (β := fun _ => c (k - i)) continuousOn_const
  have h3 : PolyFam T fun s => ∑ i ∈ Finset.range k,
      c (k - i) • ∑ a, (Dop ν (P s) a)^[2 * (k - i) + 2] (p i s) :=
    PolyFam.sum _ fun i hi => h1 i (Finset.mem_range.mp hi)
  exact h3.const_smul (β := fun _ => ν) continuousOn_const

/-- The statement of `exists_hierarchy`: `p 0, …, p K` solve the polynomial hierarchy. -/
structure IsHierarchy (T ν : ℝ) (P : ℝ → Matrix ι ι ℝ) (K : ℕ)
    (p : ℕ → ℝ → MvPolynomial ι ℝ) : Prop where
  zero : ∀ s, p 0 s = 1
  init : ∀ k, 1 ≤ k → k ≤ K → p k 0 = 0
  fam : ∀ k ≤ K, PolyFam T (p k)
  deriv : ∀ k, 1 ≤ k → k ≤ K → ∀ x, ∀ s ∈ Icc 0 T, HasDerivWithinAt (fun s => eval x (p k s))
    (eval x (polyGen ν (P s) (p k s)) + eval x (source ν (P s) (fun i => p i s) k)) (Icc 0 T) s

/-- Replace the entry after `K`. -/
def extendSeq (p : ℕ → ℝ → MvPolynomial ι ℝ) (K : ℕ) (q : ℝ → MvPolynomial ι ℝ) :
    ℕ → ℝ → MvPolynomial ι ℝ := Function.update p (K + 1) q

/-- **The hierarchy has polynomial solutions** up to any order. -/
theorem exists_hierarchy {T ν : ℝ} (hT : 0 ≤ T) {P : ℝ → Matrix ι ι ℝ}
    (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T)) (K : ℕ) :
    ∃ p, IsHierarchy T ν P K p := by
  induction K with
  | zero =>
    exact ⟨fun _ _ => 1, fun _ => rfl, fun k h1 h2 => absurd h1 (by omega),
      fun k _ => PolyFam.const _, fun k h1 h2 => absurd h1 (by omega)⟩
  | succ K ih =>
    obtain ⟨p, hp⟩ := ih
    obtain ⟨q, hq0, hqf, hqd⟩ := exists_poly_solution (ν := ν) hT hP
      (PolyFam.source (ν := ν) hP (k := K + 1) fun i hi => hp.fam i (Nat.lt_succ_iff.mp hi))
    have hnew : extendSeq p K q (K + 1) = q := Function.update_self _ _ _
    have hold : ∀ k, k ≤ K → extendSeq p K q k = p k := fun k hk =>
      Function.update_of_ne (by omega) _ _
    refine ⟨extendSeq p K q, fun s => ?_, fun k h1 h2 => ?_, fun k hk => ?_,
      fun k h1 h2 x s hs => ?_⟩
    · rw [hold 0 (Nat.zero_le K), hp.zero]
    · rcases Nat.lt_or_ge K k with h | h
      · obtain rfl : k = K + 1 := by omega
        rw [hnew, hq0]
      · rw [hold k h, hp.init k h1 h]
    · rcases Nat.lt_or_ge K k with h | h
      · obtain rfl : k = K + 1 := by omega
        rw [hnew]
        exact hqf
      · rw [hold k h]
        exact hp.fam k h
    · have hsrc : source ν (P s) (fun i => extendSeq p K q i s) k
          = source ν (P s) (fun i => p i s) k :=
        source_congr fun i hi => by rw [hold i (by omega)]
      rw [hsrc]
      rcases Nat.lt_or_ge K k with h | h
      · obtain rfl : k = K + 1 := by omega
        rw [hnew]
        exact hqd x s hs
      · rw [hold k h]
        exact hp.deriv k h1 h x s hs

end ExecutionProofs
