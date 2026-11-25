import ExecutionProofs.SignTheorem

/-!
# Uniqueness of `p₁`, homogeneity, and the radial sign change (paper, Theorem 3 (iv))
-/

open Set Matrix MvPolynomial

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [DecidableEq ι] in
lemma IsHierarchy.mono {T ν : ℝ} {P : ℝ → Matrix ι ι ℝ} {K K' : ℕ}
    {p : ℕ → ℝ → MvPolynomial ι ℝ} (h : IsHierarchy T ν P K p) (hK : K' ≤ K) :
    IsHierarchy T ν P K' p :=
  ⟨h.zero, fun k h1 h2 => h.init k h1 (h2.trans hK), fun k hk => h.fam k (hk.trans hK),
    fun k h1 h2 => h.deriv k h1 (h2.trans hK)⟩

lemma const_of_hasDeriv_zero {f : ℝ → ℝ} {T s : ℝ} (hs : s ∈ Icc 0 T)
    (hf : ∀ r ∈ Icc 0 T, HasDerivWithinAt f 0 (Icc 0 T) r) : f s = f 0 :=
  constant_of_derivWithin_zero (fun r hr => (hf r hr).differentiableWithinAt)
    (fun r hr => (hf r ⟨hr.1, hr.2.le⟩).derivWithin
      (uniqueDiffOn_Icc (hr.1.trans_lt hr.2) r ⟨hr.1, hr.2.le⟩)) s hs

namespace Model

variable {C Φ : Matrix ι ι ℝ} {lam T φ ν : ℝ} {Y Z : ℝ → Matrix ι ι ℝ}
  {A A' B B' D D' : ℝ → MvPolynomial ι ℝ}

/-- **Uniqueness** of the leading coefficient: it is determined by the lattice solutions. -/
theorem p_one_unique (h : Model C Φ lam T φ Y Z) (hν : 0 < ν) {p q : ℕ → ℝ → MvPolynomial ι ℝ}
    (hp : IsHierarchy T ν (RiccatiFlow.P Y Z) 1 p) (hq : IsHierarchy T ν (RiccatiFlow.P Y Z) 1 q)
    {s : ℝ} (hs : s ∈ Icc 0 T) : p 1 s = q 1 s := by
  refine MvPolynomial.funext fun x => ?_
  have hPc : ∀ a b, ContinuousOn (fun s => RiccatiFlow.P Y Z s a b) (Icc 0 T) := fun a b s hs =>
    (h.deriv hs a b).continuousWithinAt
  obtain ⟨θ, hθ0, hθ⟩ := exists_theta h.hT.le hPc
  obtain ⟨B₁, hB₁, h₁⟩ := grid_theorem hν h.hlam h.hT h.hCs h.hΦs h.hC h.hφ h.hΦ h.flow hθ 1 hp
  obtain ⟨B₂, hB₂, h₂⟩ := grid_theorem hν h.hlam h.hT h.hCs h.hΦs h.hC h.hφ h.hΦ h.flow hθ 1 hq
  set g := gauss ν (RiccatiFlow.P Y Z s) (θ s) x with hg
  have hg0 : 0 < g := gauss_pos _ _ _ _
  have hV : ∀ x, 0 ≤ pot ν lam C x := fun x =>
    mul_nonneg (div_nonneg h.hlam (by positivity)) (h.hC x)
  have key : ∀ Δ, 0 < Δ → g * |eval x (p 1 s) - eval x (q 1 s)| ≤ Δ ^ 2 * (B₁ + B₂) := by
    intro Δ hΔ
    obtain ⟨u, hu, hu0⟩ := exists_lattice_solution (Δ := Δ) hν.le h.hT.le x hV
      (u₀ := fun n => gauss ν (RiccatiFlow.P Y Z 0) (θ 0) (pos x Δ n)) (A := 1) (fun n => by
        rw [abs_of_pos (gauss_pos _ _ _ _), gauss_zero h.flow hθ0, Real.exp_le_one_iff]
        have hq : 0 ≤ pos x Δ n ⬝ᵥ (Φ *ᵥ pos x Δ n) :=
          le_trans (mul_nonneg h.hφ.le (Finset.sum_nonneg fun i _ => mul_self_nonneg _)) (h.hΦ _)
        exact div_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr hq) (by positivity))
    have e₁ := h₁ Δ hΔ x u hu hu0 0 s hs
    have e₂ := h₂ Δ hΔ x u hu hu0 0 s hs
    have hpos0 : pos x Δ (0 : ι → ℤ) = x := by
      funext a
      simp [pos]
    simp only [hpos0, Finset.sum_range_succ, Finset.sum_range_zero, hp.zero, hq.zero, map_one,
      ← hg] at e₁ e₂
    have h3 : Δ ^ 2 * (g * |eval x (p 1 s) - eval x (q 1 s)|) ≤ Δ ^ 2 * (Δ ^ 2 * (B₁ + B₂)) := by
      have e : Δ ^ 2 * (g * |eval x (p 1 s) - eval x (q 1 s)|)
          = |Δ ^ 2 * (eval x (p 1 s) * g) - Δ ^ 2 * (eval x (q 1 s) * g)| := by
        rw [← mul_sub, ← sub_mul, abs_mul, abs_mul, abs_of_pos hg0,
          abs_of_pos (by positivity : 0 < Δ ^ 2)]
        ring
      have e₁' : |u 0 s - (g + Δ ^ 2 * (eval x (p 1 s) * g))| ≤ Δ ^ 2 * (Δ ^ 2 * B₁) := by
        refine le_trans (le_of_eq ?_) (e₁.trans (le_of_eq ?_))
        · congr 1
          ring
        · ring
      have e₂' : |u 0 s - (g + Δ ^ 2 * (eval x (q 1 s) * g))| ≤ Δ ^ 2 * (Δ ^ 2 * B₂) := by
        refine le_trans (le_of_eq ?_) (e₂.trans (le_of_eq ?_))
        · congr 1
          ring
        · ring
      rw [e]
      rw [abs_le] at e₁' e₂' ⊢
      constructor <;> linarith [e₁'.1, e₁'.2, e₂'.1, e₂'.2]
    exact le_of_mul_le_mul_left h3 (by positivity)
  by_contra hne
  have hd : 0 < |eval x (p 1 s) - eval x (q 1 s)| := abs_pos.mpr (sub_ne_zero.mpr hne)
  set e := g * |eval x (p 1 s) - eval x (q 1 s)| with he
  have he0 : 0 < e := mul_pos hg0 hd
  set Δ := min 1 (e / (2 * (B₁ + B₂ + 1))) with hΔ
  have hΔ0 : 0 < Δ := lt_min one_pos (by positivity)
  have h4 := key Δ hΔ0
  have h5 : Δ ^ 2 ≤ Δ := by
    have : Δ ≤ 1 := min_le_left _ _
    nlinarith
  have h6 : Δ * (2 * (B₁ + B₂ + 1)) ≤ e := by
    have : Δ ≤ e / (2 * (B₁ + B₂ + 1)) := min_le_right _ _
    rwa [le_div_iff₀ (by positivity)] at this
  nlinarith

/-- `p₁ = ν⁻³ A + ν⁻² B + ν⁻¹ D`. -/
theorem p_one_eq (h : Model C Φ lam T φ Y Z) (hν : 0 < ν)
    (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D') {p : ℕ → ℝ → MvPolynomial ι ℝ}
    (hp : IsHierarchy T ν (RiccatiFlow.P Y Z) 1 p) {s : ℝ} (hs : s ∈ Icc 0 T) :
    p 1 s = (ν⁻¹ ^ 3) • A s + (ν⁻¹ ^ 2) • B s + ν⁻¹ • D s := by
  have := h.p_one_unique hν hp (g.isHierarchy hν.ne') hs
  simpa using this

/-- `A` is homogeneous of degree four. -/
lemma A_homog (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D') {s : ℝ}
    (hs : s ∈ Icc 0 T) (r : ℝ) (x : ι → ℝ) : eval (r • x) (A s) = r ^ 4 * eval x (A s) := by
  have hunit := isUnit_iff_ne_zero.mpr (h.detY hs)
  obtain ⟨ξ, hξ⟩ : ∃ ξ, Y s *ᵥ ξ = x := ⟨(Y s)⁻¹ *ᵥ x, by
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hunit, Matrix.one_mulVec]⟩
  have hd : ∀ t ∈ Icc 0 T, HasDerivWithinAt
      (fun t => eval (Y t *ᵥ (r • ξ)) (A t) - r ^ 4 * eval (Y t *ᵥ ξ) (A t)) 0 (Icc 0 T) t := by
    intro t ht
    have h1 := fun η => transport_along h.flow ht (h.detY ht)
      (f := fun r => srcA (RiccatiFlow.P Y Z r)) g.A_diff (g.A_eq t ht) η
    refine ((h1 (r • ξ)).fun_sub ((h1 ξ).const_mul (r ^ 4))).congr_deriv ?_
    simp only [eval_srcA, h.P_mulVec ht, Matrix.mulVec_smul, Pi.smul_apply, smul_eq_mul, mul_pow,
      ← Finset.mul_sum]
    ring
  have hc := const_of_hasDeriv_zero hs hd
  simp only [g.A_zero, map_zero, mul_zero, sub_zero, Matrix.mulVec_smul, hξ] at hc
  linarith

omit [Fintype ι] in
/-- Homogeneity passes to derivatives. -/
lemma pderiv_homog [Fintype ι] {q : MvPolynomial ι ℝ} {r : ℝ} (hr : r ≠ 0) {n : ℕ}
    (hq : ∀ x, eval (r • x) q = r ^ (n + 1) * eval x q) (a : ι) (x : ι → ℝ) :
    eval (r • x) (pderiv a q) = r ^ n * eval x (pderiv a q) := by
  have h1 : HasDerivAt (fun t : ℝ => eval (r • (x + Pi.single a t)) q)
      (r * eval (r • x) (pderiv a q)) 0 := by
    have e : (fun t : ℝ => eval (r • (x + Pi.single a t)) q)
        = fun t : ℝ => eval (r • x + Pi.single a (r * t)) q := by
      funext t
      rw [smul_add, ← Pi.single_smul, smul_eq_mul]
    rw [e]
    have h2 := (hasDerivAt_eval_line q (r • x) a (r * 0)).comp 0
      ((hasDerivAt_id (0 : ℝ)).const_mul r)
    refine h2.congr_deriv ?_
    simp [mul_comm]
  have h2 : HasDerivAt (fun t : ℝ => eval (r • (x + Pi.single a t)) q)
      (r ^ (n + 1) * eval x (pderiv a q)) 0 := by
    simp only [hq]
    refine ((hasDerivAt_eval_line q x a 0).const_mul _).congr_deriv ?_
    simp
  have h3 := h1.unique h2
  have h4 : r * eval (r • x) (pderiv a q) = r * (r ^ n * eval x (pderiv a q)) := by
    rw [h3, pow_succ]
    ring
  exact mul_left_cancel₀ hr h4

/-- `B` is homogeneous of degree two. -/
lemma B_homog (h : Model C Φ lam T φ Y Z) (g : ABD T (RiccatiFlow.P Y Z) A A' B B' D D') {s : ℝ}
    (hs : s ∈ Icc 0 T) {r : ℝ} (hr : r ≠ 0) (x : ι → ℝ) :
    eval (r • x) (B s) = r ^ 2 * eval x (B s) := by
  have hunit := isUnit_iff_ne_zero.mpr (h.detY hs)
  obtain ⟨ξ, hξ⟩ : ∃ ξ, Y s *ᵥ ξ = x := ⟨(Y s)⁻¹ *ᵥ x, by
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hunit, Matrix.one_mulVec]⟩
  have hd : ∀ t ∈ Icc 0 T, HasDerivWithinAt
      (fun t => eval (Y t *ᵥ (r • ξ)) (B t) - r ^ 2 * eval (Y t *ᵥ ξ) (B t)) 0 (Icc 0 T) t := by
    intro t ht
    have h1 := fun η => transport_along h.flow ht (h.detY ht)
      (f := fun r => srcB (RiccatiFlow.P Y Z r) (A r)) g.B_diff (g.B_eq t ht) η
    refine ((h1 (r • ξ)).fun_sub ((h1 ξ).const_mul (r ^ 2))).congr_deriv ?_
    have hdd : ∀ a x, eval (r • x) (pderiv a (pderiv a (A t)))
        = r ^ 2 * eval x (pderiv a (pderiv a (A t))) := fun a x =>
      pderiv_homog hr (n := 2) (pderiv_homog hr (n := 3) (h.A_homog g ht r) a) a x
    simp only [eval_srcB, h.P_mulVec ht, Matrix.mulVec_smul, Pi.smul_apply, smul_eq_mul, hdd,
      ← Finset.mul_sum]
    have e : ∑ a, RiccatiFlow.P Y Z t a a * (r * (Z t *ᵥ ξ) a) ^ 2
        = r ^ 2 * ∑ a, RiccatiFlow.P Y Z t a a * ((Z t *ᵥ ξ) a) ^ 2 := by
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun a _ => by ring
    rw [e]
    ring
  have hc := const_of_hasDeriv_zero hs hd
  simp only [g.B_zero, map_zero, mul_zero, sub_zero, Matrix.mulVec_smul, hξ] at hc
  linarith

omit [DecidableEq ι] in
/-- Euler's identity from homogeneity. -/
lemma euler {q : MvPolynomial ι ℝ} {n : ℕ}
    (hq : ∀ r : ℝ, r ≠ 0 → ∀ x, eval (r • x) q = r ^ n * eval x q) (x : ι → ℝ) :
    ∑ a, x a * eval x (pderiv a q) = n * eval x q := by
  have h1 : HasDerivWithinAt (fun r : ℝ => eval (r • x) q)
      (∑ a, x a * eval ((1 : ℝ) • x) (pderiv a q)) (Icc 0 2) 1 :=
    hasDeriv_eval_curve (y := fun r : ℝ => r • x) (fun a => by
      simpa using ((hasDerivAt_id (1 : ℝ)).mul_const (x a)).hasDerivWithinAt) q
  have h2 : HasDerivWithinAt (fun r : ℝ => eval (r • x) q) (n * eval x q) (Icc 0 2) 1 := by
    have h3 : HasDerivAt (fun r : ℝ => r ^ n * eval x q) (n * eval x q) 1 := by
      simpa using (hasDerivAt_pow n (1 : ℝ)).mul_const (eval x q)
    refine (h3.congr_of_eventuallyEq ?_).hasDerivWithinAt
    filter_upwards [isOpen_ne.mem_nhds (one_ne_zero : (1 : ℝ) ≠ 0)] with r hr
    exact hq r hr x
  have hu : UniqueDiffWithinAt ℝ (Icc (0 : ℝ) 2) 1 :=
    uniqueDiffOn_Icc (by norm_num) 1 ⟨by norm_num, by norm_num⟩
  have := hu.eq_deriv _ h1 h2
  simpa using this

/-- **Radial sign change.** The radial component `x · (-2 ν ∇p₁)` of the leading rate error is negative
below exactly one value of `ν`, zero there, and positive above it. -/
theorem radial_error_sign (h : Model C Φ lam T φ Y Z) {x : ι → ℝ} (hx : x ≠ 0) :
    ∃ νs, 0 < νs ∧ ∀ ν, 0 < ν → ∀ p, IsHierarchy T ν (RiccatiFlow.P Y Z) 1 p →
      ((∑ a, x a * (-2 * ν * eval x (pderiv a (p 1 T))) < 0 ↔ ν < νs) ∧
       (∑ a, x a * (-2 * ν * eval x (pderiv a (p 1 T))) = 0 ↔ ν = νs) ∧
       (0 < ∑ a, x a * (-2 * ν * eval x (pderiv a (p 1 T))) ↔ νs < ν)) := by
  have hPc : ∀ a b, ContinuousOn (fun s => RiccatiFlow.P Y Z s a b) (Icc 0 T) := fun a b s hs =>
    (h.deriv hs a b).continuousWithinAt
  obtain ⟨A, A', B, B', D, D', g⟩ := exists_ABD h.hT.le hPc
  have hTm : T ∈ Icc 0 T := ⟨h.hT.le, le_rfl⟩
  have hA := h.A_pos g hTm h.hT hx
  have hB := h.B_neg g hTm h.hT hx
  refine ⟨-2 * eval x (A T) / eval x (B T), div_pos_of_neg_of_neg (by linarith) hB,
    fun ν hν p hp => ?_⟩
  have hDz : ∀ a, pderiv a (D T) = 0 := fun a => pderiv_eq_zero_of_mem_zero a (g.D_mem T)
  have eA := euler (q := A T) (n := 4) (fun r _ => h.A_homog g hTm r) x
  have eB := euler (q := B T) (n := 2) (fun r hr => h.B_homog g hTm hr) x
  have hval : ∑ a, x a * (-2 * ν * eval x (pderiv a (p 1 T)))
      = -4 * ν⁻¹ ^ 2 * (2 * eval x (A T) + ν * eval x (B T)) := by
    rw [h.p_one_eq hν g hp hTm]
    simp only [map_add, Derivation.map_smul, smul_eval, hDz, map_zero, mul_zero, add_zero]
    have e : ∑ a, x a * (-2 * ν * (ν⁻¹ ^ 3 * eval x (pderiv a (A T))
        + ν⁻¹ ^ 2 * eval x (pderiv a (B T))))
        = -2 * ν * ν⁻¹ ^ 3 * ∑ a, x a * eval x (pderiv a (A T))
          + -2 * ν * ν⁻¹ ^ 2 * ∑ a, x a * eval x (pderiv a (B T)) := by
      rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun a _ => by ring
    rw [e, eA, eB]
    push_cast
    field_simp
    ring
  rw [hval]
  have hc : 0 < 4 * ν⁻¹ ^ 2 := by positivity
  have hneg : -4 * ν⁻¹ ^ 2 * (2 * eval x (A T) + ν * eval x (B T))
      = 4 * ν⁻¹ ^ 2 * (-(2 * eval x (A T) + ν * eval x (B T))) := by ring
  rw [hneg]
  refine ⟨?_, ?_, ?_⟩
  · have hm : ∀ t, 4 * ν⁻¹ ^ 2 * t < 0 ↔ t < 0 := fun t =>
      ⟨fun h' => by
        by_contra hcon
        exact absurd h' (not_lt.mpr (mul_nonneg hc.le (not_lt.mp hcon))),
       fun h' => mul_neg_of_pos_of_neg hc h'⟩
    rw [hm, lt_div_iff_of_neg hB]
    constructor <;> intro h' <;> linarith
  · rw [mul_eq_zero, eq_div_iff hB.ne]
    constructor
    · rintro (h' | h')
      · exact absurd h' hc.ne'
      · linarith
    · intro h'
      right
      linarith
  · rw [mul_pos_iff_of_pos_left hc, div_lt_iff_of_neg hB]
    constructor <;> intro h' <;> linarith

end Model

end ExecutionProofs
