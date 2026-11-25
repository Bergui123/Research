import ExecutionProofs.GridTheorem

/-!
# The leading coefficient and the rate (paper, Theorem 3 (ii)–(iii))
-/

open Set Matrix MvPolynomial

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] in
/-- Partial derivative of a polynomial from its values along the coordinate line. -/
lemma eval_pderiv_of_eval {r : MvPolynomial ι ℝ} {f : (ι → ℝ) → ℝ} (hf : ∀ x, eval x r = f x)
    {a : ι} {x : ι → ℝ} {d : ℝ} (hd : HasDerivAt (fun t : ℝ => f (x + Pi.single a t)) d 0) :
    eval x (pderiv a r) = d := by
  have h := hasDerivAt_eval_line r x a 0
  simp only [hf, Pi.single_zero, add_zero] at h
  exact h.unique hd

lemma mulVec_line (P : Matrix ι ι ℝ) (x : ι → ℝ) (a : ι) (t : ℝ) :
    (P *ᵥ (x + Pi.single a t)) a = (P *ᵥ x) a + P a a * t := by
  rw [Matrix.mulVec_add, Pi.add_apply]
  congr 1
  simp [Matrix.mulVec, dotProduct, Pi.single_apply]

/-- `∂ₐ⁴` of the Gaussian: `ℓ⁴ - 6 ℓ² γ + 3 γ²` with `ℓ = (P x)_a / 4ν`, `γ = P_aa / 4ν`. -/
lemma eval_dop_four (ν : ℝ) (P : Matrix ι ι ℝ) (a : ι) (x : ι → ℝ) :
    eval x ((Dop ν P a)^[4] 1) = ((P *ᵥ x) a / (4 * ν)) ^ 4
      - 6 * ((P *ᵥ x) a / (4 * ν)) ^ 2 * (P a a / (4 * ν)) + 3 * (P a a / (4 * ν)) ^ 2 := by
  set r : MvPolynomial ι ℝ := Dop ν P a (Dop ν P a 1) with hr
  have hr0 : ∀ x, eval x r = ((P *ᵥ x) a / (4 * ν)) ^ 2 - P a a / (4 * ν) := by
    intro x
    rw [hr, eval_dop_dop]
    simp
  have hline : ∀ x : ι → ℝ, HasDerivAt (fun t : ℝ => (P *ᵥ (x + Pi.single a t)) a / (4 * ν))
      (P a a / (4 * ν)) 0 := by
    intro x
    simp only [mulVec_line]
    simpa using (((hasDerivAt_id (0 : ℝ)).const_mul (P a a)).const_add ((P *ᵥ x) a)).div_const
      (4 * ν)
  have hr1 : ∀ x, eval x (pderiv a r) = 2 * ((P *ᵥ x) a / (4 * ν)) * (P a a / (4 * ν)) := by
    intro x
    refine eval_pderiv_of_eval hr0 ?_
    refine (((hline x).fun_pow 2).sub_const _).congr_deriv ?_
    simp
  have hr2 : ∀ x, eval x (pderiv a (pderiv a r)) = 2 * (P a a / (4 * ν)) ^ 2 := by
    intro x
    refine eval_pderiv_of_eval hr1 ?_
    refine (((hline x).const_mul 2).mul_const (P a a / (4 * ν))).congr_deriv ?_
    ring
  have h4 : (Dop ν P a)^[4] 1 = Dop ν P a (Dop ν P a r) := rfl
  rw [h4, eval_dop_dop, hr0, hr1, hr2]
  ring

/-- **The equation for `χ = p₁`.** With `φ_a = -(P x)_a / 4ν` and `φ_aa = -P_aa / 4ν` the source is
`(ν/12) ∑ₐ (φ_a⁴ + 6 φ_a² φ_aa + 3 φ_aa²)`. -/
theorem eval_source_one (ν : ℝ) (P : Matrix ι ι ℝ) {p : ℕ → MvPolynomial ι ℝ} (hp : p 0 = 1)
    (x : ι → ℝ) :
    eval x (source ν P p 1) = ν / 12 * ∑ a, ((-(P *ᵥ x) a / (4 * ν)) ^ 4
      + 6 * (-(P *ᵥ x) a / (4 * ν)) ^ 2 * (-P a a / (4 * ν)) + 3 * (-P a a / (4 * ν)) ^ 2) := by
  have h : source ν P p 1 = ν • (c 1 • ∑ a, (Dop ν P a)^[4] 1) := by
    rw [source, Finset.sum_range_one, hp]
  rw [h]
  simp only [smul_eval, map_sum, eval_dop_four, c_one]
  rw [← mul_assoc, mul_one_div, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  ring

/-! ### From the solution to its logarithm -/

/-- Elementary step: relative expansion to logarithmic expansion, with explicit constant. -/
lemma log_expansion {g a b u Δ m A B : ℝ} (hm : 0 < m) (hg : m ≤ g) (ha : |a| ≤ A) (hb : |b| ≤ A)
    (hB : 0 ≤ B) (hΔ0 : 0 ≤ Δ) (hΔ1 : Δ ≤ 1) (hsmall : Δ ^ 2 * (2 * A + B / m) ≤ 1 / 2)
    (hu : |u - g * (1 + Δ ^ 2 * a + Δ ^ 4 * b)| ≤ B * Δ ^ 6) :
    0 < u ∧ |Real.log u - Real.log g - Δ ^ 2 * a - Δ ^ 4 * (b - a ^ 2 / 2)|
      ≤ (B / m + (A + B / m) * (3 * A + B / m) / 2 + 2 * (2 * A + B / m) ^ 3) * Δ ^ 6 := by
  have hg0 : 0 < g := hm.trans_le hg
  have hA : 0 ≤ A := (abs_nonneg a).trans ha
  set E : ℝ := 2 * A + B / m with hE
  have hBm : 0 ≤ B / m := by positivity
  have hE0 : 0 ≤ E := by positivity
  set z : ℝ := u / g - 1 with hz
  have hug : u = g * (1 + z) := by
    rw [hz]
    field_simp
    ring
  set r : ℝ := z - (Δ ^ 2 * a + Δ ^ 4 * b) with hr
  have hr' : |r| ≤ B / m * Δ ^ 6 := by
    have h1 : g * r = u - g * (1 + Δ ^ 2 * a + Δ ^ 4 * b) := by
      rw [hug, hr]
      ring
    have h2 : g * |r| ≤ B * Δ ^ 6 := by
      calc g * |r| = |g * r| := by rw [abs_mul, abs_of_pos hg0]
        _ = |u - g * (1 + Δ ^ 2 * a + Δ ^ 4 * b)| := by rw [h1]
        _ ≤ B * Δ ^ 6 := hu
    have h3 : m * |r| ≤ B * Δ ^ 6 := (mul_le_mul_of_nonneg_right hg (abs_nonneg _)).trans h2
    rw [div_mul_eq_mul_div, le_div_iff₀ hm]
    linarith
  have hΔ2 : 0 ≤ Δ ^ 2 := by positivity
  have hΔ4 : Δ ^ 4 ≤ Δ ^ 2 := pow_le_pow_of_le_one hΔ0 hΔ1 (by norm_num)
  have hΔ64 : Δ ^ 6 ≤ Δ ^ 4 := pow_le_pow_of_le_one hΔ0 hΔ1 (by norm_num)
  have hΔ62 : Δ ^ 6 ≤ Δ ^ 2 := pow_le_pow_of_le_one hΔ0 hΔ1 (by norm_num)
  have hza : |z - Δ ^ 2 * a| ≤ (A + B / m) * Δ ^ 4 := by
    have h1 : z - Δ ^ 2 * a = Δ ^ 4 * b + r := by
      rw [hr]
      ring
    rw [h1]
    refine (abs_add_le _ _).trans ?_
    have h2 : |Δ ^ 4 * b| ≤ A * Δ ^ 4 := by
      rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ Δ ^ 4), mul_comm]
      exact mul_le_mul_of_nonneg_right hb (by positivity)
    have h3 : B / m * Δ ^ 6 ≤ B / m * Δ ^ 4 := mul_le_mul_of_nonneg_left hΔ64 hBm
    linarith
  have hzE : |z| ≤ E * Δ ^ 2 := by
    have h1 : |z| ≤ |z - Δ ^ 2 * a| + |Δ ^ 2 * a| := by
      have := abs_add_le (z - Δ ^ 2 * a) (Δ ^ 2 * a)
      simpa using this
    have h2 : |Δ ^ 2 * a| ≤ A * Δ ^ 2 := by
      rw [abs_mul, abs_of_nonneg hΔ2, mul_comm]
      exact mul_le_mul_of_nonneg_right ha hΔ2
    have h3 : (A + B / m) * Δ ^ 4 ≤ (A + B / m) * Δ ^ 2 :=
      mul_le_mul_of_nonneg_left hΔ4 (by positivity)
    rw [hE]
    linarith
  have hz12 : |z| ≤ 1 / 2 := hzE.trans (by linarith)
  have hz1 : |-z| < 1 := by
    rw [abs_neg]
    linarith
  have h1z : 0 < 1 + z := by
    have := neg_abs_le z
    linarith
  have hupos : 0 < u := by
    rw [hug]
    positivity
  refine ⟨hupos, ?_⟩
  have hlog : |Real.log (1 + z) - z + z ^ 2 / 2| ≤ 2 * |z| ^ 3 := by
    have h := Real.abs_log_sub_add_sum_range_le hz1 2
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, abs_neg, sub_neg_eq_add] at h
    have e1 : Real.log (1 + z) - z + z ^ 2 / 2 = 0 + (-z) ^ (0 + 1) / ((0 : ℕ) + 1 : ℝ)
        + (-z) ^ (1 + 1) / ((1 : ℕ) + 1 : ℝ) + Real.log (1 + z) := by
      push_cast
      ring
    rw [e1]
    refine h.trans ?_
    rw [div_le_iff₀ (by linarith)]
    have h5 : 0 ≤ |z| ^ 3 * (1 - 2 * |z|) :=
      mul_nonneg (by positivity) (by linarith)
    have h6 : |z| ^ (2 + 1) = |z| ^ 3 := rfl
    rw [h6]
    nlinarith
  have hz3 : 2 * |z| ^ 3 ≤ 2 * E ^ 3 * Δ ^ 6 := by
    have h1 : |z| ^ 3 ≤ (E * Δ ^ 2) ^ 3 := pow_le_pow_left₀ (abs_nonneg _) hzE 3
    have h2 : (E * Δ ^ 2) ^ 3 = E ^ 3 * Δ ^ 6 := by ring
    rw [h2] at h1
    linarith
  have hsq : |z ^ 2 - Δ ^ 4 * a ^ 2| ≤ (A + B / m) * (3 * A + B / m) * Δ ^ 6 := by
    have h1 : z ^ 2 - Δ ^ 4 * a ^ 2 = (z - Δ ^ 2 * a) * (z + Δ ^ 2 * a) := by ring
    have h2 : |z + Δ ^ 2 * a| ≤ (3 * A + B / m) * Δ ^ 2 := by
      refine (abs_add_le _ _).trans ?_
      have h3 : |Δ ^ 2 * a| ≤ A * Δ ^ 2 := by
        rw [abs_mul, abs_of_nonneg hΔ2, mul_comm]
        exact mul_le_mul_of_nonneg_right ha hΔ2
      rw [hE] at hzE
      linarith
    rw [h1, abs_mul]
    calc |z - Δ ^ 2 * a| * |z + Δ ^ 2 * a|
        ≤ ((A + B / m) * Δ ^ 4) * ((3 * A + B / m) * Δ ^ 2) :=
          mul_le_mul hza h2 (abs_nonneg _) (by positivity)
      _ = (A + B / m) * (3 * A + B / m) * Δ ^ 6 := by ring
  have hlu : Real.log u - Real.log g = Real.log (1 + z) := by
    rw [hug, Real.log_mul hg0.ne' h1z.ne']
    ring
  have hsplit : Real.log u - Real.log g - Δ ^ 2 * a - Δ ^ 4 * (b - a ^ 2 / 2)
      = (Real.log (1 + z) - z + z ^ 2 / 2) + r - (z ^ 2 - Δ ^ 4 * a ^ 2) / 2 := by
    rw [hlu, hr]
    ring
  rw [hsplit]
  have h1 := abs_sub ((Real.log (1 + z) - z + z ^ 2 / 2) + r) ((z ^ 2 - Δ ^ 4 * a ^ 2) / 2)
  have h2 := abs_add_le (Real.log (1 + z) - z + z ^ 2 / 2) r
  have h3 : |(z ^ 2 - Δ ^ 4 * a ^ 2) / 2| = |z ^ 2 - Δ ^ 4 * a ^ 2| / 2 := by
    rw [abs_div, abs_of_pos (by norm_num : (0 : ℝ) < 2)]
  rw [h3] at h1
  nlinarith

/-- The exact part: the central difference of `log gauss` is the optimal rate, for every span. -/
lemma rate_exact {ν : ℝ} (hν : ν ≠ 0) {P : Matrix ι ι ℝ} (hP : P.IsSymm) (θ : ℝ) (x : ι → ℝ) (a : ι)
    {d : ℝ} (hd : d ≠ 0) :
    -2 * ν * (Real.log (gauss ν P θ (x + Pi.single a d))
      - Real.log (gauss ν P θ (x + Pi.single a (-d)))) / (2 * d) = 1 / 2 * (P *ᵥ x) a := by
  simp only [gauss, Real.log_exp, form_line hP]
  field_simp
  ring

omit [Fintype ι] in
/-- Uniform bound of a polynomial on the coordinate segments through `x₀`. -/
lemma exists_bound_line [Finite ι] (q : MvPolynomial ι ℝ) (x₀ : ι → ℝ) (d₀ : ℝ) :
    ∃ A, 0 ≤ A ∧ ∀ a, ∀ t ∈ Icc (-d₀) d₀, |eval (x₀ + Pi.single a t) q| ≤ A := by
  have := Fintype.ofFinite ι
  have h : ∀ a, ∃ A, ∀ t ∈ Icc (-d₀) d₀, ‖eval (x₀ + Pi.single a t) q‖ ≤ A := fun a =>
    isCompact_Icc.exists_bound_of_continuousOn fun t _ =>
      (hasDerivAt_eval_line q x₀ a t).continuousAt.continuousWithinAt
  choose A hA using h
  refine ⟨∑ a, |A a|, Finset.sum_nonneg fun _ _ => abs_nonneg _, fun a t ht => ?_⟩
  exact ((hA a t ht).trans (le_abs_self _)).trans
    (Finset.single_le_sum (f := fun a => |A a|) (fun _ _ => abs_nonneg _) (Finset.mem_univ a))

/-- Lower bound of the Gaussian on the coordinate segments through `x₀`. -/
lemma exists_gauss_lower {ν : ℝ} (hν : 0 < ν) {P : Matrix ι ι ℝ} (hP : P.IsSymm) (θ : ℝ)
    (x₀ : ι → ℝ) (d₀ : ℝ) :
    ∃ m, 0 < m ∧ ∀ a, ∀ t ∈ Icc (-d₀) d₀, m ≤ gauss ν P θ (x₀ + Pi.single a t) := by
  have h : ∀ a, ∃ M, ∀ t ∈ Icc (-d₀) d₀,
      ‖x₀ ⬝ᵥ (P *ᵥ x₀) + 2 * t * (P *ᵥ x₀) a + t ^ 2 * P a a‖ ≤ M := by
    intro a
    have hc : Continuous fun t : ℝ => x₀ ⬝ᵥ (P *ᵥ x₀) + 2 * t * (P *ᵥ x₀) a + t ^ 2 * P a a := by
      fun_prop
    exact isCompact_Icc.exists_bound_of_continuousOn hc.continuousOn
  choose M hM using h
  refine ⟨Real.exp (-(∑ a, |M a|) / (8 * ν) - θ), Real.exp_pos _, fun a t ht => ?_⟩
  rw [gauss, form_line hP, Real.exp_le_exp]
  have h1 : x₀ ⬝ᵥ (P *ᵥ x₀) + 2 * t * (P *ᵥ x₀) a + t ^ 2 * P a a ≤ ∑ a, |M a| :=
    (Real.le_norm_self _).trans ((hM a t ht).trans ((le_abs_self _).trans
      (Finset.single_le_sum (f := fun a => |M a|) (fun _ _ => abs_nonneg _) (Finset.mem_univ a))))
  have h2 : -(∑ a, |M a|) / (8 * ν)
      ≤ -(x₀ ⬝ᵥ (P *ᵥ x₀) + 2 * t * (P *ᵥ x₀) a + t ^ 2 * P a a) / (8 * ν) :=
    div_le_div_of_nonneg_right (neg_le_neg h1) (by positivity)
  linarith

/-- A function on a symmetric segment with derivative bounded by `L`. -/
lemma abs_sub_le_of_deriv {f f' : ℝ → ℝ} {d₀ L : ℝ} (hf : ∀ t, HasDerivAt f (f' t) t)
    (hL : ∀ t ∈ Icc (-d₀) d₀, |f' t| ≤ L) {d : ℝ} (hd : 0 ≤ d) (hd₀ : d ≤ d₀) :
    |f d - f (-d)| ≤ L * (2 * d) := by
  have hmv := (convex_Icc (-d₀) d₀).norm_image_sub_le_of_norm_hasDerivWithin_le
    (f := f) (f' := f') (C := L) (fun t _ => (hf t).hasDerivWithinAt)
    (fun t ht => hL t ht)
    (⟨by linarith, by linarith⟩ : -d ∈ Icc (-d₀) d₀) (⟨by linarith, hd₀⟩ : d ∈ Icc (-d₀) d₀)
  simp only [Real.norm_eq_abs] at hmv
  have h2 : |d - -d| = 2 * d := by
    rw [abs_of_nonneg (by linarith)]
    ring
  rwa [h2] at hmv

/-- Lipschitz bound of a polynomial along the coordinate segments through `x₀`. -/
lemma exists_lipschitz_line (q : MvPolynomial ι ℝ) (x₀ : ι → ℝ) (d₀ : ℝ) :
    ∃ L : ℝ, 0 ≤ L ∧ ∀ (a : ι) (d : ℝ), 0 ≤ d → d ≤ d₀ →
      |eval (x₀ + Pi.single a d) q - eval (x₀ + Pi.single a (-d)) q| ≤ L * (2 * d) := by
  have h : ∃ L : ι → ℝ, ∀ a : ι, 0 ≤ L a ∧ ∀ t ∈ Icc (-d₀) d₀,
      |eval (x₀ + Pi.single a t) (pderiv a q)| ≤ L a := by
    refine ⟨fun a => (exists_bound_line (pderiv a q) x₀ d₀).choose, fun a => ?_⟩
    exact ⟨(exists_bound_line (pderiv a q) x₀ d₀).choose_spec.1,
      (exists_bound_line (pderiv a q) x₀ d₀).choose_spec.2 a⟩
  obtain ⟨L, hL⟩ := h
  refine ⟨∑ a, L a, Finset.sum_nonneg fun a _ => (hL a).1, ?_⟩
  intro a d hd hd₀
  have hle : L a ≤ ∑ a, L a :=
    Finset.single_le_sum (f := L) (fun a _ => (hL a).1) (Finset.mem_univ a)
  obtain ⟨f, hf⟩ : ∃ f : ℝ → ℝ, f = fun t : ℝ => eval (x₀ + Pi.single a t) q := ⟨_, rfl⟩
  have hfd : ∀ t, HasDerivAt f (eval (x₀ + Pi.single a t) (pderiv a q)) t := by
    intro t
    rw [hf]
    exact hasDerivAt_eval_line q x₀ a t
  have key : |f d - f (-d)| ≤ L a * (2 * d) := abs_sub_le_of_deriv hfd (hL a).2 hd hd₀
  have key' : |f d - f (-d)| ≤ (∑ a, L a) * (2 * d) :=
    key.trans (mul_le_mul_of_nonneg_right hle (by linarith))
  subst hf
  exact key'

omit [Fintype ι] in
lemma pos_single (x₀ : ι → ℝ) (Δ : ℝ) (a : ι) (e : ℤ) :
    pos x₀ Δ (Pi.single a e) = x₀ + Pi.single a (e * Δ) := by
  funext b
  by_cases hb : b = a
  · subst hb
    simp [pos, mul_comm]
  · simp [pos, hb]

/-- Final algebra of the rate expansion. -/
lemma rate_algebra {ν d Δ C₀ L lp lm gp gm cp cm kp km v : ℝ} (hν : 0 < ν) (hΔ : 0 < Δ)
    (hΔ1 : Δ ≤ 1) (hd : Δ ≤ d) (hC₀ : 0 ≤ C₀) (hex : -2 * ν * (gp - gm) / (2 * d) = v)
    (hp : |lp - gp - Δ ^ 2 * cp - Δ ^ 4 * kp| ≤ C₀ * Δ ^ 6)
    (hm : |lm - gm - Δ ^ 2 * cm - Δ ^ 4 * km| ≤ C₀ * Δ ^ 6) (hk : |kp - km| ≤ L * (2 * d)) :
    |-2 * ν * (lp - lm) / (2 * d) - (v + -2 * ν * Δ ^ 2 * (cp - cm) / (2 * d))|
      ≤ (ν * (2 * L) + 2 * ν * C₀) * Δ ^ 4 := by
  have hd0 : 0 < d := hΔ.trans_le hd
  have hid : -2 * ν * (lp - lm) / (2 * d) - (v + -2 * ν * Δ ^ 2 * (cp - cm) / (2 * d))
      = -(ν / d) * (Δ ^ 4 * (kp - km) + ((lp - gp - Δ ^ 2 * cp - Δ ^ 4 * kp)
        - (lm - gm - Δ ^ 2 * cm - Δ ^ 4 * km))) := by
    rw [← hex]
    field_simp
    ring
  rw [hid, abs_mul, abs_neg, abs_of_pos (by positivity : 0 < ν / d)]
  have hΔ4 : 0 ≤ Δ ^ 4 := by positivity
  have h1 : |Δ ^ 4 * (kp - km) + ((lp - gp - Δ ^ 2 * cp - Δ ^ 4 * kp)
      - (lm - gm - Δ ^ 2 * cm - Δ ^ 4 * km))| ≤ Δ ^ 4 * (L * (2 * d)) + 2 * (C₀ * Δ ^ 6) := by
    refine (abs_add_le _ _).trans (add_le_add ?_ ?_)
    · rw [abs_mul, abs_of_nonneg hΔ4]
      exact mul_le_mul_of_nonneg_left hk hΔ4
    · refine (abs_sub _ _).trans ?_
      linarith
  have h2 : Δ ^ 6 ≤ Δ ^ 4 * d := by
    have : Δ ^ 2 ≤ d := by nlinarith
    calc Δ ^ 6 = Δ ^ 4 * Δ ^ 2 := by ring
      _ ≤ Δ ^ 4 * d := mul_le_mul_of_nonneg_left this hΔ4
  calc ν / d * _ ≤ ν / d * (Δ ^ 4 * (L * (2 * d)) + 2 * (C₀ * Δ ^ 6)) :=
        mul_le_mul_of_nonneg_left h1 (by positivity)
    _ ≤ ν / d * (Δ ^ 4 * (L * (2 * d)) + 2 * (C₀ * (Δ ^ 4 * d))) := by
        have := mul_le_mul_of_nonneg_left h2 hC₀
        refine mul_le_mul_of_nonneg_left ?_ (by positivity)
        linarith
    _ = (ν * (2 * L) + 2 * ν * C₀) * Δ ^ 4 := by
        field_simp

/-- **Rate expansion.** The rate read from the lattice solution with span `h` (so `d = h Δ`) equals the
optimal rate plus the `χ`-shift, up to `O(Δ⁴)` uniformly for `Δ ≤ d ≤ d₀`. -/
theorem rate_expansion {ν lam T φ : ℝ} (hν : 0 < ν) (hlam : 0 ≤ lam) (hT : 0 < T)
    {C Φ : Matrix ι ι ℝ} (hCs : C.IsSymm) (hΦs : Φ.IsSymm) (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v))
    (hφ : 0 < φ) (hΦ : ∀ v, φ * (v ⬝ᵥ v) ≤ v ⬝ᵥ (Φ *ᵥ v)) {Y Z : ℝ → Matrix ι ι ℝ}
    (hYZ : RiccatiFlow C Φ lam T Y Z) {θ : ℝ → ℝ}
    (hθ : ∀ s ∈ Icc 0 T,
      HasDerivWithinAt θ (1 / 4 * (RiccatiFlow.P Y Z s).trace) (Icc 0 T) s)
    {p : ℕ → ℝ → MvPolynomial ι ℝ} (hp : IsHierarchy T ν (RiccatiFlow.P Y Z) 2 p)
    (x₀ : ι → ℝ) (d₀ : ℝ) :
    ∃ Δ₀ K₀, 0 < Δ₀ ∧ ∀ Δ, 0 < Δ → Δ ≤ Δ₀ → ∀ h : ℕ, 1 ≤ h → h * Δ ≤ d₀ →
      ∀ u : (ι → ℤ) → ℝ → ℝ, IsLatticeSolution ν Δ T x₀ (pot ν lam C) u →
      (∀ n, u n 0 = gauss ν (RiccatiFlow.P Y Z 0) (θ 0) (pos x₀ Δ n)) → ∀ a,
        |-2 * ν * (Real.log (u (Pi.single a (h : ℤ)) T) - Real.log (u (Pi.single a (-(h : ℤ))) T))
            / (2 * (h * Δ))
          - (1 / 2 * (RiccatiFlow.P Y Z T *ᵥ x₀) a
            + -2 * ν * Δ ^ 2 * (eval (x₀ + Pi.single a (h * Δ)) (p 1 T)
              - eval (x₀ + Pi.single a (-(h * Δ))) (p 1 T)) / (2 * (h * Δ)))| ≤ K₀ * Δ ^ 4 := by
  obtain ⟨B, hB, hmain⟩ := grid_theorem hν hlam hT hCs hΦs hC hφ hΦ hYZ hθ 2 hp
  have hTmem : T ∈ Icc 0 T := ⟨hT.le, le_rfl⟩
  have hPs : (RiccatiFlow.P Y Z T).IsSymm := hYZ.P_symm hlam hCs hΦs hC hφ hΦ hTmem
  generalize hPdef : RiccatiFlow.P Y Z T = P at *
  obtain ⟨κ, hκ⟩ : ∃ κ : MvPolynomial ι ℝ,
      ∀ x, eval x κ = eval x (p 2 T) - (eval x (p 1 T)) ^ 2 / 2 :=
    ⟨p 2 T - (1 / 2 : ℝ) • (p 1 T) ^ 2, fun x => by
      rw [map_sub, smul_eval, map_pow]
      ring⟩
  obtain ⟨A₁, hA₁, hb₁⟩ := exists_bound_line (p 1 T) x₀ d₀
  obtain ⟨A₂, hA₂, hb₂⟩ := exists_bound_line (p 2 T) x₀ d₀
  obtain ⟨L, hL, hLip⟩ := exists_lipschitz_line κ x₀ d₀
  obtain ⟨m, hm, hgm⟩ := exists_gauss_lower hν hPs (θ T) x₀ d₀
  obtain ⟨A, hAeq⟩ : ∃ A, A = A₁ + A₂ := ⟨_, rfl⟩
  have hA0 : 0 ≤ A := by linarith
  have hBm : 0 ≤ B / m := by positivity
  obtain ⟨C₀, hC₀⟩ : ∃ C₀, C₀ = B / m + (A + B / m) * (3 * A + B / m) / 2
    + 2 * (2 * A + B / m) ^ 3 := ⟨_, rfl⟩
  have hC₀nn : 0 ≤ C₀ := by
    rw [hC₀]
    positivity
  have hE0 : 0 ≤ 2 * A + B / m := by positivity
  refine ⟨min 1 (1 / (2 * (2 * A + B / m + 1))), ν * (2 * L) + 2 * ν * C₀,
    lt_min one_pos (by positivity), fun Δ hΔ hΔle h hh hd u hu hu0 a => ?_⟩
  have hΔ1 : Δ ≤ 1 := hΔle.trans (min_le_left _ _)
  have hΔE : Δ ≤ 1 / (2 * (2 * A + B / m + 1)) := hΔle.trans (min_le_right _ _)
  have hsmall : Δ ^ 2 * (2 * A + B / m) ≤ 1 / 2 := by
    have h1 : Δ ^ 2 ≤ Δ := by nlinarith
    rw [le_div_iff₀ (by positivity)] at hΔE
    nlinarith
  have hh' : (1 : ℝ) ≤ h := by exact_mod_cast hh
  have hdΔ : Δ ≤ h * Δ := by nlinarith
  have hdpos : 0 < (h : ℝ) * Δ := hΔ.trans_le hdΔ
  have hlogs : ∀ (e : ℤ), ((e : ℝ) * Δ) ∈ Icc (-d₀) d₀ →
      |Real.log (u (Pi.single a e) T) - Real.log (gauss ν P (θ T) (x₀ + Pi.single a (e * Δ)))
        - Δ ^ 2 * eval (x₀ + Pi.single a (e * Δ)) (p 1 T)
        - Δ ^ 4 * eval (x₀ + Pi.single a (e * Δ)) κ| ≤ C₀ * Δ ^ 6 := by
    intro e he
    rw [hκ, hC₀]
    refine (log_expansion hm (hgm a _ he) ((hb₁ a _ he).trans (by linarith))
      ((hb₂ a _ he).trans (by linarith)) hB hΔ.le hΔ1 hsmall ?_).2
    have hk := hmain Δ hΔ x₀ u hu hu0 (Pi.single a e) T hTmem
    rw [pos_single, hPdef] at hk
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, hp.zero, map_one] at hk
    rw [mul_comm B]
    refine le_trans (le_of_eq ?_) hk
    congr 1
    ring
  have hp1 := hlogs (h : ℤ) (by
    push_cast
    exact ⟨by linarith, hd⟩)
  have hm1 := hlogs (-(h : ℤ)) (by
    push_cast
    rw [neg_mul]
    exact ⟨by linarith, by linarith⟩)
  push_cast at hp1 hm1
  rw [neg_mul] at hm1
  exact rate_algebra hν hΔ hΔ1 hdΔ hC₀nn (rate_exact hν.ne' hPs (θ T) x₀ a hdpos.ne') hp1 hm1
    (hLip a _ hdpos.le hd)

end ExecutionProofs
