import ExecutionProofs.Feedback

/-!
# Recentring on a classical characteristic (paper, Proposition 10)

`ψ(x, s) = exp (α(s)·x + γ(s)) φ(x - ξ(s), s)` turns `ψ_s = ν ∇·(A∇ψ) - V ψ` into the same equation
for `φ` with the potential `V(ξ + y) - V(ξ) - G(ξ)·y`, provided `ξ' = -2ν A α`, `α' = -G(ξ)` and
`γ' = ν α·Aα + G(ξ)·ξ - V(ξ)`. The generator stays symmetric. For `G = ∇V` this is the imaginary-time
displaced frame; the identity itself holds for any vector field `G`.
-/

open Set Matrix

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

lemma hasDerivAt_add_single (x : ι → ℝ) (b : ι) :
    HasDerivAt (fun t : ℝ => x + Pi.single b t) (Pi.single b 1) 0 := by
  refine (hasDerivAt_const (0 : ℝ) x).add ?_ |>.congr_deriv (zero_add _)
  rw [hasDerivAt_pi]
  intro i
  by_cases h : i = b
  · subst h
    simpa using hasDerivAt_id' (0 : ℝ)
  · simpa [Pi.single_apply, h] using hasDerivAt_const (0 : ℝ) (0 : ℝ)

/-- The partial derivative of a Fréchet differentiable function. -/
lemma hasDerivAt_single {f : (ι → ℝ) → ℝ} {L : (ι → ℝ) →L[ℝ] ℝ} {x : ι → ℝ}
    (h : HasFDerivAt f L x) (b : ι) :
    HasDerivAt (fun t : ℝ => f (x + Pi.single b t)) (L (Pi.single b 1)) 0 := by
  have hx : x + Pi.single b (0 : ℝ) = x := by simp
  exact HasFDerivAt.comp_hasDerivAt (0 : ℝ) (by rwa [hx]) (hasDerivAt_add_single x b)

lemma dd_eq {f : (ι → ℝ) → ℝ} {L : (ι → ℝ) →L[ℝ] ℝ} {x : ι → ℝ} (h : HasFDerivAt f L x) (b : ι) :
    dd b f x = L (Pi.single b 1) :=
  (hasDerivAt_single h b).deriv

lemma dd_hasDerivAt {f : (ι → ℝ) → ℝ} {x : ι → ℝ} (h : DifferentiableAt ℝ f x) (b : ι) :
    HasDerivAt (fun t : ℝ => f (x + Pi.single b t)) (dd b f x) 0 := by
  rw [dd_eq h.hasFDerivAt]
  exact hasDerivAt_single h.hasFDerivAt b

/-- A linear functional is determined by the partial derivatives. -/
lemma clm_apply_eq_sum {f : (ι → ℝ) → ℝ} {L : (ι → ℝ) →L[ℝ] ℝ} {x : ι → ℝ}
    (h : HasFDerivAt f L x) (v : ι → ℝ) : L v = ∑ b, v b * dd b f x := by
  conv_lhs => rw [pi_eq_sum_univ' v]
  rw [map_sum]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [dd_eq h, map_smul, smul_eq_mul]

lemma dotProduct_add_single (α x : ι → ℝ) (b : ι) (t : ℝ) :
    α ⬝ᵥ (x + Pi.single b t) = α b * t + α ⬝ᵥ x := by
  rw [dotProduct_add, dotProduct_single, add_comm]

/-- First derivative of the exponential gauge. -/
lemma hasDerivAt_gauge (α x : ι → ℝ) (γ : ℝ) (b : ι) :
    HasDerivAt (fun t : ℝ => Real.exp (α ⬝ᵥ (x + Pi.single b t) + γ))
      (α b * Real.exp (α ⬝ᵥ x + γ)) 0 := by
  have hlin : ∀ t : ℝ, α ⬝ᵥ (x + Pi.single b t) + γ = α b * t + (α ⬝ᵥ x + γ) := fun t => by
    rw [dotProduct_add_single, add_assoc]
  simp only [hlin]
  refine ((((hasDerivAt_id (0 : ℝ)).const_mul (α b)).add_const (α ⬝ᵥ x + γ)).exp).congr_deriv ?_
  simp
  ring

/-- First partial derivatives of the recentred function. -/
lemma dd_gauge {f : (ι → ℝ) → ℝ} (hf : Differentiable ℝ f) (α ξ : ι → ℝ) (γ : ℝ) (b : ι) :
    dd b (fun x => Real.exp (α ⬝ᵥ x + γ) * f (x - ξ))
      = fun x => Real.exp (α ⬝ᵥ x + γ) * (α b * f (x - ξ) + dd b f (x - ξ)) := by
  funext x
  have h := (hasDerivAt_gauge α x γ b).fun_mul (dd_hasDerivAt (hf (x - ξ)) b)
  have hshift : ∀ t : ℝ, x + Pi.single b t - ξ = x - ξ + Pi.single b t := fun t => by abel
  unfold dd
  simp only [hshift]
  rw [HasDerivAt.deriv h]
  simp only [Pi.single_zero, add_zero]
  unfold dd
  ring

lemma differentiable_dd {f : (ι → ℝ) → ℝ} (hf : ContDiff ℝ 2 f) (b : ι) :
    Differentiable ℝ (dd b f) := by
  have hd : Differentiable ℝ f := hf.differentiable (by norm_num)
  have heq : dd b f = fun x => fderiv ℝ f x (Pi.single b 1) := by
    funext x
    exact dd_eq (hd x).hasFDerivAt b
  rw [heq]
  have h1 : ContDiff ℝ 1 (fderiv ℝ f) := hf.fderiv_right (by norm_num)
  exact (h1.differentiable (by norm_num)).clm_apply (differentiable_const _)

/-- Second partial derivatives of the recentred function. -/
lemma dd_dd_gauge {f : (ι → ℝ) → ℝ} (hf : ContDiff ℝ 2 f) (α ξ : ι → ℝ) (γ : ℝ) (a b : ι)
    (x : ι → ℝ) :
    dd a (dd b fun x => Real.exp (α ⬝ᵥ x + γ) * f (x - ξ)) x
      = Real.exp (α ⬝ᵥ x + γ) * (α a * α b * f (x - ξ) + α b * dd a f (x - ξ)
        + α a * dd b f (x - ξ) + dd a (dd b f) (x - ξ)) := by
  have hd : Differentiable ℝ f := hf.differentiable (by norm_num)
  have hdb := differentiable_dd hf b
  have hg : Differentiable ℝ fun y => α b * f y + dd b f y := (hd.const_mul _).add hdb
  have hsum : dd a (fun y => α b * f y + dd b f y) (x - ξ)
      = α b * dd a f (x - ξ) + dd a (dd b f) (x - ξ) :=
    (((dd_hasDerivAt (hd (x - ξ)) a).const_mul (α b)).add (dd_hasDerivAt (hdb (x - ξ)) a)).deriv
  rw [dd_gauge hd, dd_gauge (f := fun y => α b * f y + dd b f y) hg]
  simp only [hsum]
  ring

/-- **Recentring.** -/
theorem recentre {ν s : ℝ} {I : Set ℝ} (hs : s ∈ I) {A : Matrix ι ι ℝ} (hA : A.IsSymm)
    {V : (ι → ℝ) → ℝ} {G : (ι → ℝ) → ι → ℝ} {Φ : (ι → ℝ) × ℝ → ℝ} {Φ' : (ι → ℝ) × ℝ →L[ℝ] ℝ}
    {ξ α : ℝ → ι → ℝ} {γ : ℝ → ℝ} (hξ : HasDerivWithinAt ξ (-(2 * ν) • (A *ᵥ α s)) I s)
    (hα : HasDerivWithinAt α (-G (ξ s)) I s)
    (hγ : HasDerivWithinAt γ (ν * (α s ⬝ᵥ (A *ᵥ α s)) + G (ξ s) ⬝ᵥ ξ s - V (ξ s)) I s)
    (hC : ContDiff ℝ 2 fun y => Φ (y, s)) (x : ι → ℝ)
    (hΦ : HasFDerivWithinAt Φ Φ' (univ ×ˢ I) (x - ξ s, s))
    (heq : Φ' (0, 1) = ν * ∑ a, ∑ b, A a b * dd a (dd b fun y => Φ (y, s)) (x - ξ s)
      - (V x - V (ξ s) - G (ξ s) ⬝ᵥ (x - ξ s)) * Φ (x - ξ s, s)) :
    HasDerivWithinAt (fun r => Real.exp (α r ⬝ᵥ x + γ r) * Φ (x - ξ r, r))
      (ν * ∑ a, ∑ b, A a b * dd a (dd b fun x => Real.exp (α s ⬝ᵥ x + γ s) * Φ (x - ξ s, s)) x
        - V x * (Real.exp (α s ⬝ᵥ x + γ s) * Φ (x - ξ s, s))) I s := by
  set f : (ι → ℝ) → ℝ := fun y => Φ (y, s) with hfdef
  set E := Real.exp (α s ⬝ᵥ x + γ s) with hE
  -- the spatial derivative of `Φ` at fixed `s`
  have hsp : HasFDerivAt f (Φ'.comp (ContinuousLinearMap.inl ℝ (ι → ℝ) ℝ)) (x - ξ s) := by
    have hin : HasFDerivWithinAt (fun y : ι → ℝ => (y, s))
        (ContinuousLinearMap.inl ℝ (ι → ℝ) ℝ) univ (x - ξ s) :=
      (hasFDerivAt_prodMk_left (x - ξ s) s).hasFDerivWithinAt
    have hmaps : MapsTo (fun y : ι → ℝ => (y, s)) univ (univ ×ˢ I) := fun y _ => ⟨mem_univ _, hs⟩
    exact hasFDerivWithinAt_univ.mp (hΦ.comp (x - ξ s) hin hmaps)
  have hgrad : ∀ v : ι → ℝ, Φ' (v, 0) = ∑ b, v b * dd b f (x - ξ s) := fun v => by
    have h := clm_apply_eq_sum hsp v
    simpa using h
  -- time derivative of the moving argument
  have hpath : HasDerivWithinAt (fun r => (x - ξ r, r)) ((2 * ν) • (A *ᵥ α s), (1 : ℝ)) I s := by
    have h1 : HasDerivWithinAt (fun r => x - ξ r) ((2 * ν) • (A *ᵥ α s)) I s :=
      (hξ.const_sub x).congr_deriv (by simp)
    exact h1.prodMk (hasDerivWithinAt_id s I)
  have hmove : HasDerivWithinAt (fun r => Φ (x - ξ r, r))
      (Φ' ((2 * ν) • (A *ᵥ α s), 1)) I s :=
    HasFDerivWithinAt.comp_hasDerivWithinAt s hΦ hpath fun r hr => ⟨mem_univ _, hr⟩
  have hsplit : Φ' ((2 * ν) • (A *ᵥ α s), 1)
      = 2 * ν * ∑ b, (A *ᵥ α s) b * dd b f (x - ξ s) + Φ' (0, 1) := by
    have hsum : (((2 * ν) • (A *ᵥ α s), 1) : (ι → ℝ) × ℝ)
        = ((2 * ν) • (A *ᵥ α s), 0) + (0, 1) := by simp
    rw [hsum, map_add, hgrad, Finset.mul_sum]
    congr 1
    refine Finset.sum_congr rfl fun b _ => ?_
    simp only [Pi.smul_apply, smul_eq_mul]
    ring
  -- time derivative of the gauge
  have hdot : HasDerivWithinAt (fun r => α r ⬝ᵥ x) (-G (ξ s) ⬝ᵥ x) I s := by
    unfold dotProduct
    refine HasDerivWithinAt.fun_sum fun i _ => ?_
    exact ((hasDerivWithinAt_pi.mp hα) i).mul_const (x i)
  have hexp : HasDerivWithinAt (fun r => Real.exp (α r ⬝ᵥ x + γ r))
      (E * (-G (ξ s) ⬝ᵥ x + (ν * (α s ⬝ᵥ (A *ᵥ α s)) + G (ξ s) ⬝ᵥ ξ s - V (ξ s)))) I s :=
    (hdot.add hγ).exp
  refine (hexp.mul hmove).congr_deriv ?_
  -- the space derivatives
  have hAα : ∀ b, (A *ᵥ α s) b = ∑ a, A a b * α s a := fun b => by
    simp only [Matrix.mulVec, dotProduct]
    exact Finset.sum_congr rfl fun a _ => by rw [hA.apply a b]
  have hquad : ∑ a, ∑ b, A a b * (α s a * α s b) = α s ⬝ᵥ (A *ᵥ α s) := by
    simp only [Matrix.mulVec, dotProduct, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  have hcross1 : ∑ a, ∑ b, A a b * (α s b * dd a f (x - ξ s))
      = ∑ b, (A *ᵥ α s) b * dd b f (x - ξ s) := by
    refine Finset.sum_congr rfl fun a _ => ?_
    simp only [Matrix.mulVec, dotProduct, Finset.sum_mul]
    exact Finset.sum_congr rfl fun b _ => by ring
  have hcross2 : ∑ a, ∑ b, A a b * (α s a * dd b f (x - ξ s))
      = ∑ b, (A *ᵥ α s) b * dd b f (x - ξ s) := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [hAα, Finset.sum_mul]
    exact Finset.sum_congr rfl fun a _ => by ring
  have hspace : ∑ a, ∑ b, A a b * dd a (dd b fun x => Real.exp (α s ⬝ᵥ x + γ s) * f (x - ξ s)) x
      = E * ((α s ⬝ᵥ (A *ᵥ α s)) * f (x - ξ s) + 2 * ∑ b, (A *ᵥ α s) b * dd b f (x - ξ s)
        + ∑ a, ∑ b, A a b * dd a (dd b f) (x - ξ s)) := by
    rw [← hquad, two_mul]
    nth_rewrite 1 [← hcross1]
    rw [← hcross2, Finset.sum_mul, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib,
      ← Finset.sum_add_distrib, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Finset.sum_mul, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib,
      ← Finset.sum_add_distrib, Finset.mul_sum]
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [dd_dd_gauge hC]
    ring
  change _ = ν * ∑ a, ∑ b, A a b * dd a (dd b fun x => Real.exp (α s ⬝ᵥ x + γ s) * f (x - ξ s)) x
    - V x * (E * f (x - ξ s))
  rw [hspace, hsplit, heq, dotProduct_sub, neg_dotProduct]
  ring

/-- **Recentred feedback**: the reference control plus the logarithmic gradient of `φ`. -/
theorem recentred_feedback {ν : ℝ} (A : Matrix ι ι ℝ) {f : (ι → ℝ) → ℝ} (hf : Differentiable ℝ f)
    (hpos : ∀ y, 0 < f y) (α ξ : ι → ℝ) (γ : ℝ) (x : ι → ℝ) (a : ι) :
    -2 * ν * ∑ b, A a b * dd b (fun x => Real.log (Real.exp (α ⬝ᵥ x + γ) * f (x - ξ))) x
      = -2 * ν * (A *ᵥ α) a + -2 * ν * ∑ b, A a b * dd b (fun y => Real.log (f y)) (x - ξ) := by
  have hlog : ∀ b, dd b (fun x => Real.log (Real.exp (α ⬝ᵥ x + γ) * f (x - ξ))) x
      = α b + dd b (fun y => Real.log (f y)) (x - ξ) := by
    intro b
    have hsplit : (fun x => Real.log (Real.exp (α ⬝ᵥ x + γ) * f (x - ξ)))
        = fun x => (α ⬝ᵥ x + γ) + Real.log (f (x - ξ)) := by
      funext x
      rw [Real.log_mul (Real.exp_pos _).ne' (hpos _).ne', Real.log_exp]
    have hlin : HasDerivAt (fun t : ℝ => α ⬝ᵥ (x + Pi.single b t) + γ) (α b) 0 := by
      have h : ∀ t : ℝ, α ⬝ᵥ (x + Pi.single b t) + γ = α b * t + (α ⬝ᵥ x + γ) := fun t => by
        rw [dotProduct_add_single, add_assoc]
      simp only [h]
      simpa using ((hasDerivAt_id (0 : ℝ)).const_mul (α b)).add_const (α ⬝ᵥ x + γ)
    have hd : Differentiable ℝ fun y => Real.log (f y) := hf.log fun y => (hpos y).ne'
    have h := hlin.fun_add (dd_hasDerivAt (hd (x - ξ)) b)
    have hshift : ∀ t : ℝ, x + Pi.single b t - ξ = x - ξ + Pi.single b t := fun t => by abel
    rw [hsplit, ← h.deriv]
    conv_lhs => unfold dd
    simp only [hshift]
  simp only [hlog, mul_add, Finset.sum_add_distrib]
  rfl

end ExecutionProofs
