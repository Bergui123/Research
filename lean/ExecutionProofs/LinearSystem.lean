import Mathlib

/-!
# Existence of bounded solutions of linear systems

`u n' = ∑ d, κ d n s * u (nb d n) - V n * u n + b n s` on a possibly infinite index set `S`, with
bounded rates and an arbitrary (unbounded) nonnegative damping `V`. The mild formulation with an
exponential weight is a contraction. This covers the infinite lattice and, for finite `S`,
linear ordinary differential equations.
-/

open Set MeasureTheory intervalIntegral BoundedContinuousFunction

namespace ExecutionProofs

/-- Damped integral of a function bounded by `C`. -/
lemma abs_damped_integral_le {a s C : ℝ} (ha : 0 < a) (hs : 0 ≤ s) {f : ℝ → ℝ}
    (hf : Continuous f) (hC : ∀ r ∈ Icc 0 s, |f r| ≤ C) :
    |Real.exp (-a * s) * ∫ r in (0 : ℝ)..s, Real.exp (a * r) * f r| ≤ C / a := by
  have hC0 : 0 ≤ C := (abs_nonneg _).trans (hC 0 ⟨le_rfl, hs⟩)
  have hcont : Continuous fun r : ℝ => Real.exp (a * r) := by fun_prop
  have hint : ∫ r in (0 : ℝ)..s, Real.exp (a * r) = (Real.exp (a * s) - 1) / a := by
    have hd : ∀ r ∈ uIcc 0 s, HasDerivAt (fun r => Real.exp (a * r) / a) (Real.exp (a * r)) r := by
      intro r _
      have h := (((hasDerivAt_id r).const_mul a).exp).div_const a
      refine h.congr_deriv ?_
      simp only [id, mul_one]
      field_simp
    rw [integral_eq_sub_of_hasDerivAt hd (hcont.intervalIntegrable _ _)]
    simp only [mul_zero, Real.exp_zero]
    ring
  have h1 : |∫ r in (0 : ℝ)..s, Real.exp (a * r) * f r|
      ≤ ∫ r in (0 : ℝ)..s, Real.exp (a * r) * C := by
    refine (abs_integral_le_integral_abs hs).trans ?_
    refine integral_mono_on hs ?_ ?_ ?_
    · exact ((hcont.mul hf).abs).intervalIntegrable _ _
    · exact (hcont.mul continuous_const).intervalIntegrable _ _
    · intro r hr
      rw [abs_mul, abs_of_pos (Real.exp_pos _)]
      exact mul_le_mul_of_nonneg_left (hC r hr) (Real.exp_pos _).le
  rw [intervalIntegral.integral_mul_const, hint] at h1
  rw [abs_mul, abs_of_pos (Real.exp_pos _)]
  have hexp : Real.exp (-a * s) * Real.exp (a * s) = 1 := by
    rw [← Real.exp_add]
    simp
  calc Real.exp (-a * s) * |∫ r in (0 : ℝ)..s, Real.exp (a * r) * f r|
      ≤ Real.exp (-a * s) * ((Real.exp (a * s) - 1) / a * C) :=
        mul_le_mul_of_nonneg_left h1 (Real.exp_pos _).le
    _ = (1 - Real.exp (-a * s)) * (C / a) := by
        linear_combination (C / a) * hexp
    _ ≤ 1 * (C / a) := by
        refine mul_le_mul_of_nonneg_right ?_ (by positivity)
        linarith [Real.exp_pos (-a * s)]
    _ = C / a := one_mul _

/-- **Existence** for continuous coefficients on the whole line. -/
theorem exists_linear_solution_of_continuous {S D : Type*} [Fintype D] (nb : D → S → S)
    {κ : D → S → ℝ → ℝ} (hκ : ∀ d n, Continuous (κ d n)) {Kb : ℝ}
    (hκb : ∀ d n s, |κ d n s| ≤ Kb) {V : S → ℝ} (hV : ∀ n, 0 ≤ V n) {b : S → ℝ → ℝ}
    (hb : ∀ n, Continuous (b n)) {β : ℝ} (hbb : ∀ n s, |b n s| ≤ β) {u₀ : S → ℝ} {A : ℝ}
    (hA : ∀ n, |u₀ n| ≤ A) {T : ℝ} (hT : 0 ≤ T) :
    ∃ u : S → ℝ → ℝ, (∀ n, u n 0 = u₀ n) ∧ (∃ B, ∀ n, ∀ s ∈ Icc 0 T, |u n s| ≤ B) ∧
      ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n)
        (∑ d, κ d n s * u (nb d n) s - V n * u n s + b n s) (Icc 0 T) s := by
  classical
  rcases isEmpty_or_nonempty S with hS | ⟨⟨n₀⟩⟩
  · exact ⟨fun n => (hS.false n).elim, fun n => (hS.false n).elim,
      ⟨0, fun n => (hS.false n).elim⟩, fun n => (hS.false n).elim⟩
  have hA0 : 0 ≤ A := (abs_nonneg _).trans (hA n₀)
  have hβ0 : 0 ≤ β := (abs_nonneg _).trans (hbb n₀ 0)
  set Q : ℝ := Fintype.card D * |Kb| with hQ
  have hQ0 : 0 ≤ Q := by positivity
  set μ : ℝ := Q + 1 with hμ
  have hμ0 : 0 < μ := by linarith
  let _ : TopologicalSpace S := ⊥
  have : DiscreteTopology S := ⟨rfl⟩
  let X := BoundedContinuousFunction (S × Icc (0 : ℝ) T) ℝ
  let ext : X → S → ℝ → ℝ := fun v n r => v (n, projIcc 0 T hT r)
  have hext : ∀ v n, Continuous (ext v n) := fun v n =>
    v.continuous.comp (Continuous.prodMk continuous_const continuous_projIcc)
  have hext_le : ∀ (v : X) n r, |ext v n r| ≤ ‖v‖ := fun v n r => v.norm_coe_le_norm _
  let F : X → S → ℝ → ℝ := fun v n r => ∑ d, κ d n r * ext v (nb d n) r + Real.exp (-μ * r) * b n r
  have hF : ∀ v n, Continuous (F v n) := fun v n =>
    (continuous_finsetSum _ fun d _ => (hκ d n).mul (hext v _)).add
      ((by fun_prop : Continuous fun r : ℝ => Real.exp (-μ * r)).mul (hb n))
  have hsum_le : ∀ (w : S → ℝ) (R : ℝ) n r, (∀ m, |w m| ≤ R) →
      |∑ d, κ d n r * w (nb d n)| ≤ Q * R := by
    intro w R n r hw
    have hR : 0 ≤ R := (abs_nonneg _).trans (hw n)
    calc |∑ d, κ d n r * w (nb d n)| ≤ ∑ d, |κ d n r * w (nb d n)| :=
          Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ _d : D, |Kb| * R := by
          refine Finset.sum_le_sum fun d _ => ?_
          rw [abs_mul]
          exact mul_le_mul ((hκb d n r).trans (le_abs_self _)) (hw _) (abs_nonneg _) (abs_nonneg _)
      _ = Q * R := by
          rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, hQ]
          ring
  have hF_le : ∀ (v : X) n, ∀ r, 0 ≤ r → |F v n r| ≤ Q * ‖v‖ + β := by
    intro v n r hr
    have h1 := hsum_le (fun m => ext v m r) ‖v‖ n r (fun m => hext_le v m r)
    have h2 : |Real.exp (-μ * r) * b n r| ≤ β := by
      rw [abs_mul, abs_of_pos (Real.exp_pos _)]
      have : Real.exp (-μ * r) ≤ 1 := Real.exp_le_one_iff.mpr (by nlinarith)
      calc Real.exp (-μ * r) * |b n r| ≤ 1 * β :=
            mul_le_mul this (hbb n r) (abs_nonneg _) zero_le_one
        _ = β := one_mul β
    exact (abs_add_le _ _).trans (add_le_add h1 h2)
  let Γ₀ : X → S → ℝ → ℝ := fun v n s => Real.exp (-(V n + μ) * s) *
    (u₀ n + ∫ r in (0 : ℝ)..s, Real.exp ((V n + μ) * r) * F v n r)
  have hVμ : ∀ n, 0 < V n + μ := fun n => by linarith [hV n]
  have hInt : ∀ v n, Continuous fun s : ℝ =>
      ∫ r in (0 : ℝ)..s, Real.exp ((V n + μ) * r) * F v n r := fun v n =>
    intervalIntegral.continuous_primitive (fun a b =>
      (((by fun_prop : Continuous fun r : ℝ => Real.exp ((V n + μ) * r)).mul
        (hF v n)).intervalIntegrable a b)) 0
  have hΓcont : ∀ v n, Continuous (Γ₀ v n) := fun v n =>
    (by fun_prop : Continuous fun s : ℝ => Real.exp (-(V n + μ) * s)).mul
      (continuous_const.add (hInt v n))
  have hΓbd : ∀ (v : X) n, ∀ s ∈ Icc (0 : ℝ) T, |Γ₀ v n s| ≤ A + (Q * ‖v‖ + β) := by
    intro v n s hs
    have h1 : |Real.exp (-(V n + μ) * s) * u₀ n| ≤ A := by
      rw [abs_mul, abs_of_pos (Real.exp_pos _)]
      have : Real.exp (-(V n + μ) * s) ≤ 1 :=
        Real.exp_le_one_iff.mpr (by nlinarith [hVμ n, hs.1])
      calc Real.exp (-(V n + μ) * s) * |u₀ n| ≤ 1 * A :=
            mul_le_mul this (hA n) (abs_nonneg _) zero_le_one
        _ = A := one_mul A
    have h2 := abs_damped_integral_le (hVμ n) hs.1 (hF v n) (C := Q * ‖v‖ + β)
      (fun r hr => hF_le v n r hr.1)
    have h3 : (Q * ‖v‖ + β) / (V n + μ) ≤ Q * ‖v‖ + β := by
      have : 0 ≤ Q * ‖v‖ + β := by positivity
      exact div_le_self this (by linarith [hV n])
    have h4 : Γ₀ v n s = Real.exp (-(V n + μ) * s) * u₀ n + Real.exp (-(V n + μ) * s) *
        ∫ r in (0 : ℝ)..s, Real.exp ((V n + μ) * r) * F v n r := mul_add _ _ _
    rw [h4]
    exact (abs_add_le _ _).trans (add_le_add h1 (h2.trans h3))
  let Γ : X → X := fun v => BoundedContinuousFunction.ofNormedAddCommGroup
    (fun q : S × Icc (0 : ℝ) T => Γ₀ v q.1 q.2)
    (continuous_prod_of_discrete_left.mpr fun n => (hΓcont v n).comp continuous_subtype_val)
    (A + (Q * ‖v‖ + β)) (fun q => hΓbd v q.1 q.2 q.2.2)
  have hΓapp : ∀ v n (s : Icc (0 : ℝ) T), Γ v (n, s) = Γ₀ v n s := fun v n s => rfl
  have hlip : ∀ v w : X, dist (Γ v) (Γ w) ≤ Q / μ * dist v w := by
    intro v w
    rw [BoundedContinuousFunction.dist_le (by positivity)]
    rintro ⟨n, s⟩
    rw [Real.dist_eq, hΓapp, hΓapp]
    have hdiff : Γ₀ v n s - Γ₀ w n s = Real.exp (-(V n + μ) * s) *
        ∫ r in (0 : ℝ)..s, Real.exp ((V n + μ) * r) * (F v n r - F w n r) := by
      have hi : ∀ z : X, IntervalIntegrable
          (fun r => Real.exp ((V n + μ) * r) * F z n r) volume 0 s := fun z =>
        ((by fun_prop : Continuous fun r : ℝ => Real.exp ((V n + μ) * r)).mul
          (hF z n)).intervalIntegrable _ _
      simp only [mul_sub]
      rw [integral_sub (hi v) (hi w)]
      change Real.exp (-(V n + μ) * s) * _ - Real.exp (-(V n + μ) * s) * _ = _
      ring
    rw [hdiff]
    have hFd : ∀ r, |F v n r - F w n r| ≤ Q * dist v w := by
      intro r
      have h1 : F v n r - F w n r = ∑ d, κ d n r * (ext v (nb d n) r - ext w (nb d n) r) := by
        simp only [F, mul_sub, Finset.sum_sub_distrib]
        ring
      rw [h1]
      refine hsum_le (fun m => ext v m r - ext w m r) _ n r fun m => ?_
      have := BoundedContinuousFunction.dist_coe_le_dist (f := v) (g := w) (m, projIcc 0 T hT r)
      rwa [Real.dist_eq] at this
    refine (abs_damped_integral_le (hVμ n) s.2.1 ((hF v n).sub (hF w n))
      (fun r _ => hFd r)).trans ?_
    have h0 : 0 ≤ Q * dist v w := by positivity
    calc Q * dist v w / (V n + μ) ≤ Q * dist v w / μ :=
          div_le_div_of_nonneg_left h0 hμ0 (by linarith [hV n])
      _ = Q / μ * dist v w := by ring
  have hratio : Q / μ < 1 := by
    rw [div_lt_one hμ0]
    linarith
  have hcontr : ContractingWith ⟨Q / μ, by positivity⟩ Γ :=
    ⟨by exact_mod_cast hratio, LipschitzWith.of_dist_le_mul hlip⟩
  set v : X := hcontr.fixedPoint Γ with hv
  have hfix : Γ v = v := hcontr.fixedPoint_isFixedPt
  let I : S → ℝ → ℝ := fun n s => ∫ r in (0 : ℝ)..s, Real.exp ((V n + μ) * r) * F v n r
  let u : S → ℝ → ℝ := fun n s => Real.exp (-V n * s) * (u₀ n + I n s)
  have hu_ext : ∀ n, ∀ s ∈ Icc (0 : ℝ) T, u n s = Real.exp (μ * s) * ext v n s := by
    intro n s hs
    have h1 : ext v n s = Γ₀ v n s := by
      change v (n, projIcc 0 T hT s) = Γ₀ v n s
      rw [projIcc_of_mem hT hs]
      exact (congrArg (fun z : X => z (n, ⟨s, hs⟩)) hfix).symm
    rw [h1]
    change Real.exp (-V n * s) * (u₀ n + I n s)
      = Real.exp (μ * s) * (Real.exp (-(V n + μ) * s) * (u₀ n + I n s))
    rw [← mul_assoc, ← Real.exp_add]
    congr 2
    ring
  refine ⟨u, fun n => ?_, ⟨Real.exp (μ * T) * ‖v‖, fun n s hs => ?_⟩, fun n s hs => ?_⟩
  · change Real.exp (-V n * 0) * (u₀ n + ∫ r in (0 : ℝ)..0, _) = u₀ n
    simp
  · rw [hu_ext n s hs, abs_mul, abs_of_pos (Real.exp_pos _)]
    refine mul_le_mul (Real.exp_le_exp.mpr ?_) (hext_le v n s) (abs_nonneg _) (Real.exp_pos _).le
    exact mul_le_mul_of_nonneg_left hs.2 hμ0.le
  · have hexpd : HasDerivAt (fun s : ℝ => Real.exp (-V n * s)) (-V n * Real.exp (-V n * s)) s := by
      have := ((hasDerivAt_id s).const_mul (-V n)).exp
      simpa [mul_comm] using this
    have hIcont : Continuous fun r : ℝ => Real.exp ((V n + μ) * r) * F v n r :=
      (by fun_prop : Continuous fun r : ℝ => Real.exp ((V n + μ) * r)).mul (hF v n)
    have hId : HasDerivAt (I n) (Real.exp ((V n + μ) * s) * F v n s) s :=
      intervalIntegral.integral_hasDerivAt_right (hIcont.intervalIntegrable _ _)
        (hIcont.stronglyMeasurableAtFilter _ _) hIcont.continuousAt
    have hd := hexpd.mul (hId.const_add (u₀ n))
    refine (hd.hasDerivWithinAt).congr_deriv ?_
    have hnb : ∀ d, u (nb d n) s = Real.exp (μ * s) * ext v (nb d n) s := fun d =>
      hu_ext (nb d n) s hs
    have hprod : Real.exp (-V n * s) * (Real.exp ((V n + μ) * s) * F v n s)
        = ∑ d, κ d n s * u (nb d n) s + b n s := by
      rw [← mul_assoc, ← Real.exp_add]
      have he : -V n * s + (V n + μ) * s = μ * s := by ring
      rw [he]
      simp only [F, hnb, mul_add, Finset.mul_sum]
      congr 1
      · refine Finset.sum_congr rfl fun d _ => ?_
        ring
      · rw [← mul_assoc, ← Real.exp_add]
        have : μ * s + -μ * s = 0 := by ring
        rw [this, Real.exp_zero, one_mul]
    rw [hprod]
    change -V n * Real.exp (-V n * s) * (u₀ n + I n s) + (∑ d, κ d n s * u (nb d n) s + b n s)
      = ∑ d, κ d n s * u (nb d n) s - V n * (Real.exp (-V n * s) * (u₀ n + I n s)) + b n s
    ring

/-- **Existence** for coefficients that are continuous on `[0, T]`. -/
theorem exists_linear_solution {S D : Type*} [Fintype D] (nb : D → S → S)
    {κ : D → S → ℝ → ℝ} {T : ℝ} (hT : 0 ≤ T) (hκ : ∀ d n, ContinuousOn (κ d n) (Icc 0 T))
    {Kb : ℝ} (hκb : ∀ d n, ∀ s ∈ Icc 0 T, |κ d n s| ≤ Kb) {V : S → ℝ} (hV : ∀ n, 0 ≤ V n)
    {b : S → ℝ → ℝ} (hb : ∀ n, ContinuousOn (b n) (Icc 0 T)) {β : ℝ}
    (hbb : ∀ n, ∀ s ∈ Icc 0 T, |b n s| ≤ β) {u₀ : S → ℝ} {A : ℝ} (hA : ∀ n, |u₀ n| ≤ A) :
    ∃ u : S → ℝ → ℝ, (∀ n, u n 0 = u₀ n) ∧ (∃ B, ∀ n, ∀ s ∈ Icc 0 T, |u n s| ≤ B) ∧
      ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n)
        (∑ d, κ d n s * u (nb d n) s - V n * u n s + b n s) (Icc 0 T) s := by
  have hproj : ∀ (f : ℝ → ℝ), ContinuousOn f (Icc 0 T) →
      Continuous fun s => f (projIcc 0 T hT s) := fun f hf =>
    (continuousOn_iff_continuous_domRestrict.mp hf).comp continuous_projIcc
  obtain ⟨u, h0, hB, hu⟩ := exists_linear_solution_of_continuous nb
    (κ := fun d n s => κ d n (projIcc 0 T hT s)) (fun d n => hproj _ (hκ d n))
    (fun d n s => hκb d n _ (projIcc 0 T hT s).2) hV
    (b := fun n s => b n (projIcc 0 T hT s)) (fun n => hproj _ (hb n))
    (fun n s => hbb n _ (projIcc 0 T hT s).2) hA hT
  refine ⟨u, h0, hB, fun n s hs => ?_⟩
  have h := hu n s hs
  simpa only [projIcc_of_mem hT hs] using h

end ExecutionProofs
