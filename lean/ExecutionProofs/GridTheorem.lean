import ExecutionProofs.Hierarchy

/-!
# The grid theorem (paper, Theorem 3 (i))

For the linear-quadratic model with `H = 1`, `C ⪰ 0`, `Φ ⪰ φ > 0`: the lattice problem has exactly
one bounded solution, and it has a complete expansion in even powers of the spacing with explicit
Gaussian-times-polynomial coefficients, uniformly on the lattice and in the lattice shift.
-/

open Set Matrix MvPolynomial

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The potential `λ / (4 ν) x · C x`. -/
noncomputable def pot (ν lam : ℝ) (C : Matrix ι ι ℝ) (x : ι → ℝ) : ℝ :=
  lam / (4 * ν) * (x ⬝ᵥ (C *ᵥ x))

/-- Bounded solutions of the semidiscrete equation on the lattice `X + Δ ℤ^ι`. -/
structure IsLatticeSolution (ν Δ T : ℝ) (X : ι → ℝ) (V : (ι → ℝ) → ℝ)
    (u : (ι → ℤ) → ℝ → ℝ) : Prop where
  deriv : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n)
    (ν * ∑ a, (u (n + Pi.single a 1) s + u (n - Pi.single a 1) s - 2 * u n s) / Δ ^ 2
      - V (pos X Δ n) * u n s) (Icc 0 T) s
  bounded : ∃ B, ∀ n, ∀ s ∈ Icc 0 T, |u n s| ≤ B

lemma laplacian_eq_gen (ν Δ : ℝ) (u : (ι → ℤ) → ℝ) (n : ι → ℤ) :
    ν * ∑ a, (u (n + Pi.single a 1) + u (n - Pi.single a 1) - 2 * u n) / Δ ^ 2
      = gen step (fun (_ : ι × Bool) _ => ν / Δ ^ 2) u n := by
  rw [gen_step, Finset.mul_sum]
  exact Finset.sum_congr rfl fun a _ => by ring

/-- **Existence** of the bounded lattice solution. -/
theorem exists_lattice_solution {ν Δ T : ℝ} (hν : 0 ≤ ν) (hT : 0 ≤ T) (X : ι → ℝ)
    {V : (ι → ℝ) → ℝ} (hV : ∀ x, 0 ≤ V x) {u₀ : (ι → ℤ) → ℝ} {A : ℝ} (hA : ∀ n, |u₀ n| ≤ A) :
    ∃ u, IsLatticeSolution ν Δ T X V u ∧ ∀ n, u n 0 = u₀ n := by
  have hk : 0 ≤ ν / Δ ^ 2 := by positivity
  obtain ⟨u, h0, hB, hu⟩ := exists_linear_solution (nb := step) hT
    (κ := fun (_ : ι × Bool) (_ : ι → ℤ) _ => ν / Δ ^ 2) (fun _ _ => continuousOn_const)
    (Kb := ν / Δ ^ 2) (fun _ _ _ _ => (abs_of_nonneg hk).le)
    (V := fun n => V (pos X Δ n) + 2 * (ν / Δ ^ 2) * Fintype.card ι)
    (fun n => add_nonneg (hV _) (by positivity)) (b := fun _ _ => 0)
    (fun _ => continuousOn_const) (β := 0) (fun _ _ _ => by simp) hA
  refine ⟨u, ⟨fun n s hs => ?_, hB⟩, h0⟩
  refine (hu n s hs).congr_deriv ?_
  rw [laplacian_eq_gen ν Δ (fun m => u m s) n]
  simp only [gen, mul_sub, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
    Fintype.card_prod, Fintype.card_bool, nsmul_eq_mul]
  push_cast
  ring

/-- **Uniqueness** of the bounded lattice solution. -/
theorem lattice_solution_unique {ν Δ T : ℝ} (hν : 0 ≤ ν) {X : ι → ℝ} {V : (ι → ℝ) → ℝ}
    (hV : ∀ x, 0 ≤ V x) {u v : (ι → ℤ) → ℝ → ℝ} (hu : IsLatticeSolution ν Δ T X V u)
    (hv : IsLatticeSolution ν Δ T X V v) (h0 : ∀ n, u n 0 = v n 0) :
    ∀ n, ∀ s ∈ Icc 0 T, u n s = v n s := by
  intro n s hs
  obtain ⟨B, hB⟩ := hu.bounded
  obtain ⟨B', hB'⟩ := hv.bounded
  have hk : 0 ≤ ν / Δ ^ 2 := by positivity
  have key := abs_le_of_residual (nb := step) (κ := fun (_ : ι × Bool) (_ : ι → ℤ) => ν / Δ ^ 2)
    (V := fun n => V (pos X Δ n)) (ξ := barrier) (c := 2 * (ν / Δ ^ 2) * Fintype.card ι) (T := T)
    (u := fun n s => u n s - v n s) (r := fun _ _ => 0) (fun _ _ => hk) (fun n => hV _)
    one_le_barrier barrier_finite (fun n => (gen_barrier _ n).le)
    (fun n s hs => by
      refine ((hu.deriv n s hs).sub (hv.deriv n s hs)).congr_deriv ?_
      rw [laplacian_eq_gen ν Δ (fun m => u m s) n, laplacian_eq_gen ν Δ (fun m => v m s) n,
        gen_sub (fun m => u m s) (fun m => v m s)]
      ring)
    (B := B + B') (fun n s hs => (abs_sub _ _).trans (add_le_add (hB n s hs) (hB' n s hs)))
    (G := fun _ => 0) (g := fun _ => 0) (fun _ _ => hasDerivWithinAt_const _ _ _)
    (fun _ _ => le_rfl) (fun n => sub_eq_zero.mpr (h0 n)) (fun _ _ _ => by simp) n s hs
  exact sub_eq_zero.mp (abs_nonpos_iff.mp key)

omit [DecidableEq ι] in
/-- Scalar normalisation `θ' = tr P / 4`, `θ(0) = 0` (that is, `θ = ½ log det Y`). -/
theorem exists_theta {T : ℝ} (hT : 0 ≤ T) {P : ℝ → Matrix ι ι ℝ}
    (hP : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T)) :
    ∃ θ : ℝ → ℝ, θ 0 = 0 ∧
      ∀ s ∈ Icc 0 T, HasDerivWithinAt θ (1 / 4 * (P s).trace) (Icc 0 T) s := by
  have hc : ContinuousOn (fun s => 1 / 4 * (P s).trace) (Icc 0 T) := by
    simp only [Matrix.trace, Matrix.diag_apply]
    exact continuousOn_const.mul (continuousOn_finsetSum _ fun a _ => hP a a)
  obtain ⟨β, hβ⟩ := isCompact_Icc.exists_bound_of_continuousOn hc
  obtain ⟨u, h0, -, hu⟩ := exists_linear_solution (S := Unit) (D := Empty)
    (nb := fun d _ => d.elim) hT (κ := fun d _ _ => d.elim) (fun d _ => d.elim) (Kb := 0)
    (fun d _ _ _ => d.elim) (V := fun _ => 0) (fun _ => le_rfl)
    (b := fun _ s => 1 / 4 * (P s).trace) (fun _ => hc) (fun _ => hβ) (u₀ := fun _ => 0) (A := 0)
    (fun _ => by simp)
  refine ⟨u (), h0 (), fun s hs => (hu () s hs).congr_deriv ?_⟩
  simp

/-- **Grid theorem, part (i).** Complete expansion of the lattice solution with a remainder constant
that does not depend on the spacing or on the lattice shift. -/
theorem grid_theorem {ν lam T φ : ℝ} (hν : 0 < ν) (hlam : 0 ≤ lam) (hT : 0 < T)
    {C Φ : Matrix ι ι ℝ} (hCs : C.IsSymm) (hΦs : Φ.IsSymm) (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v))
    (hφ : 0 < φ) (hΦ : ∀ v, φ * (v ⬝ᵥ v) ≤ v ⬝ᵥ (Φ *ᵥ v)) {Y Z : ℝ → Matrix ι ι ℝ}
    (hYZ : RiccatiFlow C Φ lam T Y Z) {θ : ℝ → ℝ}
    (hθ : ∀ s ∈ Icc 0 T,
      HasDerivWithinAt θ (1 / 4 * (RiccatiFlow.P Y Z s).trace) (Icc 0 T) s) (K : ℕ)
    {p : ℕ → ℝ → MvPolynomial ι ℝ} (hp : IsHierarchy T ν (RiccatiFlow.P Y Z) K p) :
    ∃ B, 0 ≤ B ∧ ∀ Δ, 0 < Δ → ∀ (X : ι → ℝ) (u : (ι → ℤ) → ℝ → ℝ),
      IsLatticeSolution ν Δ T X (pot ν lam C) u →
      (∀ n, u n 0 = gauss ν (RiccatiFlow.P Y Z 0) (θ 0) (pos X Δ n)) →
      ∀ n, ∀ s ∈ Icc 0 T,
        |u n s - ∑ k ∈ Finset.range (K + 1), Δ ^ (2 * k) * (eval (pos X Δ n) (p k s)
          * gauss ν (RiccatiFlow.P Y Z s) (θ s) (pos X Δ n))| ≤ Δ ^ (2 * K + 2) * B := by
  classical
  set P := RiccatiFlow.P Y Z with hPdef
  have hPd : ∀ s ∈ Icc 0 T, ∀ i j, HasDerivWithinAt (fun s => P s i j)
      (2 * lam * C i j - 1 / 2 * (P s * P s) i j) (Icc 0 T) s := fun s hs i j =>
    hYZ.P_deriv hT hlam hC hφ hΦ hs i j
  have hPc : ∀ a b, ContinuousOn (fun s => P s a b) (Icc 0 T) := fun a b s hs =>
    (hPd s hs a b).continuousWithinAt
  have hPs : ∀ s ∈ Icc 0 T, (P s).IsSymm := fun s hs => hYZ.P_symm hlam hCs hΦs hC hφ hΦ hs
  have hθc : ContinuousOn θ (Icc 0 T) := fun s hs => (hθ s hs).continuousWithinAt
  obtain ⟨m, hm, hcoer⟩ := hYZ.exists_coercive hlam hC hφ hΦ
  have hV : ∀ x, 0 ≤ pot ν lam C x := fun x => mul_nonneg (by positivity) (hC x)
  -- derivative bounds
  have hN : ∀ k, ∃ N, 0 ≤ N ∧ (k ≤ K → ∀ a, ∀ s ∈ Icc 0 T, ∀ x,
      |pd a (2 * (K - k) + 4) (fun x => eval x (p k s) * gauss ν (P s) (θ s) x) x| ≤ N) := by
    intro k
    by_cases hk : k ≤ K
    · have h : ∀ a, ∃ N, 0 ≤ N ∧ ∀ s ∈ Icc 0 T, ∀ x,
          |eval x ((Dop ν (P s) a)^[2 * (K - k) + 4] (p k s)) * gauss ν (P s) (θ s) x| ≤ N :=
        fun a => ((hp.fam k hk).dop_iterate hPc a _).bounded hν hθc hm hcoer
      choose N hN0 hNb using h
      refine ⟨∑ a, N a, Finset.sum_nonneg fun a _ => hN0 a, fun _ a s hs x => ?_⟩
      rw [pd_polyGauss (hPs s hs)]
      exact (hNb a s hs x).trans
        (Finset.single_le_sum (f := N) (fun a _ => hN0 a) (Finset.mem_univ a))
    · exact ⟨0, le_rfl, fun h => absurd h hk⟩
  choose N hN0 hNb using hN
  set g₀ : ℝ := ν * ∑ k ∈ Finset.range (K + 1), c (K - k + 1) * ∑ _a : ι, N k with hg₀
  have hg₀nn : 0 ≤ g₀ := by
    refine mul_nonneg hν.le (Finset.sum_nonneg fun k _ => mul_nonneg (c_pos _).le ?_)
    exact Finset.sum_nonneg fun a _ => hN0 k
  refine ⟨g₀ * T, by positivity, fun Δ hΔ X u hu hu0 n s hs => ?_⟩
  have hpder : ∀ k ≤ K, ∀ x, ∀ s ∈ Icc 0 T, HasDerivWithinAt (fun s => eval x (p k s))
      (eval x (polyGen ν (P s) (p k s)) + eval x (source ν (P s) (fun i => p i s) k))
      (Icc 0 T) s := by
    intro k hk x s hs
    rcases Nat.eq_zero_or_pos k with rfl | hpos
    · have h1 : (fun s => eval x (p 0 s)) = fun _ => 1 := by
        funext s
        rw [hp.zero, map_one]
      have h2 : eval x (polyGen ν (P s) (p 0 s))
          + eval x (source ν (P s) (fun i => p i s) 0) = 0 := by
        rw [hp.zero]
        simp [polyGen, source]
      rw [h1, h2]
      exact hasDerivWithinAt_const _ _ _
    · exact hp.deriv k hpos hk x s hs
  have key := grid_expansion (ν := ν) (Δ := Δ) (T := T) hν.le hΔ X hV K
    (fun k s x => eval x (p k s) * gauss ν (P s) (θ s) x)
    (fun k _ s hs a x => contDiff_polyGauss (hPs s hs) (θ s) (p k s) x a _)
    (fun k hk => by
      obtain ⟨N, -, hN⟩ := (hp.fam k hk).bounded hν hθc hm hcoer
      exact ⟨N, hN⟩)
    (fun k hk x s hs => by
      refine ((hpder k hk x s hs).fun_mul (hasDeriv_gauss hν.ne' (hPs s hs) (hPd s hs)
        (hθ s hs) x)).congr_deriv ?_
      have hpd : ∀ i a, pd a (2 * (k - i) + 2) (fun x => eval x (p i s) * gauss ν (P s) (θ s) x) x
          = eval x ((Dop ν (P s) a)^[2 * (k - i) + 2] (p i s)) * gauss ν (P s) (θ s) x :=
        fun i a => pd_polyGauss (hPs s hs) _ _ _ _ _
      simp only [hpd, ← Finset.sum_mul, ← mul_assoc]
      rw [Finset.sum_range_succ]
      have hlast : c (k - k) * ∑ a, eval x ((Dop ν (P s) a)^[2 * (k - k) + 2] (p k s))
          = ∑ a, eval x (Dop ν (P s) a (Dop ν (P s) a (p k s))) := by
        simp only [Nat.sub_self, mul_zero, zero_add, c_zero, one_mul]
        rfl
      have hsrc : eval x (source ν (P s) (fun i => p i s) k) = ν * ∑ i ∈ Finset.range k,
          c (k - i) * ∑ a, eval x ((Dop ν (P s) a)^[2 * (k - i) + 2] (p i s)) := by
        simp only [source, smul_eval, map_sum]
      have hdd := sum_dop_dop hν.ne' (P s) x (p k s)
      rw [hlast, hsrc, pot]
      linear_combination (-gauss ν (P s) (θ s) x) * hdd)
    (fun k h1 hk x => by
      simp only [hp.init k h1 hk, map_zero, zero_mul])
    (N := fun k _ _ => N k) (fun k hk a s hs x => hNb k hk a s hs x)
    (G := fun s => g₀ * s)
    (fun s _ => (((hasDerivAt_id s).const_mul g₀).hasDerivWithinAt).congr_deriv (mul_one g₀))
    (fun s hs => mul_nonneg hg₀nn hs.1) (u := u) hu.deriv hu.bounded
    (fun n => by
      rw [hu0 n, hp.zero, map_one, one_mul]) n s hs
  refine key.trans (mul_le_mul_of_nonneg_left ?_ (by positivity))
  exact mul_le_mul_of_nonneg_left hs.2 hg₀nn

/-- At `s = 0` the Gaussian is the prescribed initial datum. -/
lemma gauss_zero {C Φ : Matrix ι ι ℝ} {lam T ν : ℝ} {Y Z : ℝ → Matrix ι ι ℝ}
    (hYZ : RiccatiFlow C Φ lam T Y Z) {θ : ℝ → ℝ} (hθ : θ 0 = 0) (x : ι → ℝ) :
    gauss ν (RiccatiFlow.P Y Z 0) (θ 0) x = Real.exp (-(x ⬝ᵥ (Φ *ᵥ x)) / (8 * ν)) := by
  rw [gauss, RiccatiFlow.P, hYZ.Y_zero, hYZ.Z_zero, hθ, inv_one, Matrix.mul_one, sub_zero]

/-- **Grid theorem, unconditional form.** Riccati flow, normalisation and polynomials exist; for every
spacing and shift the lattice problem has a unique bounded solution, and it is within `B Δ^(2K+2)` of
the partial sum. -/
theorem grid_theorem_main {ν lam T φ : ℝ} (hν : 0 < ν) (hlam : 0 ≤ lam) (hT : 0 < T)
    {C Φ : Matrix ι ι ℝ} (hCs : C.IsSymm) (hΦs : Φ.IsSymm) (hC : ∀ v, 0 ≤ v ⬝ᵥ (C *ᵥ v))
    (hφ : 0 < φ) (hΦ : ∀ v, φ * (v ⬝ᵥ v) ≤ v ⬝ᵥ (Φ *ᵥ v)) (K : ℕ) :
    ∃ (Y Z : ℝ → Matrix ι ι ℝ) (θ : ℝ → ℝ) (p : ℕ → ℝ → MvPolynomial ι ℝ) (B : ℝ),
      RiccatiFlow C Φ lam T Y Z ∧ θ 0 = 0 ∧
      (∀ s ∈ Icc 0 T, HasDerivWithinAt θ (1 / 4 * (RiccatiFlow.P Y Z s).trace) (Icc 0 T) s) ∧
      IsHierarchy T ν (RiccatiFlow.P Y Z) K p ∧ 0 ≤ B ∧
      ∀ Δ, 0 < Δ → ∀ X : ι → ℝ, ∃ u, IsLatticeSolution ν Δ T X (pot ν lam C) u ∧
        (∀ n, u n 0 = Real.exp (-(pos X Δ n ⬝ᵥ (Φ *ᵥ pos X Δ n)) / (8 * ν))) ∧
        (∀ v, IsLatticeSolution ν Δ T X (pot ν lam C) v → (∀ n, v n 0 = u n 0) →
          ∀ n, ∀ s ∈ Icc 0 T, v n s = u n s) ∧
        ∀ n, ∀ s ∈ Icc 0 T,
          |u n s - ∑ k ∈ Finset.range (K + 1), Δ ^ (2 * k) * (eval (pos X Δ n) (p k s)
            * gauss ν (RiccatiFlow.P Y Z s) (θ s) (pos X Δ n))| ≤ Δ ^ (2 * K + 2) * B := by
  obtain ⟨Y, Z, hYZ⟩ := exists_riccatiFlow C Φ lam hT.le
  have hPc : ∀ a b, ContinuousOn (fun s => RiccatiFlow.P Y Z s a b) (Icc 0 T) := fun a b s hs =>
    (hYZ.P_deriv hT hlam hC hφ hΦ hs a b).continuousWithinAt
  obtain ⟨θ, hθ0, hθ⟩ := exists_theta hT.le hPc
  obtain ⟨p, hp⟩ := exists_hierarchy (ν := ν) hT.le hPc K
  obtain ⟨B, hB, hmain⟩ := grid_theorem hν hlam hT hCs hΦs hC hφ hΦ hYZ hθ K hp
  refine ⟨Y, Z, θ, p, B, hYZ, hθ0, hθ, hp, hB, fun Δ hΔ X => ?_⟩
  have hV : ∀ x, 0 ≤ pot ν lam C x := fun x => mul_nonneg (by positivity) (hC x)
  obtain ⟨u, hu, hu0⟩ := exists_lattice_solution (Δ := Δ) hν.le hT.le X hV
    (u₀ := fun n => Real.exp (-(pos X Δ n ⬝ᵥ (Φ *ᵥ pos X Δ n)) / (8 * ν))) (A := 1) (fun n => by
      have hq : 0 ≤ pos X Δ n ⬝ᵥ (Φ *ᵥ pos X Δ n) :=
        le_trans (mul_nonneg hφ.le (Finset.sum_nonneg fun i _ => mul_self_nonneg _)) (hΦ _)
      rw [abs_of_pos (Real.exp_pos _), Real.exp_le_one_iff]
      exact div_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr hq) (by positivity))
  refine ⟨u, hu, hu0, fun v hv hv0 => lattice_solution_unique hν.le hV hv hu hv0, ?_⟩
  exact hmain Δ hΔ X u hu fun n => by rw [hu0 n, gauss_zero hYZ hθ0]

end ExecutionProofs
