import ExecutionProofs.MaxPrinciple
import ExecutionProofs.SecondDifference

/-!
# Expansion of the lattice solution in powers of the grid spacing

The lattice `X + Δ ℤ^ι` carries the semidiscrete heat equation with a nonnegative potential.
`grid_expansion` is part (i) of the grid theorem: stability (maximum principle) plus consistency
(Taylor) give an explicit remainder of order `Δ ^ (2K + 2)`, uniformly on the lattice.
-/

open Set Finset

namespace ExecutionProofs

variable {ι : Type*} [DecidableEq ι]

/-- Position of lattice node `n`. -/
def pos (X : ι → ℝ) (Δ : ℝ) (n : ι → ℤ) : ι → ℝ := fun a => X a + Δ * n a

/-- `m`-th derivative of `f` in coordinate direction `a`. -/
noncomputable def pd (a : ι) (m : ℕ) (f : (ι → ℝ) → ℝ) (x : ι → ℝ) : ℝ :=
  iteratedDeriv m (fun t : ℝ => f (x + Pi.single a t)) 0

/-- Lattice neighbours: one step up or down along each axis. -/
def step (d : ι × Bool) (n : ι → ℤ) : ι → ℤ :=
  if d.2 then n + Pi.single d.1 1 else n - Pi.single d.1 1

lemma pos_add_single (X : ι → ℝ) (Δ : ℝ) (n : ι → ℤ) (a : ι) :
    pos X Δ (n + Pi.single a 1) = pos X Δ n + Pi.single a Δ := by
  funext b
  by_cases h : b = a
  · subst h
    simp only [pos, Pi.add_apply, Pi.single_eq_same]
    push_cast
    ring
  · simp [pos, h]

lemma pos_sub_single (X : ι → ℝ) (Δ : ℝ) (n : ι → ℤ) (a : ι) :
    pos X Δ (n - Pi.single a 1) = pos X Δ n + Pi.single a (-Δ) := by
  funext b
  by_cases h : b = a
  · subst h
    simp only [pos, Pi.add_apply, Pi.sub_apply, Pi.single_eq_same]
    push_cast
    ring
  · simp [pos, h]

variable [Fintype ι]

/-- The lattice Laplacian written as a jump generator. -/
lemma gen_step (k : ℝ) (u : (ι → ℤ) → ℝ) (n : ι → ℤ) :
    gen step (fun (_ : ι × Bool) _ => k) u n
      = ∑ a, k * (u (n + Pi.single a 1) + u (n - Pi.single a 1) - 2 * u n) := by
  simp only [gen, Fintype.sum_prod_type, Fintype.sum_bool, step]
  refine sum_congr rfl fun a _ => ?_
  simp only [ite_true, Bool.false_eq_true, ite_false]
  ring

/-- Quadratic barrier `1 + ∑ n_a²`. -/
def barrier (n : ι → ℤ) : ℝ := 1 + ∑ a, (n a : ℝ) ^ 2

omit [DecidableEq ι] in
lemma one_le_barrier (n : ι → ℤ) : 1 ≤ barrier n := by
  unfold barrier
  have : 0 ≤ ∑ a, (n a : ℝ) ^ 2 := sum_nonneg fun a _ => sq_nonneg _
  linarith

lemma barrier_finite (R : ℝ) : {n : ι → ℤ | barrier n ≤ R}.Finite := by
  obtain ⟨m, hm⟩ := exists_nat_ge R
  refine (Set.finite_Icc (fun _ : ι => -(m : ℤ)) (fun _ => (m : ℤ))).subset ?_
  intro n hn
  have hn' : barrier n ≤ R := hn
  have hbound : ∀ a, (n a : ℝ) ^ 2 ≤ m := by
    intro a
    have h1 : (n a : ℝ) ^ 2 ≤ ∑ b, (n b : ℝ) ^ 2 :=
      single_le_sum (f := fun b => (n b : ℝ) ^ 2) (fun b _ => sq_nonneg _) (Finset.mem_univ a)
    unfold barrier at hn'
    linarith
  have habs : ∀ a, |n a| ≤ (m : ℤ) := by
    intro a
    have h1 : (n a) ^ 2 ≤ (m : ℤ) := by exact_mod_cast hbound a
    have h2 : |n a| ≤ (n a) ^ 2 := by
      rcases eq_or_ne (n a) 0 with h | h
      · simp [h]
      · have : 1 ≤ |n a| := Int.one_le_abs h
        rw [← sq_abs]
        nlinarith
    exact h2.trans h1
  exact ⟨fun a => (abs_le.mp (habs a)).1, fun a => (abs_le.mp (habs a)).2⟩

lemma gen_barrier (k : ℝ) (n : ι → ℤ) :
    gen step (fun (_ : ι × Bool) _ => k) barrier n = 2 * k * Fintype.card ι := by
  rw [gen_step]
  have h : ∀ a : ι, barrier (n + Pi.single a 1) + barrier (n - Pi.single a 1) - 2 * barrier n
      = 2 := by
    intro a
    have hplus : ∀ e : ℤ, barrier (n + Pi.single a e)
        = barrier n + (2 * (n a : ℝ) * e + (e : ℝ) ^ 2) := by
      intro e
      have hsum : ∑ b, (((n + Pi.single a e : ι → ℤ) b : ℤ) : ℝ) ^ 2
          = ∑ b, ((n b : ℝ) ^ 2 + if b = a then 2 * (n a : ℝ) * e + (e : ℝ) ^ 2 else 0) := by
        refine sum_congr rfl fun b _ => ?_
        by_cases hb : b = a
        · subst hb
          simp only [Pi.add_apply, Pi.single_eq_same, ite_true]
          push_cast
          ring
        · simp [hb]
      unfold barrier
      rw [hsum, sum_add_distrib, sum_ite_eq' univ a, ite_eq_left (Finset.mem_univ a)]
      ring
    have hminus : n - Pi.single a 1 = n + Pi.single a (-1) := by
      rw [sub_eq_add_neg, Pi.single_neg]
    rw [hminus, hplus, hplus]
    push_cast
    ring
  simp only [h, sum_const, card_univ, nsmul_eq_mul]
  ring

omit [Fintype ι] in
/-- Derivatives along a coordinate line are translates of the partial derivatives. -/
lemma iteratedDeriv_line (f : (ι → ℝ) → ℝ) (x : ι → ℝ) (a : ι) (m : ℕ) (t : ℝ) :
    iteratedDeriv m (fun t : ℝ => f (x + Pi.single a t)) t = pd a m f (x + Pi.single a t) := by
  unfold pd
  have h : (fun r : ℝ => f (x + Pi.single a t + Pi.single a r))
      = fun r => (fun t : ℝ => f (x + Pi.single a t)) (t + r) := by
    funext r
    simp only [Pi.single_add, add_assoc]
  have h' := congrFun (iteratedDeriv_comp_const_add m (fun t : ℝ => f (x + Pi.single a t)) t) 0
  rw [h, h', add_zero]

/-- **Grid expansion (stability + consistency).** `ψ 0, …, ψ K` solve the hierarchy
`∂ₛ ψ_k = -V ψ_k + ν ∑_{i ≤ k} c (k-i) ∑ₐ ∂ₐ^{2(k-i)+2} ψ_i` with `ψ_k(0) = 0` for `k ≥ 1`, and
`u` is a bounded solution of the lattice equation started from `ψ 0`. Then
`|u - ∑ Δ^{2k} ψ_k| ≤ Δ^{2K+2} G` on the whole lattice, for any primitive `G` of the
derivative bounds. -/
theorem grid_expansion {ν Δ T : ℝ} (hν : 0 ≤ ν) (hΔ : 0 < Δ) (X : ι → ℝ) {V : (ι → ℝ) → ℝ}
    (hV : ∀ x, 0 ≤ V x) (K : ℕ) (ψ : ℕ → ℝ → (ι → ℝ) → ℝ)
    (hsmooth : ∀ k ≤ K, ∀ s ∈ Icc 0 T, ∀ (a : ι) (x : ι → ℝ),
      ContDiff ℝ (2 * (K - k) + 4 : ℕ) (fun t : ℝ => ψ k s (x + Pi.single a t)))
    (hbound : ∀ k ≤ K, ∃ C, ∀ s ∈ Icc 0 T, ∀ x, |ψ k s x| ≤ C)
    (hode : ∀ k ≤ K, ∀ x, ∀ s ∈ Icc 0 T, HasDerivWithinAt (fun s => ψ k s x)
      (-V x * ψ k s x + ν * ∑ i ∈ range (k + 1), c (k - i) * ∑ a, pd a (2 * (k - i) + 2) (ψ i s) x)
      (Icc 0 T) s)
    (hinit : ∀ k, 1 ≤ k → k ≤ K → ∀ x, ψ k 0 x = 0)
    {N : ℕ → ι → ℝ → ℝ}
    (hN : ∀ k ≤ K, ∀ (a : ι), ∀ s ∈ Icc 0 T, ∀ x, |pd a (2 * (K - k) + 4) (ψ k s) x| ≤ N k a s)
    {G : ℝ → ℝ} (hG : ∀ s ∈ Icc 0 T, HasDerivWithinAt G
      (ν * ∑ k ∈ range (K + 1), c (K - k + 1) * ∑ a, N k a s) (Icc 0 T) s)
    (hGnn : ∀ s ∈ Icc 0 T, 0 ≤ G s)
    {u : (ι → ℤ) → ℝ → ℝ}
    (hu : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n)
      (ν * ∑ a, (u (n + Pi.single a 1) s + u (n - Pi.single a 1) s - 2 * u n s) / Δ ^ 2
        - V (pos X Δ n) * u n s) (Icc 0 T) s)
    (hub : ∃ B, ∀ n, ∀ s ∈ Icc 0 T, |u n s| ≤ B) (hu0 : ∀ n, u n 0 = ψ 0 0 (pos X Δ n)) :
    ∀ n, ∀ s ∈ Icc 0 T,
      |u n s - ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ψ k s (pos X Δ n)| ≤ Δ ^ (2 * K + 2) * G s := by
  set κ : ℝ := ν / Δ ^ 2 with hκ
  have hκnn : 0 ≤ κ := by positivity
  set sm : (ι → ℤ) → ℝ → ℝ := fun n s => ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ψ k s (pos X Δ n)
    with hsm
  -- second difference of `ψ k` along axis `a`
  set sd : ℕ → ι → ℝ → (ι → ℤ) → ℝ := fun k a s n =>
    (ψ k s (pos X Δ n + Pi.single a Δ) + ψ k s (pos X Δ n + Pi.single a (-Δ))
      - 2 * ψ k s (pos X Δ n)) / Δ ^ 2 with hsd
  set r : (ι → ℤ) → ℝ → ℝ := fun n s => ν * ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ a, (sd k a s n
    - ∑ j ∈ range (K + 1 - k), c j * Δ ^ (2 * j) * pd a (2 * j + 2) (ψ k s) (pos X Δ n)) with hr
  -- derivative of the partial sum
  have hsmderiv : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (sm n)
      (∑ k ∈ range (K + 1), Δ ^ (2 * k) * (-V (pos X Δ n) * ψ k s (pos X Δ n)
        + ν * ∑ i ∈ range (k + 1), c (k - i) * ∑ a, pd a (2 * (k - i) + 2) (ψ i s) (pos X Δ n)))
      (Icc 0 T) s := by
    intro n s hs
    refine HasDerivWithinAt.fun_sum fun k hk => ?_
    exact (hode k (Nat.lt_succ_iff.mp (mem_range.mp hk)) _ s hs).const_mul _
  -- generator of the partial sum
  have hgensm : ∀ n s, gen step (fun (_ : ι × Bool) _ => κ) (fun m => sm m s) n
      = ν * ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ a, sd k a s n := by
    intro n s
    rw [gen_step]
    simp only [hsm, hsd, pos_add_single, pos_sub_single, hκ]
    rw [mul_sum]
    simp only [mul_sum, ← sum_add_distrib, ← sum_sub_distrib]
    rw [sum_comm]
    refine sum_congr rfl fun k _ => sum_congr rfl fun a _ => ?_
    ring
  -- flip the double sum of the hierarchy
  have hflip : ∀ n s, ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ i ∈ range (k + 1),
        c (k - i) * ∑ a, pd a (2 * (k - i) + 2) (ψ i s) (pos X Δ n)
      = ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ a, ∑ j ∈ range (K + 1 - k),
          c j * Δ ^ (2 * j) * pd a (2 * j + 2) (ψ k s) (pos X Δ n) := by
    intro n s
    have h := sum_range_diag_flip (K + 1) (fun i j =>
      Δ ^ (2 * (i + j)) * (c j * ∑ a, pd a (2 * j + 2) (ψ i s) (pos X Δ n)))
    calc _ = ∑ k ∈ range (K + 1), ∑ i ∈ range (k + 1), Δ ^ (2 * (i + (k - i)))
          * (c (k - i) * ∑ a, pd a (2 * (k - i) + 2) (ψ i s) (pos X Δ n)) := by
          refine sum_congr rfl fun k _ => ?_
          rw [mul_sum]
          refine sum_congr rfl fun i hi => ?_
          rw [Nat.add_sub_cancel' (Nat.lt_succ_iff.mp (mem_range.mp hi))]
      _ = ∑ k ∈ range (K + 1), ∑ j ∈ range (K + 1 - k), Δ ^ (2 * (k + j))
          * (c j * ∑ a, pd a (2 * j + 2) (ψ k s) (pos X Δ n)) := h
      _ = ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ a, ∑ j ∈ range (K + 1 - k),
          c j * Δ ^ (2 * j) * pd a (2 * j + 2) (ψ k s) (pos X Δ n) := by
          refine sum_congr rfl fun k _ => ?_
          simp only [mul_sum]
          rw [sum_comm]
          refine sum_congr rfl fun a _ => sum_congr rfl fun j _ => ?_
          rw [mul_add, pow_add]
          ring
  have key := abs_le_of_residual (nb := step) (κ := fun (_ : ι × Bool) (_ : ι → ℤ) => κ)
    (V := fun n => V (pos X Δ n)) (ξ := barrier) (c := 2 * κ * Fintype.card ι) (T := T)
    (u := fun n s => u n s - sm n s) (r := r) (fun _ _ => hκnn) (fun n => hV _)
    one_le_barrier barrier_finite (fun n => (gen_barrier κ n).le)
    (fun n s hs => by
      refine ((hu n s hs).sub (hsmderiv n s hs)).congr_deriv ?_
      have h1 : gen step (fun (_ : ι × Bool) _ => κ) (fun m => u m s - sm m s) n
          = ν * ∑ a, (u (n + Pi.single a 1) s + u (n - Pi.single a 1) s - 2 * u n s) / Δ ^ 2
            - ν * ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ a, sd k a s n := by
        rw [gen_sub (fun m => u m s) (fun m => sm m s), hgensm, gen_step]
        have hfac : ∑ a, κ * (u (n + Pi.single a 1) s + u (n - Pi.single a 1) s - 2 * u n s)
            = ν * ∑ a, (u (n + Pi.single a 1) s + u (n - Pi.single a 1) s - 2 * u n s) / Δ ^ 2 := by
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun a _ => ?_
          rw [hκ]
          ring
        rw [hfac]
      have h2 : r n s = ν * ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ a, sd k a s n
          - ν * ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ i ∈ range (k + 1),
            c (k - i) * ∑ a, pd a (2 * (k - i) + 2) (ψ i s) (pos X Δ n) := by
        rw [hflip]
        simp only [hr, sum_sub_distrib, mul_sub]
      have h3 : ∑ k ∈ range (K + 1), Δ ^ (2 * k) * (-V (pos X Δ n) * ψ k s (pos X Δ n)
            + ν * ∑ i ∈ range (k + 1), c (k - i) * ∑ a, pd a (2 * (k - i) + 2) (ψ i s) (pos X Δ n))
          = -V (pos X Δ n) * sm n s + ν * ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ i ∈ range (k + 1),
            c (k - i) * ∑ a, pd a (2 * (k - i) + 2) (ψ i s) (pos X Δ n) := by
        simp only [hsm, mul_sum, ← sum_add_distrib]
        refine sum_congr rfl fun k _ => ?_
        rw [mul_add, mul_sum]
        congr 1
        · ring
        · refine sum_congr rfl fun i _ => ?_
          rw [mul_sum]
          refine sum_congr rfl fun a _ => ?_
          ring
      rw [h1, h2, h3]
      ring)
    (B := (hub.choose + ∑ k ∈ range (K + 1), Δ ^ (2 * k) *
      (if h : k ≤ K then (hbound k h).choose else 0)))
    (fun n s hs => by
      have h1 := hub.choose_spec n s hs
      have h2 : |sm n s| ≤ ∑ k ∈ range (K + 1), Δ ^ (2 * k) *
          (if h : k ≤ K then (hbound k h).choose else 0) := by
        refine (abs_sum_le_sum_abs _ _).trans (sum_le_sum fun k hk => ?_)
        have hk' : k ≤ K := Nat.lt_succ_iff.mp (mem_range.mp hk)
        have hp : 0 ≤ Δ ^ (2 * k) := by positivity
        rw [abs_mul, abs_of_nonneg hp, dite_eq_left hk']
        exact mul_le_mul_of_nonneg_left ((hbound k hk').choose_spec s hs _) hp
      exact (abs_sub _ _).trans (add_le_add h1 h2))
    (G := fun s => Δ ^ (2 * K + 2) * G s)
    (g := fun s => Δ ^ (2 * K + 2) * (ν * ∑ k ∈ range (K + 1), c (K - k + 1) * ∑ a, N k a s))
    (fun s hs => (hG s hs).const_mul _)
    (fun s hs => mul_nonneg (by positivity) (hGnn s hs))
    (fun n => by
      have h : sm n 0 = ψ 0 0 (pos X Δ n) := by
        simp only [hsm]
        rw [sum_range_succ']
        simp only [mul_zero, pow_zero, one_mul]
        rw [sum_eq_zero, zero_add]
        intro k hk
        rw [hinit (k + 1) (Nat.succ_pos k) (Nat.succ_le_of_lt (mem_range.mp hk)), mul_zero]
      change u n 0 - sm n 0 = 0
      rw [hu0, h, sub_self])
    (fun n s hs => by
      change |ν * ∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ a, (sd k a s n - ∑ j ∈ range (K + 1 - k),
        c j * Δ ^ (2 * j) * pd a (2 * j + 2) (ψ k s) (pos X Δ n))| ≤ _
      have hterm : ∀ k ∈ range (K + 1), ∀ a, |sd k a s n - ∑ j ∈ range (K + 1 - k),
          c j * Δ ^ (2 * j) * pd a (2 * j + 2) (ψ k s) (pos X Δ n)|
          ≤ c (K - k + 1) * Δ ^ (2 * (K - k + 1)) * N k a s := by
        intro k hk a
        have hk' : k ≤ K := Nat.lt_succ_iff.mp (mem_range.mp hk)
        have hb := second_difference_bound (L := fun t : ℝ => ψ k s (pos X Δ n + Pi.single a t))
          (J := K - k) (hsmooth k hk' s hs a _) (M := N k a s)
          (fun t => by
            rw [iteratedDeriv_line]
            exact hN k hk' a s hs _) hΔ.ne'
        have hKk : K - k + 1 = K + 1 - k := by omega
        simpa [hsd, pd, hKk] using hb
      rw [abs_mul, abs_of_nonneg hν]
      have hsum : |∑ k ∈ range (K + 1), Δ ^ (2 * k) * ∑ a, (sd k a s n - ∑ j ∈ range (K + 1 - k),
          c j * Δ ^ (2 * j) * pd a (2 * j + 2) (ψ k s) (pos X Δ n))|
          ≤ ∑ k ∈ range (K + 1), Δ ^ (2 * K + 2) * (c (K - k + 1) * ∑ a, N k a s) := by
        refine (abs_sum_le_sum_abs _ _).trans (sum_le_sum fun k hk => ?_)
        have hk' : k ≤ K := Nat.lt_succ_iff.mp (mem_range.mp hk)
        have hp : 0 ≤ Δ ^ (2 * k) := by positivity
        rw [abs_mul, abs_of_nonneg hp]
        calc Δ ^ (2 * k) * |∑ a, (sd k a s n - ∑ j ∈ range (K + 1 - k),
              c j * Δ ^ (2 * j) * pd a (2 * j + 2) (ψ k s) (pos X Δ n))|
            ≤ Δ ^ (2 * k) * ∑ a, c (K - k + 1) * Δ ^ (2 * (K - k + 1)) * N k a s :=
              mul_le_mul_of_nonneg_left
                ((abs_sum_le_sum_abs _ _).trans (sum_le_sum fun a _ => hterm k hk a)) hp
          _ = Δ ^ (2 * K + 2) * (c (K - k + 1) * ∑ a, N k a s) := by
              have hexp : Δ ^ (2 * K + 2) = Δ ^ (2 * k) * Δ ^ (2 * (K - k + 1)) := by
                rw [← pow_add]
                congr 1
                omega
              rw [hexp, mul_sum, mul_sum, mul_sum]
              refine sum_congr rfl fun a _ => ?_
              ring
      calc ν * _ ≤ ν * ∑ k ∈ range (K + 1), Δ ^ (2 * K + 2) * (c (K - k + 1) * ∑ a, N k a s) :=
            mul_le_mul_of_nonneg_left hsum hν
        _ = _ := by
            rw [← mul_sum]
            ring)
  exact key

end ExecutionProofs
