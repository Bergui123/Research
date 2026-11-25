import ExecutionProofs.Rounding

/-!
# Box truncation bound (paper, Proposition 4)

The circuit integrates a periodic box, the grid theorem speaks about the infinite lattice. On the box
the two solutions differ only through the neighbours that leave the box: the periodic problem wraps
them, the lattice problem does not. If both kinds of neighbour values stay in `[0, g]`, the solution
`U` of the box problem with zero data and ghost value `g` bounds the difference.
-/

open Set

namespace ExecutionProofs

variable {S D : Type*} [Fintype D] [Finite S] {nb : D → S → S} {κ : D → S → ℝ}

/-- **Box truncation.** `u` solves the periodic problem on the finite box `S`; `v` is the restriction
of the lattice solution, with outside values `w`; `inside d n` says that the neighbour of `n` in
direction `d` belongs to the box. -/
theorem box_truncation (hκ : ∀ d n, 0 ≤ κ d n) {V : S → ℝ} (hV : ∀ n, 0 ≤ V n)
    (inside : D → S → Prop) [∀ d n, Decidable (inside d n)] {T g : ℝ}
    {u v U : S → ℝ → ℝ} {w : D → S → ℝ → ℝ}
    (hu : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n)
      (∑ d, κ d n * (u (nb d n) s - u n s) - V n * u n s) (Icc 0 T) s)
    (hv : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (v n)
      (∑ d, κ d n * ((if inside d n then v (nb d n) s else w d n s) - v n s) - V n * v n s)
      (Icc 0 T) s)
    (hU : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (U n)
      (∑ d, κ d n * ((if inside d n then U (nb d n) s else g) - U n s) - V n * U n s) (Icc 0 T) s)
    (hwrap : ∀ d n, ¬ inside d n → ∀ s ∈ Icc 0 T, u (nb d n) s ∈ Icc 0 g)
    (hout : ∀ d n, ¬ inside d n → ∀ s ∈ Icc 0 T, w d n s ∈ Icc 0 g)
    (h0 : ∀ n, u n 0 = v n 0) (hU0 : ∀ n, U n 0 = 0) :
    ∀ n, ∀ s ∈ Icc 0 T, |u n s - v n s| ≤ U n s := by
  classical
  cases nonempty_fintype S
  -- rates inside the box, and leaving the box as killing
  set κ' : D → S → ℝ := fun d n => if inside d n then κ d n else 0 with hκ'
  set V' : S → ℝ := fun n => V n + ∑ d, if inside d n then 0 else κ d n with hV'
  set r : S → ℝ → ℝ := fun n s => ∑ d, if inside d n then 0 else κ d n * (u (nb d n) s - w d n s)
    with hr
  set q : S → ℝ → ℝ := fun n _ => ∑ d, if inside d n then 0 else κ d n * g with hq
  have hcont : ∀ f f' : S → ℝ → ℝ, (∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (f n) (f' n s) (Icc 0 T) s) →
      ∃ B, ∀ n, ∀ s ∈ Icc 0 T, |f n s| ≤ B := by
    intro f f' hf
    have h1 : ∀ n, ∃ B, ∀ s ∈ Icc 0 T, |f n s| ≤ B := fun n =>
      isCompact_Icc.exists_bound_of_continuousOn fun s hs => (hf n s hs).continuousWithinAt
    choose B hB using h1
    exact ⟨∑ n, |B n|, fun n s hs => ((hB n s hs).trans (le_abs_self _)).trans
      (Finset.single_le_sum (f := fun n => |B n|) (fun _ _ => abs_nonneg _) (Finset.mem_univ n))⟩
  have he : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (fun s => u n s - v n s)
      (gen nb κ' (fun m => u m s - v m s) n - V' n * (u n s - v n s) + r n s) (Icc 0 T) s := by
    intro n s hs
    refine ((hu n s hs).fun_sub (hv n s hs)).congr_deriv ?_
    simp only [gen, hκ', hV', hr, add_mul, Finset.sum_mul]
    rw [← sub_eq_zero]
    have : ∑ d, κ d n * (u (nb d n) s - u n s)
        - ∑ d, κ d n * ((if inside d n then v (nb d n) s else w d n s) - v n s)
        - ∑ d, (if inside d n then κ d n else 0) * (u (nb d n) s - v (nb d n) s - (u n s - v n s))
        + ∑ d, (if inside d n then 0 else κ d n) * (u n s - v n s)
        - ∑ d, (if inside d n then 0 else κ d n * (u (nb d n) s - w d n s)) = 0 := by
      rw [← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib,
        ← Finset.sum_sub_distrib]
      refine Finset.sum_eq_zero fun d _ => ?_
      split_ifs <;> ring
    linarith
  have hU' : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (U n)
      (gen nb κ' (fun m => U m s) n - V' n * U n s + q n s) (Icc 0 T) s := by
    intro n s hs
    refine (hU n s hs).congr_deriv ?_
    simp only [gen, hκ', hV', hq, add_mul, Finset.sum_mul]
    rw [← sub_eq_zero]
    have : ∑ d, κ d n * ((if inside d n then U (nb d n) s else g) - U n s)
        - ∑ d, (if inside d n then κ d n else 0) * (U (nb d n) s - U n s)
        + ∑ d, (if inside d n then 0 else κ d n) * U n s
        - ∑ d, (if inside d n then 0 else κ d n * g) = 0 := by
      rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
      refine Finset.sum_eq_zero fun d _ => ?_
      split_ifs <;> ring
    linarith
  obtain ⟨B₁, hB₁⟩ := hcont _ _ he
  obtain ⟨B₂, hB₂⟩ := hcont _ _ hU'
  refine abs_le_of_dominated (nb := nb) (κ := κ') (V := V') (ξ := fun _ => 1) (c := 0) (r := r) (q := q)
    (e := fun n s => u n s - v n s) (B := max B₁ B₂) ?_ ?_ (fun _ => le_rfl)
    (fun _ => Set.toFinite _) (fun n => (gen_const 1 n).le) he hU' ?_
    (fun n s hs => (hB₁ n s hs).trans (le_max_left _ _))
    (fun n s hs => (hB₂ n s hs).trans (le_max_right _ _)) ?_
  · intro d n
    simp only [hκ']
    split_ifs
    · exact hκ d n
    · exact le_rfl
  · intro n
    refine add_nonneg (hV n) (Finset.sum_nonneg fun d _ => ?_)
    split_ifs
    · exact le_rfl
    · exact hκ d n
  · intro n s hs
    refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun d _ => ?_)
    split_ifs with hd
    · simp
    · rw [abs_mul, abs_of_nonneg (hκ d n)]
      refine mul_le_mul_of_nonneg_left ?_ (hκ d n)
      have h1 := hwrap d n hd s hs
      have h2 := hout d n hd s hs
      rw [abs_le]
      constructor <;> linarith [h1.1, h1.2, h2.1, h2.2]
  · intro n
    simp [h0 n, hU0 n]

end ExecutionProofs
