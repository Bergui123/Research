import ExecutionProofs.MaxPrinciple

/-!
# Potential rounding (paper, Lemma 8) and the comparison behind Proposition 4

Consequences of the maximum principle, valid on the infinite lattice and on finite boxes.
-/

open Set

namespace ExecutionProofs

variable {S D : Type*} [Fintype D] {nb : D → S → S} {κ : D → S → ℝ}

/-- Nonnegative data stay nonnegative. -/
theorem nonneg_of_nonneg (hκ : ∀ d n, 0 ≤ κ d n) {V ξ : S → ℝ} (hV : ∀ n, 0 ≤ V n)
    (hξ1 : ∀ n, 1 ≤ ξ n) (hξfin : ∀ R : ℝ, {n | ξ n ≤ R}.Finite) {c : ℝ}
    (hξ : ∀ n, gen nb κ ξ n ≤ c) {T : ℝ} {u : S → ℝ → ℝ}
    (hu : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n)
      (gen nb κ (fun m => u m s) n - V n * u n s) (Icc 0 T) s)
    {B : ℝ} (hB : ∀ n, ∀ s ∈ Icc 0 T, |u n s| ≤ B) (h0 : ∀ n, 0 ≤ u n 0) :
    ∀ n, ∀ s ∈ Icc 0 T, 0 ≤ u n s := by
  intro n s hs
  have key := maximum_principle (u := fun n s => -u n s) hκ hV hξ1 hξfin hξ
    (fun n s hs => (hu n s hs).neg) (B := B) (fun n s hs => (neg_le_abs _).trans (hB n s hs))
    (fun n => neg_nonpos.mpr (h0 n))
    (fun n s hs => by
      rw [gen_neg (fun m => u m s)]
      linarith) n s hs
  exact neg_nonpos.mp key

/-- **Effect of rounding.** If `0 ≤ V - V_m ≤ δ` and the data are nonnegative and equal, then
`ψ ≤ ψ_m ≤ exp (δ s) ψ`. -/
theorem rounding (hκ : ∀ d n, 0 ≤ κ d n) {V Vm ξ : S → ℝ} (hVm : ∀ n, 0 ≤ Vm n) {δ : ℝ}
    (hlow : ∀ n, 0 ≤ V n - Vm n) (hhigh : ∀ n, V n - Vm n ≤ δ)
    (hξ1 : ∀ n, 1 ≤ ξ n) (hξfin : ∀ R : ℝ, {n | ξ n ≤ R}.Finite) {c : ℝ}
    (hξ : ∀ n, gen nb κ ξ n ≤ c) {T : ℝ} {u um : S → ℝ → ℝ}
    (hu : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n)
      (gen nb κ (fun m => u m s) n - V n * u n s) (Icc 0 T) s)
    (hum : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (um n)
      (gen nb κ (fun m => um m s) n - Vm n * um n s) (Icc 0 T) s)
    {B : ℝ} (hB : ∀ n, ∀ s ∈ Icc 0 T, |u n s| ≤ B) (hBm : ∀ n, ∀ s ∈ Icc 0 T, |um n s| ≤ B)
    (h0 : ∀ n, 0 ≤ u n 0) (heq : ∀ n, um n 0 = u n 0) :
    ∀ n, ∀ s ∈ Icc 0 T, u n s ≤ um n s ∧ um n s ≤ Real.exp (δ * s) * u n s := by
  have hV : ∀ n, 0 ≤ V n := fun n => by
    have h1 := hlow n
    have h2 := hVm n
    linarith
  have hnn := nonneg_of_nonneg hκ hV hξ1 hξfin hξ hu hB h0
  intro n s hs
  have hδ : 0 ≤ δ := (hlow n).trans (hhigh n)
  constructor
  · have key := maximum_principle (u := fun n s => u n s - um n s) hκ hVm hξ1 hξfin hξ
      (fun n s hs => (hu n s hs).sub (hum n s hs)) (B := B + B)
      (fun n s hs => by
        have h1 := (le_abs_self _).trans (hB n s hs)
        have h2 := (neg_le_abs _).trans (hBm n s hs)
        linarith)
      (fun n => by
        have := heq n
        linarith)
      (fun n s hs => by
        have h1 : 0 ≤ (V n - Vm n) * u n s := mul_nonneg (hlow n) (hnn n s hs)
        rw [gen_sub (fun m => u m s) (fun m => um m s)]
        nlinarith) n s hs
    linarith
  · have hexp : ∀ s : ℝ, HasDerivAt (fun s : ℝ => Real.exp (δ * s)) (δ * Real.exp (δ * s)) s := by
      intro s
      have := ((hasDerivAt_id s).const_mul δ).exp
      simpa [mul_comm] using this
    have key := maximum_principle (u := fun n s => um n s - Real.exp (δ * s) * u n s)
      hκ hVm hξ1 hξfin hξ
      (fun n s hs => (hum n s hs).sub ((hexp s).hasDerivWithinAt.mul (hu n s hs))) (B := B)
      (fun n s hs => by
        have h1 := (le_abs_self _).trans (hBm n s hs)
        have h2 : 0 ≤ Real.exp (δ * s) * u n s :=
          mul_nonneg (Real.exp_pos _).le (hnn n s hs)
        linarith)
      (fun n => by simp [heq n])
      (fun n s hs => by
        have hE := (Real.exp_pos (δ * s)).le
        have h1 : 0 ≤ (δ - (V n - Vm n)) * (Real.exp (δ * s) * u n s) :=
          mul_nonneg (sub_nonneg.mpr (hhigh n)) (mul_nonneg hE (hnn n s hs))
        rw [gen_sub (fun m => um m s) (fun m => Real.exp (δ * s) * u m s),
          gen_const_mul (Real.exp (δ * s)) (fun m => u m s)]
        nlinarith) n s hs
    linarith

/-- **A posteriori comparison** (used for the box truncation bound): if the residual of `e` is
dominated node by node by the source of `U`, then `|e| ≤ U`. -/
theorem abs_le_of_dominated (hκ : ∀ d n, 0 ≤ κ d n) {V ξ : S → ℝ} (hV : ∀ n, 0 ≤ V n)
    (hξ1 : ∀ n, 1 ≤ ξ n) (hξfin : ∀ R : ℝ, {n | ξ n ≤ R}.Finite) {c : ℝ}
    (hξ : ∀ n, gen nb κ ξ n ≤ c) {T : ℝ} {e r U q : S → ℝ → ℝ}
    (he : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (e n)
      (gen nb κ (fun m => e m s) n - V n * e n s + r n s) (Icc 0 T) s)
    (hU : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (U n)
      (gen nb κ (fun m => U m s) n - V n * U n s + q n s) (Icc 0 T) s)
    (hr : ∀ n, ∀ s ∈ Icc 0 T, |r n s| ≤ q n s)
    {B : ℝ} (hB : ∀ n, ∀ s ∈ Icc 0 T, |e n s| ≤ B) (hBU : ∀ n, ∀ s ∈ Icc 0 T, |U n s| ≤ B)
    (h0 : ∀ n, |e n 0| ≤ U n 0) :
    ∀ n, ∀ s ∈ Icc 0 T, |e n s| ≤ U n s := by
  intro n s hs
  rw [abs_le]
  constructor
  · have key := maximum_principle (u := fun n s => -e n s - U n s) hκ hV hξ1 hξfin hξ
      (fun n s hs => (he n s hs).neg.sub (hU n s hs)) (B := B + B)
      (fun n s hs => by
        have h1 := (neg_le_abs _).trans (hB n s hs)
        have h2 := (neg_le_abs _).trans (hBU n s hs)
        linarith)
      (fun n => by
        have := (neg_le_abs _).trans (h0 n)
        linarith)
      (fun n s hs => by
        have h1 := (neg_le_abs _).trans (hr n s hs)
        rw [gen_sub (fun m => -e m s) (fun m => U m s), gen_neg (fun m => e m s)]
        linarith) n s hs
    linarith
  · have key := maximum_principle (u := fun n s => e n s - U n s) hκ hV hξ1 hξfin hξ
      (fun n s hs => (he n s hs).sub (hU n s hs)) (B := B + B)
      (fun n s hs => by
        have h1 := (le_abs_self _).trans (hB n s hs)
        have h2 := (neg_le_abs _).trans (hBU n s hs)
        linarith)
      (fun n => by
        have := (le_abs_self _).trans (h0 n)
        linarith)
      (fun n s hs => by
        have h1 := (le_abs_self _).trans (hr n s hs)
        rw [gen_sub (fun m => e m s) (fun m => U m s)]
        linarith) n s hs
    linarith

end ExecutionProofs
