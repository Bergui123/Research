import Mathlib

/-!
# The comparator potential (paper, Lemma 7)

`B = ∑ r |r⟩⟨r| ⊗ (-1)^[r < q]` on an `m`-qubit register, averaged over the uniform state,
encodes `1 - 2 q / 2^m`; with `q` the truncation of `2^m V / V'` the encoded potential is
below `V` by less than `V' / 2^m`.
-/

open Finset

namespace ExecutionProofs

/-- Sum of the comparator signs over the register. -/
lemma sum_sign {N q : ℕ} (hq : q ≤ N) :
    ∑ r ∈ range N, (if r < q then (-1 : ℝ) else 1) = N - 2 * q := by
  have h : ∀ r ∈ range N, (if r < q then (-1 : ℝ) else 1) = 1 - 2 * (if r < q then 1 else 0) := by
    intro r _
    split_ifs <;> norm_num
  have hcard : ((range N).filter (· < q)).card = q := by
    have : (range N).filter (· < q) = range q := by
      ext r
      simp only [mem_filter, mem_range]
      omega
    rw [this, card_range]
  rw [sum_congr rfl h, sum_sub_distrib, ← mul_sum, sum_boole, hcard]
  simp

/-- **Comparator potential.** Expectation of the comparator phase in the uniform state. -/
theorem comparator_expectation {m q : ℕ} (hq : q ≤ 2 ^ m) :
    ((2 : ℝ) ^ m)⁻¹ * ∑ r ∈ range (2 ^ m), (if r < q then (-1 : ℝ) else 1)
      = 1 - 2 * q / 2 ^ m := by
  rw [sum_sign hq]
  push_cast
  field_simp

/-- The encoded potential `V' / 2 * (1 - ⟨B⟩)` equals `V' q / 2^m`. -/
theorem comparator_potential {m q : ℕ} (hq : q ≤ 2 ^ m) (V' : ℝ) :
    V' / 2 * (1 - ((2 : ℝ) ^ m)⁻¹ * ∑ r ∈ range (2 ^ m), (if r < q then (-1 : ℝ) else 1))
      = V' * q / 2 ^ m := by
  rw [comparator_expectation hq]
  ring

/-- **Truncation error.** With `q = ⌊2^m V / V'⌋` the encoded potential satisfies
`0 ≤ V - V_m < V' / 2^m`. -/
theorem comparator_truncation {m : ℕ} {V V' : ℝ} (hV : 0 ≤ V) (hV' : 0 < V') :
    0 ≤ V - V' * (⌊2 ^ m * V / V'⌋₊ : ℝ) / 2 ^ m ∧
      V - V' * (⌊2 ^ m * V / V'⌋₊ : ℝ) / 2 ^ m < V' / 2 ^ m := by
  have hp : (0 : ℝ) < 2 ^ m := by positivity
  have hx : (0 : ℝ) ≤ 2 ^ m * V / V' := by positivity
  have h1 := Nat.floor_le hx
  have h2 := Nat.lt_floor_add_one (2 ^ m * V / V')
  rw [le_div_iff₀ hV'] at h1
  rw [div_lt_iff₀ hV'] at h2
  constructor
  · rw [sub_nonneg, div_le_iff₀ hp]
    linarith
  · rw [sub_lt_iff_lt_add, ← add_div, lt_div_iff₀ hp]
    linarith

end ExecutionProofs
