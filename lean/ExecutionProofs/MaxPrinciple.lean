import Mathlib

/-!
# A maximum principle for jump generators

The semidiscrete heat equation of the paper lives on an infinite lattice (or on a finite box).
Both are covered by one abstract setting: nodes `S`, finitely many jump directions `D`,
neighbour maps `nb d : S → S` and nonnegative rates `κ d n`. A *barrier* `ξ ≥ 1` with finite
sublevel sets and `gen ξ ≤ c` replaces compactness.
-/

open Set

namespace ExecutionProofs

variable {S D : Type*} [Fintype D]

/-- Generator of a continuous-time jump process: `∑ d, κ d n * (u (nb d n) - u n)`. -/
def gen (nb : D → S → S) (κ : D → S → ℝ) (u : S → ℝ) (n : S) : ℝ :=
  ∑ d, κ d n * (u (nb d n) - u n)

variable {nb : D → S → S} {κ : D → S → ℝ}

lemma gen_sub (u v : S → ℝ) (n : S) :
    gen nb κ (fun m => u m - v m) n = gen nb κ u n - gen nb κ v n := by
  simp only [gen, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun d _ => ?_
  ring

lemma gen_add (u v : S → ℝ) (n : S) :
    gen nb κ (fun m => u m + v m) n = gen nb κ u n + gen nb κ v n := by
  simp only [gen, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun d _ => ?_
  ring

lemma gen_const_mul (a : ℝ) (u : S → ℝ) (n : S) :
    gen nb κ (fun m => a * u m) n = a * gen nb κ u n := by
  simp only [gen, Finset.mul_sum]
  refine Finset.sum_congr rfl fun d _ => ?_
  ring

lemma gen_const (a : ℝ) (n : S) : gen nb κ (fun _ => a) n = 0 := by
  simp [gen]

lemma gen_neg (u : S → ℝ) (n : S) : gen nb κ (fun m => -u m) n = -gen nb κ u n := by
  simpa using gen_const_mul (nb := nb) (κ := κ) (-1) u n

/-- At a global maximum the generator is nonpositive. -/
lemma gen_nonpos_of_max (hκ : ∀ d n, 0 ≤ κ d n) {u : S → ℝ} {n : S} (hmax : ∀ m, u m ≤ u n) :
    gen nb κ u n ≤ 0 :=
  Finset.sum_nonpos fun d _ =>
    mul_nonpos_of_nonneg_of_nonpos (hκ d n) (sub_nonpos.mpr (hmax (nb d n)))

/-- A function on `[0, T]` with a maximum at `s > 0` has nonnegative one-sided derivative there. -/
lemma deriv_nonneg_of_isMaxOn {f : ℝ → ℝ} {f' s T : ℝ} (hs : s ∈ Icc 0 T) (hpos : 0 < s)
    (hf : HasDerivWithinAt f f' (Icc 0 T) s) (hmax : IsMaxOn f (Icc 0 T) s) : 0 ≤ f' := by
  have hmem : (0 : ℝ) - s ∈ posTangentConeAt (Icc 0 T) s :=
    sub_mem_posTangentConeAt_of_segment_subset
      ((convex_Icc 0 T).segment_subset hs ⟨le_rfl, hs.1.trans hs.2⟩)
  have h := hmax.isLocalMaxOn.hasFDerivWithinAt_nonpos hf.hasFDerivWithinAt hmem
  have h' : (0 - s) * f' ≤ 0 := by simpa using h
  nlinarith

/-- **Maximum principle.** A bounded subsolution of `u' = gen u - V u` with `u(0) ≤ 0` stays
nonpositive. -/
theorem maximum_principle (hκ : ∀ d n, 0 ≤ κ d n) {V ξ : S → ℝ} (hV : ∀ n, 0 ≤ V n)
    (hξ1 : ∀ n, 1 ≤ ξ n) (hξfin : ∀ R : ℝ, {n | ξ n ≤ R}.Finite) {c : ℝ}
    (hξ : ∀ n, gen nb κ ξ n ≤ c) {T : ℝ} {u u' : S → ℝ → ℝ}
    (hu : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n) (u' n s) (Icc 0 T) s)
    {B : ℝ} (hB : ∀ n, ∀ s ∈ Icc 0 T, u n s ≤ B) (h0 : ∀ n, u n 0 ≤ 0)
    (hsub : ∀ n, ∀ s ∈ Icc 0 T, u' n s ≤ gen nb κ (fun m => u m s) n - V n * u n s) :
    ∀ n, ∀ s ∈ Icc 0 T, u n s ≤ 0 := by
  by_contra hcon
  push Not at hcon
  obtain ⟨n₀, s₀, hs₀, hpos⟩ := hcon
  set K : ℝ := |c| + 1 with hK
  have hKpos : 0 < K := by positivity
  have hT : (0 : ℝ) ≤ T := hs₀.1.trans hs₀.2
  set ε : ℝ := u n₀ s₀ / (2 * (Real.exp (K * s₀) * ξ n₀)) with hε
  have hξpos : ∀ n, 0 < ξ n := fun n => zero_lt_one.trans_le (hξ1 n)
  have hεpos : 0 < ε := by
    have := hξpos n₀
    positivity
  set w : S → ℝ → ℝ := fun n s => u n s - ε * (Real.exp (K * s) * ξ n) with hw
  have hw₀ : 0 < w n₀ s₀ := by
    have h1 : ε * (Real.exp (K * s₀) * ξ n₀) = u n₀ s₀ / 2 := by
      have := hξpos n₀
      rw [hε]
      field_simp
    have h2 : w n₀ s₀ = u n₀ s₀ - ε * (Real.exp (K * s₀) * ξ n₀) := rfl
    linarith
  have hwneg : ∀ n, ∀ s ∈ Icc 0 T, B / ε < ξ n → w n s < 0 := by
    intro n s hs hn
    have h1 : B < ε * ξ n := by
      rw [div_lt_iff₀ hεpos] at hn
      linarith
    have h2 : ε * ξ n ≤ ε * (Real.exp (K * s) * ξ n) := by
      have : 1 ≤ Real.exp (K * s) := Real.one_le_exp (mul_nonneg hKpos.le hs.1)
      have := hξpos n
      gcongr
      nlinarith
    have h3 := hB n s hs
    change u n s - ε * (Real.exp (K * s) * ξ n) < 0
    linarith
  have hwderiv : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (w n)
      (u' n s - ε * (K * Real.exp (K * s) * ξ n)) (Icc 0 T) s := by
    intro n s hs
    have hexp : HasDerivAt (fun s : ℝ => Real.exp (K * s)) (K * Real.exp (K * s)) s := by
      have := ((hasDerivAt_id s).const_mul K).exp
      simpa [mul_comm] using this
    exact (hu n s hs).sub ((hexp.mul_const (ξ n)).const_mul ε).hasDerivWithinAt
  have hmaxn : ∀ n, ∃ s ∈ Icc 0 T, IsMaxOn (w n) (Icc 0 T) s := fun n =>
    isCompact_Icc.exists_isMaxOn (nonempty_Icc.mpr hT)
      fun s hs => (hwderiv n s hs).continuousWithinAt
  choose σ hσ hσmax using hmaxn
  have hfin := hξfin (B / ε)
  have hn₀ : n₀ ∈ hfin.toFinset := by
    rw [hfin.mem_toFinset, mem_ofPred_eq]
    by_contra h
    exact absurd (hwneg n₀ s₀ hs₀ (not_le.mp h)) (not_lt.mpr hw₀.le)
  obtain ⟨n, -, hnmax⟩ := hfin.toFinset.exists_max_image (fun n => w n (σ n)) ⟨n₀, hn₀⟩
  set s := σ n with hsdef
  have hs : s ∈ Icc 0 T := hσ n
  have hwpos : 0 < w n s :=
    hw₀.trans_le ((hσmax n₀ hs₀).trans (hnmax n₀ hn₀))
  have hglobal : ∀ m, w m s ≤ w n s := by
    intro m
    by_cases hm : ξ m ≤ B / ε
    · exact (hσmax m hs).trans (hnmax m (hfin.mem_toFinset.mpr hm))
    · exact (hwneg m s hs (not_le.mp hm)).le.trans hwpos.le
  have hspos : 0 < s := by
    rcases hs.1.eq_or_lt with h | h
    · exfalso
      have h1 : w n 0 ≤ 0 := by
        have := h0 n
        have := hξpos n
        change u n 0 - ε * (Real.exp (K * 0) * ξ n) ≤ 0
        have : 0 < ε * (Real.exp (K * 0) * ξ n) := by positivity
        linarith
      rw [← h] at hwpos
      linarith
    · exact h
  have hder := deriv_nonneg_of_isMaxOn hs hspos (hwderiv n s hs) (hσmax n)
  have hgenw : gen nb κ (fun m => w m s) n ≤ 0 := gen_nonpos_of_max hκ hglobal
  have hgenu : gen nb κ (fun m => u m s) n
      = gen nb κ (fun m => w m s) n + ε * Real.exp (K * s) * gen nb κ ξ n := by
    have : (fun m => u m s) = fun m => w m s + (ε * Real.exp (K * s)) * ξ m := by
      funext m
      change u m s = u m s - ε * (Real.exp (K * s) * ξ m) + ε * Real.exp (K * s) * ξ m
      ring
    rw [this, gen_add, gen_const_mul]
  have hE : 0 < ε * Real.exp (K * s) := by positivity
  have hupos : 0 < u n s := by
    have h1 : 0 < ε * (Real.exp (K * s) * ξ n) := by
      have := hξpos n
      positivity
    have h2 : w n s = u n s - ε * (Real.exp (K * s) * ξ n) := rfl
    linarith
  have hVu : 0 ≤ V n * u n s := mul_nonneg (hV n) hupos.le
  have hcξ : ε * Real.exp (K * s) * gen nb κ ξ n ≤ ε * Real.exp (K * s) * c :=
    mul_le_mul_of_nonneg_left (hξ n) hE.le
  have hcK : c < K * ξ n := by
    have h1 : c ≤ |c| := le_abs_self c
    have h2 : K * 1 ≤ K * ξ n := mul_le_mul_of_nonneg_left (hξ1 n) hKpos.le
    linarith
  have hstrict : ε * Real.exp (K * s) * c < ε * Real.exp (K * s) * (K * ξ n) :=
    mul_lt_mul_of_pos_left hcK hE
  have hmain := hsub n s hs
  have hrewrite : ε * (K * Real.exp (K * s) * ξ n) = ε * Real.exp (K * s) * (K * ξ n) := by ring
  linarith

/-- **Comparison with a supersolution in time.** If `u' ≤ gen u - V u + g`, `G' = g`, `G ≥ 0` and
`u(0) ≤ G(0)`, then `u ≤ G`. -/
theorem le_of_subsolution (hκ : ∀ d n, 0 ≤ κ d n) {V ξ : S → ℝ} (hV : ∀ n, 0 ≤ V n)
    (hξ1 : ∀ n, 1 ≤ ξ n) (hξfin : ∀ R : ℝ, {n | ξ n ≤ R}.Finite) {c : ℝ}
    (hξ : ∀ n, gen nb κ ξ n ≤ c) {T : ℝ} {u u' : S → ℝ → ℝ}
    (hu : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n) (u' n s) (Icc 0 T) s)
    {B : ℝ} (hB : ∀ n, ∀ s ∈ Icc 0 T, u n s ≤ B) {G g : ℝ → ℝ}
    (hG : ∀ s ∈ Icc 0 T, HasDerivWithinAt G (g s) (Icc 0 T) s) (hGnn : ∀ s ∈ Icc 0 T, 0 ≤ G s)
    (h0 : ∀ n, u n 0 ≤ G 0)
    (hsub : ∀ n, ∀ s ∈ Icc 0 T, u' n s ≤ gen nb κ (fun m => u m s) n - V n * u n s + g s) :
    ∀ n, ∀ s ∈ Icc 0 T, u n s ≤ G s := by
  have key := maximum_principle (u := fun n s => u n s - G s) (u' := fun n s => u' n s - g s)
    hκ hV hξ1 hξfin hξ (fun n s hs => (hu n s hs).sub (hG s hs)) (B := B)
    (fun n s hs => by
      have h1 := hB n s hs
      have h2 := hGnn s hs
      linarith)
    (fun n => by
      have := h0 n
      linarith)
    (fun n s hs => by
      have h1 := hsub n s hs
      have h2 : 0 ≤ V n * G s := mul_nonneg (hV n) (hGnn s hs)
      have h3 : gen nb κ (fun m => u m s - G s) n = gen nb κ (fun m => u m s) n := by
        rw [gen_sub (fun m => u m s) (fun _ => G s), gen_const, sub_zero]
      rw [h3]
      nlinarith)
  intro n s hs
  have := key n s hs
  linarith

/-- Two-sided version: `|u| ≤ G` when `u' = gen u - V u + r` with `|r| ≤ g`. -/
theorem abs_le_of_residual (hκ : ∀ d n, 0 ≤ κ d n) {V ξ : S → ℝ} (hV : ∀ n, 0 ≤ V n)
    (hξ1 : ∀ n, 1 ≤ ξ n) (hξfin : ∀ R : ℝ, {n | ξ n ≤ R}.Finite) {c : ℝ}
    (hξ : ∀ n, gen nb κ ξ n ≤ c) {T : ℝ} {u r : S → ℝ → ℝ}
    (hu : ∀ n, ∀ s ∈ Icc 0 T, HasDerivWithinAt (u n)
      (gen nb κ (fun m => u m s) n - V n * u n s + r n s) (Icc 0 T) s)
    {B : ℝ} (hB : ∀ n, ∀ s ∈ Icc 0 T, |u n s| ≤ B) {G g : ℝ → ℝ}
    (hG : ∀ s ∈ Icc 0 T, HasDerivWithinAt G (g s) (Icc 0 T) s) (hGnn : ∀ s ∈ Icc 0 T, 0 ≤ G s)
    (h0 : ∀ n, u n 0 = 0) (hr : ∀ n, ∀ s ∈ Icc 0 T, |r n s| ≤ g s) :
    ∀ n, ∀ s ∈ Icc 0 T, |u n s| ≤ G s := by
  intro n s hs
  have hG0 : 0 ≤ G 0 := hGnn 0 ⟨le_rfl, hs.1.trans hs.2⟩
  rw [abs_le]
  constructor
  · have key := le_of_subsolution (u := fun n s => -u n s) hκ hV hξ1 hξfin hξ
      (fun n s hs => (hu n s hs).neg) (B := B) (fun n s hs => (neg_le_abs _).trans (hB n s hs))
      hG hGnn (fun n => by simpa [h0 n] using hG0)
      (fun n s hs => by
        have h1 := neg_le_abs (r n s)
        have h2 := hr n s hs
        rw [gen_neg (fun m => u m s)]
        linarith) n s hs
    linarith
  · exact le_of_subsolution hκ hV hξ1 hξfin hξ hu (fun n s hs => (le_abs_self _).trans (hB n s hs))
      hG hGnn (fun n => by simpa [h0 n] using hG0)
      (fun n s hs => by
        have h1 := le_abs_self (r n s)
        have h2 := hr n s hs
        linarith) n s hs

end ExecutionProofs
