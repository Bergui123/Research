import Mathlib

/-!
# Positive semidefinite forms: two trace inequalities, and the derivative of an inverse
-/

open Set Matrix

open scoped MatrixOrder

namespace ExecutionProofs

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Induction for real positive semidefinite matrices: sums of `sᵀ s`. -/
lemma psd_induction {A : Matrix ι ι ℝ} (hA : A.IsSymm) (hA0 : ∀ v, 0 ≤ v ⬝ᵥ (A *ᵥ v))
    {p : Matrix ι ι ℝ → Prop} (h0 : p 0) (hadd : ∀ M N, p M → p N → p (M + N))
    (hsq : ∀ s : Matrix ι ι ℝ, p (sᵀ * s)) : p A := by
  have hherm : A.IsHermitian := by
    rw [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial]
    exact hA
  have hpsd : A.PosSemidef := Matrix.PosSemidef.of_dotProduct_mulVec_nonneg hherm
    (fun x => by simpa using hA0 x)
  have hnn : 0 ≤ A := Matrix.nonneg_iff_posSemidef.mpr hpsd
  have hmem := StarOrderedRing.nonneg_iff.mp hnn
  refine AddSubmonoid.closure_induction (motive := fun M _ => p M) ?_ h0
    (fun M N _ _ hM hN => hadd M N hM hN) hmem
  rintro _ ⟨s, rfl⟩
  have : star s = sᵀ := by
    rw [Matrix.star_eq_conjTranspose, Matrix.conjTranspose_eq_transpose_of_trivial]
  simp only [this]
  exact hsq s

omit [DecidableEq ι] in
lemma sum_gram (s B : Matrix ι ι ℝ) :
    ∑ i, ∑ j, (sᵀ * s) i j * B i j = ∑ k, (s k) ⬝ᵥ (B *ᵥ (s k)) := by
  simp only [Matrix.mul_apply, Matrix.transpose_apply, dotProduct, Matrix.mulVec, Finset.sum_mul,
    Finset.mul_sum]
  calc ∑ i, ∑ j, ∑ k, s k i * s k j * B i j = ∑ i, ∑ k, ∑ j, s k i * s k j * B i j :=
        Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ k, ∑ i, ∑ j, s k i * s k j * B i j := Finset.sum_comm
    _ = _ := Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun i _ =>
        Finset.sum_congr rfl fun j _ => by ring

/-- The entrywise product of two positive semidefinite forms has nonnegative sum. -/
lemma sum_mul_nonneg {A B : Matrix ι ι ℝ} (hA : A.IsSymm) (hA0 : ∀ v, 0 ≤ v ⬝ᵥ (A *ᵥ v))
    (hB0 : ∀ v, 0 ≤ v ⬝ᵥ (B *ᵥ v)) : 0 ≤ ∑ i, ∑ j, A i j * B i j := by
  refine psd_induction hA hA0 (p := fun A => 0 ≤ ∑ i, ∑ j, A i j * B i j) (by simp) ?_ ?_
  · intro M N hM hN
    simp only [Matrix.add_apply, add_mul, Finset.sum_add_distrib]
    exact add_nonneg hM hN
  · intro s
    rw [sum_gram]
    exact Finset.sum_nonneg fun k _ => hB0 _

/-- Cauchy–Schwarz for a form and its inverse. -/
lemma cauchy_form {P R : Matrix ι ι ℝ} (hPs : P.IsSymm) (hP0 : ∀ v, 0 ≤ v ⬝ᵥ (P *ᵥ v))
    (hPR : P * R = 1) (x l : ι → ℝ) :
    (l ⬝ᵥ x) ^ 2 ≤ (x ⬝ᵥ (P *ᵥ x)) * (l ⬝ᵥ (R *ᵥ l)) := by
  have hl : P *ᵥ (R *ᵥ l) = l := by
    rw [Matrix.mulVec_mulVec, hPR, Matrix.one_mulVec]
  have hsym : (R *ᵥ l) ⬝ᵥ (P *ᵥ x) = x ⬝ᵥ l := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hPs, hl, dotProduct_comm]
  have h : ∀ t : ℝ, 0 ≤ (l ⬝ᵥ (R *ᵥ l)) * (t * t) + (2 * (l ⬝ᵥ x)) * t + x ⬝ᵥ (P *ᵥ x) := by
    intro t
    have h1 := hP0 (x + t • (R *ᵥ l))
    rw [Matrix.mulVec_add, Matrix.mulVec_smul, hl, add_dotProduct, dotProduct_add, dotProduct_add,
      dotProduct_smul, smul_dotProduct, smul_dotProduct, dotProduct_smul, hsym,
      dotProduct_comm (R *ᵥ l) l, dotProduct_comm x l] at h1
    simp only [smul_eq_mul] at h1
    linarith
  have hd := discrim_le_zero h
  rw [discrim] at hd
  nlinarith

/-- `x · M x ≤ (x · P x) ∑ R_ij M_ij` for `P R = 1` and `M` positive semidefinite. -/
lemma form_le_trace {P R M : Matrix ι ι ℝ} (hPs : P.IsSymm) (hP0 : ∀ v, 0 ≤ v ⬝ᵥ (P *ᵥ v))
    (hPR : P * R = 1) (hM : M.IsSymm) (hM0 : ∀ v, 0 ≤ v ⬝ᵥ (M *ᵥ v)) (x : ι → ℝ) :
    x ⬝ᵥ (M *ᵥ x) ≤ (x ⬝ᵥ (P *ᵥ x)) * ∑ i, ∑ j, M i j * R i j := by
  refine psd_induction hM hM0
    (p := fun M => x ⬝ᵥ (M *ᵥ x) ≤ (x ⬝ᵥ (P *ᵥ x)) * ∑ i, ∑ j, M i j * R i j) (by simp) ?_ ?_
  · intro M N hM hN
    simp only [Matrix.add_apply, add_mul, Finset.sum_add_distrib, Matrix.add_mulVec,
      dotProduct_add, mul_add]
    exact add_le_add hM hN
  · intro s
    rw [sum_gram, Finset.mul_sum]
    have h1 : x ⬝ᵥ ((sᵀ * s) *ᵥ x) = ∑ k, ((s k) ⬝ᵥ x) ^ 2 := by
      rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
        Matrix.transpose_transpose]
      simp only [dotProduct, Matrix.mulVec, sq]
    rw [h1]
    exact Finset.sum_le_sum fun k _ => cauchy_form hPs hP0 hPR x (s k)

end ExecutionProofs
