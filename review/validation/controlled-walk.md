# Small controlled walk: circuit and cost check

`check_controlled_walk.py` uses the compiled gate operations in
`simulation/select_circuit.py` for `n=2`, `m=2`, `frac=2`, with the same
benchmark potential and comparator scale as `block_encoding_check`. This is a
small consistency check, not a compilation of the whole algorithm. The
measurement-based AND uncomputation and Qualtran multi-control cost use the
same assumptions as the resource model. General rotations are ideal matrices;
GQSP phases, state preparation, the QFT and the amplitude-estimation circuit
are not compiled or checked here.

## Controlled SELECT

For the kinetic branch, a new work bit computes `control AND NOT t` and is
uncomputed after the controlled shifts. For the potential branch, the
arithmetic is computed unconditionally, a new work bit computes
`control AND t`, and only its phase is applied before both work bits and all
arithmetic are uncomputed. Thus control zero gives the identity and control
one gives the existing Hermitian SELECT. The two new logical AND computations
add 8 T under the paper's convention. The exhaustive unit test checks all
1,024 basis states (both control values), checks signs and work-register
cleanup, then applies the operation a second time and checks the identity.
The original no-control path and its resource counts are unchanged.

## Walk, inverse and modeled cost

Let `A = PREPARE`, `R = 2|0><0| - I` on the LCU ancillas, and `S = SELECT`.
The walk is `W = A R A^dagger S`; its inverse is
`W^dagger = S^dagger A R A^dagger`. The controlled walk can use unconditional
`A` and `A^dagger` around controlled `R`, followed by the controlled `S`.
When the control is zero, the unconditional preparation cancels. With one
control and `3+m` LCU ancillas, controlled `R` is a single controlled-Z on
`4+m` bits, with a Clifford Z on the control for the sign convention. The
checker constructs these matrices and compares controlled `W` and its inverse
against the independent ideal SELECT formula for both control values. Its
largest observed `W^dagger W - I` entry is below `1.6e-15`.

For this instance the compiled SELECT has 245 logical ANDs; controlled SELECT
has 247, or 988 T under the model. Qualtran charges 16 T for the controlled
reflection. Their combined 1,004 T is below the 1,020 T charged by the model
for the original SELECT and two larger multi-controlled reflections. The walk
uses two `R_y` preparations; the model charges four synthesized rotations per
step. This comparison does not show that unused allowance covers other parts
of the GQSP or QAE construction. It also does not independently establish
measurement-based feedforward, rotation synthesis or the full
controlled-state-preparation cost.

`check_controlled_select_profiles.py` separately runs the static gate builder
for the 18 saved benchmark and sweep configurations. For each, controlled
SELECT has two more ANDs than the saved uncontrolled circuit and the same
peak work-qubit count (with the external control counted as an input). The
controlled SELECT plus one Qualtran-costed reflection is at least 56 T below
the model's SELECT-plus-two-reflections allowance, before rotation costs.
These are counts of fixed gate builders on one basis input, not an exhaustive
correctness or phase-synthesis check at those sizes. The output is saved in
`controlled-profiles.json`.

Reproduce from the repository root:

```sh
PYTHONPATH=simulation OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s simulation/tests -p 'test_controlled_select.py' -v
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_controlled_walk.py > /tmp/controlled-walk.json
cmp /tmp/controlled-walk.json review/validation/controlled-walk.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python review/validation/check_controlled_select_profiles.py > /tmp/controlled-profiles.json
cmp /tmp/controlled-profiles.json review/validation/controlled-profiles.json
```
