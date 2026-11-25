import ExecutionProofs
open Lean

/-- Axioms the development may use. Anything else, including `sorryAx`, fails the check. -/
def expectedAxioms : List Name := [``propext, ``Classical.choice, ``Quot.sound]

def isOurs (env : Environment) (n : Name) : Bool :=
  match env.getModuleIdxFor? n with
  | none => false
  | some idx =>
      match env.header.moduleNames[idx.toNat]? with
      | some m => (`ExecutionProofs).isPrefixOf m
      | none => false

run_cmd Elab.Command.liftCoreM do
  let env ← getEnv
  let mut checked : Nat := 0
  for (n, _) in env.constants.toList do
    unless isOurs env n && !n.isInternal && !n.hasMacroScopes do continue
    checked := checked + 1
    let axs ← collectAxioms n
    if axs.contains ``sorryAx then throwError "unproved: {n}"
    let bad := axs.filter (!expectedAxioms.contains ·)
    unless bad.isEmpty do throwError "unexpected axioms on {n}: {bad.toList}"
  if checked == 0 then throwError "module filter matched nothing"
  logInfo s!"{checked} declarations checked, axioms within {expectedAxioms}"
