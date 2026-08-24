#import "../src/lib.typ": smiles, mol, mol-formula, rxn-arrow, reaction, cycle, step, atom, bond, lp, species, arrow, highlight

#let selected-case = sys.inputs.at("case", default: "")

#if selected-case == "invalid-smiles" {
  smiles("C1CC")
} else if selected-case == "unclosed-label" {
  smiles("C{OH")
} else if selected-case == "missing-arrow-endpoints" {
  smiles("CO", arrow())
} else if selected-case == "species-out-of-range" {
  reaction(
    mol("CCO"),
    rxn-arrow(),
    mol("CCO"),
    arrow(from: atom(2, 0), to: atom(1, 1)),
  )
} else if selected-case == "atom-out-of-range" {
  smiles("CO", arrow(from: atom(4), to: atom(0)))
} else if selected-case == "missing-bond" {
  smiles("CCC", highlight(bond(0, 2)))
} else if selected-case == "missing-lone-pair" {
  smiles("C", arrow(from: lp(0), to: atom(0)))
} else if selected-case == "pair-out-of-range" {
  smiles("CO", arrow(from: lp(1, pair: 2), to: atom(0)))
} else if selected-case == "opaque-atom-reference" {
  reaction(
    mol([opaque]),
    mol("CO"),
    arrow(from: atom(0, 0), to: atom(1, 0)),
  )
} else if selected-case == "ignored-annotation" {
  smiles("CO", [not an annotation])
} else if selected-case == "show-h-out-of-range" {
  smiles("CO", show-h: 5)
} else if selected-case == "mol-formula-wrong-type" {
  mol-formula(1)
} else if selected-case == "mol-formula-empty" {
  mol-formula("")
} else if selected-case == "mol-formula-wildcard" {
  mol-formula("*CC")
} else if selected-case == "annotation-out-of-range" {
  smiles("CO", atom-annotations: ((5, [note]),))
} else if selected-case == "customized-missing-bond" {
  smiles(
    "CCC",
    bond-customizations: ((bond(0, 2), (color: red)),),
  )
} else if selected-case == "duplicate-bond-customization" {
  smiles(
    "CC",
    bond-customizations: (
      (bond(0, 1), (color: red)),
      (bond(1, 0), (opacity: 50%)),
    ),
  )
} else if selected-case == "opacity-out-of-range" {
  smiles("CO", opacity: 140%)
} else if selected-case == "unknown-mol-option" {
  reaction(mol("CO", rotate: 30deg))
} else if selected-case == "content-molecule-options" {
  reaction(mol([opaque], color: false))
} else if selected-case == "invalid-reaction-item" {
  reaction(step())
} else if selected-case == "empty-reaction" {
  reaction()
} else if selected-case == "empty-cycle" {
  cycle()
} else if selected-case == "leading-cycle-step" {
  cycle(step(), "CO")
} else if selected-case == "duplicate-cycle-step" {
  cycle("CO", step(), step(), "CC")
} else if selected-case == "invalid-step-reagent" {
  cycle("CO", step(into: mol("C")), "CC")
} else {
  panic("unknown validation test case: " + selected-case)
}
