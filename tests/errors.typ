#import "../src/lib.typ": smiles, smiles-cetz, mol, mol-formula, rxn-arrow, reaction, cycle, step, atom, bond, lp, species, arrow, highlight, molecules, substructure-matches
#import "@preview/cetz:0.5.2"

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
} else if selected-case == "invalid-smarts" {
  smiles("CCO", highlight-smarts: "[O")
} else if selected-case == "unsupported-smarts" {
  substructure-matches("CCO", "[C@H]")
} else if selected-case == "empty-smarts" {
  smiles("CCO", highlight-smarts: "")
} else if selected-case == "smarts-wrong-type" {
  smiles("CCO", highlight-smarts: 1)
} else if selected-case == "smarts-invalid-item" {
  smiles("CCO", highlight-smarts: ("O", 1))
} else if selected-case == "unknown-functional-group" {
  smiles("CCO", highlight-groups: "acetylated unicorn")
} else if selected-case == "unmatched-smarts" {
  smiles("CCO", highlight-smarts: "N")
} else if selected-case == "unmatched-group" {
  reaction(mol("CCO", highlight-groups: "amide"))
} else if selected-case == "invalid-highlight-colors" {
  smiles("CCO", highlight-colors: ())
} else if selected-case == "invalid-highlight-color-item" {
  smiles("CCO", highlight-colors: ("yellow",))
} else if selected-case == "invalid-highlight-policy" {
  smiles("CCO", highlight-unmatched: "silent")
} else if selected-case == "query-invalid-smiles" {
  substructure-matches("C1CC", "C")
} else if selected-case == "highlight-request-missing-pattern" {
  smiles("CCO", highlight-smarts: (include-atoms: false))
} else if selected-case == "highlight-request-missing-group" {
  smiles("CCO", highlight-groups: (include-atoms: false))
} else if selected-case == "highlight-request-invalid-bool" {
  smiles("CCO", highlight-smarts: (pattern: "CC", include-atoms: "false"))
} else if selected-case == "highlight-group-invalid-bool" {
  reaction(mol("CCO", highlight-groups: (group: "alcohol", include-atoms: 1)))
} else if selected-case == "highlight-request-unknown-option" {
  smiles("CCO", highlight-smarts: (pattern: "CC", include-atom: false))
} else if selected-case == "unknown-library-molecule" {
  smiles(molecules.cafeine)
} else if selected-case == "one-sided-directional-bond" {
  smiles("F/C=CF")
} else if selected-case == "abbreviate-wrong-type" {
  smiles("COc1ccccc1", abbreviate: true)
} else if selected-case == "abbreviate-unknown-group" {
  smiles("COc1ccccc1", abbreviate: ("OMe", "Boc"))
} else if selected-case == "abbreviate-wrong-case" {
  smiles("COc1ccccc1", abbreviate: "ome")
} else if selected-case == "abbreviate-duplicate-group" {
  smiles("COc1ccccc1", abbreviate: ("OMe", "OMe"))
} else if selected-case == "abbreviate-all-in-list" {
  smiles("COc1ccccc1", abbreviate: ("all", "OMe"))
} else if selected-case == "abbreviate-absent-group" {
  smiles("CCO", abbreviate: "CF3")
} else if selected-case == "abbreviate-expanded-group" {
  smiles("[13CH3]Oc1ccccc1", abbreviate: "OMe")
} else if selected-case == "abbreviate-lower-priority-group" {
  smiles("COC(C)=O", abbreviate: ("OMe", "Ac"))
} else if selected-case == "abbreviate-hidden-arrow-atom" {
  smiles("COc1ccccc1", abbreviate: "OMe", arrow(from: atom(0), to: atom(3)))
} else if selected-case == "abbreviate-hidden-highlight-bond" {
  smiles("COc1ccccc1", abbreviate: "OMe", highlight(bond(0, 1)))
} else if selected-case == "abbreviate-hidden-smarts-match" {
  smiles("[O-][N+](=O)c1ccccc1", abbreviate: "NO2", highlight-groups: "nitro")
} else if selected-case == "abbreviate-hidden-annotation" {
  smiles("FC(F)(F)c1ccccc1", abbreviate: "CF3", atom-annotations: ((0, [F]),))
} else if selected-case == "abbreviate-hidden-bond-customization" {
  smiles("FC(F)(F)c1ccccc1", abbreviate: "CF3", bond-customizations: ((bond(0, 1), (color: red)),))
} else if selected-case == "abbreviate-hidden-show-h" {
  smiles("CCOc1ccccc1", abbreviate: "OEt", show-h: 0)
} else if selected-case == "abbreviate-label-show-h" {
  smiles("CC(C)(C)c1ccccc1", abbreviate: "tBu", show-h: (1,))
} else if selected-case == "abbreviate-label-lone-pair" {
  smiles("COc1ccccc1", abbreviate: "OMe", arrow(from: lp(1), to: atom(2)))
} else if selected-case == "mol-abbreviate-unknown-group" {
  reaction(mol("COc1ccccc1", abbreviate: "Me"))
} else if selected-case == "mol-abbreviate-hidden-atom" {
  reaction(
    mol("COc1ccccc1", abbreviate: "OMe"),
    rxn-arrow(),
    mol("Oc1ccccc1"),
    arrow(from: atom(0, 0), to: atom(1, 0)),
  )
} else if selected-case == "cetz-abbreviate-unknown-group" {
  cetz.canvas(smiles-cetz("COc1ccccc1", abbreviate: "Et"))
} else {
  panic("unknown validation test case: " + selected-case)
}
