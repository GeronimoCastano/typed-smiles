#import "../src/lib.typ": smiles, mol, mol-formula, rxn-arrow, reaction, cycle, step, atom, bond, lp, species, arrow, highlight, molecules, substructure-matches, align-molecules, molecule-grid, smiles-cetz

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
} else if selected-case == "ring-self-bond" {
  smiles("C11")
} else if selected-case == "ring-duplicate-bond" {
  smiles("C12C12")
} else if selected-case == "ring-conflicting-bonds" {
  smiles("C=1CCCCC-1")
} else if selected-case == "dangling-bond" {
  smiles("CC=")
} else if selected-case == "malformed-charge" {
  smiles("[Fe+++]")
} else if selected-case == "malformed-atom-map" {
  smiles("[CH3:]O")
} else if selected-case == "unbracketed-element" {
  smiles("NaCl")
} else if selected-case == "unclosed-branch" {
  smiles("CC(C")
} else if selected-case == "undepicted-octahedral" {
  smiles("C[Co@OH1](F)(Cl)(Br)(I)N")
} else if selected-case == "undepicted-in-reaction" {
  reaction(mol("NC(Br)=[C@AL1]=C(O)C"), rxn-arrow(), mol("CC"))
} else if selected-case == "undepicted-not-stereocenter" {
  smiles("[C@H2](F)Cl")
} else if selected-case == "undepicted-ring-trans" {
  smiles("C1CCC/C=C/CC1")
} else if selected-case == "undepicted-stereo-policy" {
  smiles("CCO", undepicted-stereo: "ignore")
} else if selected-case == "show-maps-type" {
  smiles("[CH3:1]O", show-maps: "yes")
} else if selected-case == "trans-double-bond-in-small-ring" {
  smiles("C1CCC/C=C/CC1")
} else if selected-case == "align-not-array" {
  align-molecules("CC(=O)c1ccccc1", scaffold: "c1ccccc1")
} else if selected-case == "align-single-molecule" {
  align-molecules(("CCO",), scaffold: "CC")
} else if selected-case == "align-molecule-type" {
  align-molecules(("CCO", 3), scaffold: "CC")
} else if selected-case == "align-aligned-input" {
  let aligned = align-molecules(("CCO", "CCN"), scaffold: "CC")
  align-molecules(aligned, scaffold: "CC")
} else if selected-case == "align-scaffold-type" {
  align-molecules(("CCO", "CCN"), scaffold: 5)
} else if selected-case == "align-no-correspondence" {
  align-molecules(("CCO", "CCN"))
} else if selected-case == "align-partial-atoms-without-scaffold" {
  align-molecules(("CCO", "CCN"), atoms: ((0, 1), auto))
} else if selected-case == "align-atoms-length" {
  align-molecules(("CCO", "CCN"), scaffold: "CC", atoms: (auto,))
} else if selected-case == "align-atoms-entry" {
  align-molecules(("CCO", "CCN"), scaffold: "CC", atoms: (auto, (0, "1")))
} else if selected-case == "align-reference-range" {
  align-molecules(("CCO", "CCN"), scaffold: "CC", reference: 2)
} else if selected-case == "align-rotation-type" {
  align-molecules(("CCO", "CCN"), scaffold: "CC", rotation: 30)
} else if selected-case == "align-mirror-value" {
  align-molecules(("CCO", "CCN"), scaffold: "CC", mirror: "diagonal")
} else if selected-case == "align-allow-reflection-type" {
  align-molecules(("CCO", "CCN"), scaffold: "CC", allow-reflection: "no")
} else if selected-case == "align-absent-scaffold" {
  align-molecules(("CCO", "CCN"), scaffold: "c1ccccc1")
} else if selected-case == "align-invalid-scaffold" {
  align-molecules(("CCO", "CCN"), scaffold: "[C")
} else if selected-case == "align-invalid-smiles" {
  align-molecules(("CCO", "C1CC"), scaffold: "CC")
} else if selected-case == "align-distinct-occurrences" {
  align-molecules(("c1ccccc1", "c1ccccc1-c1ccc(Cl)cc1"), scaffold: "c1ccccc1")
} else if selected-case == "align-atoms-not-scaffold-match" {
  align-molecules(("CCO", "CCO"), scaffold: "CO", atoms: (auto, (0, 1)))
} else if selected-case == "align-atom-out-of-range" {
  align-molecules(("CCO", "CCO"), atoms: ((0, 1), (0, 7)))
} else if selected-case == "align-repeated-atom" {
  align-molecules(("CCO", "CCO"), atoms: ((0, 1), (1, 1)))
} else if selected-case == "align-unequal-atom-lists" {
  align-molecules(("CCO", "CCO"), atoms: ((0, 1, 2), (0, 1)))
} else if selected-case == "align-single-atom-correspondence" {
  align-molecules(("CO", "CN"), scaffold: "C")
} else if selected-case == "aligned-rotation-conflict" {
  let aligned = align-molecules(("CCO", "OCC"), scaffold: "CCO")
  smiles(aligned.at(1), rotation: 30deg)
} else if selected-case == "aligned-mirror-conflict" {
  let aligned = align-molecules(("CCO", "OCC"), scaffold: "CCO")
  reaction(mol(aligned.at(1), mirror: "horizontal"))
} else if selected-case == "aligned-cetz-rotation-conflict" {
  let aligned = align-molecules(("CCO", "OCC"), scaffold: "CCO")
  smiles-cetz(aligned.at(1), rotation: 90deg)
} else if selected-case == "aligned-skeleton-hydrogens" {
  let aligned = align-molecules(("CCO", "OCC"), scaffold: "CCO")
  smiles(aligned.at(0), show-h: "skeleton")
} else if selected-case == "invalid-mol-spec" {
  reaction(mol(42))
} else if selected-case == "grid-empty" {
  molecule-grid()
} else if selected-case == "grid-columns" {
  molecule-grid(columns: 0, "CCO")
} else if selected-case == "grid-columns-type" {
  molecule-grid(columns: "3", "CCO")
} else if selected-case == "grid-scale" {
  molecule-grid(scale: -1, "CCO")
} else if selected-case == "grid-bond-length" {
  molecule-grid(bond-length: 0, "CCO")
} else if selected-case == "grid-sizing" {
  molecule-grid(sizing: "stretch", "CCO")
} else if selected-case == "grid-gutter" {
  molecule-grid(column-gutter: -1pt, "CCO")
} else if selected-case == "grid-breakable" {
  molecule-grid(breakable: "yes", "CCO")
} else if selected-case == "grid-unknown-option" {
  molecule-grid(gap: 1em, "CCO")
} else if selected-case == "grid-invalid-item" {
  molecule-grid("CCO", [plain content])
} else if selected-case == "grid-item-scale" {
  molecule-grid("CCO", mol("CCN", scale: 2))
} else if selected-case == "grid-item-bond-length" {
  molecule-grid("CCO", mol("CCN", bond-length: 2))
} else if selected-case == "grid-overflow" {
  set page(width: 6cm)
  molecule-grid(columns: 3, "CCCCCCCCCCCC", "CCO")
} else if selected-case == "grid-scaffold-type" {
  molecule-grid(scaffold: (), "CCO", "CCN")
} else if selected-case == "grid-scaffold-single" {
  molecule-grid(scaffold: "CC", "CCO")
} else if selected-case == "grid-scaffold-content" {
  molecule-grid(scaffold: "CC", "CCO", mol([content]))
} else if selected-case == "grid-scaffold-rotation" {
  molecule-grid(scaffold: "CC", "CCO", mol("CCN", rotation: 30deg))
} else if selected-case == "grid-scaffold-absent" {
  molecule-grid(scaffold: "c1ccccc1", "CCO", "CCN")
} else {
  panic("unknown validation test case: " + selected-case)
}
