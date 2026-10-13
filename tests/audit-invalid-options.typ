#import "../src/lib.typ": smiles, mol-formula
#let selected-case = sys.inputs.at("case")
#if selected-case == "audit-scale-zero" {
  smiles("CO", scale: 0)
} else if selected-case == "audit-bond-length-zero" {
  smiles("CO", bond-length: 0)
} else if selected-case == "audit-font-size-zero" {
  smiles("CO", font-size: 0pt)
} else if selected-case == "audit-bond-stroke-zero" {
  smiles("CO", bond-stroke: 0pt)
} else if selected-case == "audit-rotation-type" {
  smiles("CO", rotation: "90")
} else if selected-case == "audit-lone-pairs-value" {
  smiles("CO", lone-pairs: "rings")
} else if selected-case == "audit-atom-colors-type" {
  smiles("CO", atom-colors: (red,))
} else if selected-case == "audit-index-type" {
  smiles("CO", show-indices: "yes")
} else if selected-case == "audit-font-type" {
  smiles("CO", font: 7)
} else if selected-case == "audit-foreground-type" {
  smiles("CO", fg: "red")
} else if selected-case == "audit-aromatic-value" {
  smiles("c1ccccc1", aromatic: "dots")
} else if selected-case == "audit-opacity-type" {
  smiles("CO", opacity: "50%")
} else if selected-case == "audit-folded-hydrogen-overflow" {
  mol-formula("C" + range(300).map(_ => "([H])").join(""))
} else { panic("unknown audit validation case: " + selected-case) }
