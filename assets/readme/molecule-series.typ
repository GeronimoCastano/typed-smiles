#import "../../src/lib.typ": smiles, rxn-arrow, mol, reaction, align-molecules, molecule-grid

#set page(width: 16cm, height: auto, margin: 0.6cm)
#set text(font: "New Computer Modern", size: 10pt)

#let aligned = align-molecules(
  ("CC(=O)c1ccccc1", "CC(O)c1ccccc1"),
  scaffold: "c1ccccc1",
)

#align(center, reaction(
  mol(aligned.at(0)),
  rxn-arrow(above: [NaBH#sub[4]]),
  mol(aligned.at(1)),
))

#v(0.6cm)

#molecule-grid(
  columns: 4,
  scale: 0.5,
  column-gutter: 1em,
  scaffold: "c1ccccc1C(=O)O",
  mol("c1ccccc1C(=O)O", label: [*1* Benzoic acid]),
  mol("O=C(O)c1ccccc1O", label: [*2* Salicylic acid]),
  mol("Nc1ccc(C(=O)O)cc1", label: [*3* PABA]),
  mol("OC(=O)c1ccc(cc1)[N+](=O)[O-]", label: [*4* 4-Nitrobenzoic acid]),
)
