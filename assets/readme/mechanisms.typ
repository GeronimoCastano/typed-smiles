#import "../../src/lib.typ": smiles, rxn-arrow, mol, reaction, atom, bond, lp, arrow, highlight, brackets

#set page(width: 18cm, height: 6.6cm, margin: 0.7cm)
#set text(font: "New Computer Modern", size: 10pt)

#grid(
  columns: (1fr, 1fr, 1fr),
  column-gutter: 0.7cm,
  align: center + horizon,

  stack(
    spacing: 0.4cm,
    strong[Carbonyl highlight],
    smiles(
      "CC(=O)C",
      lone-pairs: "dots",
      highlight(bond(1, 2), fill: rgb("#FFE45C"), include-atoms: true),
    ),
  ),

  stack(
    spacing: 0.4cm,
    strong[Hydroxide attack],
    reaction(
      mol("[OH-]", lone-pairs: "dots", offset: (1.5, 1)),
      mol("C(I)(C)C"),
      arrow(from: lp(0, 0), to: atom(1, 0), bend: "left"),
    ),
  ),

  stack(
    spacing: 0.4cm,
    strong[Bracketed reaction],
    brackets(
      [#reaction(smiles("CC(=O)C"), rxn-arrow(), smiles("O=C=O"), scale: 0.55)],
      sup: [‡],
    ),
  ),
)
