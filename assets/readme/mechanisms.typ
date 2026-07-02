#import "../../src/lib.typ": smiles, rxn-arrow, mol, reaction, atom, bond, lp, arrow, highlight, brackets

#set page(width: 18cm, height: 6.6cm, margin: 0.7cm)
#set text(font: "New Computer Modern", size: 10pt)

#grid(
  columns: (1fr, 1fr, 1fr),
  column-gutter: 0.7cm,
  align: center + horizon,

  stack(
    spacing: 0.4cm,
    strong[Bond path highlight],
    smiles(
      "N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4C(=O)O",
      highlight((bond(0, 5), bond(5, 4), bond(4, 3), bond(3, 6), bond(6, 10), bond(10, 11), bond(11, 15), bond(15, 19), bond(19, 20), bond(20, 21), bond(21, 23), bond(23, 25)), fill: rgb("#96BF0D"), include-atoms: true),
      scale: 0.5,
      bond-stroke: 0.8pt,
      rotation: 90deg,
    ),
  ),

  stack(
    spacing: 0.4cm,
    strong[Hydroxide attack],
    reaction(
      mol("[OH-]", lone-pairs: "dots", offset: (1.5, 1)),
      mol("C(I)(C)C"),
      arrow(
        from: lp(0, 0, offset: (-0.3, -0.2)),
        to: atom(1, 0, offset: (0.1, -0.1)),
        bend: "right",
        color: black,
      ),
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
