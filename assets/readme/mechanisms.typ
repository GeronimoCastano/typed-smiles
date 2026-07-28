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
      highlight((bond(0, 5), bond(5, 4), bond(4, 3), bond(3, 6), bond(6, 10), bond(10, 11), bond(11, 15), bond(15, 19), bond(19, 20), bond(20, 21), bond(21, 23)), fill: rgb(150, 191, 13), include-atoms: true),
      highlight((bond(15, 16), bond(16, 18), bond(18, 17), bond(17, 16)), fill: rgb(242, 148, 1), include-atoms: true),
      highlight((bond(3, 2), bond(2, 1), bond(1, 0)), fill: rgb(137, 199, 168), include-atoms: true),
      highlight((bond(6, 7), bond(7, 8), bond(7, 9), bond(9, 12)), fill: rgb(201, 143, 75), include-atoms: true),
      highlight((bond(11, 12), bond(12, 13), bond(13, 20), bond(13, 14)), fill: rgb(236, 119, 137), include-atoms: true),
      highlight((bond(21, 22)), fill: rgb(0, 134, 203), include-atoms: true),
      color: false,
      scale: 0.5,
      bond-stroke: 0.8pt,
      rotation: 90deg,
    ),
  ),

  stack(
    spacing: 0.4cm,
    strong[Hydroxide attack],
    reaction(
      mol("[OH-]", lone-pairs: "dots", offset: (0, 1)),
      mol("C(I)(C)C"),
      arrow(
        from: lp(0, 0, offset: (0.1, -0.2)),
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
