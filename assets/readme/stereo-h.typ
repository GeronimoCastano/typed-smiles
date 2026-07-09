#import "../../src/lib.typ": smiles

#set page(width: 18cm, height: 5.8cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Manual wedge*],
  [*Manual hash*],
  [*Wavy*],
  [*Dashed*],
  [*Tetrahedral @@*],
  [*trans alkene*],

  [#smiles("C!wN")],
  [#smiles("C!hN")],
  [#smiles("C!sN")],
  [#smiles("C!dN")],
  [#smiles("N[C@@H](C)C(=O)O", scale: 0.7)],
  [#smiles("F/C=C/F", scale: 0.7)],
)
