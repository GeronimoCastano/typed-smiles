#import "../../src/lib.typ": smiles

#set page(width: 18cm, height: 5.8cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Manual wedge*],
  [*Manual hash*],
  [*Tetrahedral @@*],
  [*trans alkene*],

  [#smiles("C!wN", bond-length: 1.6)],
  [#smiles("C!hN", bond-length: 1.6)],
  [#smiles("N[C@@H](C)C(=O)O", bond-length: 0.9)],
  [#smiles("F/C=C/F", bond-length: 1.25)],
)
