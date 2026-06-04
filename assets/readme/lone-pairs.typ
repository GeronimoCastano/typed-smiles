#import "../../src/lib.typ": smiles

#set page(width: 16cm, height: 6.5cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Alcohol dots*],
  [*Amine lines*],
  [*Amide dots*],

  [#smiles("CCO", lone-pairs: "dots", bond-length: 1.15)],
  [#smiles("CCN", lone-pairs: "lines", bond-length: 1.15)],
  [#smiles("CC(=O)N", lone-pairs: "dots", bond-length: 1.05)],
)
