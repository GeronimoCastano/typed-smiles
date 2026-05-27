#import "../../src/lib.typ": smiles

#set page(width: 17cm, height: 5.7cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Solid wedge*],
  [*Hashed wedge*],
  [*Mixed*],
  [*With carbon H*],

  [#smiles("C/N", bond-length: 1.6)],
  [#smiles("C\\N", bond-length: 1.6)],
  [#smiles("F/C\\Cl", bond-length: 1.4)],
  [#smiles("CCO", bond-length: 1.2, show-all-h: true)],
)
