#import "../../src/lib.typ": smiles

#set page(width: 18cm, height: 6.6cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Original*],
  [*`mirror: "horizontal"`*],
  [*`mirror: "vertical"`*],

  [#smiles("CC(=O)OC1=CC=CC=C1C(=O)O", bond-length: 0.72)],
  [#smiles("CC(=O)OC1=CC=CC=C1C(=O)O", bond-length: 0.72, mirror: "horizontal")],
  [#smiles("CC(=O)OC1=CC=CC=C1C(=O)O", bond-length: 0.72, mirror: "vertical")],
)
