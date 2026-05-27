#import "../../src/lib.typ": smiles

#set page(width: 16cm, height: 6.5cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.2em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Small*],
  [*Default*],
  [*Large*],

  [#smiles("C1=CC=CC=C1", bond-length: 0.8)],
  [#smiles("C1=CC=CC=C1")],
  [#smiles("C1=CC=CC=C1", bond-length: 1.4)],
)
