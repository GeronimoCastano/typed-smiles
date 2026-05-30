#import "../../src/lib.typ": smiles

#set page(width: 17cm, height: 8cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Ethanol*],
  [*Alanine*],
  [*Chlorobenzene*],
  [*Furan*],

  [#smiles("CCO", bond-length: 1.2)],
  [#smiles("CC(N)C(=O)O", bond-length: 1.05)],
  [#smiles("ClC1=CC=CC=C1", bond-length: 1.2)],
  [#smiles("C1=CC=CO1", bond-length: 1.2)],
)
