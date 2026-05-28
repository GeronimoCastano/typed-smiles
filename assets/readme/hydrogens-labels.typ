#import "../../src/lib.typ": smiles

#set page(width: 21cm, height: 6.3cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Default hetero H*],
  [*All H*],
  [*Explicit bracket H*],
  [*Custom label*],
  [*Custom font*],

  [#smiles("CC(N)C(=O)O", bond-length: 1.0)],
  [#smiles("CCO", bond-length: 1.15, show-all-h: true)],
  [#smiles("[NH3]", bond-length: 1.25)],
  [#smiles("{PPh3}C=O", bond-length: 1.1)],
  [#smiles("CCN", atom-font: "Libertinus Serif", bond-length: 1.1)],
)
