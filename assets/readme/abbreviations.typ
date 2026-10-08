#import "../../src/lib.typ": smiles

#set page(width: 21cm, height: 8.2cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Full structure*],
  [*`abbreviate: ("OMe", "CF3")`*],
  [*`abbreviate: "all"`*],
  [*Rotated 90°*],

  [#smiles("COc1ccc(C(F)(F)F)cc1", bond-length: 0.9)],
  [#smiles("COc1ccc(C(F)(F)F)cc1", abbreviate: ("OMe", "CF3"), bond-length: 0.9)],
  [#smiles("CC(=O)Nc1ccc(cc1)[N+](=O)[O-]", abbreviate: "all", bond-length: 0.9)],
  [#smiles("COc1ccc(C(F)(F)F)cc1", abbreviate: "all", rotation: 90deg, bond-length: 0.75)],
)
