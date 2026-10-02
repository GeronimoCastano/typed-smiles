#import "../../src/lib.typ": smiles, molecules

#set page(width: 17cm, height: auto, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),
  inset: 8pt,

  [*Caffeine*], [*Alanine*], [*Aspirin*], [*Serotonin*],

  [#smiles(molecules.caffeine, scale: 0.7)],
  [#smiles(molecules.alanine, scale: 0.7)],
  [#smiles(molecules.aspirin, scale: 0.7)],
  [#smiles(molecules.serotonin, scale: 0.7)],
)
