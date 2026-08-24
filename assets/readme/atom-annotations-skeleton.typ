#import "../../src/lib.typ": smiles

#set page(width: 23cm, height: 6.5cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1.15fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Atom annotations*],
  [*Selected C-H*],
  [*All implicit H*],
  [*Ethanol skeleton*],

  [#smiles(
    "N[C@@H](C)C(=O)O",
    bond-length: 1.05,
    atom-annotations: (
      (1, [$alpha$], (-0.4, -0.05)),
      (2, [$beta$]),
      (3, [$gamma$], (-0.05, -0.3)),
    ),
  )],
  [#smiles("CC(N)C(=O)O", bond-length: 1.05, show-h: 1)],
  [#smiles("CC(N)C(=O)O", bond-length: 1.05, show-h: "all")],
  [#smiles("CCO", bond-length: 1.05, show-h: "skeleton")],
)
