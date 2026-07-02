#import "../../src/lib.typ": smiles

#set page(width: 18cm, height: 6.2cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#let themed-smiles = smiles.with(
  bond-length: 0.9,
  atom-colors: (O: rgb("#8B4513"), N: rgb("#008080")),
)

#table(
  columns: (1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Per-call overrides*],
  [*Project defaults*],
  [*Inline label colors*],

  [#smiles(
    "{>PPh3}C({OEt})=O",
    atom-colors: (O: rgb("#8B4513"), "{PPh3}": rgb("#7B2D8B")),
    bond-length: 1.0,
  )],
  [#themed-smiles("CC(=O)N", lone-pairs: "dots")],
  [#smiles("{Cat|teal}C(=O){Nuc|#E040FB}", bond-length: 1.0)],
)
