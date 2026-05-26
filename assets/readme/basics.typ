#import "../../src/lib.typ": smiles

#set page(width: 17cm, height: 5.7cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#let example(title, body) = [
  #stack(
    spacing: 0.5cm,
    align(center, strong(title)),
    align(center, body),
  )
]

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.2em,
  align: center,

  example("Ethanol", smiles("CCO", bond-length: 1.2)),
  example("Alanine", smiles("CC(N)C(=O)O", bond-length: 1.05)),
  example("Chlorobenzene", smiles("ClC1=CC=CC=C1", bond-length: 1.2)),
  example("Furan", smiles("C1=CC=CO1", bond-length: 1.2)),
)
