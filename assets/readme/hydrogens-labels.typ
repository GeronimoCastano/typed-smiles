#import "../../src/lib.typ": smiles

#set page(width: 17cm, height: 6.3cm, margin: 0.8cm)
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

  example("Default hetero H", smiles("CC(N)C(=O)O", bond-length: 1.0)),
  example("All H", smiles("CCO", bond-length: 1.15, show-all-h: true)),
  example("Explicit bracket H", smiles("[NH3]", bond-length: 1.25)),
  example("Custom label", smiles("{PPh3}C=O", bond-length: 1.1)),
)
