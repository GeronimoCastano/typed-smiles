#import "../../src/lib.typ": smiles

#set page(width: 16cm, height: 6.5cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#let example(title, body) = [
  #stack(
    spacing: 0.5cm,
    align(center, strong(title)),
    align(center, body),
  )
]

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.2em,
  align: center,

  example("Small", smiles("C1=CC=CC=C1", bond-length: 0.8)),
  example("Default", smiles("C1=CC=CC=C1")),
  example("Large", smiles("C1=CC=CC=C1", bond-length: 1.4)),
)
