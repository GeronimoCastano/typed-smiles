#import "../../src/lib.typ": smiles

#set page(width: 16cm, height: 6cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#let example(title, code, body) = [
  #stack(
    spacing: 4pt,
    align(center, strong(title)),
    align(center, text(size: 7pt, raw(code))),
    align(center, body),
  )
]

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.2em,
  align: center,

  example("Small", "bond-length: 0.8", smiles("C1=CC=CC=C1", bond-length: 0.8)),
  example("Default", "bond-length: 1.0", smiles("C1=CC=CC=C1")),
  example("Large", "bond-length: 1.4", smiles("C1=CC=CC=C1", bond-length: 1.4)),
)
