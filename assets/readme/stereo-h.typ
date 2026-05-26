#import "../../src/lib.typ": smiles

#set page(width: 17cm, height: 4.4cm, margin: 0.8cm)
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
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.2em,
  align: center,

  example("Solid wedge", "C/N", smiles("C/N", bond-length: 1.6)),
  example("Hashed wedge", "C\\N", smiles("C\\N", bond-length: 1.6)),
  example("Implicit H", "CCO", smiles("CCO", bond-length: 1.2, show-h: true)),
  example("Explicit H", "[NH4+]", smiles("[NH4+]", bond-length: 1.4)),
)
