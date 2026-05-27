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

  example("Solid wedge", smiles("C/N", bond-length: 1.6)),
  example("Hashed wedge", smiles("C\\N", bond-length: 1.6)),
  example("Mixed", smiles("F/C\\Cl", bond-length: 1.4)),
  example("With carbon H", smiles("CCO", bond-length: 1.2, show-all-h: true)),
)
