#import "../../src/lib.typ": smiles

#set page(width: 17cm, height: 5.2cm, margin: 0.8cm)
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

  example("Ethanol", "CCO", smiles("CCO", bond-length: 1.2)),
  example("Benzene", "C1=CC=CC=C1", smiles("C1=CC=CC=C1", bond-length: 1.2)),
  example("Acetic acid", "CC(=O)O", smiles("CC(=O)O", bond-length: 1.2)),
  example("Chlorobenzene", "ClC1=CC=CC=C1", smiles("ClC1=CC=CC=C1", bond-length: 1.2)),
)
