#import "../../src/lib.typ": smiles

#let molecule = "NC(Br)C(I)C(=O)O"

#set page(width: 18cm, height: 7.2cm, margin: 0.7cm, fill: white)
#set text(font: "New Computer Modern", size: 10pt, fill: black)

#grid(
  columns: (1fr, 1fr),
  gutter: 0.8cm,
  align: center + horizon,

  block(
    width: 100%,
    height: 5.8cm,
    fill: white,
    inset: 0.5cm,
    radius: 4pt,
    stroke: 0.4pt + rgb("#d8d8d8"),
  )[
    #align(center)[
      #text(size: 17pt, weight: "bold")[Light background]
      #v(0.6cm)
      #smiles(molecule, scale: 0.95)
    ]
  ],
  block(
    width: 100%,
    height: 5.8cm,
    fill: rgb("#1E1E24"),
    inset: 0.5cm,
    radius: 4pt,
  )[
    #set text(fill: white)
    #align(center)[
      #text(size: 17pt, weight: "bold")[Dark background]
      #v(0.6cm)
      #smiles(molecule, scale: 0.95)
    ]
  ],
)
