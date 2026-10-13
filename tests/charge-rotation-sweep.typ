#import "../src/lib.typ": smiles

#pagebreak()
= Oxonium charge placement: 15° rotation sweep

Each panel renders `C[OH+]C` with one displayed lone pair. Only the molecule
rotation changes; the labels remain upright.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  rows: (2.7cm,) * 6,
  column-gutter: 0.6em,
  row-gutter: 0.2em,
  align: center + horizon,
  ..range(24).map(index => {
    let rotation-degrees = index * 15
    [
      #text(size: 9pt, weight: "bold", str(rotation-degrees) + "°") \
      #v(0.3em)
      #smiles(
        "C[OH+]C",
        rotation: rotation-degrees * 1deg,
        lone-pairs: "dots",
        scale: 0.9,
      )
    ]
  }),
)
