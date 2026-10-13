#import "../src/lib.typ": smiles
#import "../src/chemistry.typ": _compute-layout

#let charged-atom(molecule) = _compute-layout(molecule).atoms.find(atom => atom.charge == 1)
#assert.eq(charged-atom("C[NH2+]C").hcount, 2)
#assert.eq(charged-atom("C[NH2+]C").lone_pairs, 0)
#assert.eq(charged-atom("C[ClH2+]C").hcount, 2)
#assert.eq(charged-atom("C[ClH2+]C").lone_pairs, 1)
#assert.eq(charged-atom("C[SeH2+]").hcount, 2)
#assert.eq(charged-atom("C[SeH2+]").lone_pairs, 1)

#let rotation-sweep(molecule) = grid(
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
        molecule,
        rotation: rotation-degrees * 1deg,
        lone-pairs: "dots",
        scale: 0.9,
      )
    ]
  }),
)

#pagebreak()
= Ammonium H₂ charge placement: 15° rotation sweep

`C[NH2+]C`: nitrogen has four single bonds and no lone pairs. Lone-pair
display is enabled so the absence of dots reflects the electron count.

#rotation-sweep("C[NH2+]C")

#pagebreak()
= Chlorine H₂ charge placement: hypervalent stress case

`C[ClH2+]C`: formal electron counting gives one lone pair on chlorine.
This unusual four-bond, ten-electron center tests label geometry.

#rotation-sweep("C[ClH2+]C")

#pagebreak()
= Selenium H₂ charge placement: 15° rotation sweep

`C[SeH2+]`: selenium has three single bonds and one lone pair. This is an
octet-consistent comparison with a two-letter element and an H₂ label.

#rotation-sweep("C[SeH2+]")
