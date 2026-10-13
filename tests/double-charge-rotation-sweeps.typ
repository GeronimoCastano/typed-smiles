#import "../src/lib.typ": smiles
#import "../src/chemistry.typ": _compute-layout

#let charged-atom(molecule) = _compute-layout(molecule).atoms.find(atom => atom.charge == 2)
#assert.eq(charged-atom("C[NH+2]C").hcount, 1)
#assert.eq(charged-atom("C[NH+2]C").lone_pairs, 0)
#assert.eq(charged-atom("C[ClH+2]C").hcount, 1)
#assert.eq(charged-atom("C[ClH+2]C").lone_pairs, 1)

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
= Nitrogen 2+ charge: 15° rotation sweep

`C[NH+2]C`: the label has one H and a two-character charge, 2+, with no
lone pairs. This electron-deficient structure is a charge-width stress case.

#rotation-sweep("C[NH+2]C")

#pagebreak()
= Chlorine 2+ charge: 15° rotation sweep

`C[ClH+2]C`: the label has one H, charge 2+, and one lone pair. This unusual
cation tests a two-letter element together with a two-character charge.

#rotation-sweep("C[ClH+2]C")
