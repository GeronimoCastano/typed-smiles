#import "../src/lib.typ": smiles
#import "../src/molecule/rendering.typ": _charge-label-offset, _stacked-hydrogen-charge-anchor

#let above-anchor = _stacked-hydrogen-charge-anchor(
  "above",
  (width: 8.0, height: 12.0),
  (width: 14.0, height: 10.0),
  (width: 7.0, height: 10.0),
  2,
  12.0,
)
#let below-anchor = _stacked-hydrogen-charge-anchor(
  "below",
  (width: 8.0, height: 12.0),
  (width: 14.0, height: 10.0),
  (width: 7.0, height: 10.0),
  2,
  12.0,
)
#let crowded-directions = range(16).map(index => (
  x: calc.cos(index * 22.5deg),
  y: calc.sin(index * 22.5deg),
))
#let charge-width = 4.0
#let charge-gap = 1.0
#let above-charge = _charge-label-offset(
  crowded-directions,
  (width: 8.0, height: 12.0),
  (width: charge-width, height: 4.0),
  charge-gap,
  stacked-hydrogen: above-anchor,
)
#let below-charge = _charge-label-offset(
  crowded-directions,
  (width: 8.0, height: 12.0),
  (width: charge-width, height: 4.0),
  charge-gap,
  stacked-hydrogen: below-anchor,
)
#let charge-left-edge = above-charge.x - charge-width / 2
#let h-right-edge = above-anchor.right-edge
#let above-h-bottom = above-anchor.hydrogen-center-y - 10.0 / 2
#let above-symbol-top = 12.0 / 2
#let below-h-top = below-anchor.hydrogen-center-y + 10.0 / 2
#let below-symbol-bottom = -(12.0 / 2)
#assert(calc.abs((charge-left-edge - h-right-edge) - charge-gap) < 0.0001)
#assert(above-symbol-top < above-h-bottom)
#assert(below-h-top < below-symbol-bottom)
#assert(above-charge.x == below-charge.x)

#let heavy-atom-placement = _charge-label-offset(
  ((x: 1.0, y: 0.0), (x: -1.0, y: 0.0), (x: 0.0, y: 1.0)),
  (width: 8.0, height: 12.0),
  (width: charge-width, height: 4.0),
  charge-gap,
  stacked-hydrogen: above-anchor,
)
#assert(calc.abs(heavy-atom-placement.x) < 0.0001)
#assert(heavy-atom-placement.y < 0.0)

#let lone-pairs-around-label = (
  (x: 0.707, y: 0.707),
  (x: -0.707, y: 0.707),
  (x: 0.707, y: -0.707),
  (x: -0.707, y: -0.707),
)
#let pair-crowded-placement = _charge-label-offset(
  ((x: 1.0, y: 0.0), (x: -1.0, y: 0.0), (x: 0.0, y: 1.0)),
  (width: 8.0, height: 12.0),
  (width: charge-width, height: 4.0),
  charge-gap,
  stacked-hydrogen: above-anchor,
  lone-pair-directions: lone-pairs-around-label,
)
#assert(calc.abs(pair-crowded-placement.x - above-charge.x) < 0.0001)
#assert(calc.abs(pair-crowded-placement.y - above-charge.y) < 0.0001)

#let bare-oxonium-placement = _charge-label-offset(
  ((x: 0.866, y: -0.5), (x: -0.866, y: -0.5)),
  (width: 8.56, height: 7.51),
  (width: 5.13, height: 4.4),
  charge-gap,
  uncapped-clearance: true,
  lone-pair-directions: ((x: 0.0, y: 1.0),),
)
#assert(calc.abs(bare-oxonium-placement.x) < 0.0001)
#assert(bare-oxonium-placement.y < 0.0)

#let left-anchor = _stacked-hydrogen-charge-anchor(
  "left", (width: 8.0, height: 12.0), (width: 14.0, height: 10.0),
  (width: 7.0, height: 10.0), 2, 12.0,
)
#let right-anchor = _stacked-hydrogen-charge-anchor(
  "right", (width: 8.0, height: 12.0), (width: 14.0, height: 10.0),
  (width: 7.0, height: 10.0), 2, 12.0,
)
#let left-charge = _charge-label-offset(
  crowded-directions, (width: 8.0, height: 12.0),
  (width: charge-width, height: 4.0), charge-gap, stacked-hydrogen: left-anchor,
)
#let right-charge = _charge-label-offset(
  crowded-directions, (width: 8.0, height: 12.0),
  (width: charge-width, height: 4.0), charge-gap, stacked-hydrogen: right-anchor,
)
#assert(calc.abs((left-charge.x + charge-width / 2) - (-4.0 - 14.0 - charge-gap)) < 0.0001)
#assert(calc.abs((right-charge.x - charge-width / 2) - (4.0 + 14.0 + charge-gap)) < 0.0001)

#pagebreak()
= Charged stacked hydrogen placement

#grid(columns: (1fr, 1fr, 1fr), gutter: 1.2em, row-gutter: 1.4em, align: center + horizon,
  [*Oxonium, dots* \ #smiles("C[OH+]C", lone-pairs: "dots", scale: 1.6)],
  [*Aminium, dots* \ #smiles("C[NH2+]C", lone-pairs: "dots", scale: 1.6)],
  [*Sulfonium, dots* \ #smiles("C[SH+]C", lone-pairs: "dots", scale: 1.6)],
  [*Branched aminium* \ #smiles("CC[NH2+]CC", lone-pairs: "dots", scale: 1.4)],
  [*Ring oxonium* \ #smiles("C1C[OH+]CC1", lone-pairs: "dots", scale: 1.4)],
  [*Stacked count* \ #smiles("C[NH2+]C", show-h: "all", lone-pairs: "dots", scale: 1.6)],
  [*No lone pairs* \ #smiles("C[OH+]C", lone-pairs: none, scale: 1.6)],
  [*Lone-pair lines* \ #smiles("C[OH+]C", lone-pairs: "lines", scale: 1.6)],
  [*Hydronium* \ #smiles("[OH3+]", lone-pairs: "dots", scale: 1.6)],
)

#pagebreak()
#grid(columns: (1fr, 1fr, 1fr), gutter: 1.2em, row-gutter: 1.4em, align: center + horizon,
  [*Rotation 0°* \ #smiles("C[OH+]C", rotation: 0deg, lone-pairs: "dots", scale: 1.7)],
  [*Rotation 30°* \ #smiles("C[OH+]C", rotation: 30deg, lone-pairs: "dots", scale: 1.7)],
  [*Rotation 90°* \ #smiles("C[OH+]C", rotation: 90deg, lone-pairs: "dots", scale: 1.7)],
  [*Rotation 180°* \ #smiles("C[OH+]C", rotation: 180deg, lone-pairs: "dots", scale: 1.7)],
  [*Rotation 270°* \ #smiles("C[OH+]C", rotation: 270deg, lone-pairs: "dots", scale: 1.7)],
  [*Horizontal mirror* \ #smiles("C[OH+]C", mirror: "horizontal", lone-pairs: "dots", scale: 1.7)],
  [*Vertical mirror* \ #smiles("C[OH+]C", mirror: "vertical", lone-pairs: "dots", scale: 1.7)],
  [*Tight spacing* \ #smiles("C[OH+]C", bond-length: 0.7, font-size: 16pt, lone-pairs: "dots")],
  [*Larger labels* \ #smiles("C[OH+]C", font-size: 20pt, lone-pairs: "dots")],
)
