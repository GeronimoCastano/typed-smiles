#import "../src/lib.typ": smiles, smiles-inline, smiles-cetz
#import "../src/chemistry.typ": _compute-layout
#import "../src/molecule/rendering.typ": _mirror-layout, _linearize-skeleton-layout, _skeleton-hydrogen-directions, _rendered-atom-position
#import "@preview/cetz:0.5.2"

#set text(font: "New Computer Modern", size: 11pt)

#let transform(point, rotation, mirror) = {
  let rotated = (
    x: point.x * calc.cos(rotation) - point.y * calc.sin(rotation),
    y: point.x * calc.sin(rotation) + point.y * calc.cos(rotation),
  )
  (
    x: rotated.x * if mirror == "horizontal" { -1 } else { 1 },
    y: rotated.y * if mirror == "vertical" { -1 } else { 1 },
  )
}

// Check actual H positions (parent + bond direction), including atoms without
// heavy neighbors and non-cardinal electron-domain geometry. Every H must be
// the rigidly transformed version of the unchanged 0-degree skeleton.
#for molecule in ("OCCO", "CCO", "C", "N", "O", "C[NH2]", "CC(C)O", "c1ccccc1") {
  let raw = _linearize-skeleton-layout(_compute-layout(molecule))
  let base = _mirror-layout(raw, none)
  for rotation in (0deg, 15deg, 45deg, 90deg, 135deg, 180deg) {
    for mirror in (none, "horizontal", "vertical") {
      let layout = _mirror-layout(raw, mirror, rotation: rotation)
      for (index, atom) in base.atoms.enumerate() {
        if atom.at("virtual_h", default: false) { continue }
        let count = atom.hcount + atom.at("implicit_h", default: 0)
        let original = _skeleton-hydrogen-directions(base, index, count)
        let directions = _skeleton-hydrogen-directions(layout, index, count, rotation: rotation)
        assert.eq(directions.len(), original.len())
        let parent = _rendered-atom-position(layout.atoms.at(index), rotation)
        for (before, after) in original.zip(directions) {
          let expected = transform((x: atom.pos.x + before.x * 0.88, y: atom.pos.y + before.y * 0.88), rotation, mirror)
          let actual = (x: parent.x + after.x * 0.88, y: parent.y + after.y * 0.88)
          assert(calc.abs(actual.x - expected.x) < 0.000001)
          assert(calc.abs(actual.y - expected.y) < 0.000001)
        }
      }
    }
  }
}

// Preserve the established unrotated CH2 and water geometry.
#let diol = _mirror-layout(_linearize-skeleton-layout(_compute-layout("OCCO")), none)
#assert.eq(_skeleton-hydrogen-directions(diol, 1, 2), ((x: 0.0, y: 1.0), (x: 0.0, y: -1.0)))
#let water = _mirror-layout(_compute-layout("O"), none)
#let water-directions = _skeleton-hydrogen-directions(water, 0, 2)
#assert(calc.abs(water-directions.first().x * water-directions.last().x + water-directions.first().y * water-directions.last().y - calc.cos(104.5deg)) < 0.000001)

= Skeleton hydrogen rotation

#grid(columns: (auto, auto, auto), gutter: 2em,
  [*0°* \ #smiles("OCCO", show-h: "skeleton")],
  [*45°* \ #smiles("OCCO", show-h: "skeleton", rotation: 45deg)],
  [*45° horizontal mirror* \ #smiles("OCCO", show-h: "skeleton", rotation: 45deg, mirror: "horizontal")],
  [*15° vertical mirror* \ #smiles("CCO", show-h: "skeleton", rotation: 15deg, mirror: "vertical")],
  [*Methane, 45°* \ #smiles("C", show-h: "skeleton", rotation: 45deg)],
  [*Water, 45°* \ #smiles("O", show-h: "skeleton", rotation: 45deg)],
)

Inline #smiles-inline("CCO", show-h: "skeleton", rotation: 30deg) skeleton.

#context cetz.canvas(length: 30pt, {
  smiles-cetz("OCCO", show-h: "skeleton", rotation: 45deg)
})
