#import "../src/lib.typ": smiles, atom, bond, highlight, reaction, mol, rxn-arrow, arrow
#import "../src/chemistry.typ": _compute-layout
#import "../src/molecule/rendering.typ": _mirror-layout, _linearize-skeleton-layout, _rendered-atom-position
#import "../src/mechanism/references.typ": _atom-position

#set text(font: "New Computer Modern", size: 11pt)

#context {
  for molecule in ("OCCO", "CC[OH]", "[OH-]", "O", "C[NH2]") {
    for rotation in (0deg, 15deg, 45deg, 90deg, 180deg) {
      for mirror in (none, "horizontal", "vertical") {
        let layout = _mirror-layout(_linearize-skeleton-layout(_compute-layout(molecule)), mirror, rotation: rotation)
        for scale in (0.7, 1.0, 1.5) {
          let species = (
            layout: layout, rotation: rotation, origin: (1.5, -0.5),
            mol-scale: scale, canvas-scale: 30pt, actual-font-size: 11pt * scale,
            show-h: "skeleton",
          )
          for (index, atom) in layout.atoms.enumerate() {
            if atom.at("virtual_h", default: false) { continue }
            let expected = _rendered-atom-position(atom, rotation, scale: scale)
            let actual = _atom-position(species, index)
            assert(calc.abs(actual.at(0) - expected.x - 1.5) < 0.000001)
            assert(calc.abs(actual.at(1) - expected.y + 0.5) < 0.000001)
          }
        }
      }
    }
  }
  // Inline OH still offsets the O glyph from the structural bond terminus.
  let layout = _mirror-layout(_compute-layout("OCCO"), none)
  let position = _atom-position((layout: layout, rotation: 0deg, origin: (0, 0), show-h: ()), 0)
  let base = _rendered-atom-position(layout.atoms.first(), 0deg)
  assert(calc.abs(position.at(0) - base.x) + calc.abs(position.at(1) - base.y) > 0.001)
}

= Skeleton atom highlight centers

#grid(columns: (auto, auto), gutter: 2em,
  [*0°* \ #smiles("OCCO", show-h: "skeleton", highlight((atom(0), atom(3))))],
  [*45° mirrored* \ #smiles("OCCO", show-h: "skeleton", rotation: 45deg, mirror: "horizontal", highlight((atom(0), atom(3))))],
  [*Bond endpoint shading* \ #smiles("CC[OH]", show-h: "skeleton", rotation: 30deg, highlight(bond(1, 2), include-atoms: true))],
  [*Ordinary inline OH* \ #smiles("CC[OH]", highlight(atom(2)))],
)

#reaction(
  mol("OCCO", show-h: "skeleton", rotation: 30deg, highlight((atom(0), atom(3)))),
  rxn-arrow(),
  mol("CC[OH]", show-h: "skeleton", highlight(atom(2))),
  arrow(from: atom(0, 3), to: atom(1, 2)),
)
