#import "../src/chemistry.typ": _compute-layout
#import "../src/molecule/rendering.typ": _mirror-layout, _rendered-atom-position
#import "../src/mechanism/references.typ": _atom-position

#set text(font: "New Computer Modern", size: 11pt)

#context {
  let check-center(molecule, index, mirror, rotation, terminal: false) = {
    let layout = _mirror-layout(_compute-layout(molecule), mirror)
    let species = (
      layout: layout,
      rotation: rotation,
      origin: (1.25, -0.75),
      mol-scale: 1.0,
      canvas-scale: 30pt,
      actual-font-size: 11pt,
      font: "New Computer Modern",
      show-h: (),
    )
    let expected = _rendered-atom-position(layout.atoms.at(index), rotation)
    if terminal {
      let bond = layout.bonds.find(bond => bond.from == index or bond.to == index)
      let neighbor-index = if bond.from == index { bond.to } else { bond.from }
      let neighbor = _rendered-atom-position(layout.atoms.at(neighbor-index), rotation)
      let direction-x = neighbor.x - expected.x
      let direction-y = neighbor.y - expected.y
      let element = measure(text(
        size: 11pt,
        font: "New Computer Modern",
        style: "normal",
        weight: "regular",
        layout.atoms.at(index).symbol,
      ))
      let pad = 1pt / 30pt
      expected = if calc.abs(direction-x) >= calc.abs(direction-y) {
        (
          x: expected.x + (if direction-x > 0 { -1 } else { 1 })
            * (pad + element.width / 30pt / 2),
          y: expected.y,
        )
      } else {
        (
          x: expected.x,
          y: expected.y + (if direction-y > 0 { -1 } else { 1 })
            * element.height / 30pt / 2,
        )
      }
    }
    let actual = _atom-position(species, index)
    assert(calc.abs(actual.at(0) - expected.x - 1.25) < 0.000001)
    assert(calc.abs(actual.at(1) - expected.y + 0.75) < 0.000001)
  }

  for mirror in (none, "horizontal", "vertical") {
    for rotation in (0deg, 45deg, 90deg) {
      check-center("[13CH3]C(=O)O", 0, mirror, rotation, terminal: true)
      check-center("[13C](C)(C)C", 0, mirror, rotation)
      check-center("O", 0, mirror, rotation)
      check-center("F", 0, mirror, rotation)
      check-center("Cl", 0, mirror, rotation)
      check-center("C", 0, mirror, rotation)
      check-center("[13CH4]", 0, mirror, rotation)
    }
  }
}
