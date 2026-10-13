#import "../src/lib.typ": smiles, atom, arrow
#import "../src/chemistry.typ": _compute-layout
#import "../src/mechanism/references.typ": _atom-position, _atom-arrow-attachment, _highlight-hydrogens

#context {
  for (input, index) in (("CCO", 2), ("CCN", 2), ("C[NH2+]C", 1), ("O", 0)) {
    for rotation in (0deg, 90deg, 180deg, 270deg) {
      let species = (layout: _compute-layout(input), rotation: rotation, origin: (0, 0), mol-scale: 1.0, canvas-scale: 30pt, actual-font-size: 11pt, font: "New Computer Modern", show-h: ())
      let center = _atom-position(species, index)
      let hydrogen = _highlight-hydrogens(species, index).first()
      let dx = hydrogen.position.at(0) - center.at(0)
      let dy = hydrogen.position.at(1) - center.at(1)
      let toward = (center.at(0) + 10 * dx, center.at(1) + 10 * dy)
      let endpoint = _atom-arrow-attachment((index: index), species, toward)
      if calc.abs(dx) > calc.abs(dy) {
        assert(calc.abs(endpoint.at(0) - center.at(0)) >= calc.abs(dx) + hydrogen.width / 2 - 0.00001)
      } else {
        assert(calc.abs(endpoint.at(1) - center.at(1)) >= calc.abs(dy) + hydrogen.height / 2 - 0.00001)
      }
    }
  }
  let species = (layout: _compute-layout("CCN"), rotation: 0deg, origin: (0, 0), mol-scale: 1.0, canvas-scale: 30pt, actual-font-size: 11pt, font: "New Computer Modern", show-h: ())
  let hydrogen = _highlight-hydrogens(species, 2).first()
  let endpoint = _atom-arrow-attachment((index: 3), species, (hydrogen.position.at(0) + 5, hydrogen.position.at(1)))
  assert(endpoint.at(0) >= hydrogen.position.at(0) + hydrogen.width / 2 - 0.00001)

}

= Arrowheads clear inline hydrogen labels
#grid(columns: (1fr, 1fr), gutter: 20pt,
  [*OH from the H side* \ #smiles("CCO", scale: 1.6, arrow(from: atom(0, offset: (4, 0)), to: atom(2), bend: none))],
  [*NH₂ from the H side* \ #smiles("CCN", scale: 1.6, arrow(from: atom(0, offset: (4, 0)), to: atom(2), bend: none))],
  [*Reversed OH* \ #smiles("CCO", rotation: 180deg, scale: 1.6, arrow(from: atom(0, offset: (-4, 0)), to: atom(2), bend: none))],
  [*Rotated OH* \ #smiles("CCO", rotation: 90deg, scale: 1.6, arrow(from: atom(0, offset: (0, 4)), to: atom(2), bend: none))],
  [*Direct NH₂ hydrogen reference* \ #smiles("CCN", scale: 1.6, arrow(from: atom(0, offset: (4, 0)), to: atom(3), bend: none))],
)
