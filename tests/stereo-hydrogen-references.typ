#import "../src/chemistry.typ": _compute-layout
#import "../src/molecule/rendering.typ": _mirror-layout, _rendered-atom-position
#import "../src/mechanism/references.typ": _atom-position, _highlight-hydrogens, _resolve-reference

#context {
  let smiles = "C1C[P@H]2CCC1C2"
  let layout = _compute-layout(smiles)
  let stereo-hydrogen-indices = layout.atoms.enumerate()
    .filter(((index, atom)) => atom.at("stereo_virtual_h", default: false))
    .map(((index, _)) => index)
  assert(stereo-hydrogen-indices.len() == 1)
  let hydrogen-index = stereo-hydrogen-indices.first()
  let hydrogen = layout.atoms.at(hydrogen-index)
  assert(hydrogen.virtual_h)
  assert(hydrogen.symbol == "H")
  assert(hydrogen-index == layout.atoms.len() - 1)
  assert(layout.bonds.any(bond => bond.from == 2 and bond.to == hydrogen-index and bond.virtual_bond))

  for mirror in (none, "horizontal", "vertical") {
    let transformed = _mirror-layout(layout, mirror)
    let center = transformed.atoms.at(2)
    let direction = center.stereo_h_dir
    let transformed-hydrogen = transformed.atoms.at(hydrogen-index)
    assert(calc.abs(transformed-hydrogen.pos.x - center.pos.x - direction.x * 0.82) < 0.000001)
    assert(calc.abs(transformed-hydrogen.pos.y - center.pos.y - direction.y * 0.82) < 0.000001)
    for rotation in (0deg, 45deg, 90deg) {
      let species = (
        layout: transformed,
        rotation: rotation,
        origin: (1.25, -0.75),
        mol-scale: 1.0,
        canvas-scale: 30pt,
        actual-font-size: 11pt,
        font: "New Computer Modern",
        show-h: (),
      )
      let expected = _rendered-atom-position(transformed.atoms.at(hydrogen-index), rotation)
      let actual = _atom-position(species, hydrogen-index)
      assert(calc.abs(actual.at(0) - expected.x - 1.25) < 0.000001)
      assert(calc.abs(actual.at(1) - expected.y + 0.75) < 0.000001)
      let reference-position = _resolve-reference(
        (__ref__: "atom", species: 0, index: hydrogen-index, offset: (0, 0)),
        (species,),
        0.6,
      )
      assert(calc.abs(reference-position.at(0) - actual.at(0)) < 0.000001)
      assert(calc.abs(reference-position.at(1) - actual.at(1)) < 0.000001)
      for show-h in ((), "skeleton") {
        let hydrogen-species = species + (show-h: show-h)
        let hydrogen-fragments = _highlight-hydrogens(hydrogen-species, 2)
        assert(hydrogen-fragments.len() == 1)
        assert(calc.abs(hydrogen-fragments.first().position.at(0) - actual.at(0)) < 0.000001)
        assert(calc.abs(hydrogen-fragments.first().position.at(1) - actual.at(1)) < 0.000001)
      }
    }
  }
}
