#import "../src/lib.typ": smiles, atom, highlight, arrow
#import "../src/chemistry.typ": _compute-layout
#import "../src/molecule/abbreviations.typ": _orient-abbreviation-labels
#import "../src/molecule/label-fragments.typ": _abbreviation-fragment
#import "../src/molecule/rendering.typ": _abbreviation-fragment-geometry

#context {
  let layout = _compute-layout("CC(=O)O", abbreviation-request: "CO2H")
  assert(layout.atoms.len() == 5)
  assert(layout.abbreviation_groups.first().fragments.map(fragment => fragment.index) == (1, 2, 4))
  assert(_abbreviation-fragment(layout, 2).atoms == (2, 3))
  assert(_abbreviation-fragment(layout, 3) == none)
  for rotation in (0deg, 90deg, 180deg, 270deg) {
    let oriented = _orient-abbreviation-labels(layout, rotation)
    let species = (layout: oriented, origin: (0, 0), rotation: rotation, mol-scale: 1.0, canvas-scale: 30pt, actual-font-size: 11pt, font: "New Computer Modern")
    let oxygen = _abbreviation-fragment-geometry(species, 2)
    let hydrogen = _abbreviation-fragment-geometry(species, 4)
    assert(oxygen != none and hydrogen != none)
    assert(oxygen.width > measure(text(font: "New Computer Modern", size: 11pt, "O")).width / 30pt)
    assert(calc.abs(oxygen.position.at(0) - hydrogen.position.at(0)) > oxygen.width / 2)
  }
}

#include "../research/audit-abbreviation-fragments-2026-10-12.typ"
