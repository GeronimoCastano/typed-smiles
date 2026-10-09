#import "../src/lib.typ": smiles, smiles-inline, smiles-cetz, reaction, mol, rxn-arrow, arrow, atom
#import "../src/substructure.typ": _substructure-highlights
#import "../src/chemistry.typ": _compute-layout
#import "../src/molecule/rendering.typ": _mirror-layout, _linearize-skeleton-layout, _rendered-atom-position, _skeleton-hydrogen-label-distance
#import "../src/mechanism/references.typ": _highlight-hydrogens, _atom-position
#import "@preview/cetz:0.5.2"

#context {
  // H predicates still return stable heavy-atom indices; named annotations
  // expand only the group-bearing elements into their displayed H fragments.
  for (molecule, group, element) in (
    ("CCS", "thiol", "S"), ("CC[SH]", "thiol", "S"),
    ("CCS[H]", "thiol", "S"), ("CCS[2H]", "thiol", "S"),
    ("CCO", "alcohol", "O"), ("Oc1ccccc1", "phenol", "O"),
    ("CC(=O)O", "carboxylic-acid", "O"), ("CN", "amine", "N"),
    ("CNC", "amine", "N"), ("CC(=O)N", "amide", "N"),
    ("CC=O", "aldehyde", "C"),
  ) {
    let annotations = _substructure-highlights(molecule, highlight-groups: group)
    assert.eq(annotations.first().at("hydrogen-elements"), (element,))
  }
  assert.eq(_substructure-highlights("CCS", highlight-smarts: "[SX2H1]").first().at("hydrogen-elements"), ())
  assert.eq(_substructure-highlights("CCS", highlight-smarts: "[SX2H1]").first().at("hydrogen-heteroatoms"), true)
  assert.eq(_substructure-highlights("CCS", highlight-smarts: (pattern: "[SX2H1]", include-hydrogens: false)).first().at("hydrogen-heteroatoms"), false)
  assert.eq(_substructure-highlights("CC(=O)O", highlight-groups: (group: "carboxylic-acid", include-hydrogens: false)).first().at("hydrogen-elements"), ())
  assert.eq(_substructure-highlights("CC(=O)O", highlight-groups: (group: "carboxylic-acid", include-hydrogens: true)).first().at("hydrogen-heteroatoms"), true)
  assert.eq(_substructure-highlights("CC(=O)O", highlight-groups: (group: "carboxylic-acid", include-hydrogens: true)).first().at("hydrogen-elements"), ())
  assert.eq(_substructure-highlights("COC", highlight-groups: "ether").first().at("hydrogen-elements"), ())

  // A thiol H stays at the S-H bond length after page-axis rotation/mirroring,
  // independent of notation, origin, or per-molecule scaling.
  for molecule in ("CCS", "CC[SH]", "CCS[H]") {
    let raw = _linearize-skeleton-layout(_compute-layout(molecule))
    for rotation in (0deg, 30deg, 90deg, 180deg) {
      for mirror in (none, "horizontal", "vertical") {
        for scale in (0.7, 1.0, 1.5) {
          let layout = _mirror-layout(raw, mirror, rotation: rotation)
          let species = (layout: layout, rotation: rotation, origin: (2, -1),
            mol-scale: scale, canvas-scale: 30pt, actual-font-size: 11pt * scale, show-h: "skeleton")
          let parent = _atom-position(species, 2)
          let hydrogens = _highlight-hydrogens(species, 2)
          assert.eq(hydrogens.len(), 1)
          let hydrogen = hydrogens.first().position
          let distance = calc.sqrt(calc.pow(hydrogen.at(0) - parent.at(0), 2) + calc.pow(hydrogen.at(1) - parent.at(1), 2))
          assert(calc.abs(distance - _skeleton-hydrogen-label-distance * scale) < 0.000001)
        }
      }
    }
  }

  // The inline SH highlight contains a separate visible H fragment; amine
  // hydrogens and subscripts stay grouped rather than becoming graph atoms.
  let thiol = (layout: _mirror-layout(_compute-layout("CCS"), none),
    rotation: 0deg, origin: (0, 0), show-h: ())
  let hydrogen = _highlight-hydrogens(thiol, 2).first()
  assert(calc.abs(hydrogen.position.at(0) - _atom-position(thiol, 2).at(0)) > 0.001)
  assert.eq(_highlight-hydrogens(thiol, 0), ())
  let retained = thiol + (layout: _mirror-layout(_compute-layout("CCS[2H]"), none))
  assert.eq(_highlight-hydrogens(retained, 2).len(), 1)
  let primary-amine = thiol + (layout: _mirror-layout(_compute-layout("CN"), none))
  let secondary-amine = thiol + (layout: _mirror-layout(_compute-layout("CNC"), none))
  assert(_highlight-hydrogens(primary-amine, 1).first().width > _highlight-hydrogens(secondary-amine, 1).first().width)
  let tertiary-amine = thiol + (layout: _mirror-layout(_compute-layout("CN(C)C"), none))
  assert.eq(_highlight-hydrogens(tertiary-amine, 1), ())
}

= Attached hydrogen highlights

#let blue = (rgb("#BBE1FA"),)
#grid(columns: (1fr, 1fr, 1fr), gutter: 1.5em,
  [*Inline SH* \ #smiles("CCS", highlight-groups: "thiol", highlight-colors: blue)],
  [*Bracket SH* \ #smiles("CC[SH]", highlight-groups: "thiol", highlight-colors: blue)],
  [*Folded explicit H* \ #smiles("CCS[H]", highlight-groups: "thiol", highlight-colors: blue)],
  [*Skeleton S–H* \ #smiles("CCS", show-h: "skeleton", highlight-groups: "thiol", highlight-colors: blue)],
  [*Rotated / mirrored* \ #smiles("CCS", show-h: "skeleton", rotation: 45deg, mirror: "horizontal", highlight-groups: "thiol", highlight-colors: blue)],
  [*Retained isotope H* \ #smiles("CCS[2H]", highlight-groups: "thiol", highlight-colors: blue)],
  [*Alcohol OH* \ #smiles("CCO", highlight-groups: "alcohol")],
  [*Amine NH₂* \ #smiles("CN", highlight-groups: "amine")],
  [*Amide NH₂* \ #smiles("CC(=O)N", highlight-groups: "amide")],
  [*Aldehyde C–H* \ #smiles("CC=O", show-h: "skeleton", highlight-groups: "aldehyde")],
  [*Acid OH* \ #smiles("CC(=O)O", show-h: "skeleton", highlight-groups: "carboxylic-acid")],
  [*Bond-only acid* \ #smiles("CC(=O)O", show-h: "skeleton", highlight-groups: (group: "carboxylic-acid", include-atoms: false))],
)

Inline #smiles-inline("CCS", highlight-groups: "thiol", highlight-colors: blue) thiol.

#reaction(
  mol("CCS", highlight-groups: "thiol", highlight-colors: blue), rxn-arrow(),
  mol("CCS", show-h: "skeleton", rotation: 30deg, highlight-groups: "thiol", highlight-colors: blue),
)
#reaction(
  mol("CCS", highlight-groups: "thiol", highlight-colors: blue), rxn-arrow(),
  mol("CCS", show-h: "skeleton", rotation: 30deg, scale: 0.7, highlight-groups: "thiol", highlight-colors: blue),
  arrow(from: atom(0, 2), to: atom(1, 2)),
)
#context cetz.canvas(length: 30pt, {
  smiles-cetz("CCS", show-h: "skeleton", rotation: 30deg, highlight-groups: "thiol", highlight-colors: blue)
})
