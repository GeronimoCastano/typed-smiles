#import "../src/lib.typ": smiles, atom, highlight
#import "../src/chemistry.typ": _compute-layout
#import "../src/mechanism/references.typ": _inline-hydrogen-highlight-fragment, _highlight-hydrogens, _atom-position
#import "../src/molecule/rendering.typ": _rendered-atom-position
#import "../src/mechanism/annotations.typ": _atom-highlight-pieces

#context {
  let font-size = 11pt
  let canvas-scale = 30pt
  let species = (
    layout: _compute-layout("N[C@@H](C)C(=O)O"),
    origin: (0.0, 0.0),
    rotation: 0deg,
    mol-scale: 1.0,
    canvas-scale: canvas-scale,
    actual-font-size: font-size,
    font: "New Computer Modern",
    show-h: (),
  )
  let fragment = _inline-hydrogen-highlight-fragment(species, 6)
  assert(fragment != none)
  let plain-hydrogen-width = measure(text(size: font-size, font: "New Computer Modern", "H")).width / canvas-scale
  assert(fragment.width > plain-hydrogen-width)
  assert(fragment.height > 0)
  let nitrogen-position = _atom-position(species, 0, fragment: "sym")
  let nitrogen-width = measure(text(size: font-size, font: "New Computer Modern", "N")).width / canvas-scale
  assert(calc.abs(calc.abs(fragment.position.at(0) - nitrogen-position.at(0)) - (fragment.width + nitrogen-width) / 2) < 0.00001)
  assert(calc.abs(fragment.position.at(1) - nitrogen-position.at(1)) < 0.00001)
  let highlight-piece = _atom-highlight-pieces(
    (__ref__: "atom", species: 0, index: 6, offset: (0, 0)),
    (radius: auto),
    (species,),
    (lp-offset: 0.6, atom-radius: 0.12, canvas-scale: canvas-scale),
  ).first()
  assert(highlight-piece.kind == "rounded-rectangle")
  let label-padding = font-size / canvas-scale * 0.07
  assert(highlight-piece.width >= fragment.width + 2 * label-padding)
  assert(highlight-piece.height >= fragment.height + 2 * label-padding)
  assert(highlight-piece.radius <= label-padding)
  assert(calc.abs(highlight-piece.centre.at(0) - fragment.position.at(0)) < 0.000001)
  assert(calc.abs(highlight-piece.centre.at(1) - fragment.position.at(1)) < 0.000001)
  let offset-piece = _atom-highlight-pieces(
    (__ref__: "atom", species: 0, index: 6, offset: (0.2, -0.15)),
    (radius: auto),
    (species,),
    (lp-offset: 0.6, atom-radius: 0.12, canvas-scale: canvas-scale),
  ).first()
  assert(calc.abs(offset-piece.centre.at(0) - highlight-piece.centre.at(0) - 0.2) < 0.000001)
  assert(calc.abs(offset-piece.centre.at(1) - highlight-piece.centre.at(1) + 0.15) < 0.000001)
  let large-radius-piece = _atom-highlight-pieces(
    (__ref__: "atom", species: 0, index: 6, offset: (0, 0)),
    (radius: 0.4),
    (species,),
    (lp-offset: 0.6, atom-radius: 0.12, canvas-scale: canvas-scale),
  ).first()
  assert(large-radius-piece.width >= 0.8)
  assert(large-radius-piece.height >= 0.8)

  let ammonium-layout = _compute-layout("[NH4+]")
  let ammonium = species + (layout: ammonium-layout,)
  let ammonium-fragment = _highlight-hydrogens(ammonium, 0).first()
  let nitrogen-position = _rendered-atom-position(ammonium-layout.atoms.first(), 0deg)
  let ammonium-hydrogen = text(size: font-size, font: "New Computer Modern", "H") + sub(
    text(size: font-size, font: "New Computer Modern", "4"),
  )
  let nitrogen-width = measure(text(size: font-size, font: "New Computer Modern", "N")).width / canvas-scale
  let ammonium-hydrogen-width = measure(ammonium-hydrogen).width / canvas-scale
  assert(calc.abs(
    ammonium-fragment.position.at(0)
      - nitrogen-position.x
      - nitrogen-width / 2
      - ammonium-hydrogen-width / 2
  ) < 0.000001)
}

#pagebreak()
= Whole displayed hydrogen-count highlights

#grid(columns: (1fr, 1fr), gutter: 1.2em, row-gutter: 1.4em, align: center + horizon,
  [*Alanine N-H₂, exact index 6* \ #smiles("N[C@@H](C)C(=O)O", show-indices: true, highlight(atom(6)), scale: 1.3)],
  [*Alanine H₂, no badges* \ #smiles("N[C@@H](C)C(=O)O", highlight(atom(6)), scale: 1.3)],
  [*Ammonium, first H* \ #smiles("[NH4+]", show-indices: true, highlight(atom(1)))],
  [*Aminium H₂* \ #smiles("C[NH2+]C", show-indices: true, highlight(atom(3)))],
  [*Aminium, rotation and mirror* \ #smiles("C[NH2+]C", rotation: 45deg, mirror: "horizontal", show-indices: true, highlight(atom(3)))],
  [*Aminium, vertical* \ #smiles("C[NH2+]C", rotation: 90deg, show-indices: true, highlight(atom(3)))],
  [*Aminium, inverted* \ #smiles("C[NH2+]C", rotation: 180deg, show-indices: true, highlight(atom(3)))],
  [*Plain H control* \ #smiles("C[NH]C", show-indices: true, highlight(atom(3)))],
  [*Stereo H control* \ #smiles("C1C[P@H]2CCC1C2", show-indices: true, highlight(atom(7)))],
)
