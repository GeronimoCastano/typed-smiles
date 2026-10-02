#import "../src/lib.typ": smiles, smiles-inline, smiles-cetz, substructure-matches, functional-groups, reaction, mol, rxn-arrow, highlight, atom, bond, arrow
#import "../src/substructure.typ": _substructure-highlights
#import "@preview/cetz:0.5.2"

#set text(font: "New Computer Modern", size: 11pt)
#set page(margin: 2cm)

#let aspirin = "CC(=O)OC1=CC=CC=C1C(=O)O"
#let acid = substructure-matches(aspirin, "C(=O)[OX2H1]")
#assert.eq(acid, ((atoms: (10, 11, 12), bonds: ((10, 11), (10, 12))),))
#assert.eq(substructure-matches("CCO", "N"), ())
#assert.eq(substructure-matches("CC", "CC").len(), 1)
#assert.eq(substructure-matches("C1CC1", "CCC").len(), 3)
#assert.eq(substructure-matches("c1ccccc1", "c:c").len(), 6)
#assert.eq(substructure-matches("C1=CC=CC=C1", "c"), ())
#assert.eq(substructure-matches("C[NH3+]", "[NX4H3+]").len(), 1)
#assert.eq(substructure-matches("CCO[H]", "[OX2H1]").len(), 1)
#assert.eq(substructure-matches("CC(=O)N.CN", functional-groups.amine).first().atoms, (5,))

#let highlights = _substructure-highlights("OCCO", highlight-groups: "alcohol", highlight-colors: (yellow, blue))
#assert.eq(highlights.len(), 2)
#assert.eq(highlights.map(it => it.fill), (yellow, blue))
#assert(highlights.all(it => it.include-atoms))
#assert.eq(highlights.map(it => it.ref.first().index), (0, 3))
#assert.eq(_substructure-highlights("CCO", highlight-smarts: "N", highlight-unmatched: "ignore"), ())
#assert.eq(_substructure-highlights("CCO", highlight-groups: ("amine", "alcohol"), highlight-unmatched: "ignore").len(), 1)
#assert.eq(_substructure-highlights(aspirin, highlight-groups: "Carboxylic Acid").len(), 1)
#assert.eq(_substructure-highlights("OCCOCCO", highlight-groups: "alcohol").len(), 2)

#let bonds-only = _substructure-highlights(aspirin,
  highlight-smarts: (pattern: "C=O", include-atoms: false))
#assert.eq(bonds-only.len(), 2)
#assert(bonds-only.all(it => not it.include-atoms and it.ref.all(ref => ref.__ref__ == "bond")))
#let mixed = _substructure-highlights(aspirin,
  highlight-groups: ((group: "carboxylic acid", include-atoms: false), "ester"))
#assert.eq(mixed.map(it => it.include-atoms), (false, true))
#assert.eq(_substructure-highlights("CCO", highlight-groups: (group: "alcohol")).first().include-atoms, true)
#assert.eq(_substructure-highlights("CCO", highlight-smarts: ("O", (pattern: "CC", include-atoms: false))).map(it => it.include-atoms), (true, false))
#let standalone = _substructure-highlights("N.C=O", highlight-smarts: (pattern: "N.C=O", include-atoms: false)).first()
#assert.eq(standalone.ref.map(ref => ref.__ref__), ("atom", "bond"))

= Substructure highlighting

#grid(columns: (1fr, 1fr), gutter: 2em,
  [*Aspirin acid* \ #smiles(aspirin, highlight-smarts: "C(=O)[OX2H1]")],
  [*Two carbonyls* \ #smiles(aspirin, highlight-groups: "carbonyl")],
  [*Two alcohols* \ #smiles("OCCO", highlight-groups: "alcohol", highlight-colors: (rgb("#FFE45C"), rgb("#BBE1FA")))],
  [*Aromatic ring* \ #smiles("Oc1ccccc1", aromatic: "circle", highlight-smarts: "c1ccccc1")],
  [*Ring bonds only* \ #smiles("Oc1ccccc1", aromatic: "circle", highlight-smarts: (pattern: "c1ccccc1", include-atoms: false))],
  [*Carbonyl bonds only* \ #smiles(aspirin, highlight-groups: (group: "carbonyl", include-atoms: false))],
  [*Manual + automatic* \ #smiles("CC(=O)O", highlight(atom(0), fill: rgb("#FFCAD4")), highlight-groups: "carboxylic acid", show-indices: true)],
  [*Rotated / mirrored* \ #smiles("OCCO", highlight-groups: "alcohol", rotation: 45deg, mirror: "horizontal", show-h: "skeleton")],
)

Inline ethanol #smiles-inline("CCO", highlight-groups: "alcohol") in text.

== Scheme
#reaction(
  mol("CC(=O)O", highlight-groups: (group: "carboxylic-acid", include-atoms: false)),
  rxn-arrow(),
  mol("CC(=O)OC", highlight-groups: "ester", rotation: 30deg),
)

== Mechanism canvas
#reaction(
  mol("CC(=O)O", highlight-groups: "carboxylic-acid", scale: 0.8, rotation: 15deg),
  rxn-arrow(above: mol("CCO", highlight-groups: "alcohol")),
  mol("CC(=O)OCC", highlight-groups: (group: "ester", include-atoms: false), mirror: "horizontal"),
  arrow(from: atom(0, 3), to: atom(1, 2)),
)

== Raw CeTZ
#context cetz.canvas(length: 30pt, {
  smiles-cetz("OCCO", name: "diol", highlight-groups: (group: "alcohol", include-atoms: false))
})
