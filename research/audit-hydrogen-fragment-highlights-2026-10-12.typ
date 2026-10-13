#import "../src/lib.typ": smiles, atom, highlight
#set page(width: 22cm, height: auto, margin: 1cm)
#set text(font: "New Computer Modern", size: 11pt)
= Hydrogen-count fragment highlights

#grid(columns: (1fr, 1fr), gutter: 1.5em, row-gutter: 2em, align: center + horizon,
  [*S57: indexed H₂ fragment* \ #smiles("N[C@@H](C)C(=O)O", scale: 1.5, show-indices: true, highlight(atom(6)))],
  [*S57 without indices* \ #smiles("N[C@@H](C)C(=O)O", scale: 1.5, highlight(atom(6)))],
  [*H₂ above N* \ #smiles("C[NH2+]C", scale: 1.5, highlight(atom(3)))],
  [*H₂ below N* \ #smiles("C[NH2+]C", scale: 1.5, rotation: 180deg, highlight(atom(3)))],
  [*Terminal H₃ fragment* \ #smiles("C[NH3+]", scale: 1.5, highlight(atom(2)))],
  [*Isolated H₄ fragment* \ #smiles("[NH4+]", scale: 1.5, highlight(atom(1)))],
  [*Bare N selection* \ #smiles("C[NH2+]C", scale: 1.5, highlight(atom(1), include-hydrogens: false))],
  [*Plain stereo H selection* \ #smiles("C1C[P@H]2CCC1C2", scale: 1.5, highlight(atom(7)))],
)
