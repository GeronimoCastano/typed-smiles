#import "../src/lib.typ": *
#set page(width: 24cm, height: auto, margin: 1cm)
#set text(font: "New Computer Modern", size: 11pt)
= Corrected drawing review
#grid(columns: (1fr, 1fr), gutter: 1.5em, align: center + horizon,
  [*Straight skeleton* \ #smiles("OCCO", show-h: "skeleton", scale: 1.3)],
  [*45°* \ #smiles("OCCO", show-h: "skeleton", rotation: 45deg, scale: 1.3)],
  [*45° horizontal mirror* \ #smiles("OCCO", show-h: "skeleton", rotation: 45deg, mirror: "horizontal", scale: 1.3)],
  [*15° vertical mirror* \ #smiles("OCCO", show-h: "skeleton", rotation: 15deg, mirror: "vertical", scale: 1.3)],
)
#v(1em)
#grid(columns: (1fr, 1fr), gutter: 1.5em, align: center + horizon,
  [*Sulfone* \ #smiles("CC(C)(C)c1cc(C)ccc1S(C)(=O)=O", scale: 1.3)],
  [*Centered indices* \ #smiles(molecules.alanine, show-indices: true, highlight(atom(6)), scale: 1.3)],
)
#v(1em)
#grid(columns: (1fr, 1fr, 1fr), gutter: 1.5em, align: center + horizon,
  [*F dots* \ #smiles("CF", lone-pairs: "dots", scale: 1.5)],
  [*Cl dots* \ #smiles("CCl", lone-pairs: "dots", scale: 1.5)],
  [*Br dots* \ #smiles("CBr", lone-pairs: "dots", scale: 1.5)],
  [*F lines* \ #smiles("CF", lone-pairs: "lines", scale: 1.5)],
  [*Cl lines* \ #smiles("CCl", lone-pairs: "lines", scale: 1.5)],
  [*Br lines* \ #smiles("CBr", lone-pairs: "lines", scale: 1.5)],
)
#v(1em)
#align(center, smiles("c1ccccc1", bond-stroke: 3pt, scale: 1.5))
