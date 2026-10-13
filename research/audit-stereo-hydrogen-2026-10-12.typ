#import "../src/lib.typ": *
#set page(width: 22cm, height: auto, margin: 1cm)
#set text(font: "New Computer Modern", size: 11pt)
= Addressable stereo hydrogens
#grid(columns: (1fr, 1fr), gutter: 2em, align: center + horizon,
  [*Phosphorus only* \ #smiles("C1C[P@H]2CCC1C2", scale: 1.5, show-indices: true, highlight(atom(2), include-hydrogens: false))],
  [*Default: include heteroatom H* \ #smiles("C1C[P@H]2CCC1C2", scale: 1.5, show-indices: true, highlight(atom(2)))],
  [*Stereo H only* \ #smiles("C1C[P@H]2CCC1C2", scale: 1.5, show-indices: true, highlight(atom(7)))],
  [*Gray H label and bond end* \ #smiles("C1C[P@H]2CCC1C2", scale: 1.5)],
  [*Monochrome foreground* \ #smiles("C1C[P@H]2CCC1C2", scale: 1.5, color: false, fg: maroon)],
  [*H color override* \ #smiles("C1C[P@H]2CCC1C2", scale: 1.5, atom-colors: (H: blue))],
)
