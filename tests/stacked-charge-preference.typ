#import "../src/lib.typ": smiles

#pagebreak()
= Atom-first charge placement with stacked hydrogens

#grid(columns: (1fr, 1fr, 1fr), gutter: 1.1em, row-gutter: 1.2em, align: center + horizon,
  [*Oxonium, pairs*  #smiles("C[OH+]C", lone-pairs: "dots", scale: 1.6)],
  [*Oxonium, none*  #smiles("C[OH+]C", lone-pairs: none, scale: 1.6)],
  [*Oxonium, lines*  #smiles("C[OH+]C", lone-pairs: "lines", scale: 1.6)],
  [*Aminium*  #smiles("C[NH2+]C", lone-pairs: "dots", scale: 1.6)],
  [*Sulfonium*  #smiles("C[SH+]C", lone-pairs: "dots", scale: 1.6)],
  [*Bare O+ control*  #smiles("C[O+]C", lone-pairs: "dots", scale: 1.6)],
  [*Oxonium 90°*  #smiles("C[OH+]C", rotation: 90deg, lone-pairs: "dots", scale: 1.6)],
  [*Aminium 90°*  #smiles("C[NH2+]C", rotation: 90deg, lone-pairs: none, scale: 1.6)],
  [*Bare O+ 90°*  #smiles("C[O+]C", rotation: 90deg, lone-pairs: "dots", scale: 1.6)],
  [*Oxonium mirror*  #smiles("C[OH+]C", mirror: "horizontal", lone-pairs: "dots", scale: 1.6)],
  [*Oxonium 30°*  #smiles("C[OH+]C", rotation: 30deg, lone-pairs: "dots", scale: 1.6)],
  [*Tight oxonium*  #smiles("C[OH+]C", bond-length: 0.7, font-size: 16pt, lone-pairs: "dots")],
  [*Charged stereo H*  #smiles("C1C[P@H+]2CCC1C2", rotation: 30deg, show-indices: true, lone-pairs: "dots")],
)
