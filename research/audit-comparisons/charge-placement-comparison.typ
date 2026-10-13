#import "../../src/lib.typ": *
#set page(width: 24cm, height: auto, margin: 1cm)
#set text(font: "New Computer Modern", size: 11pt)
= Charges on hydrogen-bearing labels
#grid(columns: (1fr, 1fr, 1fr), gutter: 1.5em, row-gutter: 2em, align: center + horizon,
  [*Oxonium* \ #smiles("C[OH+]C", lone-pairs: "dots", scale: 1.5)],
  [*Dimethylammonium* \ #smiles("C[NH2+]C", lone-pairs: "dots", scale: 1.5)],
  [*Sulfonium* \ #smiles("C[SH+]C", lone-pairs: "dots", scale: 1.5)],
  [*30°* \ #smiles("C[OH+]C", rotation: 30deg, lone-pairs: "dots", scale: 1.5)],
  [*90°* \ #smiles("C[OH+]C", rotation: 90deg, lone-pairs: "dots", scale: 1.5)],
  [*180°* \ #smiles("C[OH+]C", rotation: 180deg, lone-pairs: "dots", scale: 1.5)],
  [*270°* \ #smiles("C[OH+]C", rotation: 270deg, lone-pairs: "dots", scale: 1.5)],
  [*Horizontal mirror* \ #smiles("C[OH+]C", mirror: "horizontal", lone-pairs: "dots", scale: 1.5)],
  [*Vertical mirror* \ #smiles("C[OH+]C", mirror: "vertical", lone-pairs: "dots", scale: 1.5)],
  [*Compact bonds* \ #smiles("C[OH+]C", bond-length: 0.65, lone-pairs: "dots", scale: 1.5)],
  [*Line pairs* \ #smiles("C[OH+]C", lone-pairs: "lines", scale: 1.5)],
  [*Without electron pairs* \ #smiles("C[OH+]C", scale: 1.5)],
  [*Protonated ring oxygen* \ #smiles("C1C[OH+]CC1", lone-pairs: "dots", scale: 1.5)],
  [*Protonated piperidine* \ #smiles("C1CC[NH2+]CC1", lone-pairs: "dots", scale: 1.5)],
  [*Larger label font* \ #smiles("C[OH+]C", font-size: 15pt, lone-pairs: "dots", scale: 1.5)],
  [*Terminal OH₂⁺* \ #smiles("C[OH2+]", lone-pairs: "dots", scale: 1.5)],
  [*Hydronium* \ #smiles("[OH3+]", lone-pairs: "dots", scale: 1.5)],
  [*Foreground and fading* \ #smiles("C[OH+]C", lone-pairs: "dots", color: false, fg: purple, opacity: 60%, scale: 1.5)],
)
