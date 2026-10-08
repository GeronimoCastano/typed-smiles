#import "../../../src/lib.typ": *

#let tint(color) = color.transparentize(45%)
#let acetyl-color = rgb("#F5A03A")
#let carbinol-color = rgb("#5FBFA0")
#let acid-color = rgb("#F07C9B")

#set page(width: 16cm, height: auto, margin: 12pt, fill: none)
#set text(font: "New Computer Modern", size: 11pt)

#let step-label(body) = box(width: 2.9cm, align(center, text(size: 8pt, body)))

// Orients all four molecules onto the benzene ring so the rings share one
// orientation; the vertical offsets below put every ring on one horizontal line.
#let ibuprofen-series = (
  "CC(C)Cc1ccccc1",
  "CC(C)Cc1ccc(cc1)C(C)=O",
  "CC(C)Cc1ccc(cc1)C(C)O",
  "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
)
#let aligned = align-molecules(ibuprofen-series, scaffold: "c1ccccc1")
#let ring-offsets = (0.84, 0.01, 0.04, -0.175)

// Each product shades the group its step creates: acetyl, then carbinol, then acid.
#block(width: 100%, fill: white, stroke: 0.6pt + luma(215), radius: 8pt, inset: 18pt)[
  #align(center, reaction(
    gap-h: 0.5em,
    scale: 0.92,
    mol(
      aligned.at(0),
      color: false,
      scale: 0.7,
      offset: (0, ring-offsets.at(0)),
      label: move(dy: 1.1em, step-label[Isobutylbenzene]),
    ),
    rxn-arrow(above: [Ac#sub[2]O], below: ce("HF")),
    mol(
      aligned.at(1),
      color: false,
      scale: 0.7,
      offset: (0, ring-offsets.at(1)),
      highlight-smarts: "[CH3]C(=O)c",
      highlight-colors: (tint(acetyl-color),),
      label: step-label[4-Isobutyl\ acetophenone],
    ),
    rxn-arrow(above: ce("H2"), below: [Raney Ni]),
    mol(
      aligned.at(2),
      color: false,
      scale: 0.7,
      offset: (0, ring-offsets.at(2)),
      highlight-smarts: "[OH1]C(C)c",
      highlight-colors: (tint(carbinol-color),),
      label: step-label[1-(4-Isobutyl\ phenyl)ethanol],
    ),
    rxn-arrow(above: ce("CO"), below: [Pd]),
    mol(
      aligned.at(3),
      color: false,
      scale: 0.7,
      offset: (0, ring-offsets.at(3)),
      highlight-groups: "carboxylic-acid",
      highlight-colors: (tint(acid-color),),
      label: step-label[Ibuprofen],
    ),
  ))
]
