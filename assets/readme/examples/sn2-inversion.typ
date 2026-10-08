#import "../../../src/lib.typ": *

#let tint(color) = color.transparentize(45%)
#let stereocenter-color = rgb("#5BA4E6")
#let step-label(body) = box(width: 3.2cm, align(center, text(size: 8pt, body)))

#set page(width: 16cm, height: auto, margin: 12pt, fill: none)
#set text(font: "New Computer Modern", size: 11pt)

// Reactant and product are drawn in the same pose, so the leaving group's wedge
// and the incoming OH's hash sit in the same place and the inversion reads at a glance.
#block(width: 100%, height: 4.9cm, fill: white, stroke: 0.6pt + luma(215), radius: 8pt, inset: 18pt)[
  #align(center + horizon, reaction(
    gap-h: 0.6em,
    scale: 0.9,
    mol("[OH-]", color: false, lone-pairs: "dots"),
    [+],
    mol(
      "CC[C@@H](C)Br",
      color: false,
      highlight(atom(2), fill: tint(stereocenter-color), radius: 0.32),
      label: step-label[(R)-2-bromobutane],
    ),
    rxn-arrow(above: [S#sub[N]2]),
    mol(
      "CC[C@H](C)O",
      color: false,
      highlight(atom(2), fill: tint(stereocenter-color), radius: 0.32),
      label: step-label[(S)-butan-2-ol],
    ),
    [+],
    mol("[Br-]", color: false, lone-pairs: "dots"),
    arrow(from: lp(0, 0), to: atom(2, 2), bend: "left", angle: 24deg),
    arrow(from: bond(2, 2, 4), to: atom(2, 4), bend: "left", angle: 70deg),
  ))
]
