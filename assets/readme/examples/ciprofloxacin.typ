#import "../../../src/lib.typ": smiles, highlight, bond

#let ciprofloxacin = "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O"

#let tint(color) = color.transparentize(45%)
#let quinolone-ring-color = rgb("#A7B94A")
#let piperazine-color = rgb("#5FBFA0")
#let fluorine-color = rgb("#C9A06A")
#let ketone-color = rgb("#F07C9B")
#let acid-color = rgb("#5BA4E6")
#let cyclopropyl-color = rgb("#F5A03A")

#set page(width: 16cm, height: auto, margin: 12pt, fill: none)
#set text(font: "New Computer Modern", size: 11pt)

// Later highlights draw on top, so the broad ring capsules come first.
#block(width: 100%, fill: white, stroke: 0.6pt + luma(215), radius: 8pt, inset: 18pt)[
  #align(center, smiles(
    ciprofloxacin,
    color: false,
    scale: 1.1,
    rotation: 90deg,
    highlight(
      (bond(3, 4), bond(4, 5), bond(3, 13), bond(8, 13), bond(8, 9), bond(9, 10), bond(10, 11), bond(11, 12), bond(12, 13)),
      fill: tint(quinolone-ring-color),
      include-atoms: true,
    ),
    highlight(
      (bond(11, 14), bond(14, 15), bond(15, 16), bond(16, 17), bond(17, 18), bond(18, 19), bond(14, 19)),
      fill: tint(piperazine-color),
      include-atoms: true,
    ),
    highlight((bond(10, 20)), fill: tint(fluorine-color), include-atoms: true),
    highlight((bond(5, 6), bond(6, 7), bond(6, 8)), fill: tint(ketone-color), include-atoms: true),
    highlight((bond(5, 21), bond(21, 22), bond(21, 23)), fill: tint(acid-color), include-atoms: true),
    highlight((bond(0, 1), bond(1, 2), bond(0, 2), bond(2, 3)), fill: tint(cyclopropyl-color), include-atoms: true),
  ))
]
