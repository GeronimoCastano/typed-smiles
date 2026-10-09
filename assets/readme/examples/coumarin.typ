#import "../../../src/lib.typ": smiles, highlight, bond, cycle, step

#let tint(color) = color.transparentize(30%)
#let benzo-color = rgb("#5BB8E6")
#let aryl-color = rgb("#F2A38A")
#let ketone-color = rgb("#F5D23A")
#let amide-color = rgb("#F5A03A")
#let amide-reversed-color = rgb("#D9667A")
#let highlight-radius = 0.2

// Coumarin scaffold atoms: benzo ring 3, 4, 5, 6, 8, 9 (C8a, C8, C7, C6, C5, C4a).
#let benzo-bonds = (bond(3, 4), bond(4, 5), bond(5, 6), bond(6, 8), bond(8, 9), bond(9, 3))

#let benzo-highlight = highlight(
  benzo-bonds,
  fill: tint(benzo-color),
  include-atoms: true,
  radius: highlight-radius,
)

#let scaffold-pose = -135deg

#let coumarin-a = "O=C1Oc2ccc({R})cc2C=C1c1ccc({R})cc1"
#let coumarin-b = "O=C1Oc2ccc({R})cc2C=C1C(=O)c1ccc({R})cc1"
#let coumarin-c = "O=C1Oc2ccc({R})cc2C=C1C(=O)Nc1ccc({R})cc1"
#let coumarin-d = "O=C1Oc2ccc({R})cc2C=C1NC(=O)c1ccc({R})cc1"

#let aryl-highlight-a = highlight(
  (bond(12, 13), bond(13, 14), bond(14, 15), bond(15, 17), bond(17, 18), bond(18, 12)),
  fill: tint(aryl-color),
  include-atoms: true,
  radius: highlight-radius,
)
#let aryl-highlight-b = highlight(
  (bond(14, 15), bond(15, 16), bond(16, 17), bond(17, 19), bond(19, 20), bond(20, 14)),
  fill: tint(aryl-color),
  include-atoms: true,
  radius: highlight-radius,
)
#let aryl-highlight-cd = highlight(
  (bond(15, 16), bond(16, 17), bond(17, 18), bond(18, 20), bond(20, 21), bond(21, 15)),
  fill: tint(aryl-color),
  include-atoms: true,
  radius: highlight-radius,
)

#let ketone-highlight = highlight((bond(12, 13),), fill: tint(ketone-color), include-atoms: true, radius: highlight-radius)
#let amide-highlight = highlight((bond(12, 13), bond(12, 14)), fill: tint(amide-color), include-atoms: true, radius: highlight-radius)
#let reversed-amide-highlight = highlight((bond(12, 13), bond(13, 14)), fill: tint(amide-reversed-color), include-atoms: true, radius: highlight-radius)

// The cycle draws `mol()` string species without their highlights, so each
// scaffold is rendered to content first. The letter and caption are part of
// that content so the cycle's arc clearance accounts for them.
#let scaffold-species(molecule-smiles, letter, caption, annotations) = align(center, stack(
  dir: ttb,
  spacing: 4pt,
  smiles(molecule-smiles, rotation: scaffold-pose, scale: 0.55, color: false, ..annotations),
  text(weight: "bold", size: 13pt)[#letter],
  text(size: 9.5pt)[#caption],
))

#set page(width: 19cm, height: auto, margin: 12pt, fill: none)
#set text(font: "New Computer Modern", size: 11pt)

#block(width: 100%, fill: white, stroke: 0.6pt + luma(215), radius: 8pt, inset: 14pt)[
  #align(center, cycle(
    radius: 6.0,
    start: 90deg,
    scaffold-species(coumarin-a, [A], [MAO-B selective inhibitors], (benzo-highlight, aryl-highlight-a)),
    step(),
    scaffold-species(coumarin-b, [B], [MAO-A/MAO-B inhibitors], (benzo-highlight, aryl-highlight-b, ketone-highlight)),
    step(),
    scaffold-species(coumarin-c, [C], [MAO-B selective inhibitors], (benzo-highlight, aryl-highlight-cd, amide-highlight)),
    step(),
    scaffold-species(coumarin-d, [D], [MAO-B/AChE/BuChE inhibitors], (benzo-highlight, aryl-highlight-cd, reversed-amide-highlight)),
    step(),
  ))
]
