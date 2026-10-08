#import "../../../src/lib.typ": *

#let tint(color) = color.transparentize(45%)
#let core-color = rgb("#A9C7E8")

// All-carbon cyclopenta[a]phenanthrene skeleton: any carbon, any bond.
#let steroid-core = "[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6](~[#6]12)~[#6]~[#6]~[#6]4~[#6]3~[#6]~[#6]~[#6]4"

#let steroids = (
  (molecules.testosterone, [Testosterone]),
  (molecules.estradiol, [Estradiol]),
  ("CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C", [Progesterone]),
  ("C[C@@]12C[C@H](O)[C@@H]3[C@]4(CCC(=O)C=C4CC[C@H]3[C@@H]2CC[C@@]1(C(CO)=O)O)C", [Cortisol]),
)

// Orients the series into the textbook pose: ring A lower left, ring D upper right.
#let aligned-steroids = align-molecules(
  steroids.map(entry => entry.at(0)),
  scaffold: steroid-core,
  rotation: -60deg,
)

#set page(width: 16cm, height: auto, margin: 12pt, fill: none)
#set text(font: "New Computer Modern", size: 11pt)

#block(width: 100%, fill: white, stroke: 0.6pt + luma(215), radius: 8pt, inset: 18pt)[
  #align(center, molecule-grid(
    columns: 2,
    scale: 0.5,
    bond-length: 0.83,
    column-gutter: 1em,
    row-gutter: 0.6em,
    label-gap: 0.8em,
    ..steroids.zip(aligned-steroids).map(((entry, aligned)) => mol(
      aligned,
      color: false,
      highlight-smarts: steroid-core,
      highlight-colors: (tint(core-color),),
      label: text(size: 10pt, entry.at(1)),
    )),
  ))
]
