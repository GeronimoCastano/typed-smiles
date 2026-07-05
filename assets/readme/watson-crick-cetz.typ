#import "@preview/cetz:0.5.2"
#import "../../src/lib.typ": smiles-cetz

#set page(width: 13cm, height: 5cm, margin: 0.5cm)
#set text(font: "New Computer Modern", size: 10pt)

#align(center, context cetz.canvas(length: 30pt, {
  import cetz.draw: *

  smiles-cetz("Nc1ncnc2N(!s{})cnc12", name: "A")
  smiles-cetz("Cc1cN(!s{})c(=O)[nH]c1=O", name: "T", origin: (4.9, 0.42))

  let hb = (paint: rgb("#3A78C9"), thickness: 1.0pt, dash: "densely-dashed")
  let off(anchor, by) = (rel: by, to: anchor)
  line(off("A.atom-11", (0.4, -0.15)), off("T.atom-9", (-0.2, 0.06)), stroke: hb)
  line(off("A.atom-2", (0.15, 0)), off("T.atom-7", (-0.2, 0)), stroke: hb)

  content((rel: (0.2, 0.2), to: ("A.atom-11", 50%, "T.atom-9")), text(size: 7.5pt, fill: rgb("#3A78C9"))[2.9 Å])
  content((rel: (0, 0.28), to: ("A.atom-2", 50%, "T.atom-7")), text(size: 7.5pt, fill: rgb("#3A78C9"))[2.8 Å])
}))
