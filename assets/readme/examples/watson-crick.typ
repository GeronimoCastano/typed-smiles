#import "@preview/cetz:0.5.2"
#import "../../../src/lib.typ": smiles-cetz

#let hydrogen-bond-blue = rgb("#3A78C9")
#let band-yellow = rgb("#FFE45C").transparentize(55%)

#let hydrogen-bond-stroke = (paint: hydrogen-bond-blue, thickness: 1.0pt, dash: "densely-dashed")
#let distance-label(body) = text(size: 7.5pt, fill: hydrogen-bond-blue, body)
#let base-caption(body) = text(size: 8pt, style: "italic", fill: luma(110), body)

#set page(width: 17cm, height: auto, margin: 12pt, fill: none)
#set text(font: "New Computer Modern", size: 10pt)

#block(width: 100%, fill: white, stroke: 0.6pt + luma(215), radius: 8pt, inset: 18pt)[
  #align(center, context cetz.canvas(length: 40pt, {
    import cetz.draw: *

    smiles-cetz("Nc1ncnc2N(!s{})cnc12", name: "A")
    smiles-cetz("Cc1cN(!s{})c(=O)[nH]c1=O", name: "T", origin: (4.9, 0.42))

    // One soft capsule per hydrogen bond, from donor to acceptor, on a lower
    // layer so it sits behind the bases while still using their anchors.
    let band-stroke = (paint: band-yellow, thickness: 30pt, cap: "round")
    let off(anchor, by) = (rel: by, to: anchor)
    on-layer(-1, {
      line(off("A.atom-0", (0, 0.1)), "T.atom-9", stroke: band-stroke)
      line("A.atom-2", "T.atom-7", stroke: band-stroke)
    })

    line(off("A.atom-0", (0.5, 0)), off("T.atom-9", (-0.25, -0.02)), stroke: hydrogen-bond-stroke)
    line(off("A.atom-2", (0.22, -0.05)), off("T.atom-7", (-0.42, 0.02)), stroke: hydrogen-bond-stroke)

    content((2.0, 1.22), distance-label[2.9 Å])
    content((2.8, 0.42), distance-label[2.8 Å])

    content((0.45, -1.95), base-caption[adenine])
    content((4.5, -2.1), base-caption[thymine])
  }))
]
