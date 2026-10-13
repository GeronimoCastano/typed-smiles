#import "../src/lib.typ": *
#set page(width: 24cm, height: auto, margin: 1cm)
#set text(font: "New Computer Modern", size: 11pt)
= Element glyph anchors and cycle arrows
#grid(columns: (1fr, 1fr), gutter: 1.5em, align: center + horizon,
  [*Bond aligned with isotope carbon* \ #smiles("[13CH3]C(=O)O", scale: 1.6)],
  [*Oxygen highlight* \ #smiles("O", scale: 2, show-indices: true, highlight(atom(0), include-hydrogens: false))],
)
#v(1em)
*Cycle arcs clear molecule labels*
#align(center, cycle(
  mol("CCO"), step(label: [top-left], rotation: "auto"),
  mol("CC=O"), step(label: [top-right], rotation: "auto"),
  mol("CC(=O)O"), step(label: [bottom-right], rotation: "auto"),
  mol("CCOC"), step(label: [bottom-left], rotation: "auto"),
))
#v(1em)
*Reagents enter near the arc start and leave near its end*
#align(center, cycle(
  mol("CC(=O)O"), step(label: [esterify], into: ce("EtOH"), out: ce("H2O")),
  mol("CC(=O)OCC"), step(label: [hydrolysis], into: ce("H2O"), out: ce("EtOH")),
))
