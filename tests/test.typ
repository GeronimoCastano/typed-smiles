#import "../src/lib.typ": smiles, smiles-inline, smiles-cetz, ce, rxn-arrow, mol, reaction, cycle, step, atom, bond, lp, species, arrow, highlight, brackets, mol-weight
#import "@preview/cetz:0.5.2"

#set text(font: "New Computer Modern", size: 11pt)
#set page(margin: 2cm)

= SMILES Rendering Test

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Ethanol* \ #text(size: 8pt, `CCO`) \ #smiles("CCO")],
  [*Ethylene* \ #text(size: 8pt, `C=C`) \ #smiles("C=C")],
  [*Acetylene* \ #text(size: 8pt, `C#C`) \ #smiles("C#C")],

  [*Cyclohexane* \ #text(size: 8pt, `C1CCCCC1`) \ #smiles("C1CCCCC1")],
  [*Benzene* \ #text(size: 8pt, `C1=CC=CC=C1`) \ #smiles("C1=CC=CC=C1")],
  [*Naphthalene* \ #text(size: 8pt, `C1=CC2=CC=CC=C2C=C1`) \ #smiles("C1=CC2=CC=CC=C2C=C1")],

  [*Acetic acid* \ #text(size: 8pt, `CC(=O)O`) \ #smiles("CC(=O)O")],
  [*Isobutane* \ #text(size: 8pt, `CC(C)C`) \ #smiles("CC(C)C")],
  [*Alanine* \ #text(size: 8pt, `CC(N)C(=O)O`) \ #smiles("CC(N)C(=O)O")],

  [*Pyridine* \ #text(size: 8pt, `C1=CC=NC=C1`) \ #smiles("C1=CC=NC=C1")],
  [*Phenol* \ #text(size: 8pt, `OC1=CC=CC=C1`) \ #smiles("OC1=CC=CC=C1")],
  [*Chlorobenzene* \ #text(size: 8pt, `ClC1=CC=CC=C1`) \ #smiles("ClC1=CC=CC=C1")],

  [*Serine* \ #text(size: 8pt, `NC(CO)C(=O)O`) \ #smiles("NC(CO)C(=O)O")],
  [*Furan* \ #text(size: 8pt, `C1=CC=CO1`) \ #smiles("C1=CC=CO1")],
  [*Phenylalanine* \ #text(size: 8pt, `NC(CC1=CC=CC=C1)C(=O)O`) \ #smiles("NC(CC1=CC=CC=C1)C(=O)O")],
)

= Ring System Regression Tests

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Two separated rings* \
   #text(size: 7pt, `CC(N)C(=O)OCCC1=CC=CC=C1NCC1=CC=CC=C1`) \
   #smiles("CC(N)C(=O)OCCC1=CC=CC=C1NCC1=CC=CC=C1", bond-length: 0.9)],

  [*Biphenyl* \
   #text(size: 7pt, `C1=CC=CC=C1C2=CC=CC=C2`) \
   #smiles("C1=CC=CC=C1C2=CC=CC=C2", bond-length: 0.95)],

  [*Three separated rings* \
   #text(size: 7pt, `C1=CC=CC=C1CC2=CC=CC=C2CC3=CC=CC=C3`) \
   #smiles("C1=CC=CC=C1CC2=CC=CC=C2CC3=CC=CC=C3", bond-length: 0.75)],

  [*Fused bicyclic* \
   #text(size: 7pt, `C1=CC=CC(CCC2)=C12`) \
   #smiles("C1=CC=CC(CCC2)=C12", bond-length: 0.95)],
)

#v(1em)
*Steroid-like fused ring system:*

#align(center)[
  #text(size: 7pt, `C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O`) \
  #smiles("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O", bond-length: 0.75)
]

= Atom Charges

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Phenoxide* #super[-] \
   #text(size: 8pt, `[O-]C1=CC=CC=C1`) \
   #smiles("[O-]C1=CC=CC=C1")],

  [*Tetramethyl-N*#super[+] \
   #text(size: 8pt, `[N+](C)(C)(C)C`) \
   #smiles("[N+](C)(C)(C)C")],

  [*Nitro group* \
   #text(size: 8pt, `[N+](=O)[O-]`) \
   #smiles("C[N+](=O)[O-]")],
)

= Chemical Formulas (#raw("ce()") from chemformula)

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: left,

  [Sulfuric acid: #ce("H2SO4")],
  [Ammonium sulfate: #ce("(NH4)2SO4")],
  [Combustion: #ce("CH4 + 2O2 -> CO2 + 2H2O")],
  [Haber process: #ce("N2 + 3H2 <=> 2NH3")],
  [Arrow with label: #ce("->[\"KOAc\"][\"heat\"]")],
  [Disproportionation: #ce("2H2O2 -> 2H2O + O2")],
)

= Reaction Schemes (#raw("reaction()"))

*Bromination of benzene* (horizontal):

#reaction(
  mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[Benzene]),
  rxn-arrow(above: ce("Br2") + [, ] + ce("FeBr3")),
  mol(smiles("BrC1=CC=CC=C1"), label: text(size: 8pt)[Bromobenzene]),
)

#v(1.5em)
*Fischer esterification* (horizontal, with + separator):

#reaction(
  mol(smiles("CC(=O)O"), label: text(size: 8pt)[Acetic acid]),
  [+],
  mol(smiles("CCO"), label: text(size: 8pt)[Ethanol]),
  rxn-arrow(above: ce("H+"), below: [Δ]),
  mol(smiles("CCOC(=O)C"), label: text(size: 8pt)[Ethyl acetate]),
  [+],
  ce("H2O"),
)

#v(1.5em)
*Wrap-around scheme* (right → down → left):

#reaction(
  mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[*1*]),
  rxn-arrow(above: ce("Br2")),
  mol(smiles("BrC1=CC=CC=C1"), label: text(size: 8pt)[*A*]),
  rxn-arrow(dir: "down", above: ce("HNO3"), below: ce("H2SO4")),
  mol(smiles("BrC1=CC(=CC=C1)[N+](=O)[O-]"), label: text(size: 8pt)[*B*]),
  rxn-arrow(dir: "left", above: ce("Fe"), below: ce("HCl")),
  mol(smiles("BrC1=CC(=CC=C1)N"), label: text(size: 8pt)[*C*]),
)

#v(1.5em)
*Equilibrium arrows*:

#reaction(
  mol(smiles("CC(=O)O"), label: text(size: 8pt)[acid]),
  [+],
  mol(smiles("CCO"), label: text(size: 8pt)[alcohol]),
  rxn-arrow(kind: "equilibrium", above: ce("H+"), below: [heat]),
  mol(smiles("CCOC(=O)C"), label: text(size: 8pt)[ester]),
  [+],
  ce("H2O"),
)

#v(1em)
#reaction(
  mol(smiles("N#N"), label: text(size: 8pt)[#ce("N2")]),
  rxn-arrow(kind: "equilibrium", dir: "left", above: [pressure]),
  mol(smiles("[H][H]"), label: text(size: 8pt)[#ce("H2")]),
)

#v(1em)
#reaction(
  mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[A]),
  rxn-arrow(kind: "equilibrium", above: [cat.]),
  mol(smiles("C1=CC=CC=C1O"), label: text(size: 8pt)[B]),
  rxn-arrow(kind: "equilibrium", dir: "down", above: [workup]),
  mol(smiles("C1=CC=CC=C1Cl"), label: text(size: 8pt)[C]),
)

#v(1em)
*Filled equilibrium arrows*:

#reaction(
  mol(smiles("CC(=O)O"), label: text(size: 8pt)[acid]),
  [+],
  mol(smiles("CCO"), label: text(size: 8pt)[alcohol]),
  rxn-arrow(kind: "equilibrium-filled", above: ce("H+"), below: [heat]),
  mol(smiles("CCOC(=O)C"), label: text(size: 8pt)[ester]),
  [+],
  ce("H2O"),
)

#v(1em)
#reaction(
  mol(smiles("N#N"), label: text(size: 8pt)[#ce("N2")]),
  rxn-arrow(kind: "equilibrium-filled", dir: "left", above: [pressure]),
  mol(smiles("[H][H]"), label: text(size: 8pt)[#ce("H2")]),
)

#v(1em)
#reaction(
  mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[A]),
  rxn-arrow(kind: "equilibrium-filled", above: [cat.]),
  mol(smiles("C1=CC=CC=C1O"), label: text(size: 8pt)[B]),
  rxn-arrow(kind: "equilibrium-filled", dir: "down", above: [workup]),
  mol(smiles("C1=CC=CC=C1Cl"), label: text(size: 8pt)[C]),
)

= Options: color and rotation

*No color* (color: false):
#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,
  [*Pyridine* \ #smiles("C1=CC=NC=C1", color: false)],
  [*Alanine* \ #smiles("CC(N)C(=O)O", color: false)],
  [*Bromobenzene* \ #smiles("BrC1=CC=CC=C1", color: false)],
)

#v(1em)
*Rotation* (labels stay upright):
#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,
  [*0°* \ #smiles("CC(N)C(=O)O")],
  [*45°* \ #smiles("CC(N)C(=O)O", rotation: 45deg)],
  [*90°* \ #smiles("CC(N)C(=O)O", rotation: 90deg)],
  [*-90°* \ #smiles("CC(N)C(=O)O", rotation: -90deg)],
)

= Heteroatom Hydrogen Defaults

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Alcohol OH* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.15)],
  [*Carboxylic OH* \ #text(size: 8pt, `CC(=O)O`) \ #smiles("CC(=O)O", bond-length: 1.15)],
  [*Primary amine* \ #text(size: 8pt, `CN`) \ #smiles("CN", bond-length: 1.15)],
  [*Reversed amine* \ #text(size: 8pt, `NC`) \ #smiles("NC", bond-length: 1.15)],

  [*Thiol SH* \ #text(size: 8pt, `CCS`) \ #smiles("CCS", bond-length: 1.15)],
  [*Primary chain amine* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.15)],
  [*Secondary amine* \ #text(size: 8pt, `CNC`) \ #smiles("CNC", bond-length: 1.15)],
  [*Tertiary amine* \ #text(size: 8pt, `CN(C)C`) \ #smiles("CN(C)C", bond-length: 1.15)],
)

#v(1em)
*Bond meets the heteroatom, not the trailing H.* Terminal #raw("O")/#raw("N")
labels center the heavy symbol on the bond terminus and hang the H off to one
side, so the incoming bond — plain, wedge, or hash — always connects to the heavy
atom at any angle. For vertical bonds the symbol sits directly on the bond axis;
wedge and hash bonds are bicolored by atom like plain bonds.

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Vertical OH (up)*],
  [*Vertical #ce("NH2") (down)*],
  [*Bicolor wedge to OH*],
  [*Bicolor hash to #ce("NH2")*],

  [#text(size: 8pt, [`CO` (vertical)]) \ #smiles("CO", bond-length: 1.1, rotation: 120deg)],
  [#text(size: 8pt, [`CN` (steep)]) \ #smiles("CN", bond-length: 1.1, rotation: -60deg)],
  [#text(size: 8pt, `C!wO`) \ #smiles("C!wO", bond-length: 1.3, rotation: 60deg)],
  [#text(size: 8pt, `C!hN`) \ #smiles("C!hN", bond-length: 1.3, rotation: -60deg)],
)

#v(1em)
*Regression molecule with #ce("NH2") and #ce("NH") labels by default:*

#align(center)[
  #text(size: 7pt, `CC(N)C(=O)OCCC1=CC=CC=C1NCC1=CC=CC=C1`) \
  #smiles("CC(N)C(=O)OCCC1=CC=CC=C1NCC1=CC=CC=C1", bond-length: 0.85)
]

= Hydrogen Display Options

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Default hetero H* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2)],
  [*All H* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2, show-h: "all")],
  [*Explicit H* \ #text(size: 8pt, `[NH4+]`) \ #smiles("[NH4+]", bond-length: 1.2)],
  [*Literal label* \ #text(size: 8pt, `{O}`) \ #smiles("CC{O}", bond-length: 1.2)],
)

= Lone Pair Display

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Alcohol dots* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2, lone-pairs: "dots")],
  [*Amine dots* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.2, lone-pairs: "dots")],
  [*Carbonyl dots* \ #text(size: 8pt, `CC=O`) \ #smiles("CC=O", bond-length: 1.2, lone-pairs: "dots")],
  [*Halide lines* \ #text(size: 8pt, `CCl`) \ #smiles("CCl", bond-length: 1.2, lone-pairs: "lines")],

  [*Carboxylate dots* \ #text(size: 8pt, `[O-]C=O`) \ #smiles("[O-]C=O", bond-length: 1.2, lone-pairs: "dots")],
  [*Ether lines* \ #text(size: 8pt, `COC`) \ #smiles("COC", bond-length: 1.2, lone-pairs: "lines")],
  [*Ammonium, no N pair* \ #text(size: 8pt, `[NH4+]`) \ #smiles("[NH4+]", bond-length: 1.2, lone-pairs: "dots")],
  [*No lone pairs* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2)],
)

#v(1em)
*Terminal heteroatom lone-pair orientations:*

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*OH, right* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2, lone-pairs: "dots")],
  [*OH, up* \ #text(size: 8pt, `CCO rotation: 90deg`) \ #smiles("CCO", bond-length: 1.2, rotation: 90deg, lone-pairs: "dots")],
  [*OH, left* \ #text(size: 8pt, `CCO rotation: 180deg`) \ #smiles("CCO", bond-length: 1.2, rotation: 180deg, lone-pairs: "dots")],
  [*OH, down* \ #text(size: 8pt, `CCO rotation: -90deg`) \ #smiles("CCO", bond-length: 1.2, rotation: -90deg, lone-pairs: "dots")],

  [*#ce("NH2"), right* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.2, lone-pairs: "dots")],
  [*#ce("NH2"), up* \ #text(size: 8pt, `CCN rotation: 90deg`) \ #smiles("CCN", bond-length: 1.2, rotation: 90deg, lone-pairs: "dots")],
  [*#ce("NH2"), left* \ #text(size: 8pt, `CCN rotation: 180deg`) \ #smiles("CCN", bond-length: 1.2, rotation: 180deg, lone-pairs: "dots")],
  [*#ce("NH2"), down* \ #text(size: 8pt, `CCN rotation: -90deg`) \ #smiles("CCN", bond-length: 1.2, rotation: -90deg, lone-pairs: "dots")],
)

#v(1em)
*Terminal lone pairs with large labels:*

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Large OH* \ #text(size: 8pt, `font-size: 22pt`) \ #smiles("CCO", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots")],
  [*Large #ce("NH2")* \ #text(size: 8pt, `font-size: 22pt`) \ #smiles("CCN", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots")],
  [*Large OH#super[-]* \ #text(size: 8pt, `[OH-]`) \ #smiles("[OH-]", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots")],
)

#v(1em)
*Internal and multi-connected heteroatom lone pairs:*

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Secondary amine* \ #text(size: 8pt, `CNC`) \ #smiles("CNC", bond-length: 1.2, lone-pairs: "dots")],
  [*Secondary amine rotated* \ #text(size: 8pt, `CNC rotation: 90deg`) \ #smiles("CNC", bond-length: 1.2, rotation: 90deg, lone-pairs: "dots")],
  [*Tertiary amine* \ #text(size: 8pt, `CN(C)C`) \ #smiles("CN(C)C", bond-length: 1.2, lone-pairs: "dots")],
  [*Tertiary amine lines* \ #text(size: 8pt, `CN(C)C`) \ #smiles("CN(C)C", bond-length: 1.2, lone-pairs: "lines")],

  [*Ether* \ #text(size: 8pt, `COC`) \ #smiles("COC", bond-length: 1.2, lone-pairs: "dots")],
  [*Ether rotated* \ #text(size: 8pt, `COC rotation: 90deg`) \ #smiles("COC", bond-length: 1.2, rotation: 90deg, lone-pairs: "dots")],
  [*Ester oxygen* \ #text(size: 8pt, `CC(=O)OC`) \ #smiles("CC(=O)OC", bond-length: 1.2, lone-pairs: "dots")],
  [*Carbonyl plus amide* \ #text(size: 8pt, `CC(=O)N`) \ #smiles("CC(=O)N", bond-length: 1.2, lone-pairs: "dots")],
)

#v(1em)
*Lone pairs with size and font options:*

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Scale 0.8* \ #text(size: 8pt, `scale: 0.8`) \ #smiles("CNC", scale: 0.8, lone-pairs: "dots")],
  [*Scale 1.6* \ #text(size: 8pt, `scale: 1.6`) \ #smiles("CNC", scale: 1.6, lone-pairs: "dots")],
  [*Short bonds* \ #text(size: 8pt, `bond-length: 0.75`) \ #smiles("CNC", bond-length: 0.75, lone-pairs: "dots")],
  [*Long bonds* \ #text(size: 8pt, `bond-length: 1.55`) \ #smiles("CNC", bond-length: 1.55, lone-pairs: "dots")],

  [*Small labels* \ #text(size: 8pt, `font-size: 8pt`) \ #smiles("CCO", font-size: 8pt, bond-length: 1.2, lone-pairs: "dots")],
  [*Large labels* \ #text(size: 8pt, `font-size: 20pt`) \ #smiles("CCO", font-size: 20pt, bond-length: 1.2, lone-pairs: "dots")],
  [*Serif font* \ #text(size: 8pt, `Libertinus Serif`) \ #smiles("CNC", font: "Libertinus Serif", font-size: 18pt, bond-length: 1.25, lone-pairs: "dots")],
  [*Thick lines* \ #text(size: 8pt, `bond-stroke: 1.8pt`) \ #smiles("CCl", bond-stroke: 1.8pt, bond-length: 1.2, lone-pairs: "lines")],
)

#pagebreak()

= Index and Lone-Pair Alignment

Every case in this section enables both atom indices and lone pairs so the
index centers can be checked against the rendered atom glyphs while lone pairs
stay anchored to the same heavy-atom centers.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.7em,
  row-gutter: 1.7em,
  align: center + horizon,

  [*OH right* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*OH up* \ #text(size: 8pt, `CCO rotation: 90deg`) \ #smiles("CCO", bond-length: 1.15, rotation: 90deg, lone-pairs: "dots", show-indices: true)],
  [*OH left* \ #text(size: 8pt, `CCO rotation: 180deg`) \ #smiles("CCO", bond-length: 1.15, rotation: 180deg, lone-pairs: "dots", show-indices: true)],

  [*#ce("NH2") right* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*#ce("NH2") down* \ #text(size: 8pt, `CCN rotation: -90deg`) \ #smiles("CCN", bond-length: 1.15, rotation: -90deg, lone-pairs: "dots", show-indices: true)],
  [*Secondary amine* \ #text(size: 8pt, `CNC`) \ #smiles("CNC", bond-length: 1.15, lone-pairs: "dots", show-indices: true)],

  [*Charged OH#super[-]* \ #text(size: 8pt, `[OH-]`) \ #smiles("[OH-]", bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*Carboxylate* \ #text(size: 8pt, `[O-]C=O`) \ #smiles("[O-]C=O", bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*Amide* \ #text(size: 8pt, `CC(=O)N`) \ #smiles("CC(=O)N", bond-length: 1.15, lone-pairs: "dots", show-indices: true)],

  [*Ether lines* \ #text(size: 8pt, `COC`) \ #smiles("COC", bond-length: 1.15, lone-pairs: "lines", show-indices: true)],
  [*Halide lines* \ #text(size: 8pt, `CCl`) \ #smiles("CCl", bond-length: 1.15, lone-pairs: "lines", show-indices: true)],
  [*Ester rotated* \ #text(size: 8pt, `CC(=O)OC rotation: 45deg`) \ #smiles("CC(=O)OC", bond-length: 1.15, rotation: 45deg, lone-pairs: "dots", show-indices: true)],
)

#v(1em)
*Bare charged atoms:*

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Oxide* \ #text(size: 8pt, `[O-]`) \ #smiles("[O-]", font-size: 22pt, bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*Fluoride* \ #text(size: 8pt, `[F-]`) \ #smiles("[F-]", font-size: 22pt, bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*Chloride* \ #text(size: 8pt, `[Cl-]`) \ #smiles("[Cl-]", font-size: 22pt, bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*Bromide* \ #text(size: 8pt, `[Br-]`) \ #smiles("[Br-]", font-size: 22pt, bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
)

#pagebreak()

*Large labels and stereochemistry:*

#grid(
  columns: (1fr, 1fr),
  gutter: 2em,
  row-gutter: 1.8em,
  align: center + horizon,

  [*Large terminal OH* \ #text(size: 8pt, `C[OH], font-size: 22pt`) \ #smiles("C[OH]", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots", show-indices: true)],
  [*Large terminal #ce("NH2")* \ #text(size: 8pt, `C[NH2], font-size: 22pt`) \ #smiles("C[NH2]", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots", show-indices: true)],

  [*Large OH#super[-]* \ #text(size: 8pt, `[OH-], font-size: 22pt`) \ #smiles("[OH-]", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots", show-indices: true)],
  [*Large alcohol* \ #text(size: 8pt, `CCO, font-size: 22pt`) \ #smiles("CCO", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots", show-indices: true)],

  [*Chiral center* \ #text(size: 8pt, `N[C@@H](C)C(=O)O`) \ #smiles("N[C@@H](C)C(=O)O", bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*Chiral rotated* \ #text(size: 8pt, `rotation: -35deg`) \ #smiles("N[C@@H](C)C(=O)O", bond-length: 1.15, rotation: -35deg, lone-pairs: "dots", show-indices: true)],
)

#v(1em)
*Scale and bond-length variants:*

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.7em,
  row-gutter: 1.7em,
  align: center + horizon,

  [*Scale 0.8* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", scale: 0.8, lone-pairs: "dots", show-indices: true)],
  [*Scale 1.5* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", scale: 1.5, lone-pairs: "dots", show-indices: true)],
  [*Short bonds* \ #text(size: 8pt, `COC bond-length: 0.8`) \ #smiles("COC", bond-length: 0.8, lone-pairs: "lines", show-indices: true)],

  [*Long bonds* \ #text(size: 8pt, `CC(=O)N bond-length: 1.55`) \ #smiles("CC(=O)N", bond-length: 1.55, lone-pairs: "dots", show-indices: true)],
  [*Serif large* \ #text(size: 8pt, `Libertinus Serif`) \ #smiles("CCO", font: "Libertinus Serif", font-size: 20pt, bond-length: 1.2, lone-pairs: "dots", show-indices: true)],
  [*Thick halide* \ #text(size: 8pt, `bond-stroke: 1.8pt`) \ #smiles("CCl", bond-stroke: 1.8pt, bond-length: 1.2, lone-pairs: "lines", show-indices: true)],
)

#pagebreak()

#v(1em)
*Heteroatom H vs literal labels:*

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*SMILES OH* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.15)],
  [*Literal O* \ #text(size: 8pt, `CC{O|O}`) \ #smiles("CC{O|O}", bond-length: 1.15)],
  [*SMILES #ce("NH2")* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.15)],
  [*Literal N* \ #text(size: 8pt, `CC{N|N}`) \ #smiles("CC{N|N}", bond-length: 1.15)],

  [*SMILES NH* \ #text(size: 8pt, `CNC`) \ #smiles("CNC", bond-length: 1.15)],
  [*Literal N red* \ #text(size: 8pt, `C{N|red}C`) \ #smiles("C{N|red}C", bond-length: 1.15)],
  [*SMILES SH* \ #text(size: 8pt, `CCS`) \ #smiles("CCS", bond-length: 1.15)],
  [*Literal S* \ #text(size: 8pt, `CC{S|S}`) \ #smiles("CC{S|S}", bond-length: 1.15)],

  [*SMILES #ce("PH2")* \ #text(size: 8pt, `CP`) \ #smiles("CP", bond-length: 1.15)],
  [*Literal P* \ #text(size: 8pt, `C{P|P}`) \ #smiles("C{P|P}", bond-length: 1.15)],
  [*SMILES #ce("BH2")* \ #text(size: 8pt, `CB`) \ #smiles("CB", bond-length: 1.15)],
  [*Literal B* \ #text(size: 8pt, `C{B}`) \ #smiles("C{B}", bond-length: 1.15)],
)

= Font Options and Large Atom Labels

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Large default label* \
   #text(size: 8pt, `CCN`) \
   #smiles("CCN", font-size: 36pt, bond-length: 1.35)],

  [*Custom atom font* \
   #text(size: 8pt, `font: "Libertinus Serif"`) \
   #smiles("CC(=O)N", font: "Libertinus Serif", font-size: 20pt, bond-length: 1.25)],

  [*Formula font* \
   #text(size: 8pt, `ce(..., font: ..., font-size: 24pt)`) \
   #ce("NH4+ + Cl- -> NH4Cl", font: "New Computer Modern Math", font-size: 24pt)],
)

= Hidden Junction Overlap Regression

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*CC=O small* \
   #text(size: 8pt, `scale: 0.7`) \
   #smiles("CC=O", scale: 0.7)],

  [*CC=O default* \
   #text(size: 8pt, `scale: 1.0`) \
   #smiles("CC=O")],

  [*CC=O large* \
   #text(size: 8pt, `scale: 2.0`) \
   #smiles("CC=O", scale: 2.0)],

  [*Benzene small* \
   #text(size: 8pt, `scale: 0.7`) \
   #smiles("C1=CC=CC=C1", scale: 0.7)],

  [*Benzene default* \
   #text(size: 8pt, `scale: 1.0`) \
   #smiles("C1=CC=CC=C1")],

  [*Benzene large* \
   #text(size: 8pt, `scale: 2.0`) \
   #smiles("C1=CC=CC=C1", scale: 2.0)],
)

= Scale and Bond Stroke Options

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Default* \
   #text(size: 8pt, `C1=CC=CC=C1`) \
   #smiles("C1=CC=CC=C1")],

  [*Scale only* \
   #text(size: 8pt, `scale: 1.6`) \
   #smiles("C1=CC=CC=C1", scale: 1.6)],

  [*Scale + font override* \
   #text(size: 8pt, `scale: 1.6, font-size: 11pt`) \
   #smiles("CCN", scale: 1.6, font-size: 11pt)],

  [*Scale + stroke override* \
   #text(size: 8pt, `scale: 1.6, bond-stroke: 0.8pt`) \
   #smiles("C1=CC=CC=C1", scale: 1.6, bond-stroke: 0.8pt)],

  [*Thick bonds* \
   #text(size: 8pt, `bond-stroke: 2pt`) \
   #smiles("C1=CC=CC=C1", bond-stroke: 2pt)],

  [*Long thin bonds* \
   #text(size: 8pt, `bond-length: 1.6, bond-stroke: 0.8pt`) \
   #smiles("C1=CC=CC=C1", bond-length: 1.6, bond-stroke: 0.8pt)],

  [*Large labels only* \
   #text(size: 8pt, `font-size: 24pt`) \
   #smiles("CCN", font-size: 24pt)],

  [*Scale + all overrides* \
   #text(size: 8pt, `scale: 1.6, manual values`) \
   #smiles("CCN", scale: 1.6, bond-length: 1.15, font-size: 14pt, bond-stroke: 1.8pt)],
)

= Multiple Bond Rendering

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Linear C=O* \ #text(size: 8pt, `CC=O`) \ #smiles("CC=O", bond-length: 1.4)],
  [*Branched C=O* \ #text(size: 8pt, `CC(=O)C`) \ #smiles("CC(=O)C", bond-length: 1.4)],
  [*Linear O=C=O* \ #text(size: 8pt, `O=C=O`) \ #smiles("O=C=O", bond-length: 1.4)],

  [*Zig-zag CCC* \ #text(size: 8pt, `CCC`) \ #smiles("CCC", bond-length: 1.4)],
  [*Linear S=C=S* \ #text(size: 8pt, `S=C=S`) \ #smiles("S=C=S", bond-length: 1.4)],
  [*Linear Br-C#raw("#")N* \ #text(size: 8pt, `BrC#N`) \ #smiles("BrC#N", bond-length: 1.4)],
  [*Linear C-C#raw("#")N* \ #text(size: 8pt, `CC#N`) \ #smiles("CC#N", bond-length: 1.4)],
  [*P-N triple chain* \ #text(size: 8pt, `CP#N`) \ #smiles("CP#N", bond-length: 1.4)],
)

= Stereochemistry and Drawing Extensions

Standard SMILES stereo uses #raw("@")/#raw("@@") for tetrahedral centers and
#raw("/")/#raw("\\") around double bonds. typed-smiles drawing extensions use
#raw("!w") and #raw("!h") for manual wedge/hash bonds.

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Manual wedge*],
  [*Manual hash*],
  [*Tetrahedral @@*],
  [*Tetrahedral @*],

  [#text(size: 8pt, `C!wN`) \ #smiles("C!wN")],
  [#text(size: 8pt, `C!hN`) \ #smiles("C!hN")],
  [#text(size: 8pt, `N[C@@H](C)C(=O)O`) \ #smiles("N[C@@H](C)C(=O)O")],
  [#text(size: 8pt, `N[C@H](C)C(=O)O`) \ #smiles("N[C@H](C)C(=O)O")],
)

#v(1em)
*Stereo regressions* (acyclic zig-zag and ring stereochemical H labels):

#table(
  columns: (1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Chiral alcohol zig-zag*],
  [*Ring stereo H labels*],

  [#text(size: 8pt, `CC[C@@H](O)CC/C=C/CO`) \ #smiles("CC[C@@H](O)CC/C=C/CO", bond-length: 0.85)],
  [#text(size: 8pt, `C[C@]12CC[C@H]3...`) \ #smiles("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O", bond-length: 0.55)],
)

#v(1em)
*Geometric wedge parity* (estradiol). Wedge/hash directions are now derived from
the OpenSMILES neighbor ordering and the actual 2D layout, not a fixed
#raw("@")→up / #raw("@@")→down rule. Adjacent ring-fusion hydrogens at the central
junction point to opposite faces (one wedge, one hash), the 17-OH is wedged with
no redundant hydrogen, and inverting every stereocenter mirrors all wedges.

#table(
  columns: (1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Estradiol (as written)*],
  [*All stereocenters inverted*],

  [#text(size: 7pt, `C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)...`) \
   #smiles("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O", bond-length: 0.7)],
  [#text(size: 7pt, `C[C@@]12CC[C@@H]3[C@@H]([C@H]1CC[C@H]2O)...`) \
   #smiles("C[C@@]12CC[C@@H]3[C@@H]([C@H]1CC[C@H]2O)CCC4=C3C=CC(=C4)O", bond-length: 0.7)],
)

#v(1em)
*E/Z double bonds* (the directional bonds flank the #raw("=")):

#table(
  columns: (1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*trans (E)*],
  [*cis (Z)*],

  [#text(size: 8pt, `F/C=C/F`) \ #smiles("F/C=C/F")],
  [#text(size: 8pt, `F/C=C\F`) \ #smiles("F/C=C\\F")],
)

= Abbreviated Groups

Syntax: use #raw("{label}") anywhere an atom would appear. Use
#raw("{label|style}") to color a label by element symbol or named color.
Standard SMILES bonding rules apply around it.

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Simple label*],
  [*Element color*],
  [*Named color*],
  [*Bracket vs label*],

  [#text(size: 8pt, `{OEt}C=O`) \ #smiles("{OEt}C=O")],

  [#text(size: 8pt, `{>PPh3|P}C=O`) \
   #smiles("{>PPh3|P}C=O")],

  [#text(size: 8pt, `{LG|red}C=O`) \
   #smiles("{LG|red}C=O")],

  [#text(size: 8pt, `[N] vs {N}`) \
   #smiles("[N]") #h(1em) #smiles("{N}")],
)

#v(1em)
*Anchored label attachment points (`>` marker):*

#table(
  columns: (1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Anchor C, 0°*],
  [*Anchor A, 45°*],
  [*Anchor T, -45°*],

  [#text(size: 8pt, `{>CAT}C`) \
   #smiles("{>CAT}C", rotation: 0deg, show-indices: true, highlight(atom(0), fill: rgb("#D8F7C7")))],

  [#text(size: 8pt, `{C>AT}C`) \
   #smiles("{C>AT}C", rotation: 45deg, show-indices: true, highlight(atom(0), fill: rgb("#D8F7C7")))],

  [#text(size: 8pt, `{CA>T}C`) \
   #smiles("{CA>T}C", rotation: -45deg, show-indices: true, highlight(atom(0), fill: rgb("#D8F7C7")))],
)

#v(1em)
*Wedge bonds + abbreviated groups together* (transition-state style):

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Wilkinson's-like* \ #text(size: 8pt, `{PPh3}C({PPh3})=O`) \
   #smiles("{PPh3}C({PPh3})=O")],

  [*With manual wedge* \ #text(size: 8pt, `{Nu}!wC({LG|red})=O`) \
   #smiles("{Nu}!wC({LG|red})=O")],
)

#pagebreak()

= Reaction Scale and Page-Break Behaviour

*L-DOPA synthesis (default scale):*

#align(center)[
  #reaction(
    smiles("C(=O)C1=CC=C(O)C(OC)=C1", rotation:-60deg),
    [+],
    smiles("{NHAc}CC(=O)O"),
    rxn-arrow(above : [1.#ce("Ac2O") \ 2.#ce("H2O")]),
    smiles("OC(=O)C({NHAc})=CC1=CC(OC)=C({AcO})C=C1"),
    rxn-arrow(dir : "down", above : [1. Ligand, Rh#upper("(I)")-Salz / #ce("H2") \ 2. #ce("H+")]),
    mol(smiles("OC(=O)C(!hN)CC1=CC(OC)=C({AcO})C=C1"), label : [*L-DOPA Parkinsonsmittel*])
  )
]

#v(1.5em)
*L-DOPA synthesis (scale: 0.7 — entire scheme shrinks uniformly):*

#align(center)[
  #reaction(
    scale: 0.7,
    smiles("C(=O)C1=CC=C(O)C(OC)=C1", rotation:-60deg),
    [+],
    smiles("{NHAc}CC(=O)O"),
    rxn-arrow(above : [1.#ce("Ac2O") \ 2.#ce("H2O")]),
    smiles("OC(=O)C({NHAc})=CC1=CC(OC)=C({AcO})C=C1"),
    rxn-arrow(dir : "down", above : [1. Ligand, Rh#upper("(I)")-Salz / #ce("H2") \ 2. #ce("H+")]),
    mol(smiles("OC(=O)C(!hN)CC1=CC(OC)=C({AcO})C=C1"), label : [*L-DOPA Parkinsonsmittel*])
  )
]

#v(1.5em)
*Wrap-around scheme (scale: 0.85, breakable: false — moves as a unit if it spills):*

#align(center)[
  #reaction(
    scale: 0.85,
    mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[*1*]),
    rxn-arrow(above: ce("Br2")),
    mol(smiles("BrC1=CC=CC=C1"), label: text(size: 8pt)[*A*]),
    rxn-arrow(dir: "down", above: ce("HNO3"), below: ce("H2SO4")),
    mol(smiles("BrC1=CC(=CC=C1)[N+](=O)[O-]"), label: text(size: 8pt)[*B*]),
    rxn-arrow(dir: "left", above: ce("Fe"), below: ce("HCl")),
    mol(smiles("BrC1=CC(=CC=C1)N"), label: text(size: 8pt)[*C*]),
  )
]

#pagebreak()

= Color Customisation

*Extended named colors in `{label|color}` syntax:*

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*yellow*], [*pink*], [*cyan*], [*brown*],
  [#smiles("{OEt|yellow}C=O")],
  [#smiles("{Nu|pink}C=O")],
  [#smiles("{LG|cyan}C=O")],
  [#smiles("{OMe|brown}C=O")],

  [*lime*], [*teal*], [*navy*], [*silver*],
  [#smiles("{R|lime}C=O")],
  [#smiles("{X|teal}C=O")],
  [#smiles("{Y|navy}C=O")],
  [#smiles("{Z|silver}C=O")],
)

#v(1em)
*Hex RGB colors in `{label|#RRGGBB}` syntax:*

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*`#8B4513`*], [*`#E040FB`*], [*`#00BCD4`*], [*`#FF6F00`*],
  [#smiles("{OMe|#8B4513}C=O")],
  [#smiles("{Nu|#E040FB}C=O")],
  [#smiles("{LG|#00BCD4}C=O")],
  [#smiles("{R|#FF6F00}C=O")],
)

#v(1em)
*`atom-colors` parameter — per-element overrides (O: brown, N: teal):*

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Default CPK* \ #smiles("CC(N)C(=O)O")],
  [*O → brown* \ #smiles("CC(N)C(=O)O", atom-colors: (O: rgb("#8B4513")))],
  [*O → brown, N → teal* \ #smiles("CC(N)C(=O)O", atom-colors: (O: rgb("#8B4513"), N: rgb("#008080")))],
)

#v(1em)
*Preamble-wide defaults via `.with()` — rebind once, use everywhere:*

// Simulate a project-wide preamble binding with custom colors.
#let my-smiles = smiles.with(
  bond-length: 0.9,
  font: "New Computer Modern",
  atom-colors: (O: rgb("#8B4513"), N: rgb("#008080")),
)

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Alanine* \ #my-smiles("CC(N)C(=O)O")],
  [*Ethanol* \ #my-smiles("CCO")],
  [*Serine* \ #my-smiles("NC(CO)C(=O)O")],
)

#v(1em)
*Label name overrides via `"{label}"` keys — override a specific abbreviation's color:*

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Default label colors* \
   #smiles("{PPh3}C({OEt})=O")],

  [*`{PPh3}` → purple, `{OEt}` → teal* \
   #smiles("{PPh3}C({OEt})=O",
     atom-colors: ("{PPh3}": rgb("#7B2D8B"), "{OEt}": rgb("#008080")))],

  [*Mixed: atom + label overrides* \
   #smiles("{Nu}!wC({LG|red})=O",
     atom-colors: (O: rgb("#00897B"), "{Nu}": rgb("#1565C0"), "{LG}": rgb("#B71C1C")))],
)

// #pagebreak()

// #align(center)[
//   #smiles("R({H})({H})CC{Rh}({Ph3P})({Ph3P})(Cl)")
// ]


// #pagebreak()

// = Branch atoms with 3+ substituents (no overlap)

// #grid(
//   columns: (1fr, 1fr, 1fr),
//   gutter: 1.5em,
//   align: center,

//   [*Aminodichloromethane* \ #text(size: 8pt, `C(Cl)(Cl)N`) \ #smiles("C(Cl)(Cl)N")],
//   [*Chloroform* \ #text(size: 8pt, `C(Cl)(Cl)Cl`) \ #smiles("C(Cl)(Cl)Cl")],
//   [*Neopentane* \ #text(size: 8pt, `CC(C)(C)C`) \ #smiles("CC(C)(C)C")],
// )

#pagebreak()

= Reaction mechanisms

== Atom-index overlay (`show-indices`)

A development aid: stamp each atom's writing-order index, then reference it in
`atom()`, `bond()`, and `lp()`.

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*`CC(=O)C`* \ #smiles("CC(=O)C", show-indices: true)],
  [*`C(I)(C)C`* \ #smiles("C(I)(C)C", show-indices: true)],
)

== Reaction-level `show-indices`

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Scheme default on* \
   #reaction(
     show-indices: true,
     mol("CCO", label: text(size: 8pt)[ethanol]),
     rxn-arrow(above: ce("H+")),
     mol("[OH-]", lone-pairs: "dots", label: text(size: 8pt)[hydroxide]),
   )],

  [*Per-molecule opt-out* \
   #reaction(
     show-indices: true,
     mol("CC(=O)O", label: text(size: 8pt)[indexed]),
     rxn-arrow(),
     mol("CCN", show-indices: false, label: text(size: 8pt)[hidden]),
   )],

  [*Mechanism default on* \
   #reaction(
     show-indices: true,
     mol("[OH-]", lone-pairs: "dots"),
     mol("C(Br)(C)C", offset: (1.2, 0.2)),
     arrow(from: lp(0, 0), to: atom(1, 0), bend: "left", color: red),
   )],

  [*Mechanism opt-out* \
   #reaction(
     show-indices: true,
     mol("[OH-]", lone-pairs: "dots"),
     mol("C(Br)(C)C", offset: (1.2, 0.2), show-indices: false),
     arrow(from: lp(0, 0), to: atom(1, 0), bend: "left", color: red),
   )],
)

== Intramolecular arrows and highlights

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*Carbonyl: π → O, with highlight* \
   #smiles(
     "CC(=O)C",
     lone-pairs: "dots",
     highlight(bond(1, 2), fill: rgb("#FFE45C")),
     arrow(from: bond(1, 2), to: atom(2), bend: "right", color: red),
   )],
  [*Highlighted atom + bond* \
   #smiles(
     "OCC=O",
     highlight(atom(0), fill: rgb("#BBE1FA")),
     highlight(bond(2, 3), fill: rgb("#FFCAD4")),
     highlight(atom(4)),
     show-indices: true,
   )],
)

#v(1em)
*Highlight/index alignment regressions:*

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Bracket OH#super[-] atoms* \
   #text(size: 8pt, `[OH-], atom 0 and H atom 1`) \
   #smiles(
     "[OH-]",
     font-size: 22pt,
     lone-pairs: "dots",
     highlight((atom(0), atom(1)), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Terminal bracket H* \
   #text(size: 8pt, `C[OH], H atom 2`) \
   #smiles(
     "C[OH]",
     font-size: 22pt,
     lone-pairs: "dots",
     highlight(atom(2), fill: rgb("#FFE45C")),
     show-indices: true,
   )],

  [*Multiple atom refs* \
   #text(size: 8pt, `highlight((atom(0), atom(2), atom(3)))`) \
   #smiles(
     "CC(=O)O",
     highlight((atom(0), atom(2), atom(3)), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Custom label atoms* \
   #text(size: 8pt, `{PPh3}C({OEt})=O`) \
   #smiles(
     "{PPh3}C({OEt})=O",
     font-size: 20pt,
     highlight((atom(0), atom(2)), fill: rgb("#E7D7F5")),
     show-indices: true,
   )],

  [*Bond only, no atom caps* \
   #text(size: 8pt, `CC#CC, bond(1, 2)`) \
   #smiles(
     "CC#CC",
     highlight(bond(1, 2), fill: rgb("#FFCAD4")),
     show-indices: true,
   )],

  [*Multiple bond refs* \
   #text(size: 8pt, `highlight((bond(0, 1), bond(2, 3)))`) \
   #smiles(
     "CC(=O)OC",
     highlight((bond(0, 1), bond(2, 3)), fill: rgb("#FFCAD4")),
     show-indices: true,
   )],

  [*Bond with endpoint atoms* \
   #text(size: 8pt, `include-atoms: true`) \
   #smiles(
     "CC#CC",
     highlight(bond(1, 2), fill: rgb("#FFCAD4"), include-atoms: true),
     show-indices: true,
   )],

  [*Chain bonds plus atoms* \
   #text(size: 8pt, `array + include-atoms`) \
   #smiles(
     "CCOCC",
     highlight((bond(0, 1), bond(1, 2), bond(2, 3)), fill: rgb("#BBE1FA"), include-atoms: true),
     show-indices: true,
   )],

  [*Rotated OH highlight* \
   #text(size: 8pt, `CCO rotation: 90deg`) \
   #smiles(
     "CCO",
     rotation: 90deg,
     lone-pairs: "dots",
     highlight(atom(2), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Chiral H highlight* \
   #text(size: 8pt, `N[C@@H](C)C(=O)O, H atom 6`) \
   #smiles(
     "N[C@@H](C)C(=O)O",
     lone-pairs: "dots",
     highlight(atom(6), fill: rgb("#FFE45C")),
     show-indices: true,
   )],
)

== Inter-species mechanism (SN2-like)

Hydroxide attacks the central carbon; the C–I bond breaks toward the leaving group.
Indices count `mol()`/content items in written order.

#reaction(
  mol("[OH-]", lone-pairs: "dots"),
  mol("C(I)(C)C", offset: (-1, -1), angle : 10deg),
  arrow(from: lp(0, 0), to: atom(1, 0), bend: "left", label: text(size: 8pt)[attack]),
  // arrow(from: bond(1, 0, 1), to: atom(1, 1, offset: (0.8, 0)), bend: "left", color: red),
)

== String molecule options in mechanism mode

#reaction(
  mol("[OH-]", lone-pairs: "dots", font-size: 14pt, bond-stroke: 1.4pt, show-indices: false),
  mol("C(Br)(C)C", offset: (1.4, 0.35), font-size: 14pt, bond-stroke: 1.4pt, show-indices: false),
  arrow(from: lp(0, 0), to: atom(1, 0), bend: "left", color: red),
)

== Arrow to a `ce()` spectator (blob edge)

#reaction(
  mol("[OH-]", lone-pairs: "dots"),
  mol(ce("FeBr3"), offset: (2.5, 0)),
  arrow(from: lp(0, 0), to: species(1), bend: "right", color: blue, angle:15deg),
)

== Transition-state and intermediate brackets

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*Transition state (‡)* \ #brackets(smiles("CC(=O)C"), sup: [‡])],
  [*Anion bracket* \ #brackets(smiles("CC(=O)C"), sup: [#sym.minus])],
)

#pagebreak()

== Bracket-H indexing (show-indices)

The H-label group of bracket atoms (e.g. "H₄" in [NH4+]) collapses to one
addressable index placed where the label sits. Index 0 = heavy atom, index 1 = H.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 2em,
  align: center + horizon,
  [*#raw("[OH-]")* (idx 0=O, 1=H) \ #smiles("[OH-]", show-indices: true)],
  [*#raw("[NH4+]")* (idx 0=N, 1=H₄) \ #smiles("[NH4+]", show-indices: true)],
  [*Chiral center #raw("[C@@H]")* \ #smiles("N[C@@H](C)C(=O)O", show-indices: true)],
)

The bracket-H badges stay centered on the measured heavy-atom and H glyphs even
when the label size changes.

#grid(
  columns: (1fr, 1fr),
  gutter: 2em,
  align: center + horizon,
  [*Large #raw("[OH-]")* \ #smiles("[OH-]", font-size: 22pt, show-indices: true)],
  [*Large #raw("[NH4+]")* \ #smiles("[NH4+]", font-size: 22pt, show-indices: true)],
)

Terminal bracket atoms use the same measured glyph centers when the H is attached
beside a heavy atom in a larger molecule.

#grid(
  columns: (1fr, 1fr),
  gutter: 2em,
  align: center + horizon,
  [*Terminal #raw("C[OH]")* \ #smiles("C[OH]", font-size: 22pt, show-indices: true)],
  [*Terminal #raw("C[NH2]")* \ #smiles("C[NH2]", font-size: 22pt, show-indices: true)],
)

Full-size subscript and superscript labels should stay legible while index
badges remain centered on the referenced atom or H group.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 2em,
  align: center + horizon,
  [*#raw("[NH4+]")* \ #smiles("[NH4+]", font-size: 22pt, show-indices: true)],
  [*#raw("C[NH3+]")* \ #smiles("C[NH3+]", font-size: 22pt, show-indices: true)],
  [*#raw("C[N+](=O)[O-]")* \ #smiles("C[N+](=O)[O-]", font-size: 22pt, show-indices: true)],
)

Curly arrows can now target the H group by index.

#reaction(
  mol("[OH-]", lone-pairs: "dots", show-indices: true),
  mol("C(Br)(C)C", offset: (1.0, 0.0), show-indices: true),
  arrow(from: lp(0, 0), to: atom(0, 1), bend: "right", color: blue, label: text(size: 7pt)[to H]),
)

== Scheme-mode regression (must look unchanged)

#reaction(
  mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[benzene]),
  rxn-arrow(above: ce("Br2"), below: ce("FeBr3")),
  mol(smiles("BrC1=CC=CC=C1"), label: text(size: 8pt)[bromobenzene]),
)

#pagebreak()

= 0.4.0 Stress Tests

These pages intentionally combine the features most likely to disagree about
coordinates: lone pairs, bracket-H indices, charged labels, highlights, arrows,
rotations, large labels, and reaction-level index overlays.

== Glyph Centers, Lone Pairs, and Highlights

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.4em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Large charged oxide* \
   #text(size: 8pt, `[O-]`) \
   #smiles(
     "[O-]",
     font-size: 24pt,
     lone-pairs: "dots",
     highlight(atom(0), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Large hydroxide split* \
   #text(size: 8pt, `[OH-]`) \
   #smiles(
     "[OH-]",
     font-size: 24pt,
     lone-pairs: "dots",
     highlight((atom(0), atom(1)), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Ammonium H group* \
   #text(size: 8pt, `[NH4+]`) \
   #smiles(
     "[NH4+]",
     font-size: 24pt,
     highlight(atom(1), fill: rgb("#FFE45C")),
     show-indices: true,
   )],

  [*Terminal explicit OH* \
   #text(size: 8pt, `C[OH]`) \
   #smiles(
     "C[OH]",
     font-size: 22pt,
     lone-pairs: "dots",
     highlight((atom(1), atom(2)), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Terminal explicit #ce("NH2")* \
   #text(size: 8pt, `C[NH2]`) \
   #smiles(
     "C[NH2]",
     font-size: 22pt,
     lone-pairs: "dots",
     highlight((atom(1), atom(2)), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Bare halides* \
   #text(size: 8pt, `[F-] [Cl-] [Br-]`) \
   #stack(
     dir: ltr,
     spacing: 1.3em,
     smiles("[F-]", font-size: 20pt, lone-pairs: "dots", highlight(atom(0), fill: rgb("#D8F7C7")), show-indices: true),
     smiles("[Cl-]", font-size: 20pt, lone-pairs: "dots", highlight(atom(0), fill: rgb("#D8F7C7")), show-indices: true),
     smiles("[Br-]", font-size: 20pt, lone-pairs: "dots", highlight(atom(0), fill: rgb("#F6D0D0")), show-indices: true),
   )],

  [*Rotated alcohol* \
   #text(size: 8pt, `CCO rotation: 90deg`) \
   #smiles(
     "CCO",
     rotation: 90deg,
     lone-pairs: "dots",
     highlight(atom(2), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Rotated amine* \
   #text(size: 8pt, `CCN rotation: -90deg`) \
   #smiles(
     "CCN",
     rotation: -90deg,
     lone-pairs: "dots",
     highlight(atom(2), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Serif large labels* \
   #text(size: 8pt, `Libertinus Serif`) \
   #smiles(
     "CC(=O)O",
     font: "Libertinus Serif",
     font-size: 20pt,
     lone-pairs: "dots",
     highlight((atom(2), atom(3)), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],
)

#pagebreak()

== Bonds, Arrays, and Endpoint Toggles

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.4em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Triple bond only* \
   #text(size: 8pt, `CC#CC`) \
   #smiles(
     "CC#CC",
     highlight(bond(1, 2), fill: rgb("#FFCAD4")),
     show-indices: true,
   )],

  [*Triple bond + atoms* \
   #text(size: 8pt, `include-atoms: true`) \
   #smiles(
     "CC#CC",
     highlight(bond(1, 2), fill: rgb("#FFCAD4"), include-atoms: true),
     show-indices: true,
   )],

  [*Alternating chain* \
   #text(size: 8pt, `array of bonds`) \
   #smiles(
     "CCOC(=O)N",
     lone-pairs: "dots",
     highlight((bond(0, 1), bond(1, 2), bond(3, 4)), fill: rgb("#FFCAD4")),
     show-indices: true,
   )],

  [*Chain path + atoms* \
   #text(size: 8pt, `bond array + endpoint atoms`) \
   #smiles(
     "CCOCCN",
     lone-pairs: "dots",
     highlight((bond(0, 1), bond(1, 2), bond(2, 3), bond(3, 4)), fill: rgb("#BBE1FA"), include-atoms: true),
     show-indices: true,
   )],

  [*Carbonyl separated* \
   #text(size: 8pt, `atom and bond colors`) \
   #smiles(
     "CC(=O)N",
     lone-pairs: "dots",
     highlight(bond(1, 2), fill: rgb("#FFCAD4")),
     highlight(atom(3), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],

  [*Ring bond highlight* \
   #text(size: 8pt, `cyclohexene-like`) \
   #smiles(
     "C1=CCCCC1",
     highlight((bond(0, 1), bond(5, 0)), fill: rgb("#FFCAD4"), include-atoms: true),
     show-indices: true,
   )],

  [*Custom labels* \
   #text(size: 8pt, `{PPh3}C({OEt})=O`) \
   #smiles(
     "{PPh3}C({OEt})=O",
     highlight((atom(0), atom(2)), fill: rgb("#E7D7F5")),
     highlight(bond(1, 3), fill: rgb("#FFCAD4")),
     show-indices: true,
   )],

  [*Wedge + label* \
   #text(size: 8pt, `{Nu}!wC({LG|red})=O`) \
   #smiles(
     "{Nu}!wC({LG|red})=O",
     highlight((atom(0), atom(2)), fill: rgb("#BBE1FA")),
     highlight(bond(0, 1), fill: rgb("#FFCAD4"), include-atoms: true),
     show-indices: true,
   )],

  [*Hash + hetero H* \
   #text(size: 8pt, `C!hN`) \
   #smiles(
     "C!hN",
     lone-pairs: "dots",
     highlight((atom(0), atom(1), atom(2)), fill: rgb("#BBE1FA")),
     show-indices: true,
   )],
)

#pagebreak()

== Stereo, Rotation, and Arrow Targets

#grid(
  columns: (1fr, 1fr),
  gutter: 1.8em,
  row-gutter: 1.8em,
  align: center + horizon,

  [*Chiral center, default* \
   #text(size: 8pt, `N[C@@H](C)C(=O)O`) \
   #smiles(
     "N[C@@H](C)C(=O)O",
     lone-pairs: "dots",
     highlight((atom(0), atom(4), atom(5), atom(6)), fill: rgb("#BBE1FA")),
     highlight(bond(1, 0), fill: rgb("#FFCAD4")),
     show-indices: true,
   )],

  [*Chiral center, rotated* \
   #text(size: 8pt, `rotation: -35deg`) \
   #smiles(
     "N[C@@H](C)C(=O)O",
     rotation: -35deg,
     lone-pairs: "dots",
     highlight((atom(0), atom(6)), fill: rgb("#FFE45C")),
     highlight((bond(3, 4), bond(4, 5)), fill: rgb("#FFCAD4"), include-atoms: true),
     show-indices: true,
   )],

  [*Intramolecular carbonyl arrow* \
   #text(size: 8pt, `bond → O`) \
   #smiles(
     "CC(=O)C",
     lone-pairs: "dots",
     highlight(bond(1, 2), fill: rgb("#FFE45C")),
     arrow(from: bond(1, 2), to: atom(2), bend: "right", color: red),
     show-indices: true,
   )],

  [*Arrow to bracket H* \
   #text(size: 8pt, `[OH-] + tert-butyl bromide`) \
   #reaction(
     show-indices: true,
     mol("[OH-]", lone-pairs: "dots"),
     mol("C(Br)(C)C", offset: (1.2, 0.0)),
     highlight(atom(0, 1), fill: rgb("#FFE45C")),
     arrow(from: lp(0, 0), to: atom(0, 1), bend: "right", color: blue, label: text(size: 7pt)[to H]),
   )],

  [*SN2 attack + leaving group* \
   #text(size: 8pt, `reaction(show-indices: true)`) \
   #reaction(
     show-indices: true,
     mol("[OH-]", lone-pairs: "dots"),
     mol("C(Br)(C)C", offset: (1.4, 0.35)),
     highlight(bond(1, 0, 1), fill: rgb("#FFCAD4"), include-atoms: true),
     arrow(from: lp(0, 0), to: atom(1, 0), bend: "left", color: red),
     arrow(from: bond(1, 0, 1), to: atom(1, 1, offset: (0.7, 0)), bend: "left", color: red),
   )],

  [*Mixed species target* \
   #text(size: 8pt, `lp → ce() species`) \
   #reaction(
     show-indices: true,
     mol("[Cl-]", lone-pairs: "dots"),
     mol(ce("AlCl3"), offset: (2.2, 0.1)),
     highlight(atom(0, 0), fill: rgb("#D8F7C7")),
     arrow(from: lp(0, 0), to: species(1), bend: "right", color: blue),
   )],
)

#pagebreak()

== Reaction-Level Index Overlays

#grid(
  columns: (1fr),
  gutter: 1.8em,
  row-gutter: 1.8em,
  align: center + horizon,

  [*Scheme-wide indices* \
   #reaction(
     show-indices: true,
     mol("CCO", lone-pairs: "dots", label: text(size: 8pt)[ethanol]),
     rxn-arrow(above: ce("H+")),
     mol("CC(=O)O", lone-pairs: "dots", label: text(size: 8pt)[acid]),
     rxn-arrow(above: [heat]),
     mol("CCOC(=O)C", lone-pairs: "dots", label: text(size: 8pt)[ester]),
   )],

  [*Scheme opt-out* \
   #reaction(
     show-indices: true,
     mol("[OH-]", lone-pairs: "dots", label: text(size: 8pt)[indexed]),
     rxn-arrow(),
     mol("CCN", lone-pairs: "dots", show-indices: false, label: text(size: 8pt)[hidden]),
     rxn-arrow(),
     mol("[NH4+]", show-indices: true, label: text(size: 8pt)[forced]),
   )],

  [*Mechanism-wide indices* \
   #reaction(
     show-indices: true,
     mol("[O-]C=O", lone-pairs: "dots"),
     mol("CC(Br)C", offset: (1.6, 0.2)),
     highlight((atom(0, 0), bond(1, 1, 2)), fill: rgb("#BBE1FA")),
     arrow(from: lp(0, 0), to: atom(1, 1), bend: "left", color: red),
   )],

  [*Mechanism opt-out and highlights* \
   #reaction(
     show-indices: true,
     mol("[OH-]", lone-pairs: "dots"),
     mol("C(I)(C)C", offset: (1.4, 0.0), show-indices: false),
     highlight(bond(1, 0, 1), fill: rgb("#FFCAD4"), include-atoms: true),
     arrow(from: lp(0, 0), to: atom(1, 0), bend: "left", color: red),
   )],
)

#pagebreak()

= 0.4.0 Deep Stress Tests

== Charged Single-Atom Glyph Centers

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.2em,
  row-gutter: 1.4em,
  align: center + horizon,

  [*O- dots, large* \
   #smiles("[O-]", scale: 1.3, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(0), fill: rgb("#BDE0FE")),
   )],

  [*O- lines, rotated* \
   #smiles("[O-]", scale: 1.15, rotation: 35deg, show-indices: true, lone-pairs: "lines",
     highlight: highlight(atom(0), fill: rgb("#C7F9CC")),
   )],

  [*OH- dots* \
   #smiles("[OH-]", scale: 1.25, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(1)), fill: rgb("#FFD6A5")),
   )],

  [*F-, Cl-, Br-* \
   #stack(dir: ltr, spacing: 1.0em,
     smiles("[F-]", scale: 1.1, show-indices: true, lone-pairs: "dots",
       highlight: highlight(atom(0), fill: rgb("#D8F7C7"))),
     smiles("[Cl-]", scale: 1.1, show-indices: true, lone-pairs: "dots",
       highlight: highlight(atom(0), fill: rgb("#D8F7C7"))),
     smiles("[Br-]", scale: 1.1, show-indices: true, lone-pairs: "dots",
       highlight: highlight(atom(0), fill: rgb("#D8F7C7"))),
   )],

  [*NH4+ bracket Hs* \
   #smiles("[NH4+]", scale: 1.25, show-indices: true,
     highlight: highlight((atom(0), atom(1)), fill: rgb("#E7C6FF")),
   )],

  [*Small font stress* \
   #smiles("[OH-]", scale: 0.9, font-size: 8pt, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(0), fill: rgb("#FDE2E4")),
   )],
)

#pagebreak()

== Terminal Bracket Atoms in Chains

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.1em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*C[OH], default* \
   #smiles("C[OH]", scale: 1.15, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(1), atom(2)), fill: rgb("#BDE0FE")),
   )],

  [*C[OH], rotated* \
   #smiles("C[OH]", scale: 1.15, rotation: 80deg, show-indices: true, lone-pairs: "lines",
     highlight: highlight(atom(1), fill: rgb("#C7F9CC")),
   )],

  [*C[NH2], bracket Hs* \
   #smiles("C[NH2]", scale: 1.15, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(1), atom(2)), fill: rgb("#FFD6A5")),
   )],

  [*C[NH2], 180deg* \
   #smiles("C[NH2]", scale: 1.15, rotation: 180deg, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(1), fill: rgb("#E7C6FF")),
   )],

  [*CC[O-]* \
   #smiles("CC[O-]", scale: 1.15, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(2), fill: rgb("#FDE2E4")),
   )],

  [*CC[NH3+]* \
   #smiles("CC[NH3+]", scale: 1.15, show-indices: true,
     highlight: highlight((atom(2), atom(3)), fill: rgb("#D8F7C7")),
   )],
)

#pagebreak()

== Custom Labels and Label Highlights

#grid(
  columns: (1fr, 1fr),
  gutter: 1.4em,
  row-gutter: 1.6em,
  align: center + horizon,

  [*Wide terminal labels* \
   #smiles("{PPh3}C({OEt})=O", scale: 1.0, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(2), atom(3)), fill: rgb("#BDE0FE")),
   )],

  [*Colored leaving group* \
   #smiles("{Nu}!wC({LG|red})=O", scale: 1.0, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(2), bond(1, 1, 2)), fill: rgb("#FFCAD4")),
   )],

  [*Label rotation stress* \
   #smiles("{OEt}C({PPh3})F", scale: 1.0, rotation: -35deg, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(2)), fill: rgb("#C7F9CC")),
   )],

  [*Label plus bond endpoints* \
   #smiles("{MeO}C(=O){NHPh}", scale: 1.0, show-indices: true, lone-pairs: "dots",
     highlight: highlight(bond(1, 1, 3), fill: rgb("#FFD6A5"), include-atoms: true),
   )],
)

#pagebreak()

== Bond Highlight Arrays

#grid(
  columns: (1fr, 1fr),
  gutter: 1.3em,
  row-gutter: 1.6em,
  align: center + horizon,

  [*Triple bond only* \
   #smiles("CC#CC", scale: 1.1, show-indices: true,
     highlight: highlight((bond(1, 0, 1), bond(1, 1, 2), bond(1, 2, 3)), fill: rgb("#F7C6D0")),
   )],

  [*Triple bond with endpoints* \
   #smiles("CC#CC", scale: 1.1, rotation: 25deg, show-indices: true,
     highlight: highlight((bond(1, 1, 2),), fill: rgb("#BDE0FE"), include-atoms: true),
   )],

  [*Carbonyl branch bonds* \
   #smiles("CCOC(=O)N", scale: 1.05, show-indices: true, lone-pairs: "dots",
     highlight: highlight((bond(1, 2, 3), bond(1, 3, 4), bond(1, 3, 5)), fill: rgb("#C7F9CC")),
   )],

  [*Carbonyl with endpoint toggle* \
   #smiles("CCOC(=O)N", scale: 1.05, rotation: -20deg, show-indices: true, lone-pairs: "dots",
     highlight: highlight((bond(1, 3, 4), bond(1, 3, 5)), fill: rgb("#FFD6A5"), include-atoms: true),
   )],
)

#pagebreak()

== Rings, Rotations, and Overlays

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.0em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Cyclohexane atoms* \
   #smiles("C1CCCCC1", scale: 0.95, show-indices: true,
     highlight: highlight((atom(0), atom(2), atom(4)), fill: rgb("#BDE0FE")),
   )],

  [*Cyclohexane bonds* \
   #smiles("C1CCCCC1", scale: 0.95, rotation: 30deg, show-indices: true,
     highlight: highlight((bond(1, 0, 1), bond(1, 2, 3), bond(1, 4, 5)), fill: rgb("#C7F9CC")),
   )],

  [*Benzene-style uppercase* \
   #smiles("C1=CC=CC=C1", scale: 0.95, show-indices: true,
     highlight: highlight((bond(1, 0, 1), bond(1, 2, 3), bond(1, 4, 5)), fill: rgb("#FFCAD4")),
   )],

  [*Substituted ring O* \
   #smiles("OC1CCCCC1", scale: 0.95, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), bond(1, 0, 1)), fill: rgb("#FFD6A5"), include-atoms: true),
   )],

  [*Ring plus halide* \
   #smiles("ClC1CCCCC1", scale: 0.95, rotation: -45deg, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), bond(1, 0, 1)), fill: rgb("#E7C6FF")),
   )],

  [*Ring plus amine* \
   #smiles("NC1CCCCC1", scale: 0.95, rotation: 85deg, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(7)), fill: rgb("#D8F7C7")),
   )],
)

#pagebreak()

== Stereochemistry and Explicit Hydrogens

#grid(
  columns: (1fr, 1fr),
  gutter: 1.4em,
  row-gutter: 1.6em,
  align: center + horizon,

  [*Chiral center, all fragile indices* \
   #smiles("N[C@@H](C)C(=O)O", scale: 1.1, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(4), atom(5), atom(6)), fill: rgb("#BDE0FE")),
   )],

  [*Chiral center rotated* \
   #smiles("N[C@@H](C)C(=O)O", scale: 1.1, rotation: -35deg, show-indices: true, lone-pairs: "dots",
     highlight: highlight((bond(1, 1, 0), bond(1, 3, 4)), fill: rgb("#FFCAD4"), include-atoms: true),
   )],

  [*Forced hash to N* \
   #smiles("C!hN", scale: 1.25, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(1), atom(2)), fill: rgb("#C7F9CC")),
   )],

  [*Halogen stereocenter load* \
   #smiles("C(Br)(Cl)(F)I", scale: 1.05, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(1), atom(2), atom(3), atom(4)), fill: rgb("#FFD6A5")),
   )],
)

#pagebreak()

== Lone Pair Direction Stress

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.0em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Ether O dots* \
   #smiles("COC", scale: 1.1, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(1), fill: rgb("#BDE0FE")),
   )],

  [*Ether O lines* \
   #smiles("COC", scale: 1.1, rotation: 70deg, show-indices: true, lone-pairs: "lines",
     highlight: highlight(atom(1), fill: rgb("#C7F9CC")),
   )],

  [*Amine N dots* \
   #smiles("CNC", scale: 1.1, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(1), fill: rgb("#FFD6A5")),
   )],

  [*Tertiary amine* \
   #smiles("CN(C)C", scale: 1.1, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(1), fill: rgb("#E7C6FF")),
   )],

  [*Thiol sulfur* \
   #smiles("CCS", scale: 1.1, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(2), atom(3)), fill: rgb("#D8F7C7")),
   )],

  [*Phosphorus label* \
   #smiles("CP", scale: 1.1, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(1), fill: rgb("#FDE2E4")),
   )],
)

#pagebreak()

== Reaction Mechanism Highlight Load

#grid(
  columns: (1fr),
  gutter: 1.6em,
  row-gutter: 1.8em,
  align: center + horizon,

  [*Nucleophile to carbonyl* \
   #reaction(
     show-indices: true,
     mol("[O-]C=O", lone-pairs: "dots"),
     mol("{Nu}", offset: (1.8, -0.2)),
     highlight((atom(0, 0), bond(0, 1, 2), atom(1, 0)), fill: rgb("#BDE0FE")),
     arrow(from: atom(1, 0), to: atom(0, 1), bend: "left", color: red),
   )],

  [*Leaving group chain* \
   #reaction(
     show-indices: true,
     mol("{Nu}", offset: (-1.4, 0.0)),
     mol("CC(Br)C", lone-pairs: "dots"),
     highlight((bond(1, 1, 2), bond(1, 1, 3)), fill: rgb("#FFCAD4"), include-atoms: true),
     arrow(from: atom(0, 0), to: atom(1, 1), bend: "left", color: blue),
     arrow(from: bond(1, 1, 2), to: atom(1, 2), bend: "right", color: red),
   )],

  [*Lone pair to species target* \
   #reaction(
     show-indices: true,
     mol("[Cl-]", lone-pairs: "dots"),
     mol(ce("BF3"), offset: (2.0, 0.1)),
     highlight(atom(0, 0), fill: rgb("#C7F9CC")),
     arrow(from: lp(0, 0), to: species(1), bend: "right", color: blue),
   )],
)

#pagebreak()

== Dense Reaction-Level Index Stress

#grid(
  columns: (1fr),
  gutter: 1.5em,
  row-gutter: 1.8em,
  align: center + horizon,

  [*Five indexed molecules* \
   #reaction(
     show-indices: true,
     gap-h: 1.0em,
     mol("[OH-]", lone-pairs: "dots"),
     rxn-arrow(above: ce("H2O")),
     mol("CC(=O)O", lone-pairs: "dots"),
     rxn-arrow(above: [then]),
     mol("C[NH2]", lone-pairs: "dots"),
     rxn-arrow(above: ce("H+")),
     mol("[NH4+]", lone-pairs: "dots"),
     rxn-arrow(above: ce("Cl-")),
     mol("{PPh3}C({OEt})=O", lone-pairs: "dots"),
   )],

  [*Opt-out inside dense reaction* \
   #reaction(
     show-indices: true,
     gap-h: 1.0em,
     mol("[O-]", lone-pairs: "dots"),
     rxn-arrow(),
     mol("CCOC(=O)N", lone-pairs: "dots", show-indices: false),
     rxn-arrow(),
     mol("N[C@@H](C)C(=O)O", lone-pairs: "dots"),
   )],
)

#pagebreak()

== Large Scale and Tight Layout Stress

#grid(
  columns: (1fr, 1fr),
  gutter: 1.0em,
  row-gutter: 1.4em,
  align: center + horizon,

  [*Large charged oxygen* \
   #smiles("[O-]", scale: 1.8, show-indices: true, lone-pairs: "dots",
     highlight: highlight(atom(0), fill: rgb("#BDE0FE")),
   )],

  [*Large custom label* \
   #smiles("{PPh3}C=O", scale: 1.45, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(2)), fill: rgb("#FFD6A5")),
   )],

  [*Short bonds* \
   #smiles("CCOC(=O)N", scale: 1.0, bond-length: 0.7, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(2), atom(4), atom(5)), fill: rgb("#C7F9CC")),
   )],

  [*Long bonds rotated* \
   #smiles("CCOC(=O)N", scale: 1.0, bond-length: 1.4, rotation: 35deg, show-indices: true, lone-pairs: "dots",
     highlight: highlight((bond(1, 2, 3), bond(1, 3, 4)), fill: rgb("#FFCAD4")),
   )],
)

#pagebreak()

== Combined Worst-Case Sheets

#grid(
  columns: (1fr, 1fr),
  gutter: 1.1em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Charged carboxylate* \
   #smiles("[O-]C(=O)O", scale: 1.05, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(2), atom(3), bond(1, 1, 2)), fill: rgb("#BDE0FE")),
   )],

  [*Rotated charged carboxylate* \
   #smiles("[O-]C(=O)O", scale: 1.05, rotation: -55deg, show-indices: true, lone-pairs: "lines",
     highlight: highlight((atom(0), bond(1, 0, 1), bond(1, 1, 3)), fill: rgb("#C7F9CC"), include-atoms: true),
   )],

  [*Mixed labels and halides* \
   #smiles("{EtO}C(Cl)(Br)F", scale: 1.05, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(2), atom(3), atom(4)), fill: rgb("#FFD6A5")),
   )],

  [*Stereo plus labels* \
   #smiles("{Nu}!wC({LG|blue})(Cl)F", scale: 1.05, show-indices: true, lone-pairs: "dots",
     highlight: highlight((atom(0), atom(2), bond(1, 1, 3)), fill: rgb("#E7C6FF"), include-atoms: true),
   )],
)

#pagebreak()

== Crowded Branch Points

At a four-substituent branch point (e.g. the acetal carbon in `OC(C)(O)O`
fragments), the longest chain must continue away from the rest of the molecule
instead of folding back over an earlier group.

#grid(
  columns: (1fr, 1fr),
  gutter: 1.1em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Acetal diester* \ #text(size: 8pt, `BrCC(=O)OC(C)(O)OC(=O)C`) \
   #smiles("BrCC(=O)OC(C)(O)OC(=O)C", show-indices: true)],

  [*Rotated* \ #text(size: 8pt, `BrCC(=O)OC(C)(O)OC(=O)C`) \
   #smiles("BrCC(=O)OC(C)(O)OC(=O)C", rotation: 30deg)],

  [*Branched diester* \ #text(size: 8pt, `CC(=O)OC(C)(C)OC(=O)CC`) \
   #smiles("CC(=O)OC(C)(C)OC(=O)CC")],
)

#v(1.5em)

A hub whose arms are all identical is started from the hub and drawn
mirror-symmetric about the vertical axis, so equal substituents fan out evenly
instead of curling into a rotational pinwheel. The depiction no longer depends on
which atom the SMILES string starts from.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.1em,
  row-gutter: 1.5em,
  align: center + horizon,

  [*Tetraethylmethane* \ #text(size: 8pt, `CCC(CC)(CC)CC`) \
   #smiles("CCC(CC)(CC)CC")],

  [*Same, hub-first* \ #text(size: 8pt, `C(CC)(CC)(CC)CC`) \
   #smiles("C(CC)(CC)(CC)CC")],

  [*Neopentane* \ #text(size: 8pt, `CC(C)(C)C`) \
   #smiles("CC(C)(C)C")],

  [*Isobutane* \ #text(size: 8pt, `CC(C)C`) \
   #smiles("CC(C)C")],
)



= Explicit [H] atoms fold into their neighbor

Plain bracket hydrogens written as separate `[H]` atoms are folded into the
hydrogen count of the atom they bond to, just like implicit hydrogens. Carbon
hydrogens stay hidden; heteroatom hydrogens (the two amide N–H here) still show.
Charged or isotopic hydrogens such as `[2H]` are kept as drawn atoms.

#align(center)[
  #smiles("C1=C([H])C(OC([H])([H])C(=O)N([H])N([H])C(=O)C([H])([H])OC2=C([H])C(C([H])([H])[H])=C(Cl)C(C([H])([H])[H])=C2[H])=C([H])C([H])=C1OC([H])([H])[H]", scale: 0.5)
]

#grid(
  columns: 2,
  gutter: 1em,
  [*Deuterium kept* \ #text(size: 8pt, `[2H]C([H])([H])[H]`) \
   #smiles("[2H]C([H])([H])[H]")],

  [*All-H still works* \ #text(size: 8pt, `C([H])([H])([H])[H]`, ) \
   #smiles("C([H])([H])([H])[H]", show-h: "all")],
)
= Aromatic SMILES (kekulization)

Lowercase aromatic notation (`c1ccccc1`) is kekulized on input per OpenSMILES:
implicit bonds between aromatic atoms become alternating single/double bonds.
Both notations must render identically.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Benzene, aromatic* \ #text(size: 8pt, `c1ccccc1`) \ #smiles("c1ccccc1")],
  [*Benzene, Kekulé* \ #text(size: 8pt, `C1=CC=CC=C1`) \ #smiles("C1=CC=CC=C1")],
  [*Toluene* \ #text(size: 8pt, `Cc1ccccc1`) \ #smiles("Cc1ccccc1")],

  [*Pyridine* \ #text(size: 8pt, `c1ccncc1`) \ #smiles("c1ccncc1")],
  [*Pyrrole* \ #text(size: 8pt, `c1cc[nH]c1`) \ #smiles("c1cc[nH]c1")],
  [*Furan* \ #text(size: 8pt, `c1occc1`) \ #smiles("c1occc1")],

  [*Thiophene* \ #text(size: 8pt, `c1sccc1`) \ #smiles("c1sccc1")],
  [*Imidazole* \ #text(size: 8pt, `c1cnc[nH]1`) \ #smiles("c1cnc[nH]1")],
  [*N-methylpyrrole* \ #text(size: 8pt, `Cn1cccc1`) \ #smiles("Cn1cccc1")],

  [*Naphthalene* \ #text(size: 8pt, `c1ccc2ccccc2c1`) \ #smiles("c1ccc2ccccc2c1")],
  [*Indane* \ #text(size: 8pt, `c1ccc2CCCc2c1`) \ #smiles("c1ccc2CCCc2c1")],
  [*Biphenyl* \ #text(size: 8pt, `c1ccccc1-c1ccccc1`) \ #smiles("c1ccccc1-c1ccccc1", scale: 0.8)],

  [*2-Pyridinone* \ #text(size: 8pt, `O=c1cccc[nH]1`) \ #smiles("O=c1cccc[nH]1")],
  [*Pyridinium* \ #text(size: 8pt, `c1cc[nH+]cc1`) \ #smiles("c1cc[nH+]cc1")],
  [*Aniline* \ #text(size: 8pt, `Nc1ccccc1`) \ #smiles("Nc1ccccc1")],
)

#v(1em)
Atom indices follow SMILES writing order even through kekulization, so
`show-indices`, `highlight()`, and `arrow()` references address aromatic input
exactly like Kekulé input.

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Indices on aromatic input* \ #text(size: 8pt, `Cc1ccncc1`) \
   #smiles("Cc1ccncc1", show-indices: true)],

  [*Highlights on aromatic input* \ #text(size: 8pt, `c1ccccc1O` + ", " + `bond(0, 1)` + ", " + `atom(6)`) \
   #smiles(
     "c1ccccc1O",
     highlight(bond(0, 1), fill: rgb("#FFE45C"), include-atoms: true),
     highlight(atom(6), fill: rgb("#BBE1FA")),
   )],
)

= Dot-disconnected structures

The dot (`.`) means "no bond": fragments are laid out independently and
arranged left to right in SMILES writing order. Atom indices stay global
across fragments, so annotations keep working.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Sodium acetate* \ #text(size: 8pt, `CC(=O)[O-].[Na+]`) \
   #smiles("CC(=O)[O-].[Na+]")],
  [*Sodium chloride* \ #text(size: 8pt, `[Na+].[Cl-]`) \
   #smiles("[Na+].[Cl-]")],
  [*Ammonium chloride* \ #text(size: 8pt, `[NH4+].[Cl-]`) \
   #smiles("[NH4+].[Cl-]")],

  [*Sodium phenoxide* \ #text(size: 8pt, `[O-]c1ccccc1.[Na+]`) \
   #smiles("[O-]c1ccccc1.[Na+]")],
  [*Amine hydrochloride* \ #text(size: 8pt, `CC[NH3+].[Cl-]`) \
   #smiles("CC[NH3+].[Cl-]")],
  [*Hydrate* \ #text(size: 8pt, `O=C(O)C(=O)O.O.O`) \
   #smiles("O=C(O)C(=O)O.O.O", scale: 0.85)],
)

#v(1em)

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Global indices across fragments* \ #text(size: 8pt, `CC(=O)[O-].[Na+]`) \
   #smiles("CC(=O)[O-].[Na+]", show-indices: true)],

  [*Annotations across fragments* \ #text(size: 8pt, `atom(4)` + ", " + `bond(1, 3)`) \
   #smiles(
     "CC(=O)[O-].[Na+]",
     highlight(atom(4), fill: rgb("#BBE1FA")),
     highlight(bond(1, 3), fill: rgb("#FFE45C")),
   )],
)

= Extended stereo classes and quadruple bonds

Square-planar `@SP1`/`@SP2`/`@SP3` is depicted exactly: neighbors at 90°, in
the cyclic order given by the shape class (U / 4 / Z traced through the
neighbors in SMILES order). `@TB`/`@OH`/`@AL` centers are accepted and drawn
without stereo wedges. `$` renders a quadruple bond.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*cis-Platin (`@SP1`)* \ #text(size: 8pt, `N[Pt@SP1](N)(Cl)Cl`) \
   #smiles("N[Pt@SP1](N)(Cl)Cl")],
  [*trans-Platin (`@SP2`)* \ #text(size: 8pt, `N[Pt@SP2](N)(Cl)Cl`) \
   #smiles("N[Pt@SP2](N)(Cl)Cl")],
  [*Z shape (`@SP3`)* \ #text(size: 8pt, `N[Pt@SP3](N)(Cl)Cl`) \
   #smiles("N[Pt@SP3](N)(Cl)Cl")],

  [*Trigonal bipyramidal* \ #text(size: 8pt, `S[As@TB1](F)(Cl)(Br)N`) \
   #smiles("S[As@TB1](F)(Cl)(Br)N")],
  [*Octahedral* \ #text(size: 8pt, `C[Co@OH1](F)(Cl)(Br)(I)N`) \
   #smiles("C[Co@OH1](F)(Cl)(Br)(I)N")],
  [*Allene (`@AL1`)* \ #text(size: 8pt, `NC(Br)=[C@AL1]=C(O)C`) \
   #smiles("NC(Br)=[C@AL1]=C(O)C")],

  [*Re–Re quadruple bond* \ #text(size: 8pt, `[Re]$[Re]`) \
   #smiles("[Re]$[Re]")],
  [*With ligands* \ #text(size: 8pt, `Cl[Re]$[Re]Cl`) \
   #smiles("Cl[Re]$[Re]Cl")],
  [*Chromium(II) acetate core* \ #text(size: 8pt, `[Cr]$[Cr]`) \
   #smiles("[Cr]$[Cr]")],
)

= Cumulated double bonds label the sp carbon

An invisible carbon between two double bonds would make them read as one long
double bond, so cumulene centers always show an explicit `C`.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Carbon dioxide* \ #text(size: 8pt, `O=C=O`) \ #smiles("O=C=O")],
  [*Ketene* \ #text(size: 8pt, `C=C=O`) \ #smiles("C=C=O")],
  [*Allene* \ #text(size: 8pt, `NC(Br)=[C@AL1]=C(O)C`) \
   #smiles("NC(Br)=[C@AL1]=C(O)C", scale: 0.85)],
  [*Butatriene* \ #text(size: 8pt, `C=C=C=C`) \ #smiles("C=C=C=C")],
)

= Measured label trims

Bonds retreat from a bare atom label (no hydrogens, charge, or isotope) only
as far as the measured glyph extent, so bonds between two labeled atoms —
cumulenes, ethers, metal–metal bonds — no longer look disproportionately
short next to bonds ending at a bare vertex.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Ketene* \ #text(size: 8pt, `C=C=O`) \ #smiles("C=C=O")],
  [*Ether* \ #text(size: 8pt, `COC`) \ #smiles("COC")],
  [*Ester* \ #text(size: 8pt, `CC(=O)OC`) \ #smiles("CC(=O)OC")],
  [*Metal–metal* \ #text(size: 8pt, `Cl[Re]$[Re]Cl`) \ #smiles("Cl[Re]$[Re]Cl")],
)



  #smiles(
    "N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4C(=O)O",
    highlight((bond(0, 5), bond(5, 4), bond(4, 3), bond(3, 6), bond(6, 10), bond(10, 11), bond(11, 15), bond(15, 19), bond(19, 20), bond(20, 21), bond(21, 23), bond(23, 25)), fill : rgb("#96BF0D"), include-atoms:true),
    highlight((bond(15, 16), bond(16, 18), bond(18, 17), bond(17, 16)), fill : rgb("#F29401"), include-atoms:  true),
    highlight((bond(3, 2), bond(2, 1), bond(1, 0)), fill : rgb("#89C7A8"), include-atoms: true),
    highlight((bond(6, 7), bond(7, 8), bond(7, 9), bond(9, 12)), fill : rgb("#C98F4B"), include-atoms: true),
    highlight((bond(11, 12), bond(12, 13), bond(13, 20), bond(13, 14)), fill : rgb("#EC7789"), include-atoms: true),
    highlight((bond(21, 22)), fill : rgb("#0086CB"), include-atoms: true),
    rotation : 90deg,
    bond-stroke : 1.5pt,
  )

= Wavy and dashed bonds (`!s`, `!d`)

`!s` forces a wavy (squiggly) bond — the standard depiction of unspecified
stereochemistry or an attachment point — and `!d` forces a dashed bond for
hydrogen bonds, partial bonds, and coordination. Both are drawing extensions
on a single bond and bicolor at the midpoint like plain bonds.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Wavy* \ #text(size: 8pt, `C!sN`) \ #smiles("C!sN")],
  [*Dashed* \ #text(size: 8pt, `C!dN`) \ #smiles("C!dN")],
  [*Wavy in chain* \ #text(size: 8pt, `CC(!sC)C1=CC=CC=C1`) \
   #smiles("CC(!sC)C1=CC=CC=C1", scale: 0.9)],
  [*Dashed coordination* \ #text(size: 8pt, `CN!d[Cu]`) \ #smiles("CN!d[Cu]")],

  [*Wavy attachment* \ #text(size: 8pt, `OC1CCCC1!sC`) \
   #smiles("OC1CCCC1!sC", scale: 0.9)],
  [*Bicolor wavy* \ #text(size: 8pt, `O!sN`) \ #smiles("O!sN")],
  [*Bicolor dashed* \ #text(size: 8pt, `O!dN`) \ #smiles("O!dN")],
  [*Next to real stereo* \ #text(size: 8pt, `F/C=C/C!sN`) \
   #smiles("F/C=C/C!sN", scale: 0.9)],
)

= Dashed and wavy reaction arrows

`rxn-arrow(kind: "dashed")` and `rxn-arrow(kind: "wavy")`, in both
orientations, alongside the existing kinds.

#reaction(
  mol("CC=O", label: [acetaldehyde]),
  rxn-arrow(kind: "dashed", above: [formal]),
  mol("CCO", label: [ethanol]),
  rxn-arrow(kind: "wavy", above: [hv]),
  mol("C=C", label: [ethylene]),
)

#reaction(
  mol("C1CCCCC1"),
  rxn-arrow(dir: "down", kind: "dashed"),
  mol("C1=CC=CC=C1"),
)

#reaction(
  mol("CCO"),
  rxn-arrow(dir: "down", kind: "wavy"),
  mol("CC=O"),
)

= Molecular weight (`mol-weight`)

`mol-weight(smiles)` returns the molecular weight in g/mol as a float, summing
IUPAC standard atomic weights (the PubChem periodic-table values bundled with
the plugin) over all atoms and implicit/explicit hydrogens. Reference values
below are PubChem computed molecular weights. Wildcards (`*`), abbreviations
(`{...}`), and isotope labels raise a compile error; those paths are covered
by the Rust test suite.

#let mw-row(name, s, reference, digits: 2) = (
  [#name], raw(s), [#calc.round(mol-weight(s), digits: digits)], [#reference],
)

#table(
  columns: (auto, auto, auto, auto),
  align: (left, left, right, right),
  stroke: 0.4pt + rgb("#d8d8d8"),
  table.header([*Molecule*], [*SMILES*], [*mol-weight*], [*PubChem*]),

  ..mw-row([Water], "O", [18.015], digits: 3),
  ..mw-row([Ethanol], "CCO", [46.07]),
  ..mw-row([Benzene (aromatic)], "c1ccccc1", [78.11]),
  ..mw-row([Glucose], "C(C1C(C(C(C(O1)O)O)O)O)O", [180.16]),
  ..mw-row([Caffeine], "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", [194.19]),
  ..mw-row([Sodium acetate], "CC(=O)[O-].[Na+]", [82.03]),
  ..mw-row([Ammonium ion], "[NH4+]", [18.04]),
)

= Mirroring (`mirror`)

`mirror` reflects a molecule horizontally or vertically. Single-axis reflection
exchanges wedges and hashes, so a mirrored stereocenter still depicts the same
configuration.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Default* \ #text(size: 8pt, `mirror: none`) \
   #smiles("CC(=O)OC1=CC=CC=C1C(=O)O", scale: 0.85)],
  [*Horizontal* \ #text(size: 8pt, `mirror: "horizontal"`) \
   #smiles("CC(=O)OC1=CC=CC=C1C(=O)O", scale: 0.85, mirror: "horizontal")],
  [*Vertical* \ #text(size: 8pt, `mirror: "vertical"`) \
   #smiles("CC(=O)OC1=CC=CC=C1C(=O)O", scale: 0.85, mirror: "vertical")],

  [*Default stereocenter* \ #text(size: 8pt, `N[C@@H](C)C(=O)O`) \
   #smiles("N[C@@H](C)C(=O)O", scale: 0.85)],
  [*Horizontal, same config* \ #text(size: 8pt, `mirror: "horizontal"`) \
   #smiles("N[C@@H](C)C(=O)O", scale: 0.85, mirror: "horizontal")],
  [*Mirror with rotation* \ #text(size: 8pt, `mirror + rotation: 90deg`) \
   #smiles("N[C@@H](C)C(=O)O", scale: 0.85, mirror: "horizontal", rotation: 90deg)],
)

Per-molecule orientation inside a reaction:

#reaction(
  mol("OC1CCCCC1", mirror: "horizontal", label: [mirrored]),
  rxn-arrow(above: [ox.]),
  mol("O=C1CCCCC1", label: [ketone]),
)

= Branch collision avoidance

Substituent chains growing from nearby anchor points (e.g. ortho ring
positions) no longer land on top of atoms another branch already placed: the
blocking branch flips to the other side of its attachment bond, or the
colliding bond flips its zigzag turn to the opposite ideal slot. Bond angles
are never bent to arbitrary values — every resolution keeps the ideal angles
(120° at a chain carbon). In aspirin the acetyl `C=O` and the carboxyl `OH`
previously coincided exactly.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Aspirin* \ #text(size: 8pt, `CC(=O)OC1=CC=CC=C1C(=O)O`) \
   #smiles("CC(=O)OC1=CC=CC=C1C(=O)O", scale: 0.85)],
  [*Ortho chains* \ #text(size: 8pt, `CCCC1=CC=CC=C1CCC`) \
   #smiles("CCCC1=CC=CC=C1CCC", scale: 0.85)],
  [*Salicylic acid* \ #text(size: 8pt, `OC1=CC=CC=C1C(=O)O`) \
   #smiles("OC1=CC=CC=C1C(=O)O", scale: 0.85)],
)

= Dark theme and foreground (`fg`, `theme`)

`fg: auto` (the default) inherits the surrounding text color, so molecules
recolor automatically on dark slides; `theme: auto` switches to a dark CPK
variant (lifted N, O, Br, I and named label colors) when the foreground is
light. Reaction arrows inherit the text color the same way.

#block(fill: rgb("#1E1E24"), inset: 10pt, radius: 4pt, width: 100%)[
  #set text(fill: white)
  #grid(
    columns: (1fr, 1fr, 1fr),
    gutter: 1.5em,
    align: center + horizon,

    [*Caffeine, inherited* \ #text(size: 8pt, `fg: auto`) \
     #smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", scale: 0.8)],
    [*Dark N, Br, I* \ #text(size: 8pt, `NC(Br)C(I)C(=O)O`) \
     #smiles("NC(Br)C(I)C(=O)O", scale: 0.8)],
    [*Labels on dark* \ #text(size: 8pt, `{X|navy}C(=O){Y|maroon}`) \
     #smiles("{X|navy}C(=O){Y|maroon}", scale: 0.8)],
  )
  #reaction(
    mol("CCO", label: text(size: 8pt)[ethanol]),
    rxn-arrow(above: text(size: 8pt)[ox.]),
    mol("CC=O", label: text(size: 8pt)[acetaldehyde]),
    rxn-arrow(kind: "equilibrium-filled"),
    mol("CC(O)O", label: text(size: 8pt)[hydrate]),
  )
]

Explicit values still override: `fg` recolors bonds and carbon labels on any
background, and `theme: "dark"` can be forced independently.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Explicit `fg: navy`* \ #smiles("CC(N)C(=O)O", fg: navy)],
  [*`color: false` uses `fg`* \ #smiles("CC(N)C(=O)O", color: false, fg: maroon)],
  [*Forced dark palette on white* \ #text(size: 8pt, `theme: "dark"`) \
   #smiles("NC(Br)C(I)C(=O)O", theme: "dark")],
)


= Aromatic ring circles (`aromatic: "circle"`)

Rings written in aromatic (lowercase) notation can draw as single bonds with
an inscribed circle instead of alternating double bonds. Kekulé-written input
keeps its explicit bonds; in fused systems each fully aromatic ring gets its
own circle.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Benzene, Kekulé style* \ #text(size: 8pt, `c1ccccc1`) \
   #smiles("c1ccccc1")],
  [*Benzene, circle* \ #text(size: 8pt, `aromatic: "circle"`) \
   #smiles("c1ccccc1", aromatic: "circle")],
  [*Naphthalene* \ #text(size: 8pt, `c1ccc2ccccc2c1`) \
   #smiles("c1ccc2ccccc2c1", aromatic: "circle", scale: 0.9)],
  [*Pyridine* \ #text(size: 8pt, `Cc1ccncc1`) \
   #smiles("Cc1ccncc1", aromatic: "circle")],

  [*Only aromatic ring circled* \ #text(size: 8pt, `c1ccc2CCCc2c1`) \
   #smiles("c1ccc2CCCc2c1", aromatic: "circle", scale: 0.9)],
  [*Kekulé input unaffected* \ #text(size: 8pt, `C1=CC=CC=C1`) \
   #smiles("C1=CC=CC=C1", aromatic: "circle")],
  [*Circle + mirror* \ #text(size: 8pt, `mirror: "horizontal"`) \
   #smiles("Cc1ccncc1", aromatic: "circle", mirror: "horizontal")],
  [*Circle + rotation* \ #text(size: 8pt, `rotation: 30deg`) \
   #smiles("Cc1ccncc1", aromatic: "circle", rotation: 30deg)],
)

= Atom annotations (`atom-annotations`)

Small gray side labels placed on the emptiest side of the atom. Entries are
tuples: `(index, content)` or `(index, content, offset)`.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Offsets* \ #text(size: 8pt, `((1, [$alpha$], (...)), ...)`) \
   #smiles("N[C@@H](C)C(=O)O", atom-annotations: (
     (1, [$alpha$], (-0.4, -0.05)),
     (2, [$beta$], (0.2, -0.3)),
     (3, [$gamma$])
   ))],
  [*No offsets* \ #text(size: 8pt, `((1, [α]), (2, [β]))`) \
   #smiles("N[C@@H](C)C(=O)O", atom-annotations: (
     (1, [α]),
     (2, [β]),
   ))],
  [*Next to subscript labels* \ #text(size: 8pt, `NH₂ and OH untouched`) \
   #smiles(
     "NCCO",
     atom-annotations: ((0, [N1]), (3, [O1], (0.18, 0))),
   )],
  [*Custom style* \ #text(size: 8pt, `text(fill: red, size: 9pt)[..]`) \
   #smiles(
     "CC(=O)OC1=CC=CC=C1C(=O)O",
     scale: 0.8,
     atom-annotations: (
       (1, text(fill: red, size: 9pt)[Ac]),
       (10, text(fill: red, size: 9pt)[C1], (0.14, -0.10)),
     ),
   )],
)

Annotations also work per molecule inside a reaction:

#reaction(
  mol("CC(=O)O", atom-annotations: ((1, [C1]),)),
  rxn-arrow(above: text(size: 8pt)[Δ]),
  mol("CC=O", atom-annotations: ((1, [C1], (0.18, 0)),)),
)

= Per-atom hydrogen display (`show-h`)

`show-h` labels selected implicit hydrogens on carbon. It accepts a single
index, an index array, or `"all"`.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Default* \ #text(size: 8pt, `CC(N)C(=O)O`) \
   #smiles("CC(N)C(=O)O")],
  [*`show-h: 1`* \ #text(size: 8pt, [central C-H only]) \
   #smiles("CC(N)C(=O)O", show-h: 1)],
  [*`show-h: (0, 1)`* \ #text(size: 8pt, [methyl too]) \
   #smiles("CC(N)C(=O)O", show-h: (0, 1))],
  [*`show-h: "all"`* \ #text(size: 8pt, [everything]) \
   #smiles("CC(N)C(=O)O", show-h: "all")],
)

#pagebreak()

= Bond customizations (`bond-customizations`)

Per-bond style overrides keyed by `bond(i, j)` references: `color`, `stroke`
(bond width), and `opacity`.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Colored bond* \ #text(size: 8pt, `(bond(1, 2), (color: red))`) \
   #smiles("CCCC", bond-customizations: ((bond(1, 2), (color: red)),))],
  [*Thick bond* \ #text(size: 8pt, `(stroke: 2.2pt)`) \
   #smiles("CCCC", bond-customizations: ((bond(1, 2), (stroke: 2.2pt)),))],
  [*Faded bond* \ #text(size: 8pt, `(opacity: 30%)`) \
   #smiles("CCCC", bond-customizations: ((bond(1, 2), (opacity: 30%)),))],
  [*Combined* \ #text(size: 8pt, `color + stroke together`) \
   #smiles("CCCC", bond-customizations: (
     (bond(1, 2), (color: yellow, stroke: 1.8pt)),
   ))],

  [*Double bond* \ #text(size: 8pt, `both lines recolor`) \
   #smiles("CC=CC", bond-customizations: ((bond(1, 2), (color: red)),))],
  [*Ring bond* \ #text(size: 8pt, `(i, j) pair key form`) \
   #smiles("C1=CC=CC=C1", bond-customizations: (((0, 1), (color: red, stroke: 1.6pt)),))],
  [*Wedge and hash* \ #text(size: 8pt, `C!wC(!hN)O`) \
   #smiles("C!wC(!hN)O", bond-customizations: (
     (bond(0, 1), (color: red)),
     (bond(1, 2), (color: blue)),
   ))],
  [*Several bonds* \ #text(size: 8pt, [breaking bond red, forming blue]) \
   #smiles("CC(=O)OCC", bond-customizations: (
     (bond(1, 3), (color: red, opacity: 55%)),
     (bond(3, 4), (color: blue, stroke: 1.6pt)),
   ))],
)

Inside a reaction, string molecules take the same option:

#reaction(
  mol("CC(=O)OC", bond-customizations: ((bond(1, 3), (color: red)),)),
  rxn-arrow(above: ce("H2O")),
  mol("CC(=O)O"),
  [+],
  ce("CH3OH"),
)

= Molecule opacity (`opacity`)

`opacity` fades bonds, labels, charges, and lone pairs together, for ghost or
de-emphasized species.

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*100%* \ #smiles("CC(N)C(=O)O")],
  [*60%* \ #smiles("CC(N)C(=O)O", opacity: 60%)],
  [*30%* \ #smiles("CC(N)C(=O)O", opacity: 30%)],
  [*Float form, 0.15* \ #smiles("CC(N)C(=O)O", opacity: 0.15)],

  [*Monochrome ghost* \ #text(size: 8pt, `color: false`) \
   #smiles("c1ccccc1O", color: false, opacity: 35%)],
  [*Lone pairs fade too* \ #text(size: 8pt, `lone-pairs: "dots"`) \
   #smiles("CO", lone-pairs: "dots", opacity: 40%)],
  [*Charges fade too* \ #text(size: 8pt, `[NH4+]`) \
   #smiles("[NH4+]", opacity: 40%)],
  [*Ghost next to solid* \
   #reaction(
     mol("CCO", opacity: 30%),
     rxn-arrow(),
     mol("CC=O"),
   )],
)

= Inline molecules (`smiles-inline`)

`smiles-inline` scales a structure to a target height and baseline-aligns it
so it reads inline; the paragraph below keeps its normal line spacing.

#block(width: 100%)[
  The dehydration of ethanol #smiles-inline("CCO") over alumina gives ethylene
  #smiles-inline("C=C"), while oxidation gives acetaldehyde
  #smiles-inline("CC=O") and, further, acetic acid #smiles-inline("CC(=O)O").
  Aromatic solvents such as toluene #smiles-inline("Cc1ccccc1") are common;
  pyridine #smiles-inline("c1ccncc1") is a classic base. This filler sentence
  only exists so that the paragraph wraps over several lines and the line
  spacing above and below each inline structure can be inspected.
]

Heights compare as follows: default 1.4em #smiles-inline("c1ccccc1"), compact
1em #smiles-inline("c1ccccc1", height: 1em), large 2.2em
#smiles-inline("c1ccccc1", height: 2.2em) (only its own line grows), and a
lowered baseline #smiles-inline("c1ccccc1", baseline: 0.1em).

Inline molecules accept the usual drawing options: red oxygen off
#smiles-inline("CC(=O)O", color: false), rotated #smiles-inline("CCO", rotation: 90deg),
and faded #smiles-inline("CCO", opacity: 40%).

#pagebreak()

= Journal style presets (`style`)

Presets fill in bond length, label size, stroke, font, and monochrome color
from the journal's published drawing settings. `"default"` applies nothing;
explicit arguments always win.

#grid(
  columns: (1fr, 1fr, 1fr, 1fr, 1fr),
  gutter: 1em,
  align: center + horizon,

  [*default* \ #text(size: 8pt, [30 pt, 11 pt NCM]) \
   #smiles("CC(N)C(=O)O")],
  [*acs* \ #text(size: 8pt, [14.4 pt, 10 pt sans]) \
   #smiles("CC(N)C(=O)O", style: "acs")],
  [*rsc* \ #text(size: 8pt, [12.2 pt, 7 pt sans]) \
   #smiles("CC(N)C(=O)O", style: "rsc")],
  [*nature* \ #text(size: 8pt, [10.8 pt, 6 pt sans]) \
   #smiles("CC(N)C(=O)O", style: "nature")],
  [*wiley* \ #text(size: 8pt, [17 pt, 12 pt sans]) \
   #smiles("CC(N)C(=O)O", style: "wiley")],
)

Presets compose with everything else and stay overridable:

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*acs + scale 1.5* \ #text(size: 8pt, [whole preset scales]) \
   #smiles("c1ccccc1O", style: "acs", scale: 1.5)],
  [*acs + explicit font-size* \ #text(size: 8pt, `font-size: 14pt`) \
   #smiles("c1ccccc1O", style: "acs", font-size: 14pt)],
  [*acs + color true* \ #text(size: 8pt, [CPK colors opt in]) \
   #smiles("c1ccccc1O", style: "acs", color: true)],
  [*nature + circle + wedges* \ #text(size: 8pt, [options unaffected]) \
   #smiles("C!wC(!hN)c1ccccc1", style: "nature", aromatic: "circle")],
)

= CeTZ integration (`smiles-cetz`)

`smiles-cetz` draws a molecule inside an existing CeTZ canvas and registers
`atom-<i>`, `bond-<i>-<j>`, and `center` anchors, so arbitrary CeTZ drawing
attaches to real molecular positions. The Watson–Crick A–T base pair below is
two `smiles-cetz` calls plus plain CeTZ: two dashed hydrogen bonds drawn
between offset atom-anchor endpoints, distance labels, and glycosidic-bond
stubs to the sugar backbone.

#align(center, context cetz.canvas(length: 30pt, {
  import cetz.draw: *

  // Adenine on the left, pairing edge (N6, N1) facing right.
  smiles-cetz("Nc1ncnc2N(!s{})cnc12", name: "A")
  // Thymine on the right, pairing edge (N3-H, O4) facing left.
  smiles-cetz("Cc1cN(!s{})c(=O)[nH]c1=O", name: "T", origin: (4.9, 0.42))

  let hb = (paint: rgb("#3A78C9"), thickness: 1.0pt, dash: "densely-dashed")
  let off(anchor, by) = (rel: by, to: anchor)
  line(off("A.atom-11", (0.4, -0.15)), off("T.atom-9", (-0.2, 0.06)), stroke: hb)
  line(off("A.atom-2", (0.15, 0)), off("T.atom-7", (-0.2, 0)), stroke: hb)

  // Hydrogen-bond distances.
  content((rel: (0.2, 0.2), to: ("A.atom-11", 50%, "T.atom-9")), text(size: 7.5pt, fill: rgb("#3A78C9"))[2.9 Å])
  content((rel: (0, 0.28), to: ("A.atom-2", 50%, "T.atom-7")), text(size: 7.5pt, fill: rgb("#3A78C9"))[2.8 Å])

}))

Same anchor-offset pattern in a second CeTZ composition: a polarized
donor-acceptor contact with a soft colored band behind the dashed line.

#align(center, context cetz.canvas(length: 30pt, {
  import cetz.draw: *

  smiles-cetz("CC(=O)C", name: "acceptor", origin: (0, 0))
  smiles-cetz("[H]F", name: "donor", origin: (3.4, -0.05))

  let off(anchor, by) = (rel: by, to: anchor)
  let o = off("acceptor.atom-2", (0.10, 0.02))
  let h = off("donor.atom-0", (-0.10, 0.02))
  line(o, h, stroke: (paint: rgb("#6C63FF").transparentize(78%), thickness: 6pt, cap: "round"))
  line(o, h, stroke: (paint: rgb("#2D6CDF"), thickness: 0.9pt, dash: "densely-dashed"))

  content(off("acceptor.atom-2", (-0.18, 0.34)), text(size: 8pt, fill: rgb("#1565C0"))[$delta^-$])
  content(off("donor.atom-0", (0.18, 0.34)), text(size: 8pt, fill: rgb("#B71C1C"))[$delta^+$])
  content((2.05, 0.78), text(size: 7.5pt, fill: rgb("#5E35B1"))[polar contact])
}))

A larger composition — the serine-protease catalytic triad — packs four
`smiles-cetz` molecules into one scene with an active-site pocket highlight, a
nucleophile highlight disk, dashed charge-relay hydrogen bonds between anchors,
a distance label, curved electron-pushing arrows (Ser attacks the substrate
carbonyl; the C=O breaks), lone pairs, and residue labels.

#align(center, context cetz.canvas(length: 30pt, {
  import cetz.draw: *
  let off(a, by) = (rel: by, to: a)

  // Active-site pocket highlight behind the triad.
  on-layer(-1, {
    rect((-1.1, -1.7), (6.2, 1.6), radius: 0.4, stroke: none, fill: rgb("#FFF3C4"))
  })

  // The three catalytic residues, left to right: Asp — His — Ser.
  smiles-cetz("CC(=O)[O-]", name: "asp", origin: (0, 0.4), lone-pairs: "dots")
  smiles-cetz("c1cnc[nH]1", name: "his", origin: (2.7, -0.1))
  smiles-cetz("OCC", name: "ser", origin: (5.1, -0.5))
  // Substrate carbonyl, upper right, under nucleophilic attack.
  smiles-cetz("CC(=O)NC", name: "sub", origin: (7.6, 1.2))

  // Highlight the serine oxygen: the nucleophile.
  circle("ser.atom-0", radius: 0.34, stroke: none, fill: rgb("#FFB74D").transparentize(45%))

  // Charge-relay hydrogen bonds.
  let hb = (paint: rgb("#3A78C9"), thickness: 1.0pt, dash: "densely-dashed")
  line(off("asp.atom-3", (0.15, 0)), off("his.atom-4", (-0.2, 0)), stroke: hb)
  line(off("his.atom-2", (0.1, -0.1)), off("ser.atom-0", (-0.28, 0)), stroke: hb)
  content(off("asp.atom-3", (0.62, 0.42)), text(size: 7pt, fill: rgb("#3A78C9"))[2.8 Å])

  // Ser oxygen attacks the substrate carbonyl carbon (red curly arrow).
  bezier(
    off("ser.atom-0", (0.15, 0.32)), off("sub.atom-1", (-0.2, -0.34)),
    (6.1, 1.25),
    mark: (end: ">", fill: rgb("#C0392B"), size: 0.16), stroke: 1pt + rgb("#C0392B"),
  )
  // C=O pi bond breaks onto oxygen.
  bezier(
    off("sub.bond-1-2", (0.2, 0.05)), off("sub.atom-2", (0.12, -0.18)),
    (8.5, 1.7),
    mark: (end: ">", fill: rgb("#C0392B"), size: 0.14), stroke: 1pt + rgb("#C0392B"),
  )

  // Residue labels.
  content((0, -2.1), text(size: 8pt)[Asp102])
  content((2.7, -2.1), text(size: 8pt)[His57])
  content((5.1, -2.1), text(size: 8pt)[Ser195])
  content((8.0, -0.1), text(size: 8pt, style: "italic")[substrate])
  content((2.5, 2.25), text(size: 8.5pt, fill: rgb("#8B2942"))[charge-relay system])
}))

#pagebreak()

= Ring closures after branches

Ring-closure digits written after branch groups still attach to the branch
point atom. The fused example should keep the five-membered carbonyl ring at
the tail.

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*Minimal post-branch closure* \ #text(size: 8pt, `C1=CCCC(=O)1`) \
   #smiles("C1=CCCC(=O)1", scale: 0.95)],
  [*Fused tail closure* \ #text(size: 8pt, `...C5=C4CCC(=O)5`) \
   #smiles("O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5", scale: 0.58)],
)

#pagebreak()

= Multiple alkene directional markers

Directional markers can appear on both substituent bonds from the same alkene
atom when those substituents are placed on opposite sides. Shared directional
bonds in a conjugated tail should not make the next double bond fail.

#align(center)[
  #text(size: 7.2pt, `CC1=C(C(=O)C[C@@H]1OC(=O)[C@@H]2[C@H](C2(C)C)/C=C(\C)/C(=O)OC)C/C=C\C=C`) \
  #smiles(
    "CC1=C(C(=O)C[C@@H]1OC(=O)[C@@H]2[C@H](C2(C)C)/C=C(\\C)/C(=O)OC)C/C=C\\C=C",
    scale: 0.62,
  )
]

#pagebreak()

= Catalytic cycles (`cycle`)

Species on a ring with arc arrows. `step(into:/out:)` adds reagents; `merge`
fuses the side arrow with the arc; `arc-gap` tunes how close arrows sit to the
species; `label-offset`/`into-offset`/`out-offset` nudge pieces like a `mol`
offset. Wilkinson's hydrogenation (default arrows):

#let cx(b) = box(inset: 2pt, b)

#align(center, cycle(
  scale: 1.0,
  mol(cx[RhCl(PPh#sub[3])#sub[3]]),
  step(label: [−PPh#sub[3], +S]),
  mol(cx[RhCl(PPh#sub[3])#sub[2]S]),
  step(label: [oxidative\ addition], into: ce("H2")),
  mol(cx[RhH#sub[2]Cl(PPh#sub[3])#sub[2]]),
  step(label: [alkene\ insertion], into: [alkene]),
  mol(cx[RhH(R)Cl(PPh#sub[3])#sub[2]]),
  step(label: [reductive\ elimination], out: [alkane]),
))

The same cycle with `merge: true` on every step and a tight `arc-gap: 0.0`,
so reagents flow into and out of the main arrows:

#align(center, cycle(
  scale: 1.0,
  arc-gap: 0.0,
  mol(cx[RhCl(PPh#sub[3])#sub[3]]),
  step(label: [−PPh#sub[3]]),
  mol(cx[RhCl(PPh#sub[3])#sub[2]S]),
  step(label: [ox. add.], into: ce("H2"), merge: true),
  mol(cx[RhH#sub[2]Cl(PPh#sub[3])#sub[2]]),
  step(label: [insertion], into: [alkene], merge: true, into-offset: (0.4, -0.3)),
  mol(cx[RhH(R)Cl(PPh#sub[3])#sub[2]]),
  step(label: [red. elim.], out: [#reaction(mol("CCO"), [+], mol("CCCO"), flow : "up")], merge: true, label-offset: (-0.2, 0)),
))

A two-species catalytic cycle where `step(out:)` grows a branch (a nested
`reaction()`) out of the released product:

#align(center, cycle(
  scale: 1.0,
  radius: 1.6,
  mol(cx[E]),
  step(label: [binding], into: [S]),
  mol(cx[E·S]),
  step(
    label: [catalysis],
    out: reaction(
      mol("CC=O", label: [product]),
      rxn-arrow(above: [\[O\]]),
      mol("CC(=O)O"),
    ),
  ),
))

#pagebreak()

= Spiro ring systems

Two rings joined at a single shared (spiro) atom each lay out as a regular
polygon opening away from the shared atom, instead of unravelling into a chain.

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,

  [*1,6-dioxaspiro[4.4]nonane* \ #text(size: 8pt, `CC[C@H](O1)CC[C@@]12CCCO2`) \
   #smiles("CC[C@H](O1)CC[C@@]12CCCO2", scale: 0.95)],
  [*Spiro[4.5] carbocycle* \ #text(size: 8pt, `C1CCC2(CC1)CCCCC2`) \
   #smiles("C1CCC2(CC1)CCCCC2", scale: 0.95)],
)

#pagebreak()

= Reaction flow direction (`flow`)

`reaction(flow: ...)` sets the writing direction. `"left"` and `"up"` reflect
the scheme so a branch that grows from the left or bottom of a cycle reads
naturally; `rxn-arrow(dir: auto)` follows the flow.

#align(center)[*right (default)*]
#align(center, reaction(
  mol("CCO"), rxn-arrow(above: [\[O\]]), mol("CC=O"), rxn-arrow(above: [\[O\]]), mol("CC(=O)O"),
))

#v(0.4cm)
#align(center)[*left*]
#align(center, reaction(
  flow: "left",
  mol("CCO"), rxn-arrow(above: [\[O\]]), mol("CC=O"), rxn-arrow(above: [\[O\]]), mol("CC(=O)O"),
))

#v(0.4cm)
#align(center)[*left with a reactant sum (`A + B -> C` reads right-to-left)*]
#align(center, reaction(
  flow: "left",
  mol("CC(=O)O"), [+], mol("CCO"), rxn-arrow(above: ce("H+")), mol("CCOC(=O)C"),
))

#v(0.4cm)
#align(center)[*up*]
#align(center, reaction(
  flow: "up",
  mol("CCO"), rxn-arrow(), mol("CC=O"), rxn-arrow(), mol("CC(=O)O"),
))

= Nested reaction, cycle, and branch

A main reaction embeds a `cycle`, whose `step(out:)` grows a sub-`reaction()`
branch out of a released product; the main reaction then continues to the
right. Here the branch flows downward into free vertical space.

#align(center, reaction(
  mol("C=C", label: [alkene]),
  rxn-arrow(above: [cat.]),
  cycle(
    scale: 0.8,
    radius: 1.7,
    mol(box(inset: 2pt)[LnM]),
    step(label: [binding], into: [S], merge: true, label-offset: (0.5, 0)),
    mol(box(inset: 2pt)[LnM·S]),
    step(
      label: [turnover],
      merge: true,
      label-offset: (-0.5, 0),
      out: reaction(
        flow: "down",
        mol("CCO", label: [product]),
        rxn-arrow(above: [\[O\]]),
        mol("CC(=O)O"),
      ),
    ),
  ),
  rxn-arrow(above: [workup]),
  mol("CCCC", label: [alkane]),
))

#pagebreak()

= Vertical flow and cycle branch anchors

Vertical `flow` stacks ordinary reaction items even when there are no
`rxn-arrow()` separators:

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*down* \
   #reaction(flow: "down", mol("CCO"), [+], mol("CC=O"))],
  [*up* \
   #reaction(flow: "up", mol("CCO"), [+], mol("CC=O"))],
)

An outgoing cycle branch attaches at the released product side of the branch,
so the side arrow lands on the first downstream species rather than the middle
of the whole nested reaction:

#align(center, cycle(
  scale: 0.9,
  radius: 1.7,
  start: 0deg,
  mol(box(inset: 2pt)[M–R]),
  step(
    label: [release],
    merge: true,
    out: reaction(
      flow: "down",
      mol("CCO", label: [first]),
      rxn-arrow(above: [\[O\]]),
      mol("CC=O"),
    ),
  ),
  mol(box(inset: 2pt)[M]),
))

#pagebreak()

= Incoming cycle branch anchors

An incoming cycle branch attaches at the downstream side of the branch, so the
side arrow starts at the last upstream species rather than the middle of the
whole nested reaction:

#align(center, cycle(
  scale: 0.9,
  radius: 1.7,
  start: 0deg,
  mol(box(inset: 2pt)[M]),
  step(
    label: [entry],
    merge: true,
    into: reaction(
      flow: "up",
      mol("CC=O", label: [first]),
      rxn-arrow(above: [red.]),
      mol("CCO", label: [last]),
    ),
  ),
  mol(box(inset: 2pt)[M-R]),
))

#pagebreak()

= Cycle label rotation

`step(rotation:)` rotates cycle labels. `"auto"` follows the step's circle angle and
flips labels on the far side so they remain readable:

#align(center, cycle(
  scale: 0.95,
  radius: 2.4,
  start: 90deg,
  mol(box(inset: 2pt)[A]),
  step(label: [auto], rotation: "auto"),
  mol(box(inset: 2pt)[B]),
  step(label: [90deg], rotation: 90deg),
  mol(box(inset: 2pt)[C]),
  step(label: [auto], rotation: "auto"),
  mol(box(inset: 2pt)[D]),
  step(label: [straight], rotation: "straight"),
))

#pagebreak()

= Wilkinson-style cycle with vertical entry

A precursor is written first in a downward reaction flow; the catalytic cycle
then starts from the top species and continues clockwise.


#align(center)[
  #reaction(
    flow : "down",
    mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)(!w{>PPh3})"), offset: (0.5, 0)),
    rxn-arrow(above : [#ce("PPh3") solvent (S)]),
    cycle(
      radius : 6.0,
      mol(smiles("{Rh}(!w{S | S})(!h{>PPh3})(!hCl)(!w{>PPh3})")),
      step(label : [#text(size : 6pt)[oxidative addition]], label-offset : (0.3, 0.5), into : [#ce("H2")], rotation : "auto", merge : true),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({S | S})({H})(!w{H})")),
      step(into : [#smiles("C=C")], merge : true),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({})({H})({H})")),
      step(label : [#text(size : 6pt)[migratory insertion]], label-offset: (0, 0.5), rotation : "auto"),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({S | S})({H})(!wCC{H})", rotation : -90deg)),
      step(label : [#text(size : 6pt)[reductive elimination]], label-offset:  (-0.5, 0.5), out : [#smiles("{H}CC{H}")], merge:true, rotation : "auto")
    )
  )
]

#smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({S | S})({H})(!wCC{H})")
#smiles("CS")

#pagebreak()

= Arrow axes, stereobond length, and cycle radius

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*horizontal arrow axis* \
   #reaction(
     mol([A]),
     rxn-arrow(above: [very long reagent label above the arrow], below: [x]),
     mol([B]),
   )],
  [*vertical arrow axis* \
   #reaction(
     flow: "down",
     mol([A]),
     rxn-arrow(above: [right], below: [much longer left-side reagent label]),
     mol([B]),
   )],
)

#grid(
  columns: (1fr,),
  gutter: 1.5em,
  align: center + horizon,
  [*wedge/hash span* \
   #smiles("{Rh}(!w{S | S})(!h{>PPh3})(!wCl)(C)", bond-length: 1.2)],
)

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*cycle radius 1.4* \
   #cycle(
     radius: 1.4,
     scale: 0.65,
     mol([A]), step(), mol([B]), step(), mol([C]),
   )],
  [*cycle radius 2.8* \
   #cycle(
     radius: 2.8,
     scale: 0.65,
     mol([A]), step(), mol([B]), step(), mol([C]),
   )],
)

#pagebreak()

= Rotated mirror axes and forced wedge tips

`mirror: "vertical"` preserves page left/right even when the molecule is
rotated. Explicit `!w`/`!h` bonds use the written source atom as the thin end.

#let rh-alkyl = "{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({S | S})({H})(!wCC{H})"

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*rotated* \
   #smiles(rh-alkyl, rotation: -90deg, bond-length: 1.15)],
  [*rotated + vertical mirror* \
   #smiles(rh-alkyl, rotation: -90deg, mirror: "vertical", bond-length: 1.15)],
)

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*source at Rh* \
   #smiles("{Rh}(!wC)(!hC)(!wCl)(!h{>PPh3})", bond-length: 1.2)],
  [*source at carbon* \
   #smiles("C(!w{Rh})(!h{Rh})", bond-length: 1.2)],
)

#pagebreak()

= Scheme offsets stay page-relative

`mol(offset:)` in a normal scheme does not switch to mechanism layout. Positive
x moves the item right on the page even when the reaction flow is vertical or
reflected.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*right flow* \
   #reaction(
     mol([A], offset: (0.6, 0)),
     rxn-arrow(),
     mol([B]),
   )],
  [*down flow* \
   #reaction(
     flow: "down",
     mol([A], offset: (0.6, 0)),
     rxn-arrow(),
     mol([B]),
   )],
  [*left flow* \
   #reaction(
     flow: "left",
     mol([A], offset: (0.6, 0)),
     rxn-arrow(),
     mol([B]),
   )],
)

#pagebreak()

= Custom-label scripts and bracketed offsets

`_(...)` and `^(...)` script custom labels. A subscript and superscript after
the same glyph share that glyph, so the plus in `NH_4^+` belongs to H, not 4.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*Phosphine subscript* \
   #smiles("{>PPh_(3)}C=O")],
  [*Subscript and charge* \
   #smiles("{NH_4^+}") \
   #ce("NH4+")],
  [*Grouped scripts* \
   #smiles("{SO_(4)^(2-)}") \
   #ce("SO4^2-")],
)

#pagebreak()

= Custom-label script alignment with #raw("ce()")

Scripted custom labels use the same script size, spacing, and vertical
attachment as chemical formulas.

#grid(
  columns: (1fr, 1fr),
  gutter: 2em,
  row-gutter: 1em,
  align: center + horizon,
  [*Custom label*], [*Chemical formula*],
  [#smiles("{NH_4^+}", font-size: 18pt)], [#ce("NH4+", font-size: 18pt)],
  [#smiles("{SO_(4)^(2-)}", font-size: 18pt)], [#ce("SO4^2-", font-size: 18pt)],
)

#pagebreak()

= Reaction-arrow scale

`rxn-arrow(scale:)` uniformly resizes the shaft, arrowhead, spacing, and
condition labels in horizontal and vertical schemes.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center + horizon,
  [*Scale 0.7* \
   #reaction(mol([A]), rxn-arrow(scale: 0.7, above: [cat.]), mol([B]))],
  [*Scale 1.0* \
   #reaction(mol([A]), rxn-arrow(above: [cat.]), mol([B]))],
  [*Scale 1.5* \
   #reaction(mol([A]), rxn-arrow(scale: 1.5, above: [cat.]), mol([B]))],
)

#grid(
  columns: (1fr, 1fr),
  gutter: 2em,
  align: center + horizon,
  [*Vertical scale 0.7* \
   #reaction(flow: "down", mol([A]), rxn-arrow(scale: 0.7, above: [cat.]), mol([B]))],
  [*Vertical scale 1.5* \
   #reaction(flow: "down", mol([A]), rxn-arrow(scale: 1.5, above: [cat.]), mol([B]))],
)
