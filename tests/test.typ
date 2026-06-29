#import "../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction, atom, bond, lp, species, arrow, highlight, brackets

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
  rxn-arrow(above: ce("Br2") + [, FeBr#sub[3]]),
  mol(smiles("BrC1=CC=CC=C1"), label: text(size: 8pt)[Bromobenzene]),
)

#v(1.5em)
*Fischer esterification* (horizontal, with + separator):

#reaction(
  mol(smiles("CC(=O)O"), label: text(size: 8pt)[Acetic acid]),
  [+],
  mol(smiles("CCO"), label: text(size: 8pt)[Ethanol]),
  rxn-arrow(above: [H#super[+]], below: [Δ]),
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
  rxn-arrow(kind: "equilibrium", above: [H#super[+]], below: [heat]),
  mol(smiles("CCOC(=O)C"), label: text(size: 8pt)[ester]),
  [+],
  ce("H2O"),
)

#v(1em)
#reaction(
  mol(smiles("N#N"), label: text(size: 8pt)[N#sub[2]]),
  rxn-arrow(kind: "equilibrium", dir: "left", above: [pressure]),
  mol(smiles("[H][H]"), label: text(size: 8pt)[H#sub[2]]),
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
  rxn-arrow(kind: "equilibrium-filled", above: [H#super[+]], below: [heat]),
  mol(smiles("CCOC(=O)C"), label: text(size: 8pt)[ester]),
  [+],
  ce("H2O"),
)

#v(1em)
#reaction(
  mol(smiles("N#N"), label: text(size: 8pt)[N#sub[2]]),
  rxn-arrow(kind: "equilibrium-filled", dir: "left", above: [pressure]),
  mol(smiles("[H][H]"), label: text(size: 8pt)[H#sub[2]]),
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
  [*Vertical NH#sub[2] (down)*],
  [*Bicolor wedge to OH*],
  [*Bicolor hash to NH#sub[2]*],

  [#text(size: 8pt, [`CO` (vertical)]) \ #smiles("CO", bond-length: 1.1, rotation: 120deg)],
  [#text(size: 8pt, [`CN` (steep)]) \ #smiles("CN", bond-length: 1.1, rotation: -60deg)],
  [#text(size: 8pt, `C!wO`) \ #smiles("C!wO", bond-length: 1.3, rotation: 60deg)],
  [#text(size: 8pt, `C!hN`) \ #smiles("C!hN", bond-length: 1.3, rotation: -60deg)],
)

#v(1em)
*Regression molecule with NH#sub[2] and NH labels by default:*

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
  [*All H* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2, show-all-h: true)],
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

  [*NH#sub[2], right* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.2, lone-pairs: "dots")],
  [*NH#sub[2], up* \ #text(size: 8pt, `CCN rotation: 90deg`) \ #smiles("CCN", bond-length: 1.2, rotation: 90deg, lone-pairs: "dots")],
  [*NH#sub[2], left* \ #text(size: 8pt, `CCN rotation: 180deg`) \ #smiles("CCN", bond-length: 1.2, rotation: 180deg, lone-pairs: "dots")],
  [*NH#sub[2], down* \ #text(size: 8pt, `CCN rotation: -90deg`) \ #smiles("CCN", bond-length: 1.2, rotation: -90deg, lone-pairs: "dots")],
)

#v(1em)
*Terminal lone pairs with large labels:*

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Large OH* \ #text(size: 8pt, `font-size: 22pt`) \ #smiles("CCO", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots")],
  [*Large NH#sub[2]* \ #text(size: 8pt, `font-size: 22pt`) \ #smiles("CCN", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots")],
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

  [*NH#sub[2] right* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.15, lone-pairs: "dots", show-indices: true)],
  [*NH#sub[2] down* \ #text(size: 8pt, `CCN rotation: -90deg`) \ #smiles("CCN", bond-length: 1.15, rotation: -90deg, lone-pairs: "dots", show-indices: true)],
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
  [*Large terminal NH#sub[2]* \ #text(size: 8pt, `C[NH2], font-size: 22pt`) \ #smiles("C[NH2]", font-size: 22pt, bond-length: 1.2, lone-pairs: "dots", show-indices: true)],

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
  [*SMILES NH#sub[2]* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.15)],
  [*Literal N* \ #text(size: 8pt, `CC{N|N}`) \ #smiles("CC{N|N}", bond-length: 1.15)],

  [*SMILES NH* \ #text(size: 8pt, `CNC`) \ #smiles("CNC", bond-length: 1.15)],
  [*Literal N red* \ #text(size: 8pt, `C{N|red}C`) \ #smiles("C{N|red}C", bond-length: 1.15)],
  [*SMILES SH* \ #text(size: 8pt, `CCS`) \ #smiles("CCS", bond-length: 1.15)],
  [*Literal S* \ #text(size: 8pt, `CC{S|S}`) \ #smiles("CC{S|S}", bond-length: 1.15)],

  [*SMILES PH#sub[2]* \ #text(size: 8pt, `CP`) \ #smiles("CP", bond-length: 1.15)],
  [*Literal P* \ #text(size: 8pt, `C{P|P}`) \ #smiles("C{P|P}", bond-length: 1.15)],
  [*SMILES BH#sub[2]* \ #text(size: 8pt, `CB`) \ #smiles("CB", bond-length: 1.15)],
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

  [#text(size: 8pt, `{PPh3|P}C=O`) \
   #smiles("{PPh3|P}C=O")],

  [#text(size: 8pt, `{LG|red}C=O`) \
   #smiles("{LG|red}C=O")],

  [#text(size: 8pt, `[N] vs {N}`) \
   #smiles("[N]") #h(1em) #smiles("{N}")],
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

  [*Terminal explicit NH#sub[2]* \
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
