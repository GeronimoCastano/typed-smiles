#import "../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction, display-smiles

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
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*Alcohol OH* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2)],
  [*Carboxylic OH* \ #text(size: 8pt, `CC(=O)O`) \ #smiles("CC(=O)O", bond-length: 1.2)],
  [*Thiol SH* \ #text(size: 8pt, `CCS`) \ #smiles("CCS", bond-length: 1.2)],

  [*Primary amine* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.2)],
  [*Secondary amine* \ #text(size: 8pt, `CNC`) \ #smiles("CNC", bond-length: 1.2)],
  [*Tertiary amine* \ #text(size: 8pt, `CN(C)C`) \ #smiles("CN(C)C", bond-length: 1.2)],
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
  [*No hetero H* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2, show-hetero-h: false)],
  [*All H* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.2, show-all-h: true)],
  [*Explicit H* \ #text(size: 8pt, `[NH4+]`) \ #smiles("[NH4+]", bond-length: 1.2)],
)

#v(1em)
*Heteroatom H on/off comparisons:*

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  row-gutter: 1.5em,
  align: center,

  [*OH shown* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.15)],
  [*OH hidden* \ #text(size: 8pt, `CCO`) \ #smiles("CCO", bond-length: 1.15, show-hetero-h: false)],
  [*NH#sub[2] shown* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.15)],
  [*NH#sub[2] hidden* \ #text(size: 8pt, `CCN`) \ #smiles("CCN", bond-length: 1.15, show-hetero-h: false)],

  [*NH shown* \ #text(size: 8pt, `CNC`) \ #smiles("CNC", bond-length: 1.15)],
  [*NH hidden* \ #text(size: 8pt, `CNC`) \ #smiles("CNC", bond-length: 1.15, show-hetero-h: false)],
  [*SH shown* \ #text(size: 8pt, `CCS`) \ #smiles("CCS", bond-length: 1.15)],
  [*SH hidden* \ #text(size: 8pt, `CCS`) \ #smiles("CCS", bond-length: 1.15, show-hetero-h: false)],

  [*PH#sub[2] shown* \ #text(size: 8pt, `CP`) \ #smiles("CP", bond-length: 1.15)],
  [*PH#sub[2] hidden* \ #text(size: 8pt, `CP`) \ #smiles("CP", bond-length: 1.15, show-hetero-h: false)],
  [*BH#sub[2] shown* \ #text(size: 8pt, `CB`) \ #smiles("CB", bond-length: 1.15)],
  [*BH#sub[2] hidden* \ #text(size: 8pt, `CB`) \ #smiles("CB", bond-length: 1.15, show-hetero-h: false)],
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
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Linear C=O* \ #text(size: 8pt, `CC=O`) \ #smiles("CC=O", bond-length: 1.4)],
  [*Branched C=O* \ #text(size: 8pt, `CC(=O)C`) \ #smiles("CC(=O)C", bond-length: 1.4)],
  [*P-N triple branch* \ #text(size: 8pt, `CP(#N)C`) \ #smiles("CP(#N)C", bond-length: 1.4)],
  [*P-N triple chain* \ #text(size: 8pt, `CP#N`) \ #smiles("CP#N", bond-length: 1.4)],
)

= Wedge / Dash Bonds

Syntax: use #raw("/") in SMILES for a solid wedge (toward viewer) and
#raw("\\") for a hashed wedge (away from viewer).

#grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Wedge up* \ #text(size: 8pt, `C/N`) \ #smiles("C/N")],
  [*Wedge down* \ #text(size: 8pt, `C\N`) \ #smiles("C\\N")],
  [*Mixed* \ #text(size: 8pt, `F/C\Cl`) \ #smiles("F/C\\Cl")],
  [*L-Alanine* \ #text(size: 8pt, `N[C@@H](C)C(=O)O`) \
   #smiles("N/C(C)C(=O)O")],
)

#v(1em)
*E/Z double bonds* (the directional bonds flank the #raw("="):

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center,
  [*trans (E)* \ #text(size: 8pt, `F/C=C/F`) \ #smiles("F/C=C/F")],
  [*cis (Z)* \ #text(size: 8pt, `F/C=C\F`) \ #smiles("F/C=C\\F")],
)

= Abbreviated Groups

Syntax: use #raw("{label}") anywhere an atom would appear. The label is
rendered as upright text. Standard SMILES bonding rules apply around it.

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Simple abbrev* \ #text(size: 8pt, `{OEt}C=O`) \
   #smiles("{OEt}C=O")],

  [*Branch abbrev* \ #text(size: 8pt, `C({PPh3})=O`) \
   #smiles("C({PPh3})=O")],

  [*Multiple* \ #text(size: 8pt, `{L}C(=O){NHR}`) \
   #smiles("{L}C(=O){NHR}")],
)

#v(1em)
*Wedge bonds + abbreviated groups together* (transition-state style):

#grid(
  columns: (1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Wilkinson's-like* \ #text(size: 8pt, `{PPh3}C({PPh3})=O`) \
   #smiles("{PPh3}C({PPh3})=O")],

  [*With stereo* \ #text(size: 8pt, `{Nu}/C({LG})=O`) \
   #smiles("{Nu}/C({LG})=O")],
)

// #pagebreak()


// #align(center)[
//   #smiles("CC(=O)C", scale : 1)
// ]
