#import "../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction

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


#pagebreak()

#align(center)[
  #smiles("R({H})({H})CC{Rh}({Ph3P})({Ph3P})(Cl)")
]


#pagebreak()

= Branch atoms with 3+ substituents (no overlap)

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 1.5em,
  align: center,

  [*Aminodichloromethane* \ #text(size: 8pt, `C(Cl)(Cl)N`) \ #smiles("C(Cl)(Cl)N")],
  [*Chloroform* \ #text(size: 8pt, `C(Cl)(Cl)Cl`) \ #smiles("C(Cl)(Cl)Cl")],
  [*Neopentane* \ #text(size: 8pt, `CC(C)(C)C`) \ #smiles("CC(C)(C)C")],
)
