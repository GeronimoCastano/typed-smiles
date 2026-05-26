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
  rxn-arrow(dir: "down", above: [HNO#sub[3], H#sub[2]SO#sub[4]]),
  mol(smiles("BrC1=CC(=CC=C1)[N+](=O)[O-]"), label: text(size: 8pt)[*B*]),
  rxn-arrow(dir: "left", above: [Fe, HCl]),
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
rendered as italic text. Standard SMILES bonding rules apply around it.

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
