// Dark card with neon arrows and reagents: benzene -> bromobenzene -> 1-bromo-2-nitrobenzene.
// Requires the Iosevka font to be installed. Labels, reagents, captions, and
// atom labels all use it; without it Typst falls back to another sans font.

#import "../../../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction

#let neon = rgb("#59DECD")
#let card-fill = rgb("#0E1116")
#let ink = rgb("#E6E8EB")

#set page(width: 17cm, height: auto, margin: 12pt, fill: none)
#set text(font: "Iosevka", fill: ink)

#let dark-molecule(smiles-str, name) = mol(
  smiles(
    smiles-str,
    fg: ink,
    theme: "dark",
    font: "Iosevka",
  ),
  label: text(size: 8pt, fill : ink.transparentize(40%))[#name],
)

#block(
  width: 100%,
  fill: card-fill,
  radius: 10pt,
  inset: 22pt,
  stroke: 0.8pt + ink.transparentize(80%),
)[
  #align(center, reaction(
    dark-molecule("C1=CC=CC=C1", "benzene"),
    rxn-arrow(
      above: text(fill: neon, ce("Br2", font: "Iosevka")),
      below: text(fill: neon, ce("FeBr3", font: "Iosevka")),
      color: neon,
    ),
    dark-molecule("BrC1=CC=CC=C1", "bromobenzene"),
    rxn-arrow(
      above: text(fill: neon, ce("HNO3", font: "Iosevka")),
      below: text(fill: neon, ce("H2SO4", font: "Iosevka")),
      color: neon,
    ),
    dark-molecule("BrC1=CC=CC=C1[N+](=O)[O-]", "1-bromo-2-nitrobenzene"),
  ))
]
