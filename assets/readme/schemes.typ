#import "../../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction

#set page(width: 18cm, height: 13.4cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#stack(
  spacing: 1.2em,
  align(center, strong[Bromination, nitration, and reduction sequence]),
  align(center, reaction(
    mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[1]),
    rxn-arrow(above: ce("Br2"), below: ce("FeBr3")),
    mol(smiles("BrC1=CC=CC=C1"), label: text(size: 8pt)[A]),
    rxn-arrow(dir: "down", above: ce("HNO3"), below: ce("H2SO4")),
    mol(smiles("BrC1=CC(=CC=C1)[N+](=O)[O-]"), label: text(size: 8pt)[B]),
    rxn-arrow(dir: "left", above: ce("Fe"), below: ce("HCl")),
    mol(smiles("BrC1=CC(=CC=C1)N"), label: text(size: 8pt)[C]),
  )),
)
