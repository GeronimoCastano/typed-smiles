#import "../../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction

#set page(width: 18cm, height: 11cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#align(center, reaction(
  mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[1]),
  rxn-arrow(above: ce("Br2"), below: ce("FeBr3")),
  mol(smiles("BrC1=CC=CC=C1"), label: text(size: 8pt)[A]),
  rxn-arrow(dir: "down", above: [HNO#sub[3]], below: [H#sub[2]SO#sub[4]]),
  mol(smiles("BrC1=CC(=CC=C1)[N+](=O)[O-]"), label: text(size: 8pt)[B]),
  rxn-arrow(dir: "left", above: [Fe], below: [HCl]),
  mol(smiles("BrC1=CC(=CC=C1)N"), label: text(size: 8pt)[C]),
))
