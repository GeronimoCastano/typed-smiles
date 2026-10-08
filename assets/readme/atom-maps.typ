#import "../../src/lib.typ": smiles, mol, reaction, rxn-arrow

#set page(width: 18cm, height: auto, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#align(center, reaction(
  mol("[CH3:1][C:2](=[O:3])[OH:4]", show-maps: true, scale: 1.2),
  [+],
  mol("[CH3:5][OH:6]", show-maps: true, scale: 1.2),
  rxn-arrow(above: [H#super[+]]),
  mol("[CH3:1][C:2](=[O:3])[O:6][CH3:5]", show-maps: true, scale: 1.2),
  [+],
  mol("[OH2:4]", show-maps: true, scale: 1.2),
))
