#import "../../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction

#set page(width: 18cm, height: 10.2cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#stack(
  spacing: 1cm,
  stack(
    spacing: 0.4cm,
    align(center, strong[Fischer esterification]),
    align(center, [
      #reaction(
        mol(smiles("CC(=O)O"), label: text(size: 8pt)[acetic acid]),
        [+],
        mol(smiles("CCO"), label: text(size: 8pt)[ethanol]),
        rxn-arrow(above: ce("H+"), below: [heat]),
        mol(smiles("CCOC(=O)C"), label: text(size: 8pt)[ethyl acetate]),
        [+],
        ce("H2O"),
      )
    ]),
  ),
  stack(
    spacing: 0.4cm,
    align(center, strong[Electrophilic aromatic bromination]),
    align(center, [
      #reaction(
        mol(smiles("C1=CC=CC=C1"), label: text(size: 8pt)[benzene]),
        rxn-arrow(above: ce("Br2"), below: ce("FeBr3")),
        mol(smiles("BrC1=CC=CC=C1"), label: text(size: 8pt)[bromobenzene]),
      )
    ]),
  ),
)
