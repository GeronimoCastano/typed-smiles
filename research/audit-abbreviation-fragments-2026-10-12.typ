#import "../src/lib.typ": smiles, atom, highlight, arrow
#set page(width: 22cm, height: auto, margin: 12pt)
= Abbreviation fragments and retained indices
#grid(columns: (1fr, 1fr), gutter: 18pt,
  [*Expanded* \ #smiles("CC(=O)O", show-indices: true, scale: 1.6)],
  [*Abbreviated* \ #smiles("CC(=O)O", abbreviate: "CO2H", show-indices: true, scale: 1.6)],
  [*O₂ fragment* \ #smiles("CC(=O)O", abbreviate: "CO2H", show-indices: true, scale: 1.6, highlight(atom(2)))],
  [*Hydrogen fragment* \ #smiles("CC(=O)O", abbreviate: "CO2H", show-indices: true, scale: 1.6, highlight(atom(4)))],
  [*Reversed* \ #smiles("CC(=O)O", abbreviate: "CO2H", show-indices: true, rotation: 180deg, scale: 1.6, highlight(atom(2)))],
  [*Original S22* \ #smiles("OC(=O)c1ccc(O)cc1", abbreviate: "CO2H", show-indices: true, highlight(atom(10)))],
  [*CF₃* \ #smiles("CC(F)(F)F", abbreviate: "CF3", show-indices: true, highlight(atom(2)))],
  [*NO₂* \ #smiles("C[N+](=O)[O-]", abbreviate: "NO2", show-indices: true, highlight(atom(2)))],
  [*Arrow to O₂, highlight C* \ #smiles("CC(=O)O", abbreviate: "CO2H", arrow(from: atom(0), to: atom(2)), highlight(atom(1)))],
)
