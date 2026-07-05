#import "../../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction, cycle, step

#set page(width: 19cm, height: 20cm, margin: 0.5cm)
#set text(font: "New Computer Modern", size: 10pt)

#align(center)[
  #reaction(
    flow : "down",
    mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)(!w{>PPh3})"), offset: (0.5, 0)),
    rxn-arrow(above : [#ce("PPh3") solvent (S)]),
    cycle(
      radius : 6.0,
      mol(smiles("{Rh}(!w{S | S})(!h{>PPh3})(!hCl)(!w{>PPh3})")),
      step(label : [#text(size : 6pt)[oxidative addition]], label-offset : (0.3, 0.5), into : [#ce("H2")], rotation : "auto", merge : true),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({S | S})({H})(!w{H})")),
      step(into : [#smiles("C=C")], merge : true),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({})({H})({H})")),
      step(label : [#text(size : 6pt)[migratory insertion]], label-offset: (0, 0.5), rotation : "auto"),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({S | S})({H})(!wCC{H})", rotation : -90deg)),
      step(label : [#text(size : 6pt)[reductive elimination]], label-offset:  (-0.5, 0.5), out : [#smiles("{H}CC{H}")], merge:true, rotation : "auto")
    )
  )
]
