#import "../../src/lib.typ": smiles

#set page(width: 18cm, height: 7.2cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#grid(
  columns: (1fr, 1fr),
  gutter: 0.8cm,
  align: center + horizon,

  [*Automatic zigzag*],
  [*Curled with `!c`*],

  [#smiles("CC(C)=CC([O-])C({>PPh_3^+ | P})C(=O)OCC", scale: 0.82)],
  [#smiles("CC(C)=CC([O-])C({>PPh_3^+ | P})!cC(=O)OCC", scale: 0.82)],
)
