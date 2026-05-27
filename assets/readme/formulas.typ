#import "../../src/lib.typ": ce

#set page(width: 16cm, height: 8.2cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr),
  gutter: 1.2em,
  row-gutter: 0.8cm,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [*Formula*],
  [*Ions*],
  [#ce("H2SO4")],
  [#ce("(NH4)2SO4")],

  [*Combustion*],
  [*Equilibrium*],
  [#ce("CH4 + 2O2 -> CO2 + 2H2O")],
  [#ce("N2 + 3H2 <=> 2NH3")],
)
