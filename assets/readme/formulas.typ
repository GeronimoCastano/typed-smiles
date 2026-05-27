#import "../../src/lib.typ": ce

#set page(width: 16cm, height: 4.7cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#table(
  columns: (1fr, 1fr),
  gutter: 0em,
  row-gutter: 0em,
  align: center + horizon,
  stroke: 0.4pt + rgb("#d8d8d8"),

  [#stack(spacing: 0.35cm, strong[Formula], ce("H2SO4"))],
  [#stack(spacing: 0.35cm, strong[Ions], ce("(NH4)2SO4"))],
  [#stack(spacing: 0.35cm, strong[Combustion], ce("CH4 + 2O2 -> CO2 + 2H2O"))],
  [#stack(spacing: 0.35cm, strong[Equilibrium], ce("N2 + 3H2 <=> 2NH3"))],
)
