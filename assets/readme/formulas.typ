#import "../../src/lib.typ": ce

#set page(width: 16cm, height: 7.6cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#let example(title, body) = [
  #stack(
    spacing: 0.5cm,
    align(center, strong(title)),
    align(center, body),
  )
]

#grid(
  columns: (1fr, 1fr),
  gutter: 1.2em,
  row-gutter: 1.8cm,
  align: center,

  example("Formula", ce("H2SO4")),
  example("Ions", ce("(NH4)2SO4")),
  example("Combustion", ce("CH4 + 2O2 -> CO2 + 2H2O")),
  example("Equilibrium", ce("N2 + 3H2 <=> 2NH3")),
)
