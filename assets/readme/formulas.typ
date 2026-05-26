#import "../../src/lib.typ": ce

#set page(width: 16cm, height: 5cm, margin: 0.8cm)
#set text(font: "New Computer Modern", size: 10pt)

#let example(title, code, body) = [
  #stack(
    spacing: 5pt,
    align(center, strong(title)),
    align(center, text(size: 7pt, raw(code))),
    align(center, body),
  )
]

#grid(
  columns: (1fr, 1fr),
  gutter: 1.2em,
  row-gutter: 1em,
  align: center,

  example("Formula", "H2SO4", ce("H2SO4")),
  example("Ions", "(NH4)2SO4", ce("(NH4)2SO4")),
  example("Combustion", "CH4 + 2O2 -> CO2 + 2H2O", ce("CH4 + 2O2 -> CO2 + 2H2O")),
  example("Equilibrium", "N2 + 3H2 <=> 2NH3", ce("N2 + 3H2 <=> 2NH3")),
)
