#import "../src/lib.typ": smiles, atom, arrow
#set page(width: 24cm, height: auto, margin: 1cm)
#set text(font: "New Computer Modern", size: 11pt)
= Arrowheads clear inline hydrogen labels
#grid(columns: (1fr, 1fr), gutter: 20pt,
  [*OH from the H side* \ #smiles("CCO", scale: 1.6, arrow(from: atom(0, offset: (4, 0)), to: atom(2), bend: none))],
  [*NH₂ from the H side* \ #smiles("CCN", scale: 1.6, arrow(from: atom(0, offset: (4, 0)), to: atom(2), bend: none))],
  [*Reversed OH* \ #smiles("CCO", rotation: 180deg, scale: 1.6, arrow(from: atom(0, offset: (-4, 0)), to: atom(2), bend: none))],
  [*Rotated OH* \ #smiles("CCO", rotation: 90deg, scale: 1.6, arrow(from: atom(0, offset: (0, 4)), to: atom(2), bend: none))],
  [*Direct NH₂ hydrogen reference* \ #smiles("CCN", scale: 1.6, arrow(from: atom(0, offset: (4, 0)), to: atom(3), bend: none))],
)
