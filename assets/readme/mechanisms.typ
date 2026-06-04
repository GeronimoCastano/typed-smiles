#import "../../src/lib.typ": smiles, ce, mol, reaction, atom, bond, lp, species, arrow, highlight, brackets

#set page(width: 16cm, height: 7.2cm, margin: 0.7cm)
#set text(font: "New Computer Modern", size: 10pt)

#grid(
  columns: (1.15fr, 0.85fr),
  column-gutter: 0.8cm,
  align: center + horizon,

  // ── Inter-species mechanism: hydroxide attacks an alkyl halide (SN2) ──
  stack(
    spacing: 0.4cm,
    strong[Nucleophilic substitution],
    reaction(
      mol("[OH-]", lone-pairs: "dots"),
      mol("C(Br)(C)C", offset: (1.5, 0.4)),
      arrow(
        from: lp(0, 0), to: atom(1, 0),
        bend: "left", color: red, label: text(size: 8pt)[attack],
      ),
      arrow(
        from: bond(1, 0, 1), to: atom(1, 1, offset: (0.9, 0)),
        bend: "left", color: red,
      ),
    ),
  ),

  // ── Intramolecular arrows + highlight, with a transition-state bracket ──
  stack(
    spacing: 0.4cm,
    strong[Highlight + bracket],
    brackets(
      smiles(
        "CC(=O)C",
        lone-pairs: "dots",
        highlight(bond(1, 2), fill: rgb("#FFE45C")),
        arrow(from: bond(1, 2), to: atom(2), bend: "right", color: red),
      ),
      sup: [‡],
    ),
  ),
)
