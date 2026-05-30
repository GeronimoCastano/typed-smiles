// typed-smiles documentation
//
// Compile with:
//   typst compile --root . docs/documentation.typ docs/documentation.pdf

#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#import "../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction

#let version = "0.2.0"
#let accent = rgb("#239dad")
#let accent-soft = rgb("#e7f4f6")

// ── Theme ────────────────────────────────────────────────────────────────────

#set document(title: "typed-smiles User Guide", author: "Geronimo Castaño")
#set text(font: "New Computer Modern", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 1.1em)
#set heading(numbering: "1.1")
#show link: set text(fill: accent)
#show ref: set text(fill: accent)

// Code blocks never split across a page.
#show raw.where(block: true): set block(breakable: false)

#show heading.where(level: 1): it => {
  pagebreak(weak: true)
  block(width: 100%, breakable: false, sticky: true, {
    set text(fill: accent, size: 18pt, weight: "bold")
    if it.numbering != none [#counter(heading).display() #h(0.5em)]
    it.body
    v(-0.2em)
    line(length: 100%, stroke: 0.6pt + accent)
  })
}
#show heading.where(level: 2): it => block(above: 1.3em, below: 0.7em, sticky: true, {
  set text(fill: accent.darken(8%), size: 13pt, weight: "bold")
  if it.numbering != none [#counter(heading).display() #h(0.4em)]
  it.body
})
#show heading.where(level: 3): it => block(above: 1.1em, below: 0.6em, sticky: true, {
  set text(fill: accent.darken(18%), size: 11pt, weight: "bold")
  it.body
})

// ── codly ────────────────────────────────────────────────────────────────────

#show: codly-init.with()
#codly(
  languages: codly-languages,
  zebra-fill: none,
  inset: (x: 0.5em, y: 0.32em),
  radius: 0pt,
  stroke: none,
)

// ── Example helper: source on the left, live render on the right ─────────────

#let pkg-scope = (smiles: smiles, ce: ce, rxn-arrow: rxn-arrow, mol: mol, reaction: reaction)

#let example(body, side: true) = block(
  width: 100%,
  stroke: 0.6pt + luma(210),
  radius: 5pt,
  clip: true,
  breakable: false,
  {
    let rendered = block(
      width: 100%, fill: white, inset: (x: 12pt, y: 14pt),
      align(center + horizon, eval(body.text, mode: "markup", scope: pkg-scope)),
    )
    if side {
      grid(
        columns: (1.05fr, 0.95fr), column-gutter: 0pt,
        block(width: 100%, fill: luma(248), inset: (y: 4pt), body),
        block(width: 100%, stroke: (left: 0.6pt + luma(220)), rendered),
      )
    } else {
      block(width: 100%, fill: luma(248), inset: (y: 4pt), body)
      block(width: 100%, stroke: (top: 0.6pt + luma(220)), rendered)
    }
  },
)

// Keep a subsection's prose and example together on one page.
#let demo(body) = block(breakable: false, width: 100%, body)

#let callout(clr, label, body) = block(
  width: 100%, fill: clr.lighten(90%), stroke: (left: 2.5pt + clr),
  radius: (right: 4pt), inset: (x: 12pt, y: 9pt), above: 1.1em, below: 1.1em,
  { text(fill: clr, weight: "bold")[#label.#h(0.5em)]; body },
)
#let note(body) = callout(accent, "Note", body)
#let warn(body) = callout(rgb("#b54708"), "Note", body)

#let c(it) = raw(it, lang: none)

// ── Cover ────────────────────────────────────────────────────────────────────

#set page(paper: "a4", margin: (x: 2.2cm, top: 2.6cm, bottom: 2.4cm), header: none, footer: none)

#v(1.2fr)
#align(center)[
  #text(size: 40pt, weight: "bold", fill: accent)[typed-smiles]
  #v(0.2em)
  #text(size: 15pt, fill: luma(90))[Render SMILES strings as 2D molecular diagrams]
  #v(0.4em)
  #text(size: 12pt, weight: "bold")[User Guide]
  #v(1.6em)
  #block(fill: accent-soft, radius: 8pt, inset: 18pt,
    smiles("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O", scale: 1.5))
  #v(1.4em)
  #text(size: 11pt)[Version #version]
]
#v(1.6fr)

// ── Header / footer for the body ─────────────────────────────────────────────

#set page(
  header: context {
    set text(size: 8.5pt, fill: luma(130))
    grid(columns: (1fr, auto),
      align(left)[typed-smiles · User Guide],
      align(right)[v#version])
    v(-0.6em)
    line(length: 100%, stroke: 0.4pt + luma(210))
  },
  footer: context {
    set text(size: 8.5pt, fill: luma(130))
    align(center, counter(page).display("1"))
  },
)
#counter(page).update(1)

#pagebreak()
#outline(title: [Contents], indent: 1.2em, depth: 2)

// ═════════════════════════════════════════════════════════════════════════════
= Introduction
// ═════════════════════════════════════════════════════════════════════════════

#demo[
  typed-smiles renders SMILES strings as 2D skeletal molecular diagrams inside
  Typst documents. A bundled WebAssembly plugin parses the SMILES, computes atom
  coordinates, implicit hydrogens, and stereochemistry, and the Typst layer draws
  the bonds, atom labels, colors, charges, hydrogens, wedges, and reaction
  helpers.

  This guide documents every function, argument, and package-specific syntax
  extension. Each feature comes with a runnable example: the source is on the
  left, its rendered output on the right.

  #example(```typ
  #smiles("OC1=CC=CC=C1C(=O)O")
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Getting started
// ═════════════════════════════════════════════════════════════════════════════

== Installation and import

Import the package from the Typst preview namespace. A wildcard import gives you
every public symbol:

```typ
#import "@preview/typed-smiles:0.2.0": *
```

Or import only what you need:

```typ
#import "@preview/typed-smiles:0.2.0": smiles, ce, mol, rxn-arrow, reaction
```

The package exports five symbols:

#table(
  columns: (auto, 1fr), inset: 7pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Symbol*], [*Purpose*],
  [#c("smiles")], [Render a SMILES string as a molecule (the main function).],
  [#c("ce")], [Chemical formulas and equations, re-exported from #c("chemformula").],
  [#c("mol")], [Wrap a molecule with an optional caption.],
  [#c("rxn-arrow")], [A reaction arrow with conditions above and below.],
  [#c("reaction")], [Lay out a multi-step reaction scheme.],
)

== Your first molecule

#demo[
  Call #c("smiles") with a SMILES string.

  #example(```typ
  Ethanol: #smiles("CCO")
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= The #raw("smiles()") function
// ═════════════════════════════════════════════════════════════════════════════

#demo[
  #c("smiles") is the main function. Its full signature is:

  ```typ
  #smiles(
    smiles-str,
    scale: 1.0,
    bond-length: none,
    font-size: none,
    font: "New Computer Modern",
    bond-stroke: none,
    color: true,
    rotation: 0deg,
    show-all-h: false,
  )
  ```
]

== Argument reference

#table(
  columns: (auto, auto, auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Argument*], [*Type*], [*Default*], [*Description*],
  [#c("smiles-str")], [`str`], [(required)], [The SMILES string.],
  [#c("scale")], [`float`], [`1.0`], [Balanced scale for bond length, label size, and stroke at once.],
  [#c("bond-length")], [`float` / `none`], [`none`], [Bond length factor (`1.0` is 30 pt). Overrides #c("scale") for length.],
  [#c("font-size")], [`length` / `none`], [`none`], [Atom-label font size. Overrides #c("scale") for labels.],
  [#c("font")], [`str`], [`"New Computer Modern"`], [Font family for atom labels.],
  [#c("bond-stroke")], [`length` / `none`], [`none`], [Bond stroke width. Overrides #c("scale") for strokes.],
  [#c("color")], [`bool`], [`true`], [Apply Jmol CPK atom colors.],
  [#c("rotation")], [`angle`], [`0deg`], [Rotate the molecule; labels stay upright.],
  [#c("show-all-h")], [`bool`], [`false`], [Also label carbon implicit hydrogens.],
)

#note[#c("scale") sizes everything together. Use #c("bond-length"),
#c("font-size"), and #c("bond-stroke") to override one dimension on its own;
each defaults to #c("none"), meaning it follows #c("scale").]

== #raw("smiles-str")

#demo[
  The positional argument is the SMILES string.

  #example(```typ
  #smiles("CC(=O)O")
  ```)
]

== #raw("scale")

#demo[
  Scales bond length, labels, and stroke proportionally.

  #example(```typ
  #smiles("OCC(O)CO", scale: 0.7) \
  #smiles("OCC(O)CO", scale: 1.0) \
  #smiles("OCC(O)CO", scale: 1.4)
  ```)
]

== #raw("bond-length")

#demo[
  Sets bond length on its own (`1.0` is 30 pt per bond), leaving label size and
  stroke under #c("scale").

  #example(```typ
  #smiles("CCCCCC", bond-length: 0.6) \
  #smiles("CCCCCC", bond-length: 1.1)
  ```)
]

== #raw("font-size")

#demo[
  #example(```typ
  #smiles("OC(=O)CN", font-size: 8pt) \
  #smiles("OC(=O)CN", font-size: 14pt)
  ```)
]

== #raw("font")

#demo[
  Match the document typeface, or pick something with character.

  #example(```typ
  #smiles("OC(=O)CN", font: "New Computer Modern") \
  #smiles("OC(=O)CN", font: "Libertinus Serif") \
  #smiles("OC(=O)CN", font: "PT Sans") \
  #smiles("OC(=O)CN", font: "Iosevka")
  ```)
]

== #raw("bond-stroke")

#demo[
  #example(```typ
  #smiles("C1=CC=CC=C1", bond-stroke: 0.5pt) \
  #smiles("C1=CC=CC=C1", bond-stroke: 1.6pt)
  ```)
]

== #raw("color")

#demo[
  CPK colors are on by default (see @sec-colors). Set #c("color: false") for a
  monochrome diagram.

  #example(```typ
  #smiles("C[N+](=O)[O-]") \
  #smiles("C[N+](=O)[O-]", color: false)
  ```)
]

== #raw("rotation")

#demo[
  Rotates the whole molecule while keeping every label upright.

  #example(```typ
  #smiles("CC(N)C(=O)O", rotation: 0deg) \
  #smiles("CC(N)C(=O)O", rotation: 90deg)
  ```)
]

== #raw("show-all-h")

#demo[
  Hydrogens on heteroatoms are shown by default; carbon hydrogens are hidden. Set
  #c("show-all-h: true") to label every implicit hydrogen.

  #example(```typ
  #smiles("CCO") \
  #smiles("CCO", show-all-h: true)
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Atoms and charges
// ═════════════════════════════════════════════════════════════════════════════

Standard SMILES syntax (atoms, bonds, branches, ring closures) is supported. This
section covers the points specific to typed-smiles.

== Aromatic rings

#demo[
  #warn[The bundled parser does not yet read lowercase aromatic atoms
  (#c("c1ccccc1")). Write aromatic rings in Kekulé form: uppercase atoms with
  explicit alternating double bonds.]

  #example(```typ
  #smiles("C1=CC=CC=C1") \
  #smiles("C1=CC=NC=C1")
  ```)
]

== Charges

#demo[
  Formal charges inside brackets render as a raised sign after the atom.

  #example(```typ
  #smiles("[O-]C1=CC=CC=C1") \
  #smiles("C[NH3+]") \
  #smiles("C[N+](=O)[O-]")
  ```)
]

#note[A bracketed element such as #c("[N]") is a nitrogen atom. To print an
arbitrary text label instead, use the abbreviation syntax #c("{N}")
(see @sec-abbrev).]

// ═════════════════════════════════════════════════════════════════════════════
= Hydrogens <sec-hydrogens>
// ═════════════════════════════════════════════════════════════════════════════

== Default display

#demo[
  Heteroatom hydrogens are shown by default; carbon hydrogens are hidden for a
  clean skeleton. #c("show-all-h: true") labels carbon hydrogens as well.

  #example(```typ
  #smiles("OCCN") \
  #smiles("OCCN", show-all-h: true)
  ```)
]

== Label orientation

#demo[
  Terminal heteroatom labels orient so the bond meets the heavy atom rather than
  the trailing hydrogen, at any angle.

  #example(```typ
  #smiles("NCCCO")
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Stereochemistry <sec-stereo>
// ═════════════════════════════════════════════════════════════════════════════

typed-smiles supports both standard SMILES stereo notations, plus a manual
override for drawing wedges directly.

== Tetrahedral centers: #raw("@") and #raw("@@")

#demo[
  A bracket atom carrying #c("@") or #c("@@") becomes a wedge (toward the viewer)
  or hashed bond (away), computed from the 2D layout so the depiction matches the
  SMILES in the conventional orientation.

  #example(```typ
  #smiles("N[C@@H](C)C(=O)O") \
  #smiles("N[C@H](C)C(=O)O")
  ```)
]

#demo[
  The renderer wedges an exocyclic substituent (such as #ce("OH") or #ce("CH3"))
  where one exists, and an explicit hydrogen otherwise. Fused systems work too.

  #example(```typ
  #smiles(
    "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
    scale: 0.85,
  )
  ```)
]

#note[The geometry is drawn correctly, but R/S descriptors are not computed. Ring
stereochemistry between adjacent centers and bridged bicyclics may need a manual
adjustment (see @sec-limits).]

== Double bonds: #raw("/") and #raw("\\")

#demo[
  Directional bonds #c("/") and #c("\\") around a double bond set its cis/trans
  geometry. They come in pairs, one on each side of the #c("=").

  #example(```typ
  trans: #smiles("F/C=C/F")
  #h(2em)
  cis: #smiles("F/C=C\F")
  ```)
]

== Manual wedges: #raw("!w") and #raw("!h")

#demo[
  For direct control, #c("!w") draws a solid wedge and #c("!h") a hashed wedge.
  These are typed-smiles extensions, not standard SMILES. Like plain bonds, they
  are colored by the atoms at each end.

  #example(```typ
  #smiles("C!wN") \
  #smiles("C!hN") \
  #smiles("{Nu}!wC({LG|red})=O")
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Colors <sec-colors>
// ═════════════════════════════════════════════════════════════════════════════

Atoms use the Jmol CPK palette. Carbon and unlisted elements are black; common
heteroatoms get their conventional colors.

#table(
  columns: 8, inset: 6pt, align: center + horizon, stroke: 0.5pt + luma(210),
  fill: (_, y) => if y == 0 { accent-soft },
  [*N*], [*O*], [*S*], [*P*], [*F*], [*Cl*], [*Br*], [*I*],
  ..(
    rgb("#3050F8"), rgb("#FF0D0D"), rgb("#E6C800"), rgb("#FF8000"),
    rgb("#90E050"), rgb("#1FF01F"), rgb("#A62929"), rgb("#940094"),
  ).map(clr => box(width: 100%, height: 1.1em, fill: clr, radius: 2pt)),
)

#demo[
  Set #c("color: false") for a monochrome diagram.

  #example(```typ
  #smiles("OC1=CC(=CC=C1)C(=O)O") \
  #smiles("OC1=CC(=CC=C1)C(=O)O", color: false)
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Abbreviated groups <sec-abbrev>
// ═════════════════════════════════════════════════════════════════════════════

== Plain labels

#demo[
  Wrap text in braces #c("{...}") to place a labeled pseudo-atom that bonds like
  any other atom. Use it for groups you do not want to draw in full.

  #example(```typ
  #smiles("{PPh3}C=O") \
  #smiles("{OEt}C(=O){NHR}")
  ```)
]

== Colored labels

#demo[
  Add #c("|style") inside the braces to color a label. The styles are #c("red"),
  #c("blue"), #c("green"), #c("black"), #c("gray") (or #c("grey")), #c("orange"),
  and #c("purple"). Any element symbol also works and uses that element's color.

  #example(```typ
  #smiles("{Nu|blue}CC{LG|red}") \
  #smiles("{R|green}C(=O){OR|O}")
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Chemical formulas: #raw("ce()")
// ═════════════════════════════════════════════════════════════════════════════

#demo[
  For formulas and text equations, #c("ce") is re-exported from the
  #link("https://typst.app/universe/package/chemformula")[chemformula] package,
  so one import covers both structures and formulas.

  #example(```typ
  #ce("H2SO4") #h(1em) #ce("Ca^2+") #h(1em) #ce("2 H2 + O2 -> 2 H2O")
  ```)
]

== Fonts and size

#demo[
  #c("ce") accepts #c("font") and #c("font-size"); any other arguments pass
  through to chemformula.

  #example(```typ
  #ce("CuSO4 * 5 H2O") \
  #ce("CuSO4 * 5 H2O", font: "Libertinus Serif", font-size: 13pt) \
  #ce("CuSO4 * 5 H2O", font: "PT Sans", font-size: 13pt)
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Reaction schemes
// ═════════════════════════════════════════════════════════════════════════════

Three helpers compose molecules, formulas, and arrows into schemes.

== #raw("mol()")

#demo[
  #c("mol(content, label: none)") wraps any content and centers an optional
  caption beneath it.

  #example(```typ
  #mol(smiles("CCO"), label: [ethanol])
  ```)
]

== #raw("rxn-arrow()")

#c("rxn-arrow(above: none, below: none, dir: \"right\")") draws an arrow.
#c("dir") may be #c("\"right\""), #c("\"left\""), #c("\"up\""), or #c("\"down\"").
#c("above") and #c("below") carry reagents and conditions.

== #raw("reaction()")

#demo[
  #c("reaction(gap-h: 1.5em, gap-v: 1.5em, ..items)") lays out molecules and
  arrows left to right. An up or down arrow wraps the scheme onto a new row.

  #example(```typ
  #reaction(
    mol(smiles("C1=CC=CC=C1"), label: [benzene]),
    rxn-arrow(above: ce("Br2"), below: ce("FeBr3")),
    mol(smiles("BrC1=CC=CC=C1"), label: [bromobenzene]),
  )
  ```, side: false)
]

#demo[
  A wrap-around scheme using a downward arrow:

  #example(```typ
  #reaction(
    mol(smiles("C1=CC=CC=C1"), label: [1]),
    rxn-arrow(above: ce("HNO3"), below: ce("H2SO4")),
    mol(smiles("[O-][N+](=O)C1=CC=CC=C1"), label: [2]),
    rxn-arrow(dir: "down", above: [reduce]),
    mol(smiles("NC1=CC=CC=C1"), label: [3]),
  )
  ```, side: false)
]

// ═════════════════════════════════════════════════════════════════════════════
= Limitations <sec-limits>
// ═════════════════════════════════════════════════════════════════════════════

- Lowercase aromatic SMILES (#c("c1ccccc1")) are not parsed. Use Kekulé notation.
- Directional bonds (#c("/"), #c("\\")) draw the correct cis/trans geometry, but
  E/Z descriptors are not computed.
- R/S and E/Z descriptors are not calculated.
- Ring stereochemistry between adjacent centers and bridged bicyclics can overlap
  or need a manual adjustment (try #c("rotation"), or the manual #c("!w") and
  #c("!h") wedges).

// ═════════════════════════════════════════════════════════════════════════════
= Quick reference
// ═════════════════════════════════════════════════════════════════════════════

== Syntax

#table(
  columns: (auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Syntax*], [*Meaning*],
  [#c("[...]")], [Bracket atom (any element, charges, explicit H).],
  [#c("@") / #c("@@")], [Tetrahedral anticlockwise / clockwise.],
  [#c("/") #c("\\")], [Double-bond cis/trans geometry.],
  [#c("!w") / #c("!h")], [Manual solid / hashed wedge (typed-smiles extension).],
  [#c("{label}")], [Abbreviated group pseudo-atom.],
  [#c("{label|style}")], [Colored abbreviated group.],
)

== #raw("smiles()") options

#table(
  columns: (auto, auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Option*], [*Default*], [*Effect*],
  [#c("scale")], [`1.0`], [Balanced size of everything.],
  [#c("bond-length")], [`none`], [Bond length only (`1.0` is 30 pt).],
  [#c("font-size")], [`none`], [Atom-label size only.],
  [#c("font")], [`"New Computer Modern"`], [Atom-label font.],
  [#c("bond-stroke")], [`none`], [Bond width only.],
  [#c("color")], [`true`], [CPK colors on or off.],
  [#c("rotation")], [`0deg`], [Rotate, labels stay upright.],
  [#c("show-all-h")], [`false`], [Label carbon hydrogens too.],
)
