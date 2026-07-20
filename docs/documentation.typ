// typed-smiles documentation
//
// Compile with:
//   typst compile --root . docs/documentation.typ docs/documentation.pdf

#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#import "../src/lib.typ": smiles, smiles-inline, smiles-cetz, ce, rxn-arrow, mol, reaction, cycle, step, atom, bond, lp, species, arrow, highlight, brackets, mol-weight
#import "@preview/cetz:0.5.2"

#let version = "0.7.0"
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

#let pkg-scope = (
  smiles: smiles, smiles-inline: smiles-inline, smiles-cetz: smiles-cetz,
  ce: ce, rxn-arrow: rxn-arrow,
  mol: mol, reaction: reaction, cycle: cycle, step: step,
  atom: atom, bond: bond, lp: lp, species: species, arrow: arrow,
  highlight: highlight, brackets: brackets, mol-weight: mol-weight,
  cetz: cetz,
)

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
#import "@preview/typed-smiles:0.7.0": *
```

Or import only what you need:

```typ
#import "@preview/typed-smiles:0.7.0": smiles, ce, mol, rxn-arrow, reaction
```

The package exports these main symbols:

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
  [#c("atom"), #c("bond"), #c("lp"), #c("species")], [References for arrows and highlights.],
  [#c("arrow"), #c("highlight")], [Mechanism arrows and shaded atom/bond regions.],
  [#c("brackets")], [Draw square brackets with optional corner marks.],
  [#c("mol-weight")], [Compute molecular weight from a SMILES string.],
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
    style: "default",
    scale: 1.0,
    bond-length: none,
    font-size: none,
    font: auto,
    bond-stroke: none,
    color: auto,
    fg: auto,
    theme: auto,
    rotation: 0deg,
    mirror: none,
    show-h: (),
    aromatic: "kekule",
    atom-annotations: (),
    opacity: 100%,
    bond-customizations: (),
    lone-pairs: none,
    atom-colors: (:),
    show-indices: false,
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
  [#c("style")], [`str`], [`"default"`], [Journal style preset (#c("\"acs\""), #c("\"rsc\""), #c("\"nature\""), #c("\"wiley\"")) filling in bond length, label size, stroke, and font. Explicit arguments win; #c("\"default\"") applies nothing.],
  [#c("scale")], [`float`], [`1.0`], [Balanced scale for bond length, label size, and stroke at once.],
  [#c("bond-length")], [`float` / `none`], [`none`], [Bond length factor (`1.0` is 30 pt). Overrides #c("scale") for length.],
  [#c("font-size")], [`length` / `none`], [`none`], [Atom-label font size. Overrides #c("scale") for labels.],
  [#c("font")], [`auto` / `str` / `array`], [`auto`], [Font family for atom labels; `auto` is "New Computer Modern", or the preset's font when #c("style") is set.],
  [#c("bond-stroke")], [`length` / `none`], [`none`], [Bond stroke width. Overrides #c("scale") for strokes.],
  [#c("color")], [`auto` / `bool`], [`auto`], [Apply Jmol CPK atom colors. `auto` is colored for #c("\"default\"") and monochrome for journal presets.],
  [#c("fg")], [`auto` / `color`], [`auto`], [Foreground for bonds and carbon labels; `auto` inherits the surrounding text color.],
  [#c("theme")], [`auto` / `"light"` / `"dark"`], [`auto`], [CPK palette variant; `auto` picks "dark" when #c("fg") is light.],
  [#c("rotation")], [`angle`], [`0deg`], [Rotate the molecule; labels stay upright.],
  [#c("mirror")], [`none` / `"horizontal"` / `"vertical"`], [`none`], [Optional horizontal or vertical page-axis reflection.],
  [#c("show-h")], [`"all"` / `int` / `array`], [`()`], [Implicit hydrogens to label beyond the default heteroatom hydrogens. Use #c("\"all\"") for every atom.],
  [#c("aromatic")], [`"kekule"` / `"circle"`], [`"kekule"`], [Depict lowercase-aromatic rings with alternating double bonds or with an inscribed circle.],
  [#c("atom-annotations")], [`array`], [`()`], [Small gray side labels as #c("(index, content)") or #c("(index, content, offset)") tuple entries.],
  [#c("opacity")], [`ratio` / `float`], [`100%`], [Fade the whole drawing — bonds, labels, charges, lone pairs — e.g. for ghost molecules.],
  [#c("bond-customizations")], [`array`], [`()`], [Per-bond style overrides as #c("(bond(i, j), (..options..))") pairs; options are #c("color"), #c("stroke"), and #c("opacity").],
  [#c("lone-pairs")], [`none` / `"dots"` / `"lines"`], [`none`], [Draw optional non-bonding electron pairs on skeletal atom labels.],
  [#c("atom-colors")], [`dictionary`], [`(:)`], [Color overrides for elements and labels. Element-symbol keys (e.g. #c("O: red")) override CPK colors; brace-quoted label keys (e.g. #c("\"{PPh3}\": purple")) override a specific abbreviated group. See @sec-colors.],
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

== #raw("style")

#demo[
  Journal style presets approximate the corresponding ChemDraw document
  stylesheets using the journals' published drawing settings — bond length,
  atom-label size, line width, a Helvetica/Arial font stack, and monochrome
  atom labels by default:

  #table(
    columns: (auto, auto, auto, auto), inset: 5.5pt,
    align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
    fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
    [*Preset*], [*Bond*], [*Line*], [*Labels*],
    [#c("\"acs\"")], [14.4 pt], [0.6 pt], [10 pt],
    [#c("\"rsc\"")], [12.2 pt], [0.5 pt], [7 pt],
    [#c("\"nature\"")], [10.8 pt], [0.6 pt], [6 pt],
    [#c("\"wiley\"")], [17 pt], [0.7 pt], [12 pt],
  )

  #example(```typ
  #smiles("CC(N)C(=O)O") \
  #smiles("CC(N)C(=O)O", style: "acs") \
  #smiles("CC(N)C(=O)O", style: "nature") \
  #smiles("CC(N)C(=O)O", style: "acs", color: true)
  ```)
]

A preset only fills in arguments left unset: `style: "acs"` with
`font-size: 14pt` keeps the ACS bond length under your label size, and
`scale` multiplies the whole preset. `color: auto` makes journal presets
monochrome; pass `color: true` for CPK colors or `color: false` for explicit
black line art. `"default"` applies nothing, so existing documents are
unaffected. Wiley's line width and label size are derived from its 6 mm
bond-length guideline, which fixes only a minimum line width; journal-specific
geometry beyond these four settings (double-bond spacing, hash spacing) keeps
typed-smiles' defaults.

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
  CPK colors are on by default for #c("style: \"default\"") (see @sec-colors).
  Journal presets default to monochrome line art. Set #c("color: false") for an
  explicit monochrome diagram, or #c("color: true") to force CPK colors.

  #example(```typ
  #smiles("C[N+](=O)[O-]") \
  #smiles("C[N+](=O)[O-]", color: false)
  ```)
]

== #raw("fg") and #raw("theme")

#demo[
  Bond strokes and carbon labels follow #c("fg"), which defaults to #c("auto")
  and inherits the surrounding text color — on a dark slide theme (Touying,
  Polylux, or any #c("set text(fill: white)")) molecules recolor with no
  configuration. #c("theme: auto") switches to a dark CPK variant whenever the
  foreground is light, lifting only the hues that vanish on dark backgrounds
  (N, O, Br, I, and dark named label colors such as navy or maroon).

  #example(```typ
  #block(fill: rgb("#1E1E24"), inset: 8pt,
    radius: 4pt)[
    #set text(fill: white)
    #smiles("NC(Br)C(I)C(=O)O")
  ]
  ```)
]

#demo[
  Both can be forced: an explicit #c("fg") recolors the skeleton on any
  background, and #c("theme: \"dark\"") or #c("theme: \"light\"") pins the
  palette. To set a document-wide default, rebind once in the preamble:
  #c("#let smiles = smiles.with(fg: white, theme: \"dark\")").

  #example(```typ
  #smiles("CC(N)C(=O)O", fg: maroon) \
  #smiles("CC(N)C(=O)O", color: false, fg: maroon)
  ```)
]

#note[Reaction arrows inherit the surrounding text color the same way via
#c("rxn-arrow(color: auto)"), so full schemes work on dark slides too.]

== #raw("rotation")

#demo[
  Rotates the whole molecule while keeping every label upright.

  #example(```typ
  #smiles("CC(N)C(=O)O", rotation: 0deg) \
  #smiles("CC(N)C(=O)O", rotation: 90deg)
  ```)
]

== #raw("mirror")

#demo[
  #c("mirror") reflects a molecule horizontally or vertically in the final
  drawing. Useful to face a reacting group toward a reaction arrow. Inside
  #c("reaction()"), pass it per molecule:
  #c("mol(\"OC1CCCCC1\", mirror: \"horizontal\")"). With #c("rotation") set,
  #c("mirror: \"vertical\"") preserves left/right and flips top/bottom.

  #example(```typ
  #smiles("CC(=O)OC1=CC=CC=C1C(=O)O") \
  #smiles(
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    mirror: "horizontal",
  )
  ```)
]

#note[Single-axis reflection exchanges wedges and hashes: a reflected 2D
skeleton with unchanged wedges would depict the enantiomer, while reflecting
and swapping is equivalent to turning the molecule over — the depicted
configuration at every stereocenter is preserved.]

#demo[
  A stereocenter keeps its configuration when mirrored.

  #example(```typ
  #smiles("N[C@@H](C)C(=O)O") \
  #smiles(
    "N[C@@H](C)C(=O)O",
    mirror: "horizontal",
  )
  ```)
]

== #raw("show-h")

#demo[
  Hydrogens on heteroatoms are shown by default; carbon hydrogens are hidden.
  #c("show-h") selects carbon hydrogens without needing a trailing comma for one
  atom: use an integer for one atom, an array for several atoms, or #c("\"all\"")
  for every implicit hydrogen. Use #c("show-indices: true") to read off the
  indices while writing.

  #example(```typ
  #smiles("CC(N)C(=O)O") \
  #smiles("CC(N)C(=O)O", show-h: 1) \
  #smiles("CC(N)C(=O)O", show-h: (0, 1)) \
  #smiles("CC(N)C(=O)O", show-h: "all")
  ```)
]

== #raw("atom-annotations")

#demo[
  Small annotations for NMR numbering, Greek positions, or footnote marks are
  drawn small and gray on the emptiest side of an atom. #c("atom-annotations")
  is always a tuple list. Each entry is either #c("(index, content)") or
  #c("(index, content, offset)"). Values are content, so wrap them in
  #c("text()") to restyle.

  #example(```typ
  #smiles(
    "N[C@@H](C)C(=O)O",
    atom-annotations: (
      (1, [$alpha$], (0, -0.05)),
      (2, [$beta$]),
      (3, [$gamma$], (-0.4, -0.3)),
    ),
  )
  ```)
]

== #raw("opacity")

#demo[
  #c("opacity") fades every paint the drawing produces — bonds, atom labels,
  charges, lone pairs, and annotations — so a molecule can be ghosted or
  de-emphasized, e.g. spectator species in a mechanism or step-by-step slide
  builds. It accepts a ratio (#c("40%")) or a number between 0 and 1.

  #example(```typ
  #smiles("CC(N)C(=O)O")
  #h(1.5em)
  #smiles("CC(N)C(=O)O", opacity: 55%)
  #h(1.5em)
  #smiles("CC(N)C(=O)O", opacity: 0.2)
  ```)
]

#note[String molecules inside #c("reaction()") accept the same option:
#c("mol(\"CCO\", opacity: 30%)").]

== #raw("bond-customizations")

#demo[
  #c("bond-customizations") restyles individual bonds. Each entry pairs a
  #c("bond(i, j)") reference (the same atom-index form used by #c("arrow()")
  and #c("highlight()"); a plain #c("(i, j)") array also works) with an
  options dictionary:

  - #c("color") — draw the whole bond in one color instead of the bicolor halves,
  - #c("stroke") — the bond width, like a per-bond #c("bond-stroke"),
  - #c("opacity") — fade just this bond.

  Overrides apply to every part of the bond: both lines of a double bond, the
  hash lines of a stereo bond, waves and dashes.

  #example(```typ
  #smiles(
    "CCC(=O)OCC",
    bond-customizations: (
      (bond(1, 2), (color: yellow, stroke: 2pt)),
    ),
  )
  ```)
]

#note[Use #c("show-indices: true") while writing the references, exactly as for
mechanism arrows. Bond geometry (length and angles) is computed by the layout
engine and is not a per-bond style; use #c("bond-length") to size all bonds
together.]

== #raw("lone-pairs")

#demo[
  Set #c("lone-pairs") to #c("\"dots\"") or #c("\"lines\"") to annotate the
  skeletal drawing with non-bonding electron pairs on common organic heteroatoms
  and charged atoms. The default #c("none") keeps lone pairs hidden.

  #example(```typ
  #smiles("CCO", lone-pairs: "dots") \
  #smiles("CCN", lone-pairs: "lines") \
  #smiles("[O-]C=O", lone-pairs: "dots")
  ```)
]

#note[Lone pairs are inferred for a conservative organic subset such as N, O, S,
P, halogens, and simple charged forms. They are an annotation on the skeletal
structure; the carbon skeleton and implicit-hydrogen conventions do not change.
For example, ammonium (#c("[NH4+]")) has no nitrogen lone pair. On terminal
hydrogen-bearing labels such as #c("OH"), #c("OH-"), and #c("NH2"), lone pairs are
placed from the rendered heavy-atom glyph center. For charged bare atoms such as
#c("[O-]") or #c("[Br-]"), the charge mark does not move the atom center used for
lone pairs.]

// ═════════════════════════════════════════════════════════════════════════════
= Atoms and charges
// ═════════════════════════════════════════════════════════════════════════════

Standard SMILES syntax (atoms, bonds, branches, ring closures) is supported. This
section covers the points specific to typed-smiles.

#note[Ring closures may be written next to branch points. For example,
#c("C1=CCCC(=O)1") closes the ring on the same atom that carries the
#ce("C=O") branch.]

== Aromatic rings

#demo[
  Aromatic rings can be written in either OpenSMILES form: lowercase aromatic
  atoms (#c("c1ccccc1")) or Kekulé notation with explicit alternating double
  bonds (#c("C1=CC=CC=C1")). Aromatic input is kekulized on parse, so both
  render identically.

  #example(```typ
  #smiles("c1ccccc1") \
  #smiles("C1=CC=CC=C1") \
  #smiles("c1ccncc1")
  ```)
]

#demo[
  Bracket aromatic atoms carry their hydrogen count explicitly, as in pyrrole's
  #c("[nH]"). A single bond between two aromatic rings (biphenyl) may be
  written with or without the explicit #c("-").

  #example(```typ
  #smiles("c1cc[nH]c1") \
  #smiles("c1occc1") \
  #smiles("c1ccccc1-c1ccccc1", scale: 0.8)
  ```)
]

#demo[
  #c("aromatic: \"circle\"") depicts lowercase-aromatic rings as single bonds
  with an inscribed circle instead of alternating double bonds. Each fully
  aromatic ring of a fused system gets its own circle; Kekulé-written input
  always keeps its explicit bonds.

  #example(```typ
  #smiles("c1ccccc1", aromatic: "circle") \
  #smiles("c1ccc2ccccc2c1", aromatic: "circle") \
  #smiles("Cc1ccncc1", aromatic: "circle")
  ```)
]

#note[Kekulization does not renumber anything: atom indices follow SMILES
writing order for aromatic input exactly as for Kekulé input, so
#c("show-indices"), #c("highlight(...)"), and #c("arrow(...)") references work
unchanged, including on ring bonds whose double bond only appears after
kekulization. A ring that cannot be kekulized (for example pyrrole written
#c("c1ccnc1"), without the #c("[nH]") hydrogen) or an aromatic atom outside a
ring is reported as an error.]

== Charges

#demo[
  Formal charges inside brackets render as a raised sign after the atom.
  Hydrogen counts and charge marks use the atom-label size so they remain
  legible when labels are scaled.

  #example(```typ
  #smiles("[O-]C1=CC=CC=C1") \
  #smiles("C[NH3+]") \
  #smiles("C[N+](=O)[O-]")
  ```)
]

#note[A bracketed element such as #c("[N]") is a nitrogen atom. To print an
arbitrary text label instead, use the abbreviation syntax #c("{N}")
(see @sec-abbrev).]

== Salts and disconnected structures

#demo[
  A dot (#c(".")) means "no bond": it separates disconnected fragments, as in
  salts, counterions, and hydrates. Each fragment is laid out on its own and
  the fragments are arranged left to right in writing order.

  #example(```typ
  #smiles("CC(=O)[O-].[Na+]") \
  #smiles("[NH4+].[Cl-]") \
  #smiles("CC[NH3+].[Cl-]")
  ```)
]

#note[Atom indices stay global across fragments, still following SMILES
writing order, so #c("show-indices"), #c("highlight(...)"), and
#c("arrow(...)") can reference atoms in any fragment of the same
#c("smiles()") call. A ring closure written across a dot (as in the
OpenSMILES example #c("C1.C1"), ethane) still forms its bond.]

// ═════════════════════════════════════════════════════════════════════════════
= Hydrogens <sec-hydrogens>
// ═════════════════════════════════════════════════════════════════════════════

== Default display

#demo[
  Heteroatom hydrogens are shown by default; carbon hydrogens are hidden for a
  clean skeleton. #c("show-h: \"all\"") labels carbon hydrogens as well.

  #example(```typ
  #smiles("OCCN") \
  #smiles("OCCN", show-h: "all")
  ```)
]

== Explicit bracket hydrogens

#demo[
  A hydrogen written as its own bracket atom (for example the #c("[H]") atoms in
  #c("C([H])([H])[H]")) is folded into its neighbor's hydrogen count, exactly
  like an implicit hydrogen. Carbon-bound hydrogens stay hidden and heteroatom
  hydrogens are still labeled, so a fully hydrogen-suppressed SMILES depicts the
  same as its skeletal form. Hydrogen counts written inside a single bracket,
  such as #c("[OH-]") or #c("[NH4+]"), are unaffected and still shown.

  #example(```typ
  #smiles("C([H])([H])([H])[H]") \
  #smiles("N([H])([H])C(=O)O")
  ```)
]

== Isotopes

#demo[
  A bracket isotope is drawn as a leading superscript mass number. Isotopically
  labeled hydrogens such as deuterium #c("[2H]") are kept as drawn atoms rather
  than folded away.

  #example(```typ
  #smiles("[2H]OC([2H])([2H])[2H]") \
  #smiles("[13CH3]C(=O)O")
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
  geometry. They come in pairs, one on each side of the #c("="). If an alkene
  carbon has two marked substituent bonds, the markers must place those
  substituents on opposite sides.

  #example(```typ
  trans: #smiles("F/C=C/F")
  #h(2em)
  cis: #smiles("F/C=C\F")
  ```)
]

== Square-planar centers: #raw("@SP1"), #raw("@SP2"), #raw("@SP3")

#demo[
  Square-planar stereochemistry is depicted exactly: the four neighbors sit at
  90° around the center, arranged so the line traced through them in SMILES
  order forms the class's shape — 'U' for #c("@SP1"), '4' for #c("@SP2"), 'Z'
  for #c("@SP3"). The classic example is cis- versus trans-platin.

  #example(```typ
  cis: #smiles("N[Pt@SP1](N)(Cl)Cl")
  #h(2em)
  trans: #smiles("N[Pt@SP2](N)(Cl)Cl")
  ```)
]

#note[Trigonal-bipyramidal (#c("@TB1")–#c("@TB20")), octahedral
(#c("@OH1")–#c("@OH30")), and allenal (#c("@AL1"), #c("@AL2")) centers are
accepted and drawn with correct connectivity, but without stereo decoration —
their 3D arrangement is not translated into wedges.]

#note[Cumulated double bonds always label the central sp carbon (#ce("O=C=O"),
ketenes, allenes): the neighbors are collinear, so without the #c("C") the two
double bonds would read as one long double bond.]

== Manual wedges: #raw("!w") and #raw("!h")

#demo[
  For direct control, #c("!w") draws a solid wedge and #c("!h") a hashed wedge.
  These are typed-smiles extensions, not standard SMILES. Like plain bonds, they
  are colored by the atoms at each end.

  #example(```typ
  #smiles("C!wN") \
  #smiles("C!hN") \
  #smiles("C(!w{Nu|teal})({LG|green})=O")
  ```)
]

== Wavy and dashed bonds: #raw("!s") and #raw("!d")

#demo[
  Two more drawing extensions on single bonds: #c("!s") draws a wavy
  (squiggly) bond — the conventional depiction of unspecified stereochemistry
  or an attachment point — and #c("!d") draws a dashed bond for hydrogen
  bonds, partial bonds, and coordination. Like plain bonds, they are colored
  by the atoms at each end.

  #example(```typ
  #smiles("C!sN") \
  #smiles("C!dN") \
  #smiles("CC(!sC)C1=CC=CC=C1") \
  #smiles("CN!d[Cu]")
  ```)
]

#note[#c("!s") and #c("!d") carry no cis/trans meaning and never consume a
double-bond geometry slot; they can appear next to real #c("/") and #c("\\")
directional bonds.]
#note[Tip: Use empty custom labels #c("{}") or #c("!w{}") to draw bonds without a connecting atom.]

// ═════════════════════════════════════════════════════════════════════════════
= Colors <sec-colors>
// ═════════════════════════════════════════════════════════════════════════════

== Jmol CPK palette

Atoms use the Jmol CPK palette by default. Carbon and unlisted elements are
black; common heteroatoms get their conventional colors.

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
  #smiles("OC1=CC(=NC=C1)C(=S)Cl") \
  #smiles("OC1=CC(=NC=C1)C(=S)Cl", color: false)
  ```)
]

== Custom colors (`atom-colors`)

`atom-colors` is a single dictionary that overrides colors for both elements
and arbitrary abbreviated groups. Any Typst color value works — #c("rgb()"),
#c("luma()"), named constants, or anything else.

=== Element-symbol keys

#demo[
  Use an element symbol as the key to replace that element's CPK color
  everywhere it appears, including when referenced as an inline label style
  (#c("{OMe|O}") would pick up the overridden oxygen color).

  #example(```typ
  // default CPK
  #smiles("CC(N)C(=O)O")
  // oxygen → dark brown; nitrogen → teal
  #smiles("CC(N)C(=O)O",
    atom-colors: (O: rgb("#8B4513"), N: rgb("#008080")))
  ```, side: false)
]

=== Label-name keys

#demo[
  Wrap the label text in braces to target a specific abbreviated group,
  regardless of any inline style it carries. The key is written as a quoted
  string because `{` is not a valid Typst identifier character.

  #example(```typ
  // default — colors come from inline styles or element fallback
  #smiles("{>PPh3|P}C({OEt|O})=O")
  // override by label name — inline style is ignored for these
  #smiles("{>PPh3|P}C({OEt|O})=O",
    atom-colors: ("{PPh3}": rgb("#7B2D8B"), "{OEt}": rgb("#008080")))
  ```, side: false)
]

=== Mixing both key forms

#demo[
  Element keys and label keys can coexist in the same dictionary.

  #example(```typ
  // O (element) → teal; {Nu} (label) → navy; {LG} (label) → maroon
  #smiles("C(!w{Nu})({LG|red})=O",
    atom-colors: (O: rgb("#00897B"), "{Nu}": rgb("#1565C0"), "{LG}": rgb("#B71C1C")))
  ```)
]

=== Priority

Color is resolved in this order for each atom or label:

+ Label-name key in `atom-colors` (e.g. #c("\"{PPh3}\"")), if the atom is an abbreviation.
+ Element-symbol key in `atom-colors` (e.g. #c("O")), for both real atoms and element-style labels.
+ The inline `{label|style}` style from the SMILES string (named color, hex, or element symbol).
+ The default CPK palette.

#note[`color: false` is a hard override that sits above all of the above. When
it is set, every atom and label renders in black regardless of any `atom-colors`
entries or inline styles. If you want a mostly-monochrome diagram with one or
two highlighted groups, keep or set `color: true` and put everything else in
`atom-colors`.]

== Label colors (`{label|style}`)

#demo[
  In the SMILES string, write #c("{label|style}") to color an abbreviated group
  independently of its element. The style can be a named color, an element
  symbol, or a hex code.

  *Named colors:*

  #table(
    columns: (auto, auto, auto, auto), inset: 7pt,
    align: (x, y) => if x == 0 { left } else { center + horizon },
    stroke: 0.5pt + luma(210),
    fill: (_, y) => if y == 0 { accent-soft },
    [*Name*], [*Swatch*], [*Name*], [*Swatch*],
    [red],    [#box(width: 2em, height: 0.9em, fill: rgb("#FF0D0D"), radius: 2pt)],
    [orange], [#box(width: 2em, height: 0.9em, fill: rgb("#FF8000"), radius: 2pt)],
    [yellow], [#box(width: 2em, height: 0.9em, fill: rgb("#E6C800"), radius: 2pt)],
    [brown],  [#box(width: 2em, height: 0.9em, fill: rgb("#8B4513"), radius: 2pt)],
    [green],  [#box(width: 2em, height: 0.9em, fill: rgb("#1FA51F"), radius: 2pt)],
    [lime],   [#box(width: 2em, height: 0.9em, fill: rgb("#32CD32"), radius: 2pt)],
    [teal],   [#box(width: 2em, height: 0.9em, fill: rgb("#008080"), radius: 2pt)],
    [cyan],   [#box(width: 2em, height: 0.9em, fill: rgb("#00B4D8"), radius: 2pt)],
    [blue],   [#box(width: 2em, height: 0.9em, fill: rgb("#3050F8"), radius: 2pt)],
    [navy],   [#box(width: 2em, height: 0.9em, fill: rgb("#000080"), radius: 2pt)],
    [purple], [#box(width: 2em, height: 0.9em, fill: rgb("#940094"), radius: 2pt)],
    [pink],   [#box(width: 2em, height: 0.9em, fill: rgb("#FF69B4"), radius: 2pt)],
    [black],  [#box(width: 2em, height: 0.9em, fill: black, radius: 2pt)],
    [gray],   [#box(width: 2em, height: 0.9em, fill: rgb("#777777"), radius: 2pt)],
    [silver], [#box(width: 2em, height: 0.9em, fill: rgb("#C0C0C0"), radius: 2pt)],
    [maroon], [#box(width: 2em, height: 0.9em, fill: rgb("#800000"), radius: 2pt)],
    [white],  [#box(width: 2em, height: 0.9em, fill: white, stroke: 0.5pt + luma(210), radius: 2pt)],
    [],       [],
  )
]

#demo[
  *Hex colors:* write a `#RRGGBB` hex code after the pipe to use any color.
  The `#` is just a regular character inside a Typst string literal.

  #example(```typ
  #smiles("{OMe|#8B4513}C(=O){NHAc|#5B2E8C}")
  ```)
]

#demo[
  *Element symbol:* using a symbol as the style applies that element's (possibly
  overridden) CPK color to the label and its bonds.

  #example(```typ
  #smiles("{>PPh3|P}C({OEt|O})=O")
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Abbreviated groups <sec-abbrev>
// ═════════════════════════════════════════════════════════════════════════════

== Plain labels

#demo[
  Wrap text in braces #c("{...}") to place a labeled pseudo-atom that bonds like
  any other atom. Use #c(">") inside the label to choose the attachment glyph;
  the marker is not displayed. Rotate the molecule when needed so the bond
  approaches the chosen glyph roughly perpendicular to the written label.

  #example(```typ
  #smiles("{>PPh3}C=O") \
  #smiles("{Me>O}C=O") \
  #smiles("{CA>T}C=O") \
  #smiles("{OEt}C(=O){NHR}")
  ```)
]

== Colored labels

#demo[
  Add #c("|style") inside the braces to color a label. Use a named color,
  an element symbol, or a hex code — see @sec-colors for the full list.

  #example(```typ
  #smiles("{Nu|blue}CC{LG|red}") \
  #smiles("{R|green}C(=O){OR|O}") \
  #smiles("{Cat|teal}C(=O){Nuc|#E040FB}")
  ```)
]

== Scripted labels

#demo[
  Use #c("_(...)"), #c("^(...)"), or their one-character forms #c("_4")
  and #c("^+") for explicit subscripts and superscripts in a label. When both
  follow one glyph, they are paired with that glyph: #c("{NH_4^+}") places the
  plus above H rather than above 4. Script size, spacing, and vertical placement
  match #c("ce()") notation. Parenthesize multi-character scripts.

  #example(```typ
  #smiles("{>PPh_(3)}C=O") \
  #smiles("{NH_4^+}") \
  #smiles("{SO_(4)^(2-)}")
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
= Molecular weight: #raw("mol-weight()")
// ═════════════════════════════════════════════════════════════════════════════

#demo[
  #c("mol-weight(smiles)") returns the molecular weight in g/mol as a plain
  #c("float"): the sum of IUPAC standard atomic weights over every atom,
  including implicit and explicit hydrogens. Round it for display with
  #c("calc.round").

  #example(```typ
  Ethanol: #calc.round(mol-weight("CCO"), digits: 2) g/mol \
  Caffeine: #calc.round(
    mol-weight("CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    digits: 2,
  ) g/mol \
  Sodium acetate: #calc.round(
    mol-weight("CC(=O)[O-].[Na+]"),
    digits: 2,
  ) g/mol
  ```)
]

Dot-separated fragments are summed together, so salts and hydrates work as
written. Charges are accepted; as in standard formula-weight arithmetic, the
electron mass difference of ions is ignored.

#note[Compilation fails with a descriptive error when the weight is undefined:
wildcard #c("*") atoms, #c("{label}") abbreviations (no defined composition),
and isotope-labeled atoms such as #c("[2H]") (a specific nuclide needs a
nuclide mass, not a standard atomic weight).]

// ═════════════════════════════════════════════════════════════════════════════
= Inline molecules: #raw("smiles-inline()")
// ═════════════════════════════════════════════════════════════════════════════

#demo[
  #c("smiles-inline(smiles-str, height: 1.4em, baseline: auto, ..args)") scales
  a structure to a target height and baseline-aligns it so it reads inline,
  like a word in the sentence. Extra named arguments are passed straight to
  #c("smiles()").

  #example(```typ
  Dehydrating ethanol #smiles-inline("CCO")
  over alumina gives ethylene
  #smiles-inline("C=C"); toluene
  #smiles-inline("Cc1ccccc1") is a common
  solvent.
  ```)
]

Sizing works in two steps: the drawing is scaled to fit the #c("height")
(default #c("1.4em"), so ordinary line spacing stays essentially unchanged),
and the scale is capped so one bond never exceeds most of that height —
otherwise a flat, mostly horizontal molecule, which measures almost no height,
would blow up to fill it. Typst grows a line to fit tall inline content, so
neighboring lines are never overlapped; passing a larger #c("height") only
makes the molecule's own line taller.

#note[#c("baseline") sets how far the drawing's vertical center sits above the
text baseline (default centers it on the lowercase body). Atom labels scale
down with the drawing: at inline sizes, label-heavy molecules read better with
a larger #c("height") or as a display-mode #c("smiles()").]

// ═════════════════════════════════════════════════════════════════════════════
= CeTZ integration: #raw("smiles-cetz()")
// ═════════════════════════════════════════════════════════════════════════════

#demo[
  #c("smiles-cetz(smiles-str, name:, origin:, fg:, theme:, ..opts)") draws a
  molecule *inside an existing CeTZ canvas* and registers named anchors on the
  group, so arbitrary CeTZ drawing attaches to real molecular positions:

  - #c("\"<name>.atom-<i>\"") — every atom center (writing-order indices, as
    shown by #c("show-indices")),
  - #c("\"<name>.bond-<i>-<j>\"") — every bond midpoint (with $i < j$),
  - #c("\"<name>.center\"") — the molecule origin.

  Coordinates are in bond-length units; give the canvas
  #c("length: 30pt * scale") so sizes match #c("smiles(scale: ...)"), and wrap
  the canvas in #c("context") (atom-label measurement needs it — plain
  #c("smiles()") provides it automatically, a raw canvas does not).

  The Watson–Crick A–T base pair below is two #c("smiles-cetz") calls plus
  plain CeTZ: two dashed hydrogen bonds drawn between atom anchors of the two
  bases, endpoint nudges so the dashed lines do not start at atom centers,
  distance labels, and glycosidic-bond stubs to the sugar backbone.

  #example(```typ
  #align(center, context cetz.canvas(length: 30pt, {
  import cetz.draw: *

  smiles-cetz("Nc1ncnc2N(!s{})cnc12", name: "A")
  smiles-cetz("Cc1cN(!s{})c(=O)[nH]c1=O", name: "T", origin: (4.9, 0.42))

  let hb = (paint: rgb("#3A78C9"), thickness: 1.0pt, dash: "densely-dashed")
  let off(anchor, by) = (rel: by, to: anchor)
  line(off("A.atom-11", (0.4, -0.15)), off("T.atom-9", (-0.2, 0.06)), stroke: hb)
  line(off("A.atom-2", (0.15, 0)), off("T.atom-7", (-0.2, 0)), stroke: hb)

  // Hydrogen-bond distances.
  content((rel: (0.2, 0.2), to: ("A.atom-11", 50%, "T.atom-9")), text(size: 7.5pt, fill: rgb("#3A78C9"))[2.9 Å])
  content((rel: (0, 0.28), to: ("A.atom-2", 50%, "T.atom-7")), text(size: 7.5pt, fill: rgb("#3A78C9"))[2.8 Å])

  }))
  ```, side: false)
]

The remaining #c("..opts") are #c("smiles()") drawing options (#c("color"),
#c("rotation"), #c("mirror"), #c("show-h"), #c("aromatic"), #c("opacity"),
#c("bond-customizations"), …). #c("fg") and #c("theme") default to black and
light because `auto` foreground resolution is not available inside a raw
canvas — pass them explicitly on dark backgrounds.

Use regular CeTZ relative coordinates to offset any connection endpoint:
#c("(rel: (0.12, 0.04), to: \"A.atom-0\")"). This is the same idea as
#c("atom(..., offset:)") in mechanism arrows, but in raw CeTZ syntax.

#note[Turn on #c("show-indices: true") while writing anchor references, exactly
as for mechanism arrows, to read off the atom indices. Bond-midpoint anchors
(#c("bond-<i>-<j>")) attach to the middle of a bond — useful for marking a
reacting or attachment bond.]

// ═════════════════════════════════════════════════════════════════════════════
= Reaction schemes
// ═════════════════════════════════════════════════════════════════════════════

Three helpers compose molecules, formulas, and arrows into schemes.

== #raw("mol()")

#demo[
  #c("mol(spec, label: none, offset: (0,0), ..opts)") is a reaction item. #c("spec")
  is either any content (#c("smiles(...)"), #c("ce(...)"), text) or a SMILES *string*.
  Passing a string lets #c("reaction") render the molecule itself so its atoms become
  addressable by curly arrows (see #link(<sec-mech>)[Reaction mechanisms]). String
  molecules accept common drawing options such as #c("scale"), #c("font-size"),
  #c("font"), #c("bond-stroke"), #c("color"), #c("rotation"), #c("show-h"),
  #c("lone-pairs"), #c("opacity"), #c("bond-customizations"), #c("atom-colors"),
  and #c("show-indices"). #c("offset") nudges the molecule in page coordinates, in
  bond-length units: positive x moves right and positive y moves up regardless
  of #c("reaction(flow:)"). In mechanism mode, #c("reaction(scale: ...)") sets
  the shared canvas scale (and its offsets), and each #c("mol(scale: ...)")
  multiplies it for that molecule alone. A #c("mol()") item can also be placed in
  a #c("rxn-arrow(above:/below:)") slot to draw a molecule over or under the
  arrow. In ordinary schemes, an offset also changes the reaction's layout
  bounds, so #c("brackets()") encloses the shifted edge.

  #example(```typ
  #reaction(mol(smiles("CCO"), label: [ethanol]))
  ```)
]

== #raw("rxn-arrow()")

#c("rxn-arrow(above: none, below: none, dir: auto, kind: \"single\", scale: 1.0)") draws an arrow.
#c("dir") may be #c("\"right\""), #c("\"left\""), #c("\"up\""), #c("\"down\""), or
#c("auto") (the default, following the enclosing #c("reaction(flow:)")).
#c("kind") may be #c("\"single\""), #c("\"equilibrium\""), #c("\"equilibrium-filled\""),
#c("\"dashed\""), or #c("\"wavy\"").
#c("above") and #c("below") carry reagents and conditions — any content, or a
#c("mol()") item to draw a molecule in the label slot (in mechanism mode a
#c("mol()") there also becomes a referenceable species; see
#link(<sec-mech>)[Reaction mechanisms]).
#c("scale") uniformly resizes the complete arrow, including its condition labels.
#c("color") defaults to #c("auto"), inheriting the surrounding text color
(so arrows recolor on dark slides).

#table(
  columns: (auto, auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Argument*], [*Default*], [*Effect*],
  [#c("above")], [`none`], [Label above a horizontal arrow, or right of a vertical arrow. Content or a #c("mol()") item.],
  [#c("below")], [`none`], [Label below a horizontal arrow, or left of a vertical arrow. Content or a #c("mol()") item.],
  [#c("dir")], [`auto`], [Arrow direction: #c("\"right\""), #c("\"left\""), #c("\"up\""), #c("\"down\""), or #c("auto") (follows #c("reaction(flow:)")).],
  [#c("kind")], [`"single"`], [Arrow style: #c("\"single\""), #c("\"equilibrium\""), #c("\"equilibrium-filled\""), #c("\"dashed\""), or #c("\"wavy\"").],
  [#c("scale")], [`1.0`], [Uniform arrow scale, including condition labels.],
  [#c("color")], [`auto`], [Arrow color; `auto` inherits the surrounding text color.],
)

#demo[
  Scale an individual arrow without changing the rest of its reaction. The
  shaft, arrowhead, spacing, and condition labels resize together.

  #example(```typ
  #reaction(
    ce("A"),
    rxn-arrow(scale: 1.4, above: [cat.]),
    ce("B"),
  )
  ```, side: false)
]

#demo[
  An arrow label slot also accepts a #c("mol()") item, drawing a molecule above
  or below the arrow with its own per-molecule options.

  #example(```typ
  #reaction(
    mol("C=C"),
    rxn-arrow(above: mol("BrBr", scale: 0.7), below: [CCl#sub[4]]),
    mol("BrC(Br)C"),
  )
  ```, side: false)
]

#demo[
  Equilibrium arrows use paired half-heads. The filled variant uses filled
  half-heads.

  #example(```typ
  #reaction(
    ce("A"),
    rxn-arrow(kind: "equilibrium", above: ce("H+"), below: [heat]),
    ce("B"),
    rxn-arrow(kind: "equilibrium-filled", above: [cat.]),
    ce("C"),
  )
  ```, side: false)
]

#demo[
  A dashed arrow marks a hypothetical or formal transformation; a wavy arrow
  marks a distorted or non-elementary one.

  #example(```typ
  #reaction(
    ce("A"),
    rxn-arrow(kind: "dashed", above: [formal]),
    ce("B"),
    rxn-arrow(kind: "wavy", above: [hv]),
    ce("C"),
  )
  ```, side: false)
]

== #raw("reaction()")

#demo[
  `reaction(gap-h: 1.5em, gap-v: 1.5em, scale: 1.0, breakable: false, show-indices: false, flow: "right", ..items)`
  lays out molecules and arrows left to right. An up or down arrow wraps the
  scheme onto a new row.

  #example(```typ
  #reaction(
    mol(smiles("C1=CC=CC=C1"), label: [benzene]),
    rxn-arrow(above: ce("Br2"), below: ce("FeBr3")),
    mol(smiles("BrC1=CC=CC=C1"), label: [bromobenzene]),
  )
  ```, side: false)
]

#note[
  Set #c("reaction(show-indices: true)") to stamp indices on every string SMILES
  molecule rendered through #c("mol(\"...\")") in that reaction. This is useful
  while authoring large mechanisms with many arrows or highlights. A molecule can
  still opt out with #c("mol(\"...\", show-indices: false)"). Molecules already
  pre-rendered as content, such as #c("mol(smiles(\"...\"))"), keep their own
  #c("smiles()") options.
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

=== Flow direction

#demo[
  #c("flow") sets the writing direction: #c("\"right\"") (default), #c("\"left\""),
  #c("\"up\""), or #c("\"down\""). #c("\"left\"") and #c("\"up\"") reflect the
  scheme along the flow axis and flip the arrowheads, so a scheme written in
  reading order displays right-to-left or bottom-to-top. Ordinary non-arrow
  items follow the same direction, so #c("flow: \"up\"") and #c("\"down\"")
  stack molecules, plus signs, and text vertically even when no #c("rxn-arrow()")
  separators are present. This is the natural choice for a branch that grows out
  of the left or bottom of a #c("cycle()") — its released product ends up on the
  side facing the cycle. Arrows created with #c("rxn-arrow(dir: auto)") (the
  default) follow the flow; an explicit #c("dir") is kept.

  #example(```typ
  #reaction(flow: "left",
    mol(smiles("CCO")),
    rxn-arrow(above: [\[O\]]),
    mol(smiles("CC=O")),
    rxn-arrow(above: [\[O\]]),
    mol(smiles("CC(=O)O")),
  )
  ```, side: false)
]

#note[Because #c("step(out:)") of a #c("cycle()") accepts any content, a nested
#c("reaction(flow: \"left\")") (or #c("\"down\"")) grows a full sub-reaction out
of a released molecule that reads naturally away from the ring, while the outer
reaction continues in its own direction. Flow applies to the scheme layout;
mechanism mode (curly arrows) is unaffected.]

=== Uniform scale

#demo[
  Pass `scale` to shrink or enlarge the entire scheme uniformly — molecules,
  arrows, labels, and separators all change together. This is especially useful
  for multi-step or wrap-around schemes that would otherwise be too wide.

  #example(```typ
  // same scheme at three different scales
  #reaction(
    scale: 0.75,
    mol(smiles("C1=CC=CC=C1"), label: text(size: 7pt)[1]),
    rxn-arrow(above: ce("Br2"), below: ce("FeBr3")),
    mol(smiles("BrC1=CC=CC=C1"), label: text(size: 7pt)[A]),
    rxn-arrow(dir: "down", above: ce("HNO3"), below: ce("H2SO4")),
    mol(smiles("BrC1=CC(=CC=C1)[N+](=O)[O-]"), label: text(size: 7pt)[B]),
    rxn-arrow(dir: "left", above: ce("Fe"), below: ce("HCl")),
    mol(smiles("BrC1=CC(=CC=C1)N"), label: text(size: 7pt)[C]),
  )
  ```, side: false)
]

#note[
  `reaction(scale: ...)` and `smiles(scale: ...)` are independent. The reaction
  scale is applied on top of however each individual molecule is sized. To
  resize everything from a single place, use `reaction(scale: ...)` and let each
  `smiles()` call use its default.
]

=== Page-break behaviour

#demo[
  By default, #c("reaction") sets #c("breakable: false"), so the whole scheme
  moves to the next page as a unit if it does not fit on the current one. This
  prevents a molecule or vertical branch from stranding on a different page from
  the rest of the scheme. Set #c("breakable: true") only for very long schemes
  that must span pages.

  ```typ
  // Default — the scheme always stays on one page:
  #reaction( … )                      // breakable: false

  // Opt in to page splitting for very long schemes:
  #reaction(breakable: true, … )
  ```
]

// ═════════════════════════════════════════════════════════════════════════════
= Reaction mechanisms <sec-mech>
// ═════════════════════════════════════════════════════════════════════════════

#c("reaction()") also draws electron-pushing mechanisms: curly arrows that flow
from a lone pair, bond, or atom into another bond or atom — within one structure
or across separate species. It switches from grid layout to a single shared canvas
as soon as any curly #c("arrow()") or #c("highlight()") is present.
#c("mol(offset:)") by itself stays in scheme layout and nudges the species
inside its grid cell.

== Referencing atoms

Atoms are addressed by their *writing-order index* (0-based), so the SMILES string is
never modified. Inside #c("smiles()") use single-index references; inside
#c("reaction()") prefix them with the species index #c("s") — the position of the
#c("mol()")/content item in written order. #c("mol()") items inside
#c("rxn-arrow(above:/below:)") slots are counted too, at the arrow's position in
the sequence (#c("above") before #c("below")); the arrows themselves, plain
arrow-label content, and annotations are not counted.

#table(
  columns: (auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Reference*], [*Resolves to*],
  [#c("atom(i)") / #c("atom(s, i)")], [Atom center.],
  [#c("bond(i, j)") / #c("bond(s, i, j)")], [Midpoint of the bond between atoms i and j.],
  [#c("lp(i)") / #c("lp(s, i)")], [A lone pair on atom i (use #c("pair: n") to pick one).],
  [#c("species(k)")], [Bounding-box edge of a whole item (e.g. a #c("ce()") formula).],
)

#note[
  Every reference accepts an #c("offset: (dx, dy)") nudge in bond-length units — the
  escape hatch for fine-tuning an endpoint or pointing at a single-electron dot.
  Set #c("show-indices: true") on a #c("smiles()") / #c("mol()") to stamp the indices
  on the diagram while you write the arrows.

  Hydrogen atoms declared inside bracket notation (e.g. #c("[OH-]"), #c("[NH4+]"))
  are also addressable. Each H gets the next available index after all heavy atoms,
  so #c("atom(1)") in #c("[OH-]") refers to the hydrogen. With #c("show-indices"),
  the heavy-atom and H badges are centered on their rendered label glyphs. Charge
  marks are excluded from the atom center, so #c("atom(0)") in #c("[O-]") points
  to the O glyph rather than the combined #c("O-") label.
]

== #raw("arrow()") and #raw("highlight()")

#demo[
  #c("arrow(from:, to:, label: none, color: red, bend: \"left\", angle: 15deg, half: false, heads: \"end\", style: \"solid\")")
  draws a curly arrow between two references. #c("bend") is #c("\"left\""),
  #c("\"right\""), or #c("none"); #c("angle") sets how strongly it bows;
  #c("half: true") draws fishhook (single-electron) heads. #c("heads") selects
  which ends carry an arrowhead — #c("\"end\"") (default), #c("\"both\""), or
  #c("\"none\"") — and #c("style") the shaft: #c("\"solid\"") (default),
  #c("\"dashed\""), or #c("\"wavy\""); both combine freely with #c("bend"). #c("highlight(ref, fill:)")
  shades an atom (disk) or bond (capsule) behind the structure. Pass an array of
  references to shade several atoms or bonds with one call. Bond highlights are
  trimmed away from atom centers by default; set #c("include-atoms: true") to also
  shade the endpoint atoms of each highlighted bond. Abbreviated group labels such
  as #c("{PPh3}") are highlighted as measured label-width capsules rather than
  atom-sized disks.

  #example(```typ
  #smiles(
    "N1CCN(CC1)C(C(F)=C2)=CC(=C2C4=O)N(C3CC3)C=C4C(=O)O",
    highlight((bond(0, 5), bond(5, 4), bond(4, 3), bond(3, 6), bond(6, 10), bond(10, 11), bond(11, 15), bond(15, 19), bond(19, 20), bond(20, 21), bond(21, 23), bond(23, 25)), fill: rgb(150, 191, 13), include-atoms: true),
    highlight((bond(15, 16), bond(16, 18), bond(18, 17), bond(17, 16)), fill: rgb(242, 148, 1), include-atoms: true),
    highlight((bond(3, 2), bond(2, 1), bond(1, 0)), fill: rgb(137, 199, 168), include-atoms: true),
    highlight((bond(6, 7), bond(7, 8), bond(7, 9), bond(9, 12)), fill: rgb(201, 143, 75), include-atoms: true),
    highlight((bond(11, 12), bond(12, 13), bond(13, 20), bond(13, 14)), fill: rgb(236, 119, 137), include-atoms: true),
    highlight((bond(21, 22)), fill: rgb(0, 134, 203), include-atoms: true),
    color: false,
    rotation: 90deg,
    bond-stroke: 0.8pt,
    scale: 0.5,
  )
  #smiles(
    "CC(=O)O",
    highlight((atom(0), atom(2), atom(3), atom(4)), fill: rgb("#BBE1FA")),
  )
  #smiles(
    "CC(=O)OC",
    highlight((bond(2, 1)), include-atoms:true, fill: rgb("#FFCAD4")),
  )
  #smiles(
    "CCOCC",
    highlight(
      (bond(0, 1), bond(1, 2), bond(2, 3)),
      fill: rgb("#BBE1FA"),
      include-atoms: true,
    ),
  )
  ```)
]

#demo[
  #c("heads") and #c("style") tune the arrow shape — for instance a
  double-headed dashed arrow for an interaction, or a wavy arrow for a
  photochemical or homolytic step.

  #example(```typ
  #smiles("CC(=O)O",
    arrow(from: atom(3), to: atom(2), heads: "both", color: blue))
  #smiles("CC(=O)O",
    arrow(from: atom(3), to: atom(2), style: "dashed"))
  #smiles("OCCCCO",
    arrow(from: atom(0), to: atom(5), style: "wavy", heads: "both", half: true))
  ```)
]

== A mechanism across species

#demo[
  Hydroxide attacks the central carbon; the C–I bond breaks toward iodide. The
  nucleophile is offset so the curly arrow reads cleanly.

  #example(```typ
  #reaction(
    mol("[OH-]", lone-pairs: "dots", offset: (0, 1)),
    mol("C(I)(C)C"),
    arrow(from: lp(0, 0, offset:(0.1, -0.2)), to: atom(1, 0, offset : (0.1, -0.1)),
          bend: "right", color : black),
  )
  ```, side: false)
]

== Per-molecule scale in mechanisms

#demo[
  In mechanism mode #c("reaction(scale:)") sets the shared canvas scale and each
  #c("mol(scale:)") multiplies it for its own species. Curly arrows and
  highlights land on the scaled coordinates.

  #example(```typ
  #reaction(
    mol("CC(=O)O", scale: 0.7, label: [minor]),
    mol("CC(=O)O", scale: 1.2, label: [major]),
    highlight(atom(0, 3)),
    arrow(from: atom(0, 3), to: atom(1, 1), color: blue),
  )
  ```, side: false)
]

== Referencing molecules above an arrow

#demo[
  A #c("mol()") placed in a #c("rxn-arrow(above:/below:)") slot is lifted onto
  the shared canvas as a species of its own, so its atoms take part in the
  mechanism like any other species.

  #example(```typ
  #reaction(
    mol("C=C"),
    rxn-arrow(above: mol("BrBr", scale: 0.8), below: [CCl#sub[4]]),
    mol("BrC(Br)C"),
    arrow(from: bond(0, 0, 1), to: atom(1, 0), color: red),
    arrow(from: atom(1, 0), to: atom(1, 1), color: red, bend: "right"),
  )
  ```, side: false)
]

== #raw("brackets()")

#demo[
  #c("brackets(body, sup: none, sub: none)") encloses any content in square
  brackets including another reaction, with optional marks at the corners — a #c("[‡]") for a transition state
  or a charge for a reactive intermediate.

  #example(```typ
  #brackets(
    [#reaction(smiles("CC(=O)C"), rxn-arrow(), smiles("O=C=O"), scale : 0.7)],
    sup: [‡])
  ```)
]

// ═════════════════════════════════════════════════════════════════════════════
= Catalytic cycles: #raw("cycle()")
// ═════════════════════════════════════════════════════════════════════════════

  #c("cycle(..items, radius:, start:, clockwise:, scale:, reagent-bend:, ...)")
  arranges
  species on a circle with arc arrows between them — the standard depiction of a
  catalytic cycle. Items alternate species and step()s, just as #c("reaction()")
  alternates molecules and arrows, but the sequence closes into a ring (the last
  step returns to the first species). Species are mol() items or any content;
  SMILES strings are rendered by #c("smiles()").

  #c("step(label:, into:, out:, bend:, merge:, rotation:, label-offset:, into-offset:, out-offset:)")
  describes the arc that follows a species: #c("label") names the transformation,
  #c("into") adds a reagent merging into the arc from outside the ring, and
  #c("out") a product leaving it. Those side arrows are bookkeeping arrows for
  material entering or leaving the catalytic manifold, not electron-pushing
  arrows. #c("merge: true") draws a side arrow tangent to the arc so it visually
  fuses with the main cycle arrow. #c("rotation") rotates the step label:
  #c("\"straight\"") keeps it horizontal, #c("\"auto\"") follows the step's
  circle angle while keeping text upright, and an explicit angle is used as given. The three
  #c("*-offset") arguments nudge the label and the into/out reagents like a
  #c("mol") offset.


  #example(```typ
  #align(center)[
  #reaction(
    scale : 0.6,
    flow : "down",
    mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)(!w{>PPh3})"), offset: (0.5, 0)),
    rxn-arrow(above : [#ce("PPh3") solvent (S)]),
    cycle(
      radius : 6.0,
      mol(smiles("{Rh}(!w{S | S})(!h{>PPh3})(!hCl)(!w{>PPh3})")),
      step(label : [#text(size : 8pt)[oxidative addition]], label-offset : (0.3, 0.5), into : [#ce("H2")], rotation : "auto", merge : true),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({S | S})({H})(!w{H})")),
      step(into : [#smiles("C=C")], merge : true),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({})({H})({H})")),
      step(label : [#text(size : 8pt)[migratory insertion]], label-offset: (0, 0.5), rotation : "auto"),
      mol(smiles("{Rh}(!w{>PPh3})(!h{>PPh3})(!hCl)({S | S})({H})(!wCC{H})", rotation : -90deg)),
      step(label : [#text(size : 8pt)[reductive elimination]], label-offset:  (-0.5, 0.5), out : [#smiles("{H}CC{H}")], merge:true, rotation : "auto")
    )
  )
]
  ```, side: false)


Because #c("step(out:)") accepts any content, passing a #c("reaction()") there
grows a full linear branch out of a released molecule — the released product
becomes the first species of a downstream sequence. Reagents (both #c("into")
and #c("out")) are offset by their own measured size, so even a wide nested
reaction clears the ring. Incoming content is attached on the side facing the
cycle so the side arrow starts at the last upstream item; outgoing content uses
the same side so the leaving arrow lands at the first downstream item instead
of the center of the whole branch.

#note[#c("radius") defaults to #c("auto"), which fits the species without
overlap; set it explicitly for a tighter or looser ring. #c("start") is the
angle of the first species (#c("90deg") is the top) and #c("clockwise") the
travel direction. #c("scale") sets the canvas bond length as in
#c("reaction()"). #c("arc-gap") tunes how close the arc arrows sit to the
species (in bond-length units beyond each molecule's radius; smaller or even
negative brings them closer). #c("reagent-bend") sets the default curvature of
into/out side arrows; #c("step(bend:)") overrides one step, and #c("bend: 0")
makes that side arrow nearly radial. A one-species #c("cycle()") draws that
species alone.]

#demo[
  With #c("merge: true") the into/out arrows fuse tangentially with the main
  cycle arrows, and #c("arc-gap") pulls the arrows closer to the species:

  #example(```typ
  #cycle(
    scale: 0.9,
    arc-gap: 0.0,
    mol(box(inset: 2pt)[E]),
    step(label: [binding], into: [S], merge: true),
    mol(box(inset: 2pt)[E·S]),
    step(label: [turnover], out: [P], merge: true),
  )
  ```, side: false)
]

// ═════════════════════════════════════════════════════════════════════════════
= Project-wide defaults
// ═════════════════════════════════════════════════════════════════════════════

Typst functions support partial application with `.with()`. Calling
#c("smiles.with(...)") returns a new function with the given arguments
pre-filled. Rebind the name in your preamble and every subsequent
#c("#smiles(...)") call uses those defaults — including inside #c("#reaction()"),
because the content is evaluated before #c("reaction") sees it.

Any argument can be pre-filled this way: #c("bond-length"), #c("font"),
#c("scale"), #c("font-size"), #c("atom-colors"), or any other #c("smiles")
parameter.

```typ
// ── preamble ───────────────────────────────────────────────────────────
#import "@preview/typed-smiles:0.7.0": *

#let smiles = smiles.with(
  bond-length: 0.9,
  font:        "Libertinus Serif",
  atom-colors: (O: rgb("#8B4513"), N: rgb("#008080")),
)

// ── anywhere in the document ───────────────────────────────────────────
#smiles("CC(N)C(=O)O")   // uses bond-length 0.9, Libertinus, custom O/N
#reaction(
  smiles("CCO"),
  rxn-arrow(),
  smiles("CC(=O)O"),     // same custom smiles — preamble applies here too
)
```

#note[Any argument passed directly at the call site still overrides the
`.with()` default for that call alone.]

// ═════════════════════════════════════════════════════════════════════════════
= Limitations <sec-limits>
// ═════════════════════════════════════════════════════════════════════════════

- Directional bonds (#c("/"), #c("\\")) draw the correct cis/trans geometry, but
  E/Z descriptors are not computed.
- R/S and E/Z descriptors are not calculated.
- Trigonal-bipyramidal (#c("@TB")), octahedral (#c("@OH")), and allenal
  (#c("@AL")) centers are accepted but drawn without stereo wedges.
- Ring stereochemistry between adjacent centers and bridged bicyclics can overlap
  or need a manual adjustment (try #c("rotation"), or the manual #c("!w") and
  #c("!h") wedges). Ordinary substituent branches resolve crowding automatically —
  a blocking branch flips to the other side of its attachment bond (e.g. the
  ortho acetyl and carboxyl groups of aspirin), or the colliding bond flips
  its zigzag turn — always keeping ideal bond angles, so plain branches do
  not overlap.

// ═════════════════════════════════════════════════════════════════════════════
= Quick reference
// ═════════════════════════════════════════════════════════════════════════════

== Syntax

#table(
  columns: (auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Syntax*], [*Meaning*],
  [#c("c") #c("n") #c("o") …], [Aromatic atom (lowercase); kekulized on parse.],
  [#c(".")], [Disconnected fragment (no bond): salts, hydrates.],
  [#c("[...]")], [Bracket atom (any element, charges, explicit H).],
  [#c("@") / #c("@@")], [Tetrahedral anticlockwise / clockwise.],
  [#c("@SP1")–#c("@SP3")], [Square planar (U / 4 / Z shape), depicted exactly.],
  [#c("@TB") #c("@OH") #c("@AL")], [Extended stereo: accepted, drawn without wedges.],
  [#c("$")], [Quadruple bond (four parallel lines).],
  [#c("/") #c("\\")], [Double-bond cis/trans geometry.],
  [#c("!w") / #c("!h")], [Manual solid / hashed wedge (typed-smiles extension).],
  [#c("!s") / #c("!d")], [Wavy / dashed bond (typed-smiles extension).],
  [#c("{label}")], [Abbreviated group pseudo-atom.],
  [#c("{>label}")], [Abbreviated group anchored at the glyph after #c(">").],
  [#c("{label|style}")], [Colored abbreviated group.],
)

== #raw("smiles()") options

#table(
  columns: (auto, auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Option*], [*Default*], [*Effect*],
  [#c("style")], [`"default"`], [Journal preset: #c("\"acs\""), #c("\"rsc\""), #c("\"nature\""), #c("\"wiley\"").],
  [#c("scale")], [`1.0`], [Balanced size of everything.],
  [#c("bond-length")], [`none`], [Bond length only (`1.0` is 30 pt).],
  [#c("font-size")], [`none`], [Atom-label size only.],
  [#c("font")], [`"New Computer Modern"`], [Atom-label font.],
  [#c("bond-stroke")], [`none`], [Bond width only.],
  [#c("color")], [`auto`], [CPK colors for #c("\"default\""); monochrome for journal presets.],
  [#c("fg")], [`auto`], [Foreground for bonds/carbon labels; `auto` inherits the text color.],
  [#c("theme")], [`auto`], [CPK palette variant; `auto` goes dark when #c("fg") is light.],
  [#c("rotation")], [`0deg`], [Rotate, labels stay upright.],
  [#c("mirror")], [`none`], [Optional #c("\"horizontal\"") or #c("\"vertical\"") reflection.],
  [#c("show-h")], [`()`], [Label selected implicit hydrogens; use #c("\"all\"") for every atom.],
  [#c("aromatic")], [`"kekule"`], [Lowercase-aromatic rings as doubles or #c("\"circle\"").],
  [#c("atom-annotations")], [`()`], [Small gray side labels as #c("(index, content)") or #c("(index, content, offset)") tuples.],
  [#c("opacity")], [`100%`], [Fade the whole drawing (ghost molecules).],
  [#c("bond-customizations")], [`()`], [Per-bond #c("color"), #c("stroke"), #c("opacity") keyed by #c("bond(i, j)").],
  [#c("lone-pairs")], [`none`], [Draw lone pairs as #c("\"dots\"") or #c("\"lines\"").],
  [#c("atom-colors")], [`(:)`], [Color overrides: #c("O: red") for elements, #c("\"{PPh3}\": blue") for labels.],
  [#c("show-indices")], [`false`], [Stamp atom indices for writing references.],
  [#c("..annotations")], [—], [#c("arrow()") / #c("highlight()") items on this molecule.],
)

== #raw("reaction()") options

#table(
  columns: (auto, auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Option*], [*Default*], [*Effect*],
  [#c("gap-h")], [`1.5em`], [Horizontal gap between items.],
  [#c("gap-v")], [`1.5em`], [Vertical gap between rows.],
  [#c("scale")], [`1.0`], [Uniform scale applied to the whole scheme.],
  [#c("breakable")], [`false`], [Allow the scheme to split across pages.],
  [#c("flow")], [`"right"`], [Writing direction: #c("\"right\""), #c("\"left\""), #c("\"up\""), or #c("\"down\"") (left/up reflect the scheme).],
  [#c("show-indices")], [`false`], [Default atom-index overlay for string SMILES molecules in this reaction.],
)

== #raw("rxn-arrow()") options

#table(
  columns: (auto, auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Option*], [*Default*], [*Effect*],
  [#c("above")], [`none`], [Condition label (content or #c("mol()")) above a horizontal arrow, or right of a vertical arrow.],
  [#c("below")], [`none`], [Condition label (content or #c("mol()")) below a horizontal arrow, or left of a vertical arrow.],
  [#c("dir")], [`auto`], [Arrow direction: #c("\"right\""), #c("\"left\""), #c("\"up\""), #c("\"down\""), or #c("auto") (follows #c("reaction(flow:)")).],
  [#c("kind")], [`"single"`], [Arrow style: #c("\"single\""), #c("\"equilibrium\""), #c("\"equilibrium-filled\""), #c("\"dashed\""), or #c("\"wavy\"").],
  [#c("scale")], [`1.0`], [Uniform arrow scale, including condition labels.],
  [#c("color")], [`auto`], [Arrow color; `auto` inherits the surrounding text color.],
)

== Data helpers

#table(
  columns: (auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Helper*], [*Purpose*],
  [#c("mol-weight(smiles)")], [Molecular weight in g/mol as a #c("float"); errors on wildcards, abbreviations, and isotopes.],
  [#c("smiles-inline(smiles, height:, baseline:, ..args)")], [Molecule scaled and baseline-aligned for running text.],
  [#c("smiles-cetz(smiles, name:, origin:, fg:, theme:, ..opts)")], [Molecule as CeTZ elements with #c("atom-<i>") / #c("bond-<i>-<j>") / #c("center") anchors.],
)

== Mechanism helpers

#table(
  columns: (auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Helper*], [*Purpose*],
  [#c("mol(spec, label:, offset:, ..opts)")], [Reaction item; #c("spec") as a string is rendered with addressable atoms; #c("offset") nudges in page coordinates; per-molecule #c("scale") also works in mechanism mode.],
  [#c("cycle(..items, radius:, start:, clockwise:, scale:, reagent-bend:, arc-gap:)")], [Catalytic cycle: species on a ring with arc arrows.],
  [#c("step(label:, into:, out:, bend:, merge:, rotation:, *-offset:)")], [A cycle arc: transformation label, reagent in/out, side-arrow bend, tangential merge, label rotation, and per-piece offsets.],
  [#c("atom / bond / lp / species")], [Atom-index references (optional #c("offset:")).],
  [#c("arrow(from:, to:, label:, color:, bend:, angle:, half:, heads:, style:)")], [Curly electron-pushing arrow; #c("heads") is #c("\"end\"") / #c("\"both\"") / #c("\"none\""), #c("style") is #c("\"solid\"") / #c("\"dashed\"") / #c("\"wavy\"").],
  [#c("highlight(ref, fill:, stroke:, radius:, include-atoms:)")], [Shade one atom/bond reference or an array of references.],
  [#c("brackets(body, sup:, sub:)")], [Square brackets with optional corner marks.],
)

== Label color names

Named styles for #c("{label|style}") — any of these plus any `#RRGGBB` hex:

#table(
  columns: (auto, auto) * 3, inset: (x: 8pt, y: 5pt),
  align: (x, _) => if calc.even(x) { left } else { center + horizon },
  stroke: 0.5pt + luma(210),
  fill: (x, y) => if y == 0 { accent-soft },
  [*Name*], [*Swatch*], [*Name*], [*Swatch*], [*Name*], [*Swatch*],
  [red],    [#box(width: 2em, height: 0.75em, fill: rgb("#FF0D0D"), radius: 2pt)],
  [orange], [#box(width: 2em, height: 0.75em, fill: rgb("#FF8000"), radius: 2pt)],
  [yellow], [#box(width: 2em, height: 0.75em, fill: rgb("#E6C800"), radius: 2pt)],
  [brown],  [#box(width: 2em, height: 0.75em, fill: rgb("#8B4513"), radius: 2pt)],
  [green],  [#box(width: 2em, height: 0.75em, fill: rgb("#1FA51F"), radius: 2pt)],
  [lime],   [#box(width: 2em, height: 0.75em, fill: rgb("#32CD32"), radius: 2pt)],
  [teal],   [#box(width: 2em, height: 0.75em, fill: rgb("#008080"), radius: 2pt)],
  [cyan],   [#box(width: 2em, height: 0.75em, fill: rgb("#00B4D8"), radius: 2pt)],
  [blue],   [#box(width: 2em, height: 0.75em, fill: rgb("#3050F8"), radius: 2pt)],
  [navy],   [#box(width: 2em, height: 0.75em, fill: rgb("#000080"), radius: 2pt)],
  [purple], [#box(width: 2em, height: 0.75em, fill: rgb("#940094"), radius: 2pt)],
  [pink],   [#box(width: 2em, height: 0.75em, fill: rgb("#FF69B4"), radius: 2pt)],
  [black],  [#box(width: 2em, height: 0.75em, fill: black,          radius: 2pt)],
  [gray],   [#box(width: 2em, height: 0.75em, fill: rgb("#777777"), radius: 2pt)],
  [silver], [#box(width: 2em, height: 0.75em, fill: rgb("#C0C0C0"), radius: 2pt)],
  [maroon], [#box(width: 2em, height: 0.75em, fill: rgb("#800000"), radius: 2pt)],
  [white],  [#box(width: 2em, height: 0.75em, fill: white, stroke: 0.4pt + luma(200), radius: 2pt)],
  [],       [],
)
