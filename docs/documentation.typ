// typed-smiles documentation
//
// Compile with:
//   typst compile --root . docs/documentation.typ docs/documentation.pdf

#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#import "../src/lib.typ": smiles, ce, rxn-arrow, mol, reaction, atom, bond, lp, species, arrow, highlight, brackets

#let version = "0.4.2"
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
  smiles: smiles, ce: ce, rxn-arrow: rxn-arrow, mol: mol, reaction: reaction,
  atom: atom, bond: bond, lp: lp, species: species, arrow: arrow,
  highlight: highlight, brackets: brackets,
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
#import "@preview/typed-smiles:0.4.2": *
```

Or import only what you need:

```typ
#import "@preview/typed-smiles:0.4.2": smiles, ce, mol, rxn-arrow, reaction
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
    lone-pairs: none,
    atom-colors: (:),
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
  #smiles("OC1=CC(=CC=C1)C(=O)O") \
  #smiles("OC1=CC(=CC=C1)C(=O)O", color: false)
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
  #smiles("{PPh3|P}C({OEt|O})=O")
  // override by label name — inline style is ignored for these
  #smiles("{PPh3|P}C({OEt|O})=O",
    atom-colors: ("{PPh3}": rgb("#7B2D8B"), "{OEt}": rgb("#008080")))
  ```, side: false)
]

=== Mixing both key forms

#demo[
  Element keys and label keys can coexist in the same dictionary.

  #example(```typ
  // O (element) → teal; {Nu} (label) → navy; {LG} (label) → maroon
  #smiles("{Nu}!wC({LG|red})=O",
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
two highlighted groups, keep `color: true` (the default) and put everything else
in `atom-colors`.]

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
  #smiles("{PPh3|P}C({OEt|O})=O")
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
  Add #c("|style") inside the braces to color a label. Use a named color,
  an element symbol, or a hex code — see @sec-colors for the full list.

  #example(```typ
  #smiles("{Nu|blue}CC{LG|red}") \
  #smiles("{R|green}C(=O){OR|O}") \
  #smiles("{Cat|teal}C(=O){Nuc|#E040FB}")
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
  #c("mol(spec, label: none, offset: (0,0), ..opts)") is a reaction item. #c("spec")
  is either any content (#c("smiles(...)"), #c("ce(...)"), text) or a SMILES *string*.
  Passing a string lets #c("reaction") render the molecule itself so its atoms become
  addressable by curly arrows (see #link(<sec-mech>)[Reaction mechanisms]). String
  molecules accept common drawing options such as #c("font-size"), #c("font"),
  #c("bond-stroke"), #c("color"), #c("rotation"), #c("show-all-h"), #c("lone-pairs"),
  #c("atom-colors"), and #c("show-indices"). #c("offset") nudges the molecule in
  bond-length units and switches the reaction into mechanism mode; use
  #c("reaction(scale: ...)") to resize a shared mechanism canvas.

  #example(```typ
  #reaction(mol(smiles("CCO"), label: [ethanol]))
  ```)
]

== #raw("rxn-arrow()")

#c("rxn-arrow(above: none, below: none, dir: \"right\", kind: \"single\")") draws an arrow.
#c("dir") may be #c("\"right\""), #c("\"left\""), #c("\"up\""), or #c("\"down\"").
#c("kind") may be #c("\"single\""), #c("\"equilibrium\""), or #c("\"equilibrium-filled\"").
#c("above") and #c("below") carry reagents and conditions.

#table(
  columns: (auto, auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Argument*], [*Default*], [*Effect*],
  [#c("above")], [`none`], [Label above a horizontal arrow, or right of a vertical arrow.],
  [#c("below")], [`none`], [Label below a horizontal arrow, or left of a vertical arrow.],
  [#c("dir")], [`"right"`], [Arrow direction: #c("\"right\""), #c("\"left\""), #c("\"up\""), or #c("\"down\"").],
  [#c("kind")], [`"single"`], [Arrow style: #c("\"single\""), #c("\"equilibrium\""), or #c("\"equilibrium-filled\"").],
)

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

== #raw("reaction()")

#demo[
  `reaction(gap-h: 1.5em, gap-v: 1.5em, scale: 1.0, breakable: false, show-indices: false, ..items)`
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
as soon as any curly #c("arrow()") or #c("highlight()") is present, or any molecule
carries an #c("offset"). Schemes without these render exactly as before.

== Referencing atoms

Atoms are addressed by their *writing-order index* (0-based), so the SMILES string is
never modified. Inside #c("smiles()") use single-index references; inside
#c("reaction()") prefix them with the species index #c("s") — the position of the
#c("mol()")/content item in written order (#c("rxn-arrow")s and annotations are not
counted).

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
  #c("arrow(from:, to:, label: none, color: red, bend: \"left\", angle: 35deg, half: false)")
  draws a curly arrow between two references. #c("bend") is #c("\"left\""),
  #c("\"right\""), or #c("none"); #c("angle") sets how strongly it bows;
  #c("half: true") draws a fishhook (single-electron) head. #c("highlight(ref, fill:)")
  shades an atom (disk) or bond (capsule) behind the structure. Pass an array of
  references to shade several atoms or bonds with one call. Bond highlights are
  trimmed away from atom centers by default; set #c("include-atoms: true") to also
  shade the endpoint atoms of each highlighted bond. Abbreviated group labels such
  as #c("{PPh3}") are highlighted as measured label-width capsules rather than
  atom-sized disks.

  #example(```typ
  #smiles(
    "CC(=O)C",
    lone-pairs: "dots",
    highlight(bond(1, 2), fill: rgb("#FFE45C"), include-atoms:true),
  )
  #smiles(
    "CC(=O)O",
    highlight((atom(0), atom(2), atom(3), atom(4)), fill: rgb("#BBE1FA")),
  )
  #smiles(
    "CC(=O)OC",
    highlight((bond(0, 1), bond(2, 1)), include-atoms:true, fill: rgb("#FFCAD4")),
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

== A mechanism across species

#demo[
  Hydroxide attacks the central carbon; the C–I bond breaks toward iodide. The
  nucleophile is offset so the curly arrow reads cleanly.

  #example(```typ
  #reaction(
    mol("[OH-]", lone-pairs: "dots", offset: (1.5, 1)),
    mol("C(I)(C)C"),
    arrow(from: lp(0, 0), to: atom(1, 0),
          bend: "left"),
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
#import "@preview/typed-smiles:0.4.2": *

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
  [#c("show-indices")], [`false`], [Default atom-index overlay for string SMILES molecules in this reaction.],
)

== #raw("rxn-arrow()") options

#table(
  columns: (auto, auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Option*], [*Default*], [*Effect*],
  [#c("above")], [`none`], [Condition label above a horizontal arrow, or right of a vertical arrow.],
  [#c("below")], [`none`], [Condition label below a horizontal arrow, or left of a vertical arrow.],
  [#c("dir")], [`"right"`], [Arrow direction: #c("\"right\""), #c("\"left\""), #c("\"up\""), or #c("\"down\"").],
  [#c("kind")], [`"single"`], [Arrow style: #c("\"single\""), #c("\"equilibrium\""), or #c("\"equilibrium-filled\"").],
)

== Mechanism helpers

#table(
  columns: (auto, 1fr), inset: 6.5pt,
  align: (x, y) => if y == 0 { center + horizon } else { left + horizon },
  fill: (_, y) => if y == 0 { accent-soft }, stroke: 0.5pt + luma(210),
  [*Helper*], [*Purpose*],
  [#c("mol(spec, label:, offset:, ..opts)")], [Reaction item; #c("spec") as a string is rendered with addressable atoms.],
  [#c("atom / bond / lp / species")], [Atom-index references (optional #c("offset:")).],
  [#c("arrow(from:, to:, label:, color:, bend:, angle:, half:)")], [Curly electron-pushing arrow.],
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
