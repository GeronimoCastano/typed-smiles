// Typst SMILES Package
// Renders SMILES strings as 2D molecular structure diagrams via a WASM plugin.
// Also re-exports `ce` from typsium for chemical formula notation.

#import "@preview/cetz:0.3.4"
#import "@preview/chemformula:0.1.3": ch

#let smiles-plugin = plugin("../plugin/typst_smiles_plugin.wasm")

// Re-export as ce so users only need one import line.
// chemformula uses math mode internally, giving proper operator spacing.
#let ce = ch

// ── Internal helpers ──────────────────────────────────────────────────────────

#let _has-label(atom) = {
  let has-abbrev = atom.at("abbrev", default: "") != ""
  let has-hetero = (atom.symbol != "C" and atom.symbol != "c" and atom.symbol != "*") or (atom.charge != 0)
  has-abbrev or has-hetero
}

#let _atom-color(sym) = {
  if sym == "N" or sym == "n"      { rgb("#3050F8") }
  else if sym == "O" or sym == "o" { rgb("#FF0D0D") }
  else if sym == "S" or sym == "s" { rgb("#E6C800") }
  else if sym == "P"               { rgb("#FF8000") }
  else if sym == "F"               { rgb("#90E050") }
  else if sym == "Cl"              { rgb("#1FF01F") }
  else if sym == "Br"              { rgb("#A62929") }
  else if sym == "I"               { rgb("#940094") }
  else { black }
}

// ── SMILES renderer ───────────────────────────────────────────────────────────

/// Renders a SMILES string as a 2D skeletal molecular diagram.
///
/// - smiles-str (str): A valid SMILES string, e.g. "C1=CC=CC=C1" for benzene.
///   Use Kekulé notation for aromatic rings (C not c) until aromatic support lands.
/// - bond-length (float): Scale factor; 1.0 = 30 pt per bond. Default: 1.0.
/// - atom-font-size (length): Font size for heteroatom labels. Default: 11pt.
/// - color (bool): Apply Jmol CPK atom colors. Default: true.
/// - rotation (angle): Rotate the molecule by this angle. Atom labels stay upright.
///   Example: rotation: 90deg. Default: 0deg.
/// -> content
#let smiles(
  smiles-str,
  bond-length: 1.0,
  atom-font-size: 11pt,
  color: true,
  rotation: 0deg,
) = {
  let raw-bytes = smiles-plugin.layout(bytes(smiles-str))
  let layout = json(raw-bytes)

  let scale = bond-length * 30pt

  let label-margin    = 0.27
  let double-gap      = 0.065
  let ring-double-gap = 0.09
  let inner-trim      = 0.07
  let stroke-w        = 0.9pt

  let atom-clr = if color { _atom-color } else { (sym) => black }

  let cos-a = calc.cos(rotation)
  let sin-a = calc.sin(rotation)
  let rx(x, y) = x * cos-a - y * sin-a
  let ry(x, y) = x * sin-a + y * cos-a

  cetz.canvas(length: scale, {
    import cetz.draw: *

    for bond in layout.bonds {
      let af = layout.atoms.at(bond.from)
      let at = layout.atoms.at(bond.to)
      let p1x = rx(af.pos.x, af.pos.y)
      let p1y = ry(af.pos.x, af.pos.y)
      let p2x = rx(at.pos.x, at.pos.y)
      let p2y = ry(at.pos.x, at.pos.y)
      let c1 = atom-clr(af.symbol)
      let c2 = atom-clr(at.symbol)

      let dx = p2x - p1x
      let dy = p2y - p1y
      let len = calc.sqrt(dx * dx + dy * dy)
      let ux = if len > 0.001 { dx / len } else { 1.0 }
      let uy = if len > 0.001 { dy / len } else { 0.0 }

      let s1 = if _has-label(af) { label-margin } else { 0.0 }
      let s2 = if _has-label(at) { label-margin } else { 0.0 }
      let q1x = p1x + ux * s1
      let q1y = p1y + uy * s1
      let q2x = p2x - ux * s2
      let q2y = p2y - uy * s2
      let mx = (q1x + q2x) / 2
      let my = (q1y + q2y) / 2

      let stereo = bond.at("stereo", default: "none")

      if stereo == "wedge_up" {
        // Solid filled triangle: point at q1 (from-atom), wide base at q2 (to-atom).
        let half-w = 0.10
        let ox = -uy * half-w
        let oy =  ux * half-w
        line(
          (q1x, q1y), (q2x + ox, q2y + oy), (q2x - ox, q2y - oy),
          close: true, fill: c2, stroke: none,
        )

      } else if stereo == "wedge_down" {
        // Hashed wedge: parallel lines perpendicular to bond, widening toward q2.
        let n-lines = 7
        let half-w = 0.10
        for i in range(n-lines) {
          let t = (i + 1) / (n-lines + 1)
          let cx = q1x + (q2x - q1x) * t
          let cy = q1y + (q2y - q1y) * t
          let w  = half-w * t
          let ox = -uy * w
          let oy =  ux * w
          let c  = if t < 0.5 { c1 } else { c2 }
          line((cx - ox, cy - oy), (cx + ox, cy + oy), stroke: stroke-w + c)
        }

      } else if bond.order == 1 {
        line((q1x, q1y), (mx, my),   stroke: stroke-w + c1)
        line((mx, my),   (q2x, q2y), stroke: stroke-w + c2)

      } else if bond.order == 2 {
        let is-ring-bond = bond.inner_x != 0.0 or bond.inner_y != 0.0

        if is-ring-bond {
          line((q1x, q1y), (mx, my),   stroke: stroke-w + c1)
          line((mx, my),   (q2x, q2y), stroke: stroke-w + c2)
          let ix  = rx(bond.inner_x, bond.inner_y) * ring-double-gap
          let iy  = ry(bond.inner_x, bond.inner_y) * ring-double-gap
          let lx1 = q1x + ux * inner-trim + ix
          let ly1 = q1y + uy * inner-trim + iy
          let lx2 = q2x - ux * inner-trim + ix
          let ly2 = q2y - uy * inner-trim + iy
          let lmx = (lx1 + lx2) / 2
          let lmy = (ly1 + ly2) / 2
          line((lx1, ly1), (lmx, lmy), stroke: stroke-w + c1)
          line((lmx, lmy), (lx2, ly2), stroke: stroke-w + c2)
        } else {
          // At junction ends (degree > 1) extend slightly to close corner gap;
          // at terminal ends (degree 1) keep flat.
          let deg-f = layout.bonds.filter(b => b.from == bond.from or b.to == bond.from).len()
          let deg-t = layout.bonds.filter(b => b.from == bond.to   or b.to == bond.to  ).len()
          let ext = 0.04
          let ox = -uy * double-gap
          let oy =  ux * double-gap
          let e1x = q1x - ux * (if deg-f > 1 { ext } else { 0.0 })
          let e1y = q1y - uy * (if deg-f > 1 { ext } else { 0.0 })
          let e2x = q2x + ux * (if deg-t > 1 { ext } else { 0.0 })
          let e2y = q2y + uy * (if deg-t > 1 { ext } else { 0.0 })
          line((e1x + ox, e1y + oy), (mx + ox, my + oy), stroke: stroke-w + c1)
          line((mx + ox,  my + oy),  (e2x + ox, e2y + oy), stroke: stroke-w + c2)
          line((e1x - ox, e1y - oy), (mx - ox, my - oy), stroke: stroke-w + c1)
          line((mx - ox,  my - oy),  (e2x - ox, e2y - oy), stroke: stroke-w + c2)
        }

      } else if bond.order == 3 {
        let ox = -uy * double-gap * 1.3
        let oy =  ux * double-gap * 1.3
        line((q1x,      q1y     ), (mx,      my     ), stroke: stroke-w + c1)
        line((mx,       my      ), (q2x,     q2y    ), stroke: stroke-w + c2)
        line((q1x + ox, q1y + oy), (mx + ox, my + oy), stroke: stroke-w + c1)
        line((mx + ox,  my + oy ), (q2x + ox, q2y + oy), stroke: stroke-w + c2)
        line((q1x - ox, q1y - oy), (mx - ox, my - oy), stroke: stroke-w + c1)
        line((mx - ox,  my - oy ), (q2x - ox, q2y - oy), stroke: stroke-w + c2)

      } else if bond.order == 4 {
        line((q1x, q1y), (mx, my),   stroke: stroke-w + c1)
        line((mx, my),   (q2x, q2y), stroke: stroke-w + c2)
        let ox = -uy * double-gap
        let oy =  ux * double-gap
        line(
          (q1x + ox, q1y + oy), (mx + ox, my + oy),
          stroke: (paint: c1, thickness: stroke-w, dash: "dashed"),
        )
        line(
          (mx + ox, my + oy), (q2x + ox, q2y + oy),
          stroke: (paint: c2, thickness: stroke-w, dash: "dashed"),
        )
      }
    }

    // Atom labels — heteroatoms, charged atoms, and abbreviated groups.
    // Positions are rotated; text content stays upright.
    for atom in layout.atoms {
      if _has-label(atom) {
        let abbrev = atom.at("abbrev", default: "")
        let fill   = atom-clr(atom.symbol)

        let label-content = if abbrev != "" {
          // Abbreviated group: render the label text in black, italic, serif.
          text(
            size: atom-font-size,
            font: "New Computer Modern",
            style: "italic",
            weight: "regular",
            fill: black,
            abbrev,
          )
        } else {
          let charge-str = if atom.charge == 1        { "+" }
                           else if atom.charge == -1  { "\u{2212}" }
                           else if atom.charge > 1    { str(atom.charge) + "+" }
                           else if atom.charge < -1   { str(-atom.charge) + "\u{2212}" }
                           else                       { "" }

          let sym-text = text(
            size: atom-font-size,
            font: "New Computer Modern",
            weight: "regular",
            fill: fill,
            atom.symbol,
          )

          if charge-str == "" {
            sym-text
          } else {
            sym-text + super(text(
              fill: fill,
              font: "New Computer Modern",
              weight: "regular",
              charge-str,
            ))
          }
        }

        content(
          (rx(atom.pos.x, atom.pos.y), ry(atom.pos.x, atom.pos.y)),
          label-content,
          anchor: "center",
          padding: 1pt,
        )
      }
    }
  })
}

// ── Reaction scheme helpers ───────────────────────────────────────────────────

/// Creates a reaction arrow for use inside #reaction().
///
/// - above (content): Label above a horizontal arrow / to the right of a vertical one.
/// - below (content): Label below a horizontal arrow / to the left of a vertical one.
/// - dir (str): Arrow direction — "right" (default), "left", "down", or "up".
/// -> dictionary  (consumed by #reaction)
#let rxn-arrow(above: none, below: none, dir: "right") = (
  __rxn_arrow__: true,
  above: above,
  below: below,
  dir: dir,
)

// Render a horizontal arrow (→ or ←)
#let _horiz-arrow(above, below, dir) = {
  let arrow-parts = ()
  if above != none { arrow-parts.push(align(center, text(size: 8pt, above))) }
  let (sx, ex) = if dir == "left" { (52, 0) } else { (0, 52) }
  arrow-parts.push(cetz.canvas(length: 1pt, {
    import cetz.draw: *
    line((sx, 0), (ex, 0), mark: (end: ">", fill: black, size: 5), stroke: 0.8pt + black)
  }))
  if below != none { arrow-parts.push(align(center, text(size: 8pt, below))) }
  align(center + horizon, stack(spacing: 3pt, ..arrow-parts))
}

// Render a vertical arrow (↓ or ↑) with optional labels to the right
#let _vert-arrow(above, below, dir) = {
  let (from-y, to-y) = if dir == "up" { (0, 52) } else { (52, 0) }
  let arrow-canvas = cetz.canvas(length: 1pt, {
    import cetz.draw: *
    line((0, from-y), (0, to-y), mark: (end: ">", fill: black, size: 5), stroke: 0.8pt + black)
  })
  let label-parts = ()
  if above != none { label-parts.push(text(size: 8pt, above)) }
  if below != none { label-parts.push(text(size: 8pt, below)) }
  if label-parts.len() == 0 {
    align(center + horizon, arrow-canvas)
  } else {
    grid(
      columns: (auto, auto),
      column-gutter: 4pt,
      align: center + horizon,
      arrow-canvas,
      stack(spacing: 2pt, ..label-parts),
    )
  }
}

/// Wraps any content with an optional centred label below it.
///
/// - content (content): The molecule or formula to display.
/// - label (content): Optional label shown below. Default: none.
/// -> content
#let mol(content, label: none) = {
  if label == none { content }
  else { stack(spacing: 4pt, content, align(center, label)) }
}

/// Lays out a reaction scheme supporting multi-directional arrows.
///
/// Items are any mix of content (smiles(), ce(), text…) and rxn-arrow()
/// dictionaries.  rxn-arrow() accepts dir: "right"|"left"|"down"|"up" so
/// schemes can wrap across the page.
///
/// - gap-h (length): Gap between items in horizontal direction. Default: 1.5em.
/// - gap-v (length): Gap between items in vertical direction. Default: 1.5em.
/// -> content
#let reaction(gap-h: 1.5em, gap-v: 1.5em, ..items) = {
  let steps = items.pos()

  // ── Phase 1: assign each item a (grid-row, grid-col) position ────────────
  // Logical positions: molecules at (lr, lc), arrows in between.
  // Grid mapping: mol at (lr,lc) → grid (2*lr, 2*lc)
  //   horizontal arrow after (lr,lc) going right → grid (2*lr, 2*lc+1)
  //   horizontal arrow after (lr,lc) going left  → grid (2*lr, 2*lc-1)
  //   vertical   arrow after (lr,lc) going down  → grid (2*lr+1, 2*lc)
  //   vertical   arrow after (lr,lc) going up    → grid (2*lr-1, 2*lc)

  let lr = 0
  let lc = 0
  let placed = ()

  for item in steps {
    if type(item) == dictionary and item.at("__rxn_arrow__", default: false) {
      let dir = item.at("dir", default: "right")
      // Each non-arrow item increments lc, so arrows slot into the gap
      // just before the current lc position:
      //   right → gap at 2*lc-1 (between prev mol and next mol)
      //   left  → gap at 2*lc-3 (going backward two slots)
      //   down/up → gap at row ±1, same col as last mol (2*lc-2)
      let (gr, gc) = if dir == "right"     { (2 * lr,     2 * lc - 1) }
                     else if dir == "left"  { (2 * lr,     2 * lc - 3) }
                     else if dir == "down"  { (2 * lr + 1, 2 * lc - 2) }
                     else                   { (2 * lr - 1, 2 * lc - 2) }
      placed.push((gr: gr, gc: gc, kind: "arrow", data: item))
      if dir == "left"      { lc -= 2 }
      else if dir == "down" { lr += 1; lc -= 1 }
      else if dir == "up"   { lr -= 1; lc -= 1 }
      // right: lc unchanged — next mol lands exactly at 2*lc
    } else {
      placed.push((gr: 2 * lr, gc: 2 * lc, kind: "mol", data: item))
      lc += 1
    }
  }

  // ── Phase 2: normalise to 0-based ────────────────────────────────────────
  if placed.len() == 0 { return [] }

  let min-gr = placed.fold(placed.first().gr, (m, c) => calc.min(m, c.gr))
  let min-gc = placed.fold(placed.first().gc, (m, c) => calc.min(m, c.gc))
  let max-gr = placed.fold(placed.first().gr, (m, c) => calc.max(m, c.gr))
  let max-gc = placed.fold(placed.first().gc, (m, c) => calc.max(m, c.gc))

  let n-rows = max-gr - min-gr + 1
  let n-cols = max-gc - min-gc + 1

  // ── Phase 3: build sparse lookup ─────────────────────────────────────────
  let lookup = (:)
  for p in placed {
    let key = str(p.gr - min-gr) + "," + str(p.gc - min-gc)
    lookup.insert(key, p)
  }

  // ── Phase 4: render as Typst grid ────────────────────────────────────────
  let flat-cells = ()
  for gr in range(n-rows) {
    for gc in range(n-cols) {
      let key = str(gr) + "," + str(gc)
      let p = lookup.at(key, default: none)
      if p == none {
        flat-cells.push([])
      } else if p.kind == "mol" {
        flat-cells.push(align(center + horizon, p.data))
      } else {
        let item = p.data
        let d = item.at("dir", default: "right")
        flat-cells.push(
          if d == "right" or d == "left" { _horiz-arrow(item.above, item.below, d) }
          else                           { _vert-arrow(item.above, item.below, d) }
        )
      }
    }
  }

  grid(
    columns: (auto,) * n-cols,
    rows:    (auto,) * n-rows,
    column-gutter: gap-h / 2,
    row-gutter:    gap-v / 2,
    align: center + horizon,
    ..flat-cells,
  )
}

/// Alias for #smiles.
#let display-smiles = smiles
