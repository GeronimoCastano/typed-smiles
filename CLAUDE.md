# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Commit and Attribution Rules

Never include any of the following in commit messages, PR descriptions, code comments, or project files:
- `Co-Authored-By:` lines
- Generated-by or tool-attribution footers
- Any assistant/tooling attribution

## Current Status

`typed-smiles` is a Typst package that renders SMILES chemical structure strings as 2D molecular diagrams.

Architecture:

```text
SMILES string -> Rust WASM plugin -> JSON layout data -> CeTZ drawing in Typst
```

The Rust plugin parses SMILES and computes 2D coordinates. The Typst layer renders bonds, atom labels, colors, hydrogens, formulas, and reaction helpers.

The current package name is `typed-smiles`. The package submission PR for `0.1.0` lives in the `typst/packages` fork, under:

```text
packages/preview/typed-smiles/0.1.0/
```

Use `main` as the source of truth for development. The `package-submission` branch is only a staging branch for the runtime package bundle.

## Commands

```sh
# Run Rust tests
cargo test --manifest-path plugin/Cargo.toml

# Type-check Rust without running tests
cargo check --manifest-path plugin/Cargo.toml

# Build the bundled WASM plugin
./build.sh

# Compile the visual smoke test
typst compile --root . tests/test.typ tests/test.pdf

# Regenerate README assets
typst compile --root . --ppi 300 assets/readme/basics.typ assets/readme/basics.png
typst compile --root . --ppi 300 assets/readme/scaling.typ assets/readme/scaling.png
typst compile --root . --ppi 300 assets/readme/formulas.typ assets/readme/formulas.png
typst compile --root . --ppi 300 assets/readme/reactions.typ assets/readme/reactions.png
typst compile --root . --ppi 300 assets/readme/schemes.typ assets/readme/schemes.png
typst compile --root . --ppi 300 assets/readme/stereo-h.typ assets/readme/stereo-h.png

# Prepare a Typst packages preview copy for a version
scripts/package-preview.sh 0.1.1 /path/to/typst/packages
```

## Project Structure

```text
typst.toml          Package manifest
src/lib.typ         Typst API: #smiles(), #display-smiles, ce, reaction helpers
plugin/
  Cargo.toml        Rust dependencies
  src/
    lib.rs          WASM entrypoint and native test entrypoint
    graph.rs        Chain AST walker -> MoleculeGraph
    layout.rs       2D coordinate generation and implicit H calculation
    render.rs       JSON output types
    error.rs        SmilesError type
assets/readme/      Typst source and rendered PNG README examples
tests/test.typ      Visual smoke test
scripts/            Release helper scripts
```

## Typst API

Main molecule renderer:

```typst
#smiles(
  smiles-str,
  bond-length: 1.0,
  atom-font-size: 11pt,
  color: true,
  rotation: 0deg,
  show-h: false,
)
```

Also exported:

```typst
#display-smiles
#ce
#rxn-arrow(above: none, below: none, dir: "right")
#mol(content, label: none)
#reaction(gap-h: 1.5em, gap-v: 1.5em, ..items)
```

`ce` is re-exported from `chemformula`.

## Implemented Features

- SMILES molecule rendering through `#smiles(...)`
- Alias `#display-smiles`
- Jmol-like CPK atom colors
- Heteroatom labels
- Charge labels
- Abbreviation labels using `{label}` syntax, for example `{PPh3}C=O`
- Upright abbreviation labels
- Reaction helpers: `ce`, `rxn-arrow`, `mol`, `reaction`
- Horizontal and vertical reaction arrows
- Wrap-around reaction schemes
- Wedge and hashed wedge rendering from `/` and `\`
- Wedge direction heuristic with narrow tip at the likely stereocenter
- Improved acyclic double-bond rendering
- Triple-bond rendering with shortened outer lines
- `show-h: true` implicit hydrogen labels
- Explicit bracket hydrogens such as `[NH4+]`

## Key Design Decisions

**Manual Chain AST walker:** `smiles_parser::graph::MoleculeGraph::from_chain` ignores ring-closure bonds, panics on some heteroatoms/halogens, and misreads bond order. The local walker in `graph.rs` avoids those issues.

**Bond order invariant:** In the smiles-parser Chain AST, `chain.bond_or_dot` is the bond from the current atom outward to `chain.chain`. It must be passed forward as `incoming_bond` when recursing.

**Aromatic atoms:** `smiles-parser` 0.4 does not parse unbracketed aromatic atoms such as `c`, `n`, and `o`. Use Kekule SMILES such as `C1=CC=CC=C1` until preprocessing or Kekulization is added.

**2D layout:** Ring placement uses regular polygons, fused rings are placed from shared edges, and acyclic chains are placed with alternating zigzag angles. The layout is centered before rendering.

**Hydrogens:** Rust computes `implicit_h` from a small standard valence map and skips implicit hydrogen generation when explicit bracket hydrogens are present.

**JSON contract:** Rust returns JSON matching `LayoutOutput`. Typst reads it with `json(raw-bytes)`. Coordinates are in bond-length units; Typst scales them with `bond-length * 30pt`.

## Dependencies

Rust plugin:
- `smiles-parser = "0.4"`
- `ptable = { package = "periodic-table-on-an-enum" }`
- `wasm-minimal-protocol = "0.2"`
- `serde` and `serde_json`

Typst imports in `src/lib.typ`:
- `@preview/cetz:0.5.2`
- `@preview/chemformula:0.1.3`

Do not add a top-level `[dependencies]` table to `typst.toml`; the Typst package checker rejects unknown manifest sections.

## Packaging Notes

For Typst package PRs, copy only runtime/package files:

```text
typst.toml
README.md
LICENSE
src/lib.typ
plugin/typst_smiles_plugin.wasm
assets/readme/*.png
```

Do not include development-only files in the Typst package PR:

```text
plugin/src/*
plugin/Cargo.toml
plugin/Cargo.lock
build.sh
tests/*
assets/readme/*.typ
CLAUDE.md
AGENTS.md
.DS_Store
plugin/target/*
```

README PNGs should be committed to `typst/packages` for Typst Universe rendering, but excluded from the downloadable package archive with:

```toml
exclude = ["assets/readme/*.png"]
```

## Known Limitations / Next Steps

- Aromatic SMILES (`c1ccccc1`) need a Kekulization preprocessor.
- Directional `/` and `\` bonds are rendered visually as wedge/hash bonds, but full stereochemical interpretation is not implemented.
- The package does not compute or validate R/S or E/Z descriptors.
- Ring stereochemistry and stereobonds between stereocenters are not handled semantically.
- Bridged bicyclics can still have overlap and would benefit from templates.
