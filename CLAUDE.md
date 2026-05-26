# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commit and Attribution Rules

Never include any of the following in commit messages, PR descriptions, code comments, or any project files:
- `Co-Authored-By:` lines
- "Generated with Claude Code" or similar footers
- Any mention of AI tools or assistants

## What This Is

A Typst package that renders SMILES chemical structure strings as 2D molecular diagrams. Architecture:

```
SMILES string → [Rust WASM plugin] → JSON layout data → CeTZ drawing in Typst
```

The Rust plugin handles the hard parts (parsing + 2D layout). The Typst layer is a thin wrapper.

## Commands

```sh
# Run all tests (native Rust, no WASM needed)
cargo test --manifest-path plugin/Cargo.toml

# Type-check without running tests
cargo check --manifest-path plugin/Cargo.toml

# Build WASM plugin (required to use the Typst package)
rustup target add wasm32-unknown-unknown
cargo build --manifest-path plugin/Cargo.toml --target wasm32-unknown-unknown --release
# Output: plugin/target/wasm32-unknown-unknown/release/typst_smiles_plugin.wasm

# Test a single Rust test by name
cargo test --manifest-path plugin/Cargo.toml <test_name>
```

## Project Structure

```
typst.toml          — Typst package manifest
src/lib.typ         — Typst API: #smiles(), #display-smiles()
plugin/
  Cargo.toml        — Rust deps: smiles-parser, wasm-minimal-protocol, serde_json
  src/
    lib.rs          — WASM entrypoint + native test entrypoint
    graph.rs        — Chain AST walker → MoleculeGraph (atoms + bonds + adj list)
    layout.rs       — 2D coordinate generation (ring placement + chain DFS)
    render.rs       — Output types: LayoutOutput, AtomOutput, BondOutput, Vec2
    error.rs        — SmilesError type
```

## Key Design Decisions

**Why we walk the Chain AST manually (graph.rs)**: `smiles_parser::graph::MoleculeGraph::from_chain` has bugs — it ignores ring-closure bonds, panics on N/S/halogens, and misreads bond order (reads it from the wrong chain node). Our walker fixes all three issues.

**Bond order invariant**: In the smiles-parser Chain AST, `chain.bond_or_dot` is the bond from the **current** atom **outward** to `chain.chain`. It must be passed **forward** as `incoming_bond` when recursing into `chain.chain`, not read by the child. This is the non-obvious thing that was wrong in the upstream graph module.

**Aromatic atoms**: `smiles-parser` 0.4 does not parse unbracketed aromatic atoms (`c`, `n`, `o`). Use Kekulé SMILES (`C1=CC=CC=C1`) for aromatic molecules until we add a preprocessing/Kekulization step.

**2D layout**: No pure-Rust 2D layout library exists. Our algorithm (layout.rs):
1. DFS ring detection → rings as atom-index lists
2. Place each ring as a regular n-gon (circumradius = 1/(2·sin(π/n)) for unit bond length)
3. Fused rings: find the shared edge, compute center on the far side, place remaining atoms
4. Acyclic chains: DFS from the first placed atom, alternating ±120° zigzag
5. Center by subtracting centroid

**WASM entrypoint** (lib.rs): Uses `wasm-minimal-protocol` v0.2 macros. The `initiate_protocol!()` + `#[wasm_func]` macros handle the memory protocol. On non-WASM targets, `layout_native()` is exposed for tests.

**JSON contract** (render.rs → lib.typ): The plugin returns JSON matching `LayoutOutput`. Typst reads it via `json(str(raw-bytes))`. Coordinates are in "bond-length units" (1.0 = one bond length); Typst scales with `bond-length * 30pt` as the base.

## Dependencies

- `smiles-parser = "0.4"` — nom-based OpenSMILES parser
- `ptable = { package = "periodic-table-on-an-enum" }` — needed to call `.get_symbol()` on elements (smiles-parser uses the same alias)
- `wasm-minimal-protocol = "0.2"` — Typst plugin protocol macros
- `serde_json = "1"` — serialize LayoutOutput to JSON bytes

## Known Limitations / Next Steps

- Aromatic SMILES (`c1ccccc1`) need a Kekulization preprocessor
- Stereochemistry (wedge/dash bonds) is parsed but not rendered
- Bridged bicyclics (e.g., norbornane) may have atom overlap — template matching needed
- No explicit hydrogen rendering (shown only when in bracket atoms `[NH4+]`)
- WASM binary not yet built/bundled (the `.wasm` file needs to be placed at `plugin/pkg/`)
