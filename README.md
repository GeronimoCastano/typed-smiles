# smiles

A Typst package for rendering SMILES chemical structure strings as 2D molecular diagrams.

## Usage

```typst
#import "@preview/smiles:0.1.0": smiles, display-smiles

// Ethanol
#smiles("CCO")

// Benzene (Kekulé form)
#smiles("C1=CC=CC=C1")

// Caffeine
#smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")

// Custom size
#smiles("C1=CC=CC=C1", width: 150pt, height: 150pt)
```

## API

### `#smiles(smiles-str, width, height, bond-length, atom-font-size)`

Renders a SMILES string as a 2D molecular diagram.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `smiles-str` | `str` | required | OpenSMILES string |
| `width` | `length` | `200pt` | Canvas width |
| `height` | `length` | `200pt` | Canvas height |
| `bond-length` | `float` | `1.0` | Bond length scale factor |
| `atom-font-size` | `length` | `8pt` | Heteroatom label font size |

`#display-smiles` is an alias for `#smiles`.

## SMILES Support

The package uses the [smiles-parser](https://crates.io/crates/smiles-parser) crate (OpenSMILES spec) for parsing. Current limitations:

- **Aromatic atoms**: Use Kekulé notation (`C1=CC=CC=C1`) rather than lowercase aromatic notation (`c1ccccc1`). Aromatic notation support is planned.
- **Stereochemistry**: Wedge bonds are not yet rendered (stereo information is parsed but ignored in the layout).
- **Complex ring systems**: Fused bicyclics work; bridged bicyclics may have atom overlap.

## Architecture

The package is implemented as a Typst WASM plugin:

```
SMILES string → [Rust WASM plugin] → JSON layout → CeTZ drawing
```

The Rust plugin (`plugin/`) handles parsing and 2D coordinate generation. The Typst layer (`src/lib.typ`) calls the plugin and renders with [CeTZ](https://github.com/cetz-package/cetz).

## Building

```sh
# Run tests
cargo test --manifest-path plugin/Cargo.toml

# Build WASM (requires wasm-pack or wasm32-unknown-unknown target)
rustup target add wasm32-unknown-unknown
cargo build --manifest-path plugin/Cargo.toml --target wasm32-unknown-unknown --release
```

## License

MIT
