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
#smiles("C1=CC=CC=C1", bond-length: 1.4)
```

## API

### `#smiles(smiles-str, bond-length, atom-font-size, color, rotation, show-h)`

Renders a SMILES string as a 2D molecular diagram.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `smiles-str` | `str` | required | OpenSMILES string |
| `bond-length` | `float` | `1.0` | Uniform bond length scale factor; `1.0` equals 30pt per bond |
| `atom-font-size` | `length` | `11pt` | Heteroatom label font size |
| `color` | `bool` | `true` | Apply Jmol CPK atom colors |
| `rotation` | `angle` | `0deg` | Rotate the molecule while keeping atom labels upright |
| `show-h` | `bool` | `false` | Show computed implicit hydrogens |

Explicit bracket hydrogens, such as `[NH4+]`, are always shown. Computed
implicit hydrogens are hidden unless `show-h: true` is set:

```typst
#smiles("CCO", show-h: true)
```

### Sizing and aspect ratio

`#smiles` currently does not take `width` or `height` arguments. Use `bond-length` to scale the drawing up or down:

```typst
#smiles("CCO", bond-length: 0.8)
#smiles("CCO", bond-length: 1.5)
```

Scaling is uniform, so the molecule's aspect ratio is preserved. The drawing is not independently stretched along the x and y axes.

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
