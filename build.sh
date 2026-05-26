#!/usr/bin/env bash
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cargo build \
  --manifest-path "$SCRIPT_DIR/plugin/Cargo.toml" \
  --target wasm32-unknown-unknown \
  --release
cp "$SCRIPT_DIR/plugin/target/wasm32-unknown-unknown/release/typst_smiles_plugin.wasm" \
   "$SCRIPT_DIR/plugin/typst_smiles_plugin.wasm"
echo "WASM built → plugin/typst_smiles_plugin.wasm ($(du -h "$SCRIPT_DIR/plugin/typst_smiles_plugin.wasm" | cut -f1))"
