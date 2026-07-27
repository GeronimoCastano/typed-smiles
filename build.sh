#!/usr/bin/env bash
set -e
project_root="$(cd "$(dirname "$0")" && pwd)"
cargo build \
  --manifest-path "$project_root/plugin/Cargo.toml" \
  --target wasm32-unknown-unknown \
  --release
built_plugin="$project_root/plugin/target/wasm32-unknown-unknown/release/typst_smiles_plugin.wasm"
bundled_plugin="$project_root/plugin/typst_smiles_plugin.wasm"
cp "$built_plugin" "$bundled_plugin"
plugin_size="$(du -h "$bundled_plugin" | cut -f1)"
echo "WASM built → plugin/typst_smiles_plugin.wasm ($plugin_size)"
