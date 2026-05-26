# Release Checklist

Use `main` as the source of truth. The `package-submission` branch and the `typst/packages` fork are release artifacts, not development branches.

## Version Choice

- Patch release, such as `0.1.1`: bug fixes, documentation fixes, or small rendering improvements without API changes.
- Minor release, such as `0.2.0`: new public API, behavior changes, or larger rendering features.

Typst packages are immutable in practice. After a version is accepted, publish fixes as a new version instead of editing the old version directory.

## Local Release Steps

1. Make changes on `main`.
2. Update `typst.toml` to the new version.
3. Update README imports and examples to use the new version.
4. Build the WASM plugin:

```sh
./build.sh
```

5. Run Rust tests:

```sh
cargo test --manifest-path plugin/Cargo.toml
```

6. Compile the visual smoke test:

```sh
typst compile --root . tests/test.typ tests/test.pdf
```

7. Regenerate README assets if examples or rendering changed:

```sh
typst compile --root . --ppi 300 assets/readme/basics.typ assets/readme/basics.png
typst compile --root . --ppi 300 assets/readme/scaling.typ assets/readme/scaling.png
typst compile --root . --ppi 300 assets/readme/formulas.typ assets/readme/formulas.png
typst compile --root . --ppi 300 assets/readme/reactions.typ assets/readme/reactions.png
typst compile --root . --ppi 300 assets/readme/schemes.typ assets/readme/schemes.png
typst compile --root . --ppi 300 assets/readme/stereo-h.typ assets/readme/stereo-h.png
```

8. Commit and push `main`.

## Typst Packages PR Steps

1. Sync the existing fork of `typst/packages` with upstream `main`.
2. Create a new branch in that fork, for example:

```sh
git switch -c add-typed-smiles-0.1.1
```

3. Copy the package files:

```sh
scripts/package-preview.sh 0.1.1 /path/to/typst/packages
```

This creates:

```text
/path/to/typst/packages/packages/preview/typed-smiles/0.1.1/
```

4. From the `typst/packages/packages` directory, run:

```sh
typst-package-check check @preview/typed-smiles:0.1.1
```

5. Commit the new version directory in the fork:

```sh
git add packages/preview/typed-smiles/0.1.1
git commit -m "Add typed-smiles 0.1.1"
git push -u fork add-typed-smiles-0.1.1
```

6. Open a PR to `typst/packages:main` with the title format:

```text
typed-smiles:0.1.1
```

Use the official package submission checklist in the PR body.

## Package Bundle Contents

Include:

```text
typst.toml
README.md
LICENSE
src/lib.typ
plugin/typst_smiles_plugin.wasm
assets/readme/*.png
```

Do not include:

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

README PNGs are committed for Typst Universe rendering, but excluded from the downloaded package bundle with:

```toml
exclude = ["assets/readme/*.png"]
```
