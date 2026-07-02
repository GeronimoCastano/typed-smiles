# Packaging And Release Checklist

Use `main` as the source of truth. The `package-submission` branch and the
`typst/packages` fork are release artifacts, not development branches.

Agents must read this file before preparing a Typst package release, copying
files into the Typst packages checkout, or opening a Typst packages PR.

## Repositories

Source package repository:

```text
/Users/gerocastano8/Documents/Coding/Projects/typed-smiles
https://github.com/GeronimoCastano/typed-smiles
```

Local Typst packages checkout:

```text
/Users/gerocastano8/Documents/Coding/Projects/typst-packages
```

Typst packages fork and upstream:

```text
origin   https://github.com/GeronimoCastano/packages.git
upstream https://github.com/typst/packages.git
```

The package PR target is always:

```text
typst/packages:main
```

The package PR head is normally:

```text
GeronimoCastano/packages:<release-branch>
```

## Version Choice

- Patch release, such as `0.1.1`: bug fixes, documentation fixes, or small rendering improvements without API changes.
- Minor release, such as `0.2.0`: new public API, behavior changes, or larger rendering features.

Typst packages are immutable in practice. After a version is accepted, publish fixes as a new version instead of editing the old version directory.

Set the version once while following the checklist:

```sh
VERSION=0.4.1
PKG_REPO=/Users/gerocastano8/Documents/Coding/Projects/typst-packages
BRANCH=typed-smiles-package-$VERSION
```

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
typst compile --root . --ppi 300 assets/readme/hydrogens-labels.typ assets/readme/hydrogens-labels.png
typst compile --root . --ppi 300 assets/readme/formulas.typ assets/readme/formulas.png
typst compile --root . --ppi 300 assets/readme/reactions.typ assets/readme/reactions.png
typst compile --root . --ppi 300 assets/readme/schemes.typ assets/readme/schemes.png
typst compile --root . --ppi 300 assets/readme/stereo-h.typ assets/readme/stereo-h.png
typst compile --root . --ppi 300 assets/readme/mirror.typ assets/readme/mirror.png
```

8. Commit and push `main`.

## Typst Packages PR Steps

1. Go to the local Typst packages checkout:

```sh
cd "$PKG_REPO"
```

2. Make sure the checkout is clean:

```sh
git status -sb
```

3. Sync the existing fork with upstream `main`:

```sh
git remote get-url upstream || git remote add upstream https://github.com/typst/packages.git
git fetch upstream
```

4. Create a new package update branch from upstream `main`:

```sh
git switch -C "$BRANCH" upstream/main
```

5. Copy the package files from the source repository:

```sh
cd /Users/gerocastano8/Documents/Coding/Projects/typed-smiles
scripts/package-preview.sh "$VERSION" "$PKG_REPO"
```

This creates:

```text
/Users/gerocastano8/Documents/Coding/Projects/typst-packages/packages/preview/typed-smiles/<VERSION>/
```

6. From the `typst-packages/packages` directory, run:

```sh
cd "$PKG_REPO/packages"
typst-package-check check @preview/typed-smiles:$VERSION
```

7. Commit the new version directory in the fork:

```sh
cd "$PKG_REPO"
git add "packages/preview/typed-smiles/$VERSION"
git commit -m "typed-smiles:$VERSION"
git push -u origin "$BRANCH"
```

8. Open a PR to `typst/packages:main`.

PR title format:

```text
typed-smiles:<VERSION>
```

Example:

```text
typed-smiles:0.4.1
```

For package updates, use only the compact update PR body. Do not include a
validation section, a checks section, generated artifact notes, or the full
new-package checklist.

PR body format:

```markdown
I am submitting

- [ ] a new package
- [x] an update for a package

Changes:

- <meaningful user-visible change>
- <meaningful user-visible change>
```

Example update PR body:

```markdown
I am submitting

- [ ] a new package
- [x] an update for a package

Changes:

- Add equilibrium reaction arrow styles via `rxn-arrow(kind: "equilibrium")` and `rxn-arrow(kind: "equilibrium-filled")`.
- Document the new `rxn-arrow(kind:)` option.
```

Reserve the full Typst package submission checklist for the initial package
submission only.

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
