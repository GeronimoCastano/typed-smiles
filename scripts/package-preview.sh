#!/usr/bin/env sh
set -eu

usage() {
  echo "Usage: scripts/package-preview.sh <version> <typst-packages-repo>" >&2
  echo "Example: scripts/package-preview.sh 0.2.0 /path/to/typst/packages" >&2
}

if [ "$#" -ne 2 ]; then
  usage
  exit 2
fi

version="$1"
packages_repo="${2%/}"
package_name="typed-smiles"

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
package_target="$packages_repo/packages/preview/$package_name/$version"

if [ ! -d "$packages_repo/packages/preview" ]; then
  echo "error: not a typst/packages checkout: $packages_repo" >&2
  exit 1
fi

manifest_version=$(awk -F '"' '/^version = / { print $2; exit }' "$repo_root/typst.toml")
if [ "$manifest_version" != "$version" ]; then
  echo "error: typst.toml version is $manifest_version, not $version" >&2
  exit 1
fi

for required_file in \
  "$repo_root/typst.toml" \
  "$repo_root/README.md" \
  "$repo_root/LICENSE" \
  "$repo_root/src/lib.typ" \
  "$repo_root/plugin/typst_smiles_plugin.wasm"
do
  if [ ! -f "$required_file" ]; then
    echo "error: missing required file: $required_file" >&2
    exit 1
  fi
done

if [ -e "$package_target" ]; then
  echo "error: target already exists: $package_target" >&2
  echo "Remove it manually if you intentionally want to recreate it." >&2
  exit 1
fi

mkdir -p "$package_target/src" "$package_target/plugin" "$package_target/assets/readme"

cp "$repo_root/typst.toml" "$package_target/typst.toml"
cp "$repo_root/README.md" "$package_target/README.md"
cp "$repo_root/LICENSE" "$package_target/LICENSE"
cp "$repo_root/src/lib.typ" "$package_target/src/lib.typ"
cp "$repo_root/plugin/typst_smiles_plugin.wasm" "$package_target/plugin/typst_smiles_plugin.wasm"

found_readme_image=0
for image in "$repo_root"/assets/readme/*.png; do
  if [ -f "$image" ]; then
    cp "$image" "$package_target/assets/readme/"
    found_readme_image=1
  fi
done

if [ "$found_readme_image" -eq 0 ]; then
  echo "warning: no README PNG assets found in assets/readme" >&2
fi

echo "Prepared $package_name $version at:"
echo "$package_target"
echo
echo "Next:"
echo "  cd $packages_repo/packages"
echo "  typst-package-check check @preview/$package_name:$version"
