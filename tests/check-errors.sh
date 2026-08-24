#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "$0")/.." && pwd)"
temporary_directory="$(mktemp -d "${TMPDIR:-/tmp}/typed-smiles-errors.XXXXXX")"
trap 'rm -rf "$temporary_directory"' EXIT

expect_error() {
  local case_name="$1"
  local expected_message="$2"
  local output_file="$temporary_directory/$case_name.txt"
  local pdf_file="$temporary_directory/$case_name.pdf"

  if typst compile \
    --root "$project_root" \
    --input "case=$case_name" \
    "$project_root/tests/errors.typ" \
    "$pdf_file" \
    >"$output_file" 2>&1
  then
    echo "Expected validation case '$case_name' to fail, but it compiled."
    return 1
  fi

  if ! grep -Fq "$expected_message" "$output_file"; then
    echo "Validation case '$case_name' did not contain: $expected_message"
    sed -n '1,120p' "$output_file"
    return 1
  fi
}

expect_error "invalid-smiles" "invalid SMILES"
expect_error "unclosed-label" "unclosed custom label"
expect_error "missing-arrow-endpoints" "arrow from is invalid"
expect_error "species-out-of-range" "rxn-arrow() itself does not count as a species"
expect_error "atom-out-of-range" "atom index is invalid"
expect_error "missing-bond" "are not joined by a visible bond"
expect_error "missing-lone-pair" "has no addressable lone pairs"
expect_error "pair-out-of-range" "pair index is invalid"
expect_error "opaque-atom-reference" "opaque content and has no addressable atoms"
expect_error "ignored-annotation" "expected arrow() or highlight()"
expect_error "show-h-out-of-range" "show-h atom index is invalid"
expect_error "mol-formula-wrong-type" "mol-formula SMILES expression is invalid"
expect_error "mol-formula-empty" "mol-formula SMILES expression is invalid"
expect_error "mol-formula-wildcard" "wildcard atom"
expect_error "annotation-out-of-range" "atom-annotations atom index is invalid"
expect_error "customized-missing-bond" "bond-customizations bond reference is invalid"
expect_error "duplicate-bond-customization" "is customized more than once"
expect_error "opacity-out-of-range" "opacity is invalid"
expect_error "unknown-mol-option" "mol option"
expect_error "content-molecule-options" "mol content options is invalid"
expect_error "invalid-reaction-item" "reaction item is invalid"
expect_error "empty-reaction" "reaction items is invalid"
expect_error "empty-cycle" "cycle items is invalid"
expect_error "leading-cycle-step" "a step appears before the first species"
expect_error "duplicate-cycle-step" "more than one step follows the same species"
expect_error "invalid-step-reagent" "step into is invalid"

echo "All editor-visible validation cases passed."
