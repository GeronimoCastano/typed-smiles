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
expect_error "invalid-smarts" "invalid SMARTS"
expect_error "unsupported-smarts" "unsupported atom predicate"
expect_error "empty-smarts" "highlight-smarts is invalid"
expect_error "smarts-wrong-type" "highlight-smarts is invalid"
expect_error "smarts-invalid-item" "highlight-smarts is invalid"
expect_error "unknown-functional-group" "unknown functional group"
expect_error "unmatched-smarts" "no match for"
expect_error "unmatched-group" "no match for"
expect_error "invalid-highlight-colors" "highlight-colors is invalid"
expect_error "invalid-highlight-color-item" "highlight-colors is invalid"
expect_error "invalid-highlight-policy" "highlight-unmatched is invalid"
expect_error "query-invalid-smiles" "invalid SMILES"
expect_error "highlight-request-missing-pattern" "expected a non-empty pattern string"
expect_error "highlight-request-missing-group" "expected a non-empty group string"
expect_error "highlight-request-invalid-bool" "highlight-smarts include-atoms is invalid"
expect_error "highlight-group-invalid-bool" "highlight-groups include-atoms is invalid"
expect_error "highlight-request-unknown-option" "unknown request option"
expect_error "unknown-library-molecule" 'dictionary does not contain key "cafeine"'
expect_error "one-sided-directional-bond" "must mark both ends of a double bond"

expect_error "align-not-array" "align-molecules molecules is invalid"
expect_error "align-single-molecule" "expected an array of at least two SMILES strings"
expect_error "align-molecule-type" "align-molecules molecule 1 is invalid"
expect_error "align-aligned-input" "aligned molecules cannot be aligned again"
expect_error "align-scaffold-type" "align-molecules scaffold is invalid"
expect_error "align-no-correspondence" "align-molecules correspondence is invalid"
expect_error "align-partial-atoms-without-scaffold" "align-molecules correspondence is invalid"
expect_error "align-atoms-length" "one entry per molecule (2)"
expect_error "align-atoms-entry" "align-molecules atoms entry 1 is invalid"
expect_error "align-reference-range" "align-molecules reference is invalid"
expect_error "align-rotation-type" "align-molecules rotation is invalid"
expect_error "align-mirror-value" "align-molecules mirror is invalid"
expect_error "align-allow-reflection-type" "align-molecules allow-reflection is invalid"
expect_error "align-absent-scaffold" "does not occur in molecule 0"
expect_error "align-invalid-scaffold" "invalid SMARTS"
expect_error "align-invalid-smiles" "molecule 1: invalid SMILES"
expect_error "align-distinct-occurrences" "occurs at 2 chemically distinct places in molecule 1"
expect_error "align-atoms-not-scaffold-match" "is not an occurrence of scaffold"
expect_error "align-atom-out-of-range" "atom 7 does not exist in molecule 1"
expect_error "align-repeated-atom" "repeats an atom"
expect_error "align-unequal-atom-lists" "same number of corresponding atoms"
expect_error "align-single-atom-correspondence" "at least two corresponding atoms"
expect_error "aligned-rotation-conflict" "smiles rotation is invalid: an aligned molecule already carries its rotation"
expect_error "aligned-mirror-conflict" "mol mirror is invalid: an aligned molecule already carries its reflection"
expect_error "aligned-cetz-rotation-conflict" "smiles-cetz rotation is invalid"
expect_error "aligned-skeleton-hydrogens" "would undo the alignment"
expect_error "invalid-mol-spec" "expected a SMILES string, aligned molecule, or content"
expect_error "grid-empty" "molecule-grid items is invalid"
expect_error "grid-columns" "molecule-grid columns is invalid"
expect_error "grid-columns-type" "molecule-grid columns is invalid"
expect_error "grid-scale" "molecule-grid scale is invalid"
expect_error "grid-bond-length" "molecule-grid bond-length is invalid"
expect_error "grid-sizing" "molecule-grid sizing is invalid"
expect_error "grid-gutter" "molecule-grid column-gutter is invalid"
expect_error "grid-breakable" "molecule-grid breakable is invalid"
expect_error "grid-unknown-option" "the option is not supported. Use columns, scale"
expect_error "grid-invalid-item" "molecule-grid item 1 is invalid"
expect_error "grid-item-scale" "molecule-grid item 1 scale is invalid"
expect_error "grid-item-bond-length" "molecule-grid item 1 bond-length is invalid"
expect_error "grid-overflow" "molecule-grid item 0 is invalid: the molecule is"
expect_error "grid-scaffold-type" "molecule-grid scaffold is invalid"
expect_error "grid-scaffold-single" "alignment needs at least two molecules"
expect_error "grid-scaffold-content" "scaffold alignment needs a SMILES molecule"
expect_error "grid-scaffold-rotation" "the grid scaffold sets each molecule's orientation"
expect_error "grid-scaffold-absent" "does not occur in molecule 0"

echo "All editor-visible validation cases passed."
