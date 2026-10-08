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
# Typst escapes quotation marks in panic messages, so quoted values appear as \".
expect_error "abbreviate-wrong-type" 'smiles abbreviate is invalid: expected none, \"all\", a group name, or an array of group names, got true'
expect_error "abbreviate-unknown-group" 'smiles abbreviate is invalid: unknown group \"Boc\". Available groups are \"tBu\", \"CF3\", \"NO2\", \"CN\", \"OEt\", \"OMe\", \"Ac\".'
expect_error "abbreviate-wrong-case" "Group names are case-sensitive"
expect_error "abbreviate-duplicate-group" 'group \"OMe\" is listed more than once'
expect_error "abbreviate-all-in-list" '\"all\" appears inside a list of group names'
expect_error "abbreviate-absent-group" 'smiles abbreviate group \"CF3\" is invalid: no terminal CF3 group in \"CCO\" can be drawn as a label'
expect_error "abbreviate-expanded-group" "Groups with isotopes, stereo marks, unexpected charges, or bracket hydrogens stay expanded"
expect_error "abbreviate-lower-priority-group" 'smiles abbreviate group \"Ac\" is invalid'
expect_error "abbreviate-hidden-arrow-atom" "atom 0 is hidden inside the automatic abbreviation OMe (atoms 0, 1)"
expect_error "abbreviate-hidden-highlight-bond" 'Reference the labeled attachment atom 1 instead, or remove \"OMe\" from abbreviate'
expect_error "abbreviate-hidden-smarts-match" "is hidden inside the automatic abbreviation NO2"
expect_error "abbreviate-hidden-annotation" "atom-annotations atom index is invalid: atom 0 is hidden inside the automatic abbreviation CF3"
expect_error "abbreviate-hidden-bond-customization" "bond-customizations first atom index is invalid: atom 0 is hidden"
expect_error "abbreviate-hidden-show-h" "show-h atom index is invalid: atom 0 is hidden inside the automatic abbreviation OEt"
expect_error "abbreviate-label-show-h" "1 is drawn as the automatic abbreviation tBu, which shows no hydrogens"
expect_error "abbreviate-label-lone-pair" "is drawn as the automatic abbreviation OMe, which has no addressable lone pairs"
expect_error "mol-abbreviate-unknown-group" 'mol abbreviate is invalid: unknown group \"Me\"'
expect_error "mol-abbreviate-hidden-atom" "in species 0 is invalid: atom 0 is hidden inside the automatic abbreviation OMe"
expect_error "cetz-abbreviate-unknown-group" 'smiles-cetz abbreviate is invalid: unknown group \"Et\"'

echo "All editor-visible validation cases passed."
