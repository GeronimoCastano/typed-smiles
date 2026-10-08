# Functional-group verification corpus

`substructure-pubchem.json` pins 63 compounds retrieved from PubChem PUG REST
on 2026-10-02 (UTC). Each record retains its CID, source URL, IUPAC name, and
original `ConnectivitySMILES`. No network access is needed to run the tests.

The corpus includes everyday molecules, medicines, amino acids, carbohydrates,
salts, and deliberate negative controls. All 16 named groups have positive
controls. Counts check classification independently of the selected atom/bond
sets: acids versus carboxylates, alcohols versus phenols/enols, amines versus
amides/guanidines/sulfonamides/cyanamides, carboxylic esters versus
anhydrides/carbonates/carbamates, and nitro groups versus nitrate esters.

## Independent oracle

The `reviewed-counts` refer to the chemically interpreted aromatic structures.
RDKit 2026.03.1 generated the exact atom and bond sets from the exported SMARTS
registry. Atom indices refer to each variant's own SMILES string. Query edges
are recorded separately from atom sets, and automorphisms are deduplicated by
both sets. Recursive context does not contribute highlighted atoms or bonds.

An amide match is a local O=C–N motif. Urea has two overlapping matches;
caffeine has three motifs sharing two carbonyls. Alcohol, phenol, amine,
ether, and thiol return their O, N, or S query atom. Named-group drawing also
highlights displayed attached hydrogen for OH, NH, and SH groups, including
skeleton H atoms and bonds. These display fragments have no separate query
index; the pinned RDKit atom/bond sets remain unchanged. Carbon attachments
are context and stay outside the highlight. `substructure-hydrogens.typ`
checks this display expansion separately, including implicit, bracket,
folded explicit, and retained isotope H.

The package preserves explicit aromatic flags instead of inferring aromaticity
from uppercase Kekulé SMILES. Lowercase aromatic variants check chemical
classification; separate raw PubChem variants check that notation contract.
For raw variants, RDKit aromaticity perception is disabled. These can miss
phenols or classify explicit benzene C=C bonds as alkenes, so use lowercase
aromatic notation for named groups on aromatic molecules. Neither test nor
package changes protonation or tautomers.

## Repeat the checks

`cargo test --manifest-path plugin/Cargo.toml` checks all exact match sets
through the native matcher. `typst compile --root .
tests/substructure-chemistry.typ tests/substructure-chemistry.pdf` checks the
shipped WASM API and the highlight annotations with both `include-atoms`
settings, and renders representative positive and negative controls.

For an independent rerun, use `python3 tests/verify-substructures.py` in an
environment containing RDKit. `--write` regenerates the pinned match sets
after checking the reviewed counts; it does not rewrite those counts. RDKit
is not a runtime or ordinary test dependency. Changes to chemical definitions
require reviewing the counts and rerunning all three checks.

Sources: [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest),
[RDKit substructure matching](https://www.rdkit.org/docs/GettingStartedInPython.html#substructure-searching),
and [RDKit functional-group definitions](https://github.com/rdkit/rdkit/blob/master/Data/Functional_Group_Hierarchy.txt).
Individual PubChem compound URLs are recorded in the JSON fixture.

# SMILES conformance corpus

`smiles-conformance.json` holds categorized SMILES: aromatic systems, salts,
charges, isotopes, atom maps, ring closures, syntax errors, tetrahedral,
double-bond, and extended stereochemistry, long chains, deep branches, cages,
crowded fused rings, macrocycles, and typed-smiles extensions. Each case records a reviewed
expectation:

- `expect`: whether typed-smiles accepts the input, with an `error` fragment
  that every rejection's diagnostic must contain;
- `undepicted`: whether the accepted input has stereochemistry the drawing
  reports as not shown;
- `note`: why typed-smiles deliberately differs from RDKit, required for every
  such difference.

Two fields are generated: `rdkit` (`accept`, `reject`, or `syntax-only` when
RDKit parses the syntax but rejects the chemistry) and `layout`, the geometric
quality of the drawing (`clean`, or any of `distorted-bonds`,
`overlapping-atoms`, `crossing-bonds`). Parsing and drawing are judged
separately: a cage can parse correctly and still draw with crossing bonds.

## Independent checks

`python3 tests/verify-conformance.py` rebuilds every accepted molecule in RDKit
from the drawing alone (coordinates, bond orders, wedge tips, charges,
isotopes, hydrogens, atom maps, and square-planar geometry) and requires the
canonical isomeric SMILES of the input. `--random N` repeats this for N random
atom orders of each molecule; each must round-trip or report its
stereochemistry as undepicted. `--write` regenerates `rdkit` and `layout` after
the reviewed expectations pass. RDKit is a verification tool only.

`cargo test --manifest-path plugin/Cargo.toml` checks every reviewed
expectation and that `clean` layouts stay clean, without RDKit.

Where the OpenSMILES specification is silent, typed-smiles follows RDKit: a
three-coordinate stereocenter's lone pair counts as its last neighbor. Atom
class 0 means "no class", as the specification defines; RDKit keeps `:0` as a
label, so the comparison ignores it.

## Performance

`cargo test --release --manifest-path plugin/Cargo.toml --lib
measure_pipeline_stages -- --ignored --nocapture` times parsing,
kekulization, layout, and JSON serialization natively. `python3
tests/measure-typst-performance.py` times the shipped WASM plugin (first and
later calls) and Typst rendering.
