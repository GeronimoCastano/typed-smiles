#!/usr/bin/env python3
"""Independently check the pinned PubChem corpus with RDKit.

Run from any directory with a Python environment containing RDKit. --write
updates exact match sets after checking the separately reviewed group counts.
RDKit is a verification tool only, not a package or test-suite dependency.
"""

import argparse
import json
import re
from pathlib import Path

from rdkit import Chem, rdBase


def match_sets(molecule, query):
    matches = set()
    for mapping in molecule.GetSubstructMatches(query, uniquify=False, maxMatches=10000):
        atoms = tuple(sorted(mapping))
        bonds = tuple(sorted(
            tuple(sorted((mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()])))
            for bond in query.GetBonds()
        ))
        matches.add((atoms, bonds))
    return [dict(atoms=list(atoms), bonds=[list(pair) for pair in bonds])
            for atoms, bonds in sorted(matches)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    fixture = root / 'tests/fixtures/substructure-pubchem.json'
    corpus = json.loads(fixture.read_text())
    registry = (root / 'src/substructure.typ').read_text().split('#let _highlight-palette')[0]
    patterns = dict(re.findall(r'^  ([\w-]+): "(.+)",$', registry, re.MULTILINE))
    queries = {name: Chem.MolFromSmarts(pattern) for name, pattern in patterns.items()}
    assert all(query is not None for query in queries.values())
    checks = 0
    for case in corpus['cases']:
        aromatic = Chem.MolFromSmiles(case['pubchem-smiles'])
        smiles = Chem.MolToSmiles(aromatic, canonical=False)
        variants = [('aromatic', smiles)]
        if smiles != case['pubchem-smiles']:
            variants.append(('pubchem', case['pubchem-smiles']))
        generated = []
        for notation, smiles in variants:
            # The package retains explicit aromatic flags rather than perceiving
            # aromaticity in uppercase Kekule SMILES. Test that contract too.
            if notation == 'aromatic':
                molecule = Chem.MolFromSmiles(smiles)
            else:
                molecule = Chem.MolFromSmiles(smiles, sanitize=False)
                Chem.SanitizeMol(molecule, sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                ))
            found = {name: match_sets(molecule, query) for name, query in queries.items()}
            checks += len(queries)
            if notation == 'aromatic':
                counts = {name: len(hits) for name, hits in found.items() if hits}
                assert counts == case['reviewed-counts'], (
                    f"{case['name']} (CID {case['cid']}): {counts} != {case['reviewed-counts']}"
                )
            generated.append(dict(notation=notation, smiles=smiles,
                                  matches={name: hits for name, hits in found.items() if hits}))
        if args.write:
            case['variants'] = generated
        else:
            assert generated == case['variants'], f"RDKit match sets changed: {case['name']}"
    if args.write:
        corpus['patterns'] = patterns
        corpus['rdkit-version'] = rdBase.rdkitVersion
        fixture.write_text(json.dumps(corpus, indent=2) + '\n')
    else:
        assert corpus['patterns'] == patterns, 'Pattern registry changed; review counts before --write'
    print(f"RDKit {rdBase.rdkitVersion}: {len(corpus['cases'])} compounds, {checks} group/notation checks passed.")


if __name__ == '__main__':
    main()
