#!/usr/bin/env python3
"""Independently check the SMILES conformance corpus with RDKit.

Run from any directory with a Python environment containing RDKit and a Rust
toolchain. The script builds the native `layout_json` example, then checks for
every case in tests/fixtures/smiles-conformance.json:

  * that typed-smiles accepts or rejects it as the reviewed `expect` says;
  * whether RDKit accepts its syntax and chemistry (recorded as `rdkit`);
  * that RDKit, rebuilding the molecule from the drawing alone (coordinates,
    bond orders, wedge tips, charges, isotopes, hydrogens, atom maps, and
    square-planar geometry), recovers the input's canonical isomeric SMILES,
    unless the case expects undepicted stereo;
  * the geometric layout class (recorded as `layout`), kept separate from
    parsing so drawing problems are not mistaken for parser problems.

--write records `rdkit` and `layout` after the reviewed expectations pass.
--random N additionally redraws N random atom orders of every round-trip case;
each must round-trip or report its stereochemistry as undepicted.
RDKit is a verification tool only, not a package or test-suite dependency.
"""

import argparse
import json
import math
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Geometry import Point3D

RDLogger.DisableLog('rdApp.*')
ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / 'tests/fixtures/smiles-conformance.json'
TYPED_SMILES_EXTENSIONS = ('{', '!')


def typed_smiles_layouts(smiles_list):
    """Layout JSON (or the diagnostic) for each SMILES from the native build."""
    subprocess.run(['cargo', 'build', '--quiet', '--release', '--example', 'layout_json',
                    '--manifest-path', str(ROOT / 'plugin/Cargo.toml')], check=True)
    binary = ROOT / 'plugin/target/release/examples/layout_json'
    lines = ''.join(f'{index}\t{smiles}\n' for index, smiles in enumerate(smiles_list))
    output = subprocess.run([binary], input=lines, capture_output=True, text=True, check=True)
    return [json.loads(line) for line in output.stdout.splitlines()]


def rdkit_verdict(smiles):
    if Chem.MolFromSmiles(smiles, sanitize=False) is None:
        return 'reject'
    return 'accept' if Chem.MolFromSmiles(smiles) is not None else 'syntax-only'


def drawn_atoms_and_bonds(layout):
    atoms = [index for index, atom in enumerate(layout['atoms']) if not atom.get('virtual_h')]
    bonds = [bond for bond in layout['bonds'] if not bond.get('virtual_bond')]
    return atoms, bonds


def layout_class(layout):
    """Clean, or the geometric problems of a drawing in bond-length units."""
    atoms, bonds = drawn_atoms_and_bonds(layout)
    position = {index: layout['atoms'][index]['pos'] for index in atoms}
    problems = set()
    for bond in bonds:
        length = math.dist(*((position[end]['x'], position[end]['y']) for end in (bond['from'], bond['to'])))
        if not 0.75 <= length <= 1.35:
            problems.add('distorted-bonds')
    bonded = {frozenset((bond['from'], bond['to'])) for bond in bonds}
    for first_position, first in enumerate(atoms):
        for second in atoms[first_position + 1:]:
            if frozenset((first, second)) in bonded:
                continue
            if math.dist((position[first]['x'], position[first]['y']),
                         (position[second]['x'], position[second]['y'])) < 0.5:
                problems.add('overlapping-atoms')
    segments = [((position[b['from']]['x'], position[b['from']]['y']),
                 (position[b['to']]['x'], position[b['to']]['y']), {b['from'], b['to']}) for b in bonds]
    for first_index, (a, b, first_atoms) in enumerate(segments):
        for c, d, second_atoms in segments[first_index + 1:]:
            if first_atoms & second_atoms:
                continue
            if segments_cross(a, b, c, d):
                problems.add('crossing-bonds')
    return ','.join(sorted(problems)) or 'clean'


def segments_cross(a, b, c, d):
    def orientation(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    return (orientation(a, b, c) * orientation(a, b, d) < -1e-9
            and orientation(c, d, a) * orientation(c, d, b) < -1e-9)


def molecule_from_drawing(layout):
    """Rebuilds the molecule from what a reader of the drawing can see."""
    atoms, bonds = drawn_atoms_and_bonds(layout)
    atom_outputs = layout['atoms']
    square_planar = [index for index in atoms if atom_outputs[index].get('chirality') == 'square_planar']
    molecule = Chem.RWMol()
    index_map = {}
    for index in atoms:
        atom = atom_outputs[index]
        rdkit_atom = Chem.Atom(0 if atom['symbol'] == '*' else atom['symbol'])
        rdkit_atom.SetFormalCharge(atom['charge'])
        rdkit_atom.SetIsotope(atom.get('isotope', 0))
        rdkit_atom.SetAtomMapNum(atom.get('atom_map', 0))
        hydrogens = atom['hcount'] + atom.get('implicit_h', 0)
        if atom.get('stereo_h', 'none') != 'none':
            hydrogens -= 1
        rdkit_atom.SetNumExplicitHs(hydrogens)
        rdkit_atom.SetNoImplicit(True)
        index_map[index] = molecule.AddAtom(rdkit_atom)
    positions = {index_map[index]: atom_outputs[index]['pos'] for index in atoms}

    def square_planar_order(bond):
        # Bonds of a square-planar center are added in drawn angular order, so
        # RDKit's neighbor order traces the drawn square (@SP1, the U shape).
        for center in square_planar:
            if center in (bond['from'], bond['to']):
                other = bond['to'] if bond['from'] == center else bond['from']
                offset_x = atom_outputs[other]['pos']['x'] - atom_outputs[center]['pos']['x']
                offset_y = atom_outputs[other]['pos']['y'] - atom_outputs[center]['pos']['y']
                return (0, center, math.atan2(offset_y, offset_x))
        return (1, 0, 0.0)

    order_type = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
                  3: Chem.BondType.TRIPLE, 4: Chem.BondType.QUADRUPLE}
    wedges = []
    for bond in sorted(bonds, key=square_planar_order):
        begin, end = bond['from'], bond['to']
        is_wedge = bond['stereo'] in ('wedge_up', 'wedge_down')
        if is_wedge and bond.get('stereo_tip') == end:
            begin, end = end, begin
        molecule.AddBond(index_map[begin], index_map[end], order_type[bond['order']])
        if is_wedge:
            wedges.append((index_map[begin], index_map[end], bond['stereo']))
    for index in atoms:
        atom = atom_outputs[index]
        stereo = atom.get('stereo_h', 'none')
        if stereo == 'none':
            continue
        hydrogen = molecule.AddAtom(Chem.Atom(1))
        direction = atom['stereo_h_dir']
        positions[hydrogen] = {'x': atom['pos']['x'] + direction['x'], 'y': atom['pos']['y'] + direction['y']}
        molecule.AddBond(index_map[index], hydrogen, Chem.BondType.SINGLE)
        wedges.append((index_map[index], hydrogen, stereo))

    conformer = Chem.Conformer(molecule.GetNumAtoms())
    conformer.Set3D(False)
    for index, position in positions.items():
        conformer.SetAtomPosition(index, Point3D(position['x'], position['y'], 0.0))
    molecule.AddConformer(conformer, assignId=True)
    for begin, end, stereo in wedges:
        direction = Chem.BondDir.BEGINWEDGE if stereo == 'wedge_up' else Chem.BondDir.BEGINDASH
        molecule.GetBondBetweenAtoms(begin, end).SetBondDir(direction)
    # A drawn double bond only claims cis/trans when the input marked both ends.
    directed = [(bond['from'], bond['to']) for bond in bonds if bond.get('direction', 'none') != 'none']
    for bond in bonds:
        if bond['order'] != 2:
            continue
        ends = (bond['from'], bond['to'])
        both_marked = all(any(end in pair and other not in pair for pair in directed)
                          for end, other in (ends, ends[::-1]))
        if not both_marked:
            molecule.GetBondBetweenAtoms(index_map[ends[0]], index_map[ends[1]]).SetBondDir(Chem.BondDir.EITHERDOUBLE)
    for center in square_planar:
        rdkit_atom = molecule.GetAtomWithIdx(index_map[center])
        if rdkit_atom.GetDegree() == 4:
            rdkit_atom.SetChiralTag(Chem.ChiralType.CHI_SQUAREPLANAR)
            rdkit_atom.SetUnsignedProp('_chiralPermutation', 1)
    molecule = molecule.GetMol()
    Chem.SanitizeMol(molecule)
    Chem.AssignChiralTypesFromBondDirs(molecule)
    Chem.DetectBondStereoChemistry(molecule, molecule.GetConformer())
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    for bond in molecule.GetBonds():
        if bond.GetBondDir() == Chem.BondDir.EITHERDOUBLE:
            bond.SetStereo(Chem.BondStereo.STEREONONE)
    return Chem.MolToSmiles(Chem.RemoveHs(molecule))


def canonical_reference(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    # OpenSMILES defines atom class 0 as no class; RDKit keeps `:0` as a label.
    for atom in molecule.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule)


def round_trip_failure(smiles, record):
    expected = canonical_reference(smiles)
    drawn = molecule_from_drawing(record['layout'])
    return None if drawn == expected else f'RDKit reads {expected}, the drawing shows {drawn}'


def can_round_trip(case, record):
    return ('layout' in record and not case.get('undepicted')
            and not any(marker in case['smiles'] for marker in TYPED_SMILES_EXTENSIONS)
            and Chem.MolFromSmiles(case['smiles']) is not None)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--random', type=int, default=0, metavar='N')
    args = parser.parse_args()
    corpus = json.loads(FIXTURE.read_text())
    cases = corpus['cases']
    records = typed_smiles_layouts([case['smiles'] for case in cases])
    failures = []
    round_trips = 0
    for case, record in zip(cases, records):
        label = f"{case['category']}/{case['name']}"
        accepted = 'layout' in record
        if accepted != (case['expect'] == 'accept'):
            failures.append(f"{label}: expected {case['expect']}, got {record.get('error', 'accept')}")
            continue
        if not accepted and case.get('error') not in record['error']:
            failures.append(f"{label}: diagnostic lacks {case.get('error')!r}: {record['error']}")
        if accepted:
            reported = bool(record['layout']['undepicted_stereo'])
            if reported != bool(case.get('undepicted')):
                failures.append(f"{label}: undepicted stereo {record['layout']['undepicted_stereo']}")
        verdict = rdkit_verdict(case['smiles'])
        layout = layout_class(record['layout']) if accepted else None
        if verdict.startswith('accept') != accepted and verdict != 'syntax-only' and 'note' not in case:
            failures.append(f"{label}: RDKit {verdict}s it; explain the difference in a note")
        if can_round_trip(case, record):
            round_trips += 1
            failure = round_trip_failure(case['smiles'], record)
            if failure:
                failures.append(f'{label}: {failure}')
        if args.write:
            case['rdkit'] = verdict
            if layout:
                case['layout'] = layout
        else:
            if case.get('rdkit') != verdict:
                failures.append(f"{label}: RDKit verdict changed to {verdict}")
            if accepted and case.get('layout') != layout:
                failures.append(f"{label}: layout class changed to {layout}")

    random_checks = 0
    reported_limitations = 0
    if args.random:
        variants, origins = [], []
        for case, record in zip(cases, records):
            if not can_round_trip(case, record):
                continue
            molecule = Chem.RWMol(Chem.MolFromSmiles(case['smiles']))
            # RDKit writes metal-amine bonds as dative arrows, an RDKit extension.
            for bond in molecule.GetBonds():
                if bond.GetBondType() == Chem.BondType.DATIVE:
                    bond.SetBondType(Chem.BondType.SINGLE)
            for _ in range(args.random):
                variants.append(Chem.MolToSmiles(molecule, doRandom=True, canonical=False))
                origins.append(case)
        for smiles, case, record in zip(variants, origins, typed_smiles_layouts(variants)):
            random_checks += 1
            label = f"{case['category']}/{case['name']} as {smiles}"
            if 'layout' not in record:
                failures.append(f"{label}: {record['error']}")
            elif record['layout']['undepicted_stereo']:
                # The drawing reports what it cannot show rather than showing it wrongly.
                reported_limitations += 1
            else:
                failure = round_trip_failure(smiles, record)
                if failure:
                    failures.append(f'{label}: {failure}')

    if failures:
        raise SystemExit('\n'.join(failures))
    if args.write:
        corpus['rdkit-version'] = rdBase.rdkitVersion
        FIXTURE.write_text(json.dumps(corpus, indent=2, ensure_ascii=False) + '\n')
    layout_problems = sum(1 for case in cases if case.get('layout', 'clean') != 'clean')
    print(f"RDKit {rdBase.rdkitVersion}: {len(cases)} cases, {round_trips} drawing round trips, "
          f"{random_checks} random-order drawings checked ({reported_limitations} reported as "
          f"undepicted, none drawn wrongly); {layout_problems} accepted cases have recorded "
          f"layout problems.")


if __name__ == '__main__':
    main()
