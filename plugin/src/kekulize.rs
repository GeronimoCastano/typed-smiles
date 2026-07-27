/// Kekulization of aromatic SMILES input (OpenSMILES §3.5).
///
/// Aromatic (lowercase) notation leaves the single/double bond assignment to
/// the reader. This module converts the aromatic bonds of a parsed molecule
/// into an explicit Kekulé structure:
///
///   1. Implicit bonds (and `:` bonds) between aromatic atoms are aromatic;
///      explicit `-` bonds are not (e.g. the biphenyl linker).
///   2. Aromatic bonds outside any ring degrade to single bonds; an aromatic
///      atom left without an aromatic ring bond is an error, since aromatic
///      atoms must be part of a ring.
///   3. Every aromatic atom that still has room for one more bond — by its
///      normal valence, adjusted for charge and declared hydrogens — must
///      receive exactly one double bond. Atoms whose valence is already
///      saturated (pyrrole [nH], furan o, exocyclic C=O carbons, …)
///      contribute a lone pair or nothing instead.
///   4. A perfect matching over those atoms picks the double bonds; the
///      remaining aromatic bonds become single. If no perfect matching
///      exists, the aromatic system cannot be kekulized and the SMILES is
///      rejected.
///
/// Atom and bond indices are untouched, so writing-order references
/// (`show-indices`, highlights, mechanism arrows) stay valid for aromatic
/// and Kekulé input alike.
use crate::graph::{BondOrder, MoleculeGraph};

pub(crate) fn kekulize(
    molecule: &mut MoleculeGraph,
    implicit_bonds: &[bool],
) -> Result<(), String> {
    let has_aromatic_input = molecule.atoms.iter().any(|atom| atom.aromatic)
        || molecule
            .bonds
            .iter()
            .any(|bond| bond.order == BondOrder::Aromatic);
    if !has_aromatic_input {
        return Ok(());
    }

    mark_aromatic_bonds(molecule, implicit_bonds)?;
    replace_non_ring_aromatic_bonds(molecule);
    validate_aromatic_atoms_are_in_rings(molecule)?;

    // The renderer needs the input aromaticity after Kekulé bond orders replace
    // the temporary aromatic order.
    for bond in &mut molecule.bonds {
        if bond.order == BondOrder::Aromatic {
            bond.aromatic = true;
        }
    }

    assign_aromatic_double_bonds(molecule)?;

    for bond in &mut molecule.bonds {
        if bond.order == BondOrder::Aromatic {
            bond.order = BondOrder::Single;
        }
    }
    Ok(())
}

/// Wildcard atoms may sit inside an aromatic ring and bond aromatically to
/// their aromatic neighbors without being aromatic themselves.
fn can_form_aromatic_bond(molecule: &MoleculeGraph, first_atom: usize, second_atom: usize) -> bool {
    let is_aromatic_capable = |atom_index: usize| {
        molecule.atoms[atom_index].aromatic || molecule.atoms[atom_index].symbol == "*"
    };
    is_aromatic_capable(first_atom)
        && is_aromatic_capable(second_atom)
        && (molecule.atoms[first_atom].aromatic || molecule.atoms[second_atom].aromatic)
}

fn mark_aromatic_bonds(
    molecule: &mut MoleculeGraph,
    implicit_bonds: &[bool],
) -> Result<(), String> {
    for bond_index in 0..molecule.bonds.len() {
        let bond = &molecule.bonds[bond_index];
        let endpoints_are_aromatic = can_form_aromatic_bond(molecule, bond.from, bond.to);

        match bond.order {
            BondOrder::Aromatic if !endpoints_are_aromatic => {
                return Err("Aromatic bond ':' must connect two aromatic atoms".to_string());
            }
            BondOrder::Single
                if implicit_bonds.get(bond_index).copied().unwrap_or(false)
                    && endpoints_are_aromatic =>
            {
                molecule.bonds[bond_index].order = BondOrder::Aromatic;
            }
            _ => {}
        }
    }
    Ok(())
}

fn replace_non_ring_aromatic_bonds(molecule: &mut MoleculeGraph) {
    for bond_index in 0..molecule.bonds.len() {
        if molecule.bonds[bond_index].order == BondOrder::Aromatic
            && !is_bond_in_ring(molecule, bond_index)
        {
            molecule.bonds[bond_index].order = BondOrder::Single;
        }
    }
}

/// A bond lies in a ring iff its endpoints stay connected without it.
fn is_bond_in_ring(molecule: &MoleculeGraph, excluded_bond: usize) -> bool {
    let bond = &molecule.bonds[excluded_bond];
    let mut visited = vec![false; molecule.n_atoms()];
    let mut pending_atoms = vec![bond.from];
    visited[bond.from] = true;

    while let Some(atom_index) = pending_atoms.pop() {
        for &(neighbor, bond_index) in &molecule.adj[atom_index] {
            if bond_index == excluded_bond || visited[neighbor] {
                continue;
            }
            if neighbor == bond.to {
                return true;
            }
            visited[neighbor] = true;
            pending_atoms.push(neighbor);
        }
    }
    false
}

fn validate_aromatic_atoms_are_in_rings(molecule: &MoleculeGraph) -> Result<(), String> {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        if atom.aromatic
            && !molecule.adj[atom_index]
                .iter()
                .any(|&(_, bond_index)| molecule.bonds[bond_index].order == BondOrder::Aromatic)
        {
            return Err(format!(
                "Aromatic atom '{}' (atom {atom_index}) must be part of an aromatic ring",
                atom.symbol.to_lowercase()
            ));
        }
    }
    Ok(())
}

/// Normal valences of the elements that can be written aromatic.
fn aromatic_valence(symbol: &str) -> Option<i16> {
    match symbol {
        "B" => Some(3),
        "C" => Some(4),
        "N" | "P" | "As" => Some(3),
        "O" | "S" | "Se" => Some(2),
        _ => None,
    }
}

/// Whether an aromatic atom must receive one double bond during kekulization.
///
/// Counting aromatic bonds as single, the atom needs a double bond exactly
/// when its bonds plus declared hydrogens leave room under the charge-adjusted
/// normal valence (one π electron in the OpenSMILES table); a saturated atom
/// contributes a lone pair or an empty orbital instead (zero or two).
fn aromatic_atom_needs_double_bond(molecule: &MoleculeGraph, atom_index: usize) -> bool {
    let atom = &molecule.atoms[atom_index];
    if !atom.aromatic {
        return false;
    }
    let Some(valence) = aromatic_valence(&atom.symbol) else {
        return false;
    };
    let sigma_bond_order: i16 = molecule.adj[atom_index]
        .iter()
        .map(|&(_, bond_index)| match molecule.bonds[bond_index].order {
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
        })
        .sum();
    valence + atom.charge as i16 - atom.hcount as i16 - sigma_bond_order >= 1
}

fn assign_aromatic_double_bonds(molecule: &mut MoleculeGraph) -> Result<(), String> {
    let atom_count = molecule.n_atoms();
    let requires_double_bond: Vec<bool> = (0..atom_count)
        .map(|atom_index| aromatic_atom_needs_double_bond(molecule, atom_index))
        .collect();
    let can_accept_double_bond: Vec<bool> = molecule
        .atoms
        .iter()
        .map(|atom| atom.symbol == "*")
        .collect();
    let mut matched_partner = vec![None; atom_count];

    if !find_aromatic_matching(
        molecule,
        &requires_double_bond,
        &can_accept_double_bond,
        &mut matched_partner,
        0,
    ) {
        return Err(
            "Cannot kekulize aromatic system: no valid alternating double-bond assignment \
             exists (check hydrogen counts, e.g. pyrrole is c1cc[nH]c1)"
                .to_string(),
        );
    }

    for (atom_index, partner) in matched_partner.iter().copied().enumerate() {
        let Some(partner_index) = partner else {
            continue;
        };
        if atom_index >= partner_index {
            continue;
        }

        let bond_index = molecule.adj[atom_index]
            .iter()
            .find_map(|&(neighbor, bond_index)| {
                (neighbor == partner_index
                    && molecule.bonds[bond_index].order == BondOrder::Aromatic)
                    .then_some(bond_index)
            })
            .expect("matched atoms share an aromatic bond");
        molecule.bonds[bond_index].order = BondOrder::Double;
    }
    Ok(())
}

/// Backtracking search for a matching that pairs every atom in `needs` with an
/// aromatic-bonded partner. Aromatic systems are small, so exhaustive search
/// with backtracking is fast enough.
fn find_aromatic_matching(
    molecule: &MoleculeGraph,
    requires_double_bond: &[bool],
    can_accept_double_bond: &[bool],
    matched_partner: &mut [Option<usize>],
    start: usize,
) -> bool {
    let unmatched_atom = (start..requires_double_bond.len()).find(|&atom_index| {
        requires_double_bond[atom_index] && matched_partner[atom_index].is_none()
    });
    let Some(atom_index) = unmatched_atom else {
        return true;
    };

    for &(neighbor, bond_index) in &molecule.adj[atom_index] {
        if molecule.bonds[bond_index].order != BondOrder::Aromatic
            || matched_partner[neighbor].is_some()
            || !(requires_double_bond[neighbor] || can_accept_double_bond[neighbor])
        {
            continue;
        }

        matched_partner[atom_index] = Some(neighbor);
        matched_partner[neighbor] = Some(atom_index);
        if find_aromatic_matching(
            molecule,
            requires_double_bond,
            can_accept_double_bond,
            matched_partner,
            atom_index + 1,
        ) {
            return true;
        }
        matched_partner[atom_index] = None;
        matched_partner[neighbor] = None;
    }
    false
}
