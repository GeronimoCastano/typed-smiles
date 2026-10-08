//! Ring perception: the smallest set of smallest rings and the ring systems
//! they form.

use std::collections::{HashSet, VecDeque};

use crate::graph::MoleculeGraph;

pub(crate) fn ring_bond_set(molecule: &MoleculeGraph, rings: &[Vec<usize>]) -> HashSet<usize> {
    let mut ring_bonds = HashSet::new();
    for ring in rings {
        for (first_atom, second_atom) in ring_edges(ring) {
            if let Some(bond_index) = molecule.bond_between(first_atom, second_atom) {
                ring_bonds.insert(bond_index);
            }
        }
    }
    ring_bonds
}

/// Returns cycles as lists of atom indices (the ring path).
pub(crate) fn find_rings(molecule: &MoleculeGraph) -> Vec<Vec<usize>> {
    let target_count = cycle_rank(molecule);
    if target_count == 0 {
        return Vec::new();
    }

    let bit_words = molecule.bonds.len().div_ceil(64);
    let mut seen_cycles: HashSet<Vec<u64>> = HashSet::new();
    let mut candidates: Vec<(Vec<usize>, Vec<u64>)> = Vec::new();

    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        if let Some(ring) = shortest_path_excluding_bond(molecule, bond.from, bond.to, bond_index) {
            if ring.len() < 3 {
                continue;
            }
            let bits = ring_bond_bits(molecule, &ring, bit_words);
            if seen_cycles.insert(bits.clone()) {
                candidates.push((ring, bits));
            }
        }
    }

    candidates.sort_by(|(first_ring, first_bits), (second_ring, second_bits)| {
        first_ring
            .len()
            .cmp(&second_ring.len())
            .then_with(|| first_bits.cmp(second_bits))
    });

    let mut basis: Vec<Vec<u64>> = Vec::new();
    let mut rings = Vec::new();
    for (ring, bits) in candidates {
        if add_independent_cycle(&mut basis, bits) {
            rings.push(ring);
            if rings.len() == target_count {
                break;
            }
        }
    }

    rings
}

fn cycle_rank(molecule: &MoleculeGraph) -> usize {
    if molecule.n_atoms() == 0 {
        return 0;
    }

    let mut seen = vec![false; molecule.n_atoms()];
    let mut components = 0;
    for start in 0..molecule.n_atoms() {
        if seen[start] {
            continue;
        }
        components += 1;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(atom) = stack.pop() {
            for &(neighbor, _) in &molecule.adj[atom] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
    }

    molecule.bonds.len() + components - molecule.n_atoms()
}

/// BFS from `from` to `to`, intentionally skipping one bond so the path plus
/// that skipped bond is a ring candidate.
fn shortest_path_excluding_bond(
    molecule: &MoleculeGraph,
    from: usize,
    to: usize,
    excluded_bond: usize,
) -> Option<Vec<usize>> {
    let atom_count = molecule.n_atoms();
    let mut bfs_parent = vec![None; atom_count];

    bfs_parent[from] = Some(from); // root sentinel
    let mut queue = VecDeque::new();
    queue.push_back(from);

    while let Some(atom_index) = queue.pop_front() {
        for &(neighbor, bond_index) in &molecule.adj[atom_index] {
            if bond_index == excluded_bond {
                continue;
            }
            if bfs_parent[neighbor].is_some() {
                continue;
            }
            bfs_parent[neighbor] = Some(atom_index);

            if neighbor == to {
                let mut path = Vec::new();
                let mut current_atom = to;
                loop {
                    path.push(current_atom);
                    if current_atom == from {
                        break;
                    }
                    current_atom = bfs_parent[current_atom]?;
                }
                path.reverse();
                return Some(path);
            }
            queue.push_back(neighbor);
        }
    }

    None
}

fn ring_bond_bits(molecule: &MoleculeGraph, ring: &[usize], bit_words: usize) -> Vec<u64> {
    let mut bits = vec![0_u64; bit_words];
    for (first_atom, second_atom) in ring_edges(ring) {
        if let Some(bond_index) = molecule.bond_between(first_atom, second_atom) {
            bits[bond_index / 64] |= 1_u64 << (bond_index % 64);
        }
    }
    bits
}

pub(crate) fn ring_edges(ring: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    ring.iter()
        .copied()
        .zip(ring.iter().copied().cycle().skip(1))
        .take(ring.len())
}

fn add_independent_cycle(basis: &mut Vec<Vec<u64>>, bits: Vec<u64>) -> bool {
    let mut candidate = bits;
    for existing in basis.iter() {
        if let Some(pivot) = pivot_bit(existing) {
            if bit_is_set(&candidate, pivot) {
                xor_assign(&mut candidate, existing);
            }
        }
    }

    let Some(pivot) = pivot_bit(&candidate) else {
        return false;
    };

    for existing in basis.iter_mut() {
        if bit_is_set(existing, pivot) {
            xor_assign(existing, &candidate);
        }
    }
    basis.push(candidate);
    basis.sort_by_key(|bits| pivot_bit(bits).unwrap_or(usize::MAX));
    true
}

fn pivot_bit(bits: &[u64]) -> Option<usize> {
    for (word_idx, word) in bits.iter().enumerate() {
        if *word != 0 {
            return Some(word_idx * 64 + word.trailing_zeros() as usize);
        }
    }
    None
}

fn bit_is_set(bits: &[u64], bit: usize) -> bool {
    bits[bit / 64] & (1_u64 << (bit % 64)) != 0
}

fn xor_assign(lhs: &mut [u64], rhs: &[u64]) {
    for (a, b) in lhs.iter_mut().zip(rhs.iter()) {
        *a ^= *b;
    }
}

pub(crate) fn rings_share_edge(first_ring: &[usize], second_ring: &[usize]) -> bool {
    shared_atoms(first_ring, second_ring)
        .windows(2)
        .any(|atom_pair| {
            ring_has_edge(first_ring, atom_pair[0], atom_pair[1])
                && ring_has_edge(second_ring, atom_pair[0], atom_pair[1])
        })
}

pub(crate) fn shared_atoms(first_ring: &[usize], second_ring: &[usize]) -> Vec<usize> {
    let mut atoms: Vec<usize> = first_ring
        .iter()
        .copied()
        .filter(|atom| second_ring.contains(atom))
        .collect();
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

pub(crate) fn ring_has_edge(ring: &[usize], first_atom: usize, second_atom: usize) -> bool {
    ring_edges(ring).any(|(edge_start, edge_end)| {
        (edge_start == first_atom && edge_end == second_atom)
            || (edge_start == second_atom && edge_end == first_atom)
    })
}

/// Indices of every ring in the ring system that contains `ring_index`.
///
/// Rings sharing at least one atom (fused, bridged, or spiro) belong to the
/// same system, so the whole system can be laid out as one unit.
pub(crate) fn ring_system_rings(rings: &[Vec<usize>], ring_index: usize) -> Vec<usize> {
    let mut in_system = vec![false; rings.len()];
    in_system[ring_index] = true;
    let mut pending = vec![ring_index];
    while let Some(current_ring) = pending.pop() {
        for (other_index, other_ring) in rings.iter().enumerate() {
            if !in_system[other_index] && !shared_atoms(&rings[current_ring], other_ring).is_empty()
            {
                in_system[other_index] = true;
                pending.push(other_index);
            }
        }
    }
    (0..rings.len()).filter(|&index| in_system[index]).collect()
}

/// Atoms of the given rings, sorted and without repeats.
pub(crate) fn ring_system_atoms(rings: &[Vec<usize>], system_rings: &[usize]) -> Vec<usize> {
    let mut atoms: Vec<usize> = system_rings
        .iter()
        .flat_map(|&ring_index| rings[ring_index].iter().copied())
        .collect();
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rings_of(smiles: &str) -> Vec<Vec<usize>> {
        find_rings(&crate::parse_molecule(smiles).unwrap())
    }

    #[test]
    fn fused_bridged_and_spiro_rings_form_one_system() {
        for smiles in ["C1CC2CCCCC2CC1", "C1CC2CCC1C2", "C1CCC2(C1)CCCC2"] {
            let rings = rings_of(smiles);
            assert_eq!(ring_system_rings(&rings, 0), vec![0, 1], "{smiles}");
        }
    }

    #[test]
    fn rings_joined_by_a_chain_form_separate_systems() {
        let rings = rings_of("C1CCCCC1CCC1CCCCC1");
        assert_eq!(ring_system_rings(&rings, 0), vec![0]);
        assert_eq!(ring_system_rings(&rings, 1), vec![1]);
    }

    #[test]
    fn system_atoms_are_sorted_without_repeats() {
        let rings = rings_of("C1CC2CCC1C2");
        assert_eq!(
            ring_system_atoms(&rings, &[0, 1]),
            (0..7).collect::<Vec<_>>()
        );
    }
}
