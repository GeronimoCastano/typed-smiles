//! Drawing templates for ring skeletons without a clean generic depiction.
//!
//! Bicyclo[2.2.2]octane, bicyclo[1.1.1]pentane, adamantane, and cubane cannot
//! be assembled from polygons and bridges without atoms landing on each other,
//! so chemists draw them in fixed perspective views. Small bicyclics whose
//! one-atom bridge carries substituents, such as camphor and cocaine, are
//! drawn in perspective too, because the flat drawing puts that bridge inside
//! the ring where its substituents would cross ring bonds. Each template
//! stores one such view in bond-length units. Porphyrins use the same
//! mechanism for their square, fourfold-symmetric drawing. Templates match
//! the skeleton's shape alone, so heteroatoms, unsaturation, and fused rings
//! (quinuclidine, DABCO, triptycene) reuse the same drawing.

use crate::graph::MoleculeGraph;
use crate::render::Vec2;

struct CageTemplate {
    bonds: &'static [(usize, usize)],
    coordinates: &'static [(f64, f64)],
    applicability: TemplateApplicability,
}

/// When a matched template should be used instead of the generic layout.
enum TemplateApplicability {
    /// The skeleton has no clean generic drawing.
    Always,
    /// Only when the given template atom carries substituents or further
    /// rings, which the generic drawing would place inside the ring.
    SubstitutedAt(usize),
    /// Only when the given template atoms carry nothing outside the
    /// template, because they sit inside the drawing.
    UnsubstitutedAt(&'static [usize]),
}

/// Four regular pyrroles pointing their nitrogens inward, joined by meso
/// carbons: the square depiction used for porphyrins, which no hexagon
/// lattice outline can reproduce.
const PORPHINE: CageTemplate = CageTemplate {
    bonds: &[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 0),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 5),
        (10, 11),
        (11, 12),
        (12, 13),
        (13, 14),
        (14, 10),
        (15, 16),
        (16, 17),
        (17, 18),
        (18, 19),
        (19, 15),
        (20, 4),
        (20, 6),
        (21, 9),
        (21, 11),
        (22, 14),
        (22, 16),
        (23, 19),
        (23, 1),
    ],
    coordinates: &[
        (0.000, 1.474),
        (0.809, 2.062),
        (0.500, 3.013),
        (-0.500, 3.013),
        (-0.809, 2.062),
        (-1.474, 0.000),
        (-2.062, 0.809),
        (-3.013, 0.500),
        (-3.013, -0.500),
        (-2.062, -0.809),
        (0.000, -1.474),
        (-0.809, -2.062),
        (-0.500, -3.013),
        (0.500, -3.013),
        (0.809, -2.062),
        (1.474, 0.000),
        (2.062, -0.809),
        (3.013, -0.500),
        (3.013, 0.500),
        (2.062, 0.809),
        (-1.763, 1.763),
        (-1.763, -1.763),
        (1.763, -1.763),
        (1.763, 1.763),
    ],
    applicability: TemplateApplicability::Always,
};

/// Perspective view with the apex methylene on top and one bond passing in
/// front of another.
const ADAMANTANE: CageTemplate = CageTemplate {
    bonds: &[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (7, 9),
        (5, 0),
        (8, 1),
        (9, 3),
    ],
    coordinates: &[
        (0.474, -1.358),
        (-0.666, -0.995),
        (-0.561, -0.184),
        (0.101, 0.810),
        (1.236, 0.437),
        (1.142, -0.371),
        (0.561, 0.184),
        (-0.578, 0.555),
        (-1.236, -0.437),
        (-0.474, 1.358),
    ],
    applicability: TemplateApplicability::Always,
};

/// Two offset squares joined at their corners, the usual drawing of a cube.
const CUBANE: CageTemplate = CageTemplate {
    bonds: &[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (3, 0),
        (5, 0),
        (6, 1),
        (7, 2),
        (7, 4),
    ],
    coordinates: &[
        (-0.697, 0.604),
        (-0.977, -0.027),
        (0.076, 0.394),
        (0.355, 1.025),
        (0.977, 0.027),
        (-0.076, -0.394),
        (-0.355, -1.025),
        (0.697, -0.604),
    ],
    applicability: TemplateApplicability::Always,
};

/// A hexagon through both bridgeheads with the third bridge drawn as a
/// zigzag across its middle.
const BICYCLO_2_2_2_OCTANE: CageTemplate = CageTemplate {
    bonds: &[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (5, 0),
        (7, 2),
    ],
    coordinates: &[
        (0.915, -0.528),
        (0.915, 0.528),
        (0.000, 1.057),
        (-0.915, 0.528),
        (-0.915, -0.528),
        (0.000, -1.057),
        (-0.370, -0.264),
        (0.370, 0.264),
    ],
    applicability: TemplateApplicability::Always,
};

/// Bicyclo[3.2.1]octane seen from the front: the one-atom bridge on top,
/// the three-atom bridge in front, and the two-atom bridge behind.
const BICYCLO_3_2_1_OCTANE_PERSPECTIVE: CageTemplate = CageTemplate {
    bonds: &[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 0),
        (0, 7),
        (7, 4),
    ],
    coordinates: &[
        (-1.055, 0.448),
        (-0.897, -0.554),
        (0.000, -1.029),
        (0.897, -0.554),
        (1.055, 0.448),
        (0.422, 0.026),
        (-0.422, 0.026),
        (0.000, 1.187),
    ],
    applicability: TemplateApplicability::SubstitutedAt(7),
};

/// Bicyclo[2.2.1]heptane seen from the front: the one-atom bridge on top
/// and one two-atom bridge folded behind the other.
const BICYCLO_2_2_1_HEPTANE_PERSPECTIVE: CageTemplate = CageTemplate {
    bonds: &[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
        (0, 6),
        (6, 3),
    ],
    coordinates: &[
        (-1.002, 0.306),
        (-0.580, -0.665),
        (0.580, -0.665),
        (1.002, 0.306),
        (0.401, -0.137),
        (-0.401, -0.137),
        (0.000, 0.991),
    ],
    applicability: TemplateApplicability::SubstitutedAt(6),
};

/// Bicyclo[3.1.1]heptane as a hexagon with the second one-atom bridge at
/// its center, which shows the four-membered ring without any distortion.
const BICYCLO_3_1_1_HEPTANE: CageTemplate = CageTemplate {
    bonds: &[
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
        (0, 6),
        (6, 4),
    ],
    coordinates: &[
        (-0.866, 0.500),
        (-0.866, -0.500),
        (0.000, -1.000),
        (0.866, -0.500),
        (0.866, 0.500),
        (0.000, 1.000),
        (0.000, 0.000),
    ],
    applicability: TemplateApplicability::UnsubstitutedAt(&[6]),
};

/// A diamond through both bridgeheads with the third bridge folded inside.
const BICYCLO_1_1_1_PENTANE: CageTemplate = CageTemplate {
    bonds: &[(0, 1), (1, 2), (2, 3), (3, 4), (3, 0), (4, 1)],
    coordinates: &[
        (0.000, 0.771),
        (-0.732, -0.059),
        (0.000, -0.888),
        (0.732, -0.059),
        (0.000, 0.234),
    ],
    applicability: TemplateApplicability::Always,
};

/// Larger cages come first so a skeleton is drawn with the most specific
/// template that fits.
const CAGE_TEMPLATES: [&CageTemplate; 8] = [
    &PORPHINE,
    &ADAMANTANE,
    &CUBANE,
    &BICYCLO_2_2_2_OCTANE,
    &BICYCLO_3_2_1_OCTANE_PERSPECTIVE,
    &BICYCLO_2_2_1_HEPTANE_PERSPECTIVE,
    &BICYCLO_3_1_1_HEPTANE,
    &BICYCLO_1_1_1_PENTANE,
];

/// Positions for the atoms of the first cage template found in a ring
/// system, or `None` when the system contains no templated cage.
///
/// A template matches when its skeleton appears among the system's atoms with
/// exactly the template's bonds between them. Among the symmetric ways to
/// fit it, atoms that carry substituents or further rings are given the
/// template's outermost positions, where those groups have room.
pub(crate) fn cage_template_positions(
    molecule: &MoleculeGraph,
    system_atoms: &[usize],
) -> Option<Vec<(usize, Vec2)>> {
    CAGE_TEMPLATES.iter().find_map(|template| {
        let assignment = best_template_assignment(molecule, system_atoms, template)?;
        Some(
            assignment
                .iter()
                .enumerate()
                .map(|(template_atom, &molecule_atom)| {
                    let (x, y) = template.coordinates[template_atom];
                    (molecule_atom, Vec2::new(x, y))
                })
                .collect(),
        )
    })
}

fn best_template_assignment(
    molecule: &MoleculeGraph,
    system_atoms: &[usize],
    template: &CageTemplate,
) -> Option<Vec<usize>> {
    let template_size = template.coordinates.len();
    if system_atoms.len() < template_size {
        return None;
    }
    let template_neighbors = template_adjacency(template);
    let search = TemplateSearch {
        molecule,
        system_atoms,
        template,
        template_neighbors: &template_neighbors,
        search_order: connected_search_order(&template_neighbors),
    };

    let mut best: Option<(f64, Vec<usize>)> = None;
    let mut assignment = vec![usize::MAX; template_size];
    search.extend(0, &mut assignment, &mut |complete_assignment| {
        if !search.is_applicable(complete_assignment) {
            return;
        }
        let score = search.substituent_room(complete_assignment);
        if best
            .as_ref()
            .is_none_or(|(best_score, _)| score > *best_score + 1e-9)
        {
            best = Some((score, complete_assignment.to_vec()));
        }
    });
    best.map(|(_, assignment)| assignment)
}

fn template_adjacency(template: &CageTemplate) -> Vec<Vec<usize>> {
    let mut neighbors = vec![Vec::new(); template.coordinates.len()];
    for &(first, second) in template.bonds {
        neighbors[first].push(second);
        neighbors[second].push(first);
    }
    neighbors
}

/// Template atoms ordered so that each one after the first is bonded to an
/// earlier one, which lets the search extend matches along bonds.
fn connected_search_order(template_neighbors: &[Vec<usize>]) -> Vec<usize> {
    let mut order = vec![0];
    let mut visited = vec![false; template_neighbors.len()];
    visited[0] = true;
    let mut cursor = 0;
    while cursor < order.len() {
        for &neighbor in &template_neighbors[order[cursor]] {
            if !visited[neighbor] {
                visited[neighbor] = true;
                order.push(neighbor);
            }
        }
        cursor += 1;
    }
    order
}

struct TemplateSearch<'a> {
    molecule: &'a MoleculeGraph,
    system_atoms: &'a [usize],
    template: &'a CageTemplate,
    template_neighbors: &'a [Vec<usize>],
    search_order: Vec<usize>,
}

impl TemplateSearch<'_> {
    /// Assigns template atoms in search order, calling `on_match` for every
    /// complete assignment whose bonds agree exactly with the template.
    fn extend(&self, depth: usize, assignment: &mut [usize], on_match: &mut dyn FnMut(&[usize])) {
        if depth == self.search_order.len() {
            on_match(assignment);
            return;
        }
        let template_atom = self.search_order[depth];
        for candidate in self.candidates_for(template_atom, assignment) {
            if assignment.contains(&candidate)
                || !self.agrees_with_assigned(template_atom, candidate, assignment)
            {
                continue;
            }
            assignment[template_atom] = candidate;
            self.extend(depth + 1, assignment, on_match);
            assignment[template_atom] = usize::MAX;
        }
    }

    /// Molecule atoms that could take `template_atom`: system atoms with at
    /// least as many bonds inside the system as the template atom has. Every
    /// template atom after the first is bonded to an assigned one, so only
    /// the neighbors of that assigned atom need to be tried.
    fn candidates_for(&self, template_atom: usize, assignment: &[usize]) -> Vec<usize> {
        let required_bonds = self.template_neighbors[template_atom].len();
        let has_enough_bonds = |atom: &usize| self.system_bond_count(*atom) >= required_bonds;
        let assigned_neighbor = self.template_neighbors[template_atom]
            .iter()
            .map(|&neighbor| assignment[neighbor])
            .find(|&assigned| assigned != usize::MAX);
        match assigned_neighbor {
            Some(anchor) => self.molecule.adj[anchor]
                .iter()
                .map(|&(neighbor, _)| neighbor)
                .filter(|neighbor| self.system_atoms.binary_search(neighbor).is_ok())
                .filter(has_enough_bonds)
                .collect(),
            None => self
                .system_atoms
                .iter()
                .copied()
                .filter(has_enough_bonds)
                .collect(),
        }
    }

    fn system_bond_count(&self, atom: usize) -> usize {
        self.molecule.adj[atom]
            .iter()
            .filter(|(neighbor, _)| self.system_atoms.binary_search(neighbor).is_ok())
            .count()
    }

    /// True when `candidate` is bonded to exactly the molecule atoms assigned
    /// to the template neighbors of `template_atom`, and to no other assigned
    /// atom.
    fn agrees_with_assigned(
        &self,
        template_atom: usize,
        candidate: usize,
        assignment: &[usize],
    ) -> bool {
        assignment
            .iter()
            .enumerate()
            .filter(|&(_, &assigned)| assigned != usize::MAX)
            .all(|(other_template_atom, &assigned)| {
                let template_bonded =
                    self.template_neighbors[template_atom].contains(&other_template_atom);
                let molecule_bonded = self.molecule.bond_between(candidate, assigned).is_some();
                template_bonded == molecule_bonded
            })
    }

    fn is_applicable(&self, assignment: &[usize]) -> bool {
        match self.template.applicability {
            TemplateApplicability::Always => true,
            TemplateApplicability::SubstitutedAt(template_atom) => {
                self.outside_neighbor_count(assignment, template_atom) > 0
            }
            TemplateApplicability::UnsubstitutedAt(template_atoms) => template_atoms
                .iter()
                .all(|&template_atom| self.outside_neighbor_count(assignment, template_atom) == 0),
        }
    }

    fn outside_neighbor_count(&self, assignment: &[usize], template_atom: usize) -> usize {
        self.molecule.adj[assignment[template_atom]]
            .iter()
            .filter(|(neighbor, _)| !assignment.contains(neighbor))
            .count()
    }

    /// Sum over atoms with neighbors outside the template of their distance
    /// from the template's center, weighted by how many such neighbors they
    /// carry. Larger values leave more room for those neighbors.
    fn substituent_room(&self, assignment: &[usize]) -> f64 {
        let atom_count = self.template.coordinates.len() as f64;
        let (center_x, center_y) = self
            .template
            .coordinates
            .iter()
            .fold((0.0, 0.0), |(sum_x, sum_y), &(x, y)| (sum_x + x, sum_y + y));
        let center = Vec2::new(center_x / atom_count, center_y / atom_count);
        assignment
            .iter()
            .enumerate()
            .map(|(template_atom, _)| {
                let outside_neighbors = self.outside_neighbor_count(assignment, template_atom);
                let (x, y) = self.template.coordinates[template_atom];
                outside_neighbors as f64 * Vec2::new(x, y).distance_to(center)
            })
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rings::{find_rings, ring_system_atoms, ring_system_rings};

    fn template_positions_for(smiles: &str) -> Option<Vec<(usize, Vec2)>> {
        let preprocessed = crate::preprocess_smiles(smiles).unwrap();
        let molecule = MoleculeGraph::from_smiles(
            &preprocessed.smiles,
            preprocessed.forced_direction_markers,
            preprocessed.aromatic_atom_markers,
        )
        .unwrap();
        let rings = find_rings(&molecule);
        let system_atoms = ring_system_atoms(&rings, &ring_system_rings(&rings, 0));
        cage_template_positions(&molecule, &system_atoms)
    }

    #[test]
    fn cages_without_flat_drawings_use_templates() {
        for smiles in [
            "C1C2CC3CC1CC(C2)C3",
            "C12C3C4C1C5C2C3C45",
            "C1CN2CCN1CC2",
            "C1C2CC1C2",
            "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3",
        ] {
            assert!(template_positions_for(smiles).is_some(), "{smiles}");
        }
    }

    #[test]
    fn bicyclics_with_bare_one_atom_bridges_keep_the_flat_drawing() {
        for smiles in ["C1CC2CCC1C2", "C1C2CCC1C=C2", "C1CC2CCC(C1)C2"] {
            assert!(template_positions_for(smiles).is_none(), "{smiles}");
        }
    }

    #[test]
    fn substituted_one_atom_bridges_use_perspective_templates() {
        // Camphor's gem-dimethyl bridge and tropane's N-methyl bridge land on
        // the apex of the perspective view.
        for (smiles, bridge_atom) in [("CC1(C)C2CCC1(C)C(=O)C2", 1), ("CN1C2CCC1CC(O)C2", 1)] {
            let positions = template_positions_for(smiles).expect(smiles);
            let apex = positions
                .iter()
                .max_by(|first, second| first.1.y.partial_cmp(&second.1.y).unwrap())
                .unwrap();
            assert_eq!(apex.0, bridge_atom, "{smiles}");
        }
    }

    #[test]
    fn pinane_puts_the_substituted_bridge_on_the_outline() {
        // Atom 7 carries the gem-dimethyl; the hexagon's center is reserved
        // for the bare bridge.
        let positions = template_positions_for("CC1CCC2CC1C2(C)C").expect("pinane template");
        let center = positions
            .iter()
            .find(|(_, position)| position.length() < 1e-6)
            .unwrap();
        assert_ne!(center.0, 7);
    }
}
