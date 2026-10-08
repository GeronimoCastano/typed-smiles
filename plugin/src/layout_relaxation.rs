//! Geometric refinement for ring systems whose ideal shapes cannot all be
//! drawn at once.
//!
//! Ring systems are built from ideal pieces: regular polygons, lattice
//! macrocycles, cage templates, and bridges drawn as arcs. Strained
//! polycycles, helical overlaps, and crowded bridges leave some of those
//! pieces stretched or colliding. A small force field then spreads the strain
//! over the whole system: bonds pull toward unit length, each pair of bonds
//! at an atom keeps the angle it was built with, and atoms that are not
//! neighbors in the graph push apart. Pinned atoms, such as a cage template,
//! never move. The descent is deterministic, so a molecule always receives the
//! same drawing.

use crate::graph::MoleculeGraph;
use crate::render::Vec2;

/// Distance below which two atoms that share no bond or neighbor repel.
const CLASH_DISTANCE: f64 = 0.85;
const BOND_WEIGHT: f64 = 1.0;
const ANGLE_WEIGHT: f64 = 0.2;
const CLASH_WEIGHT: f64 = 1.5;
const MAX_ITERATIONS: usize = 3000;
/// Largest remaining force, in bond-length units, at which the drawing is
/// considered settled; smaller corrections are invisible.
const CONVERGED_FORCE: f64 = 1e-4;
/// Neighbors of one atom are held at least this far apart, which keeps
/// every bond angle at or above the 60° of an equilateral triangle.
const MIN_ANGLE_PAIR_DISTANCE: f64 = 1.0;
/// Largest distance one atom may move in a single descent step, which keeps
/// atoms from jumping across bonds.
const MAX_STEP_DISPLACEMENT: f64 = 0.1;

/// Tolerances that decide whether a constructed ring system needs refinement.
const BOND_LENGTH_TOLERANCE: f64 = 0.08;
const MIN_ATOM_SEPARATION: f64 = 0.5;
const MIN_ATOM_TO_BOND_SEPARATION: f64 = 0.3;

/// The atoms of one ring system together with the geometric relations the
/// refinement preserves.
struct RelaxationModel {
    atoms: Vec<usize>,
    movable: Vec<bool>,
    bonds: Vec<(usize, usize)>,
    /// Neighbor pairs of a shared atom with the distance that holds their
    /// constructed bond angle.
    angle_pairs: Vec<(usize, usize, f64)>,
    clash_pairs: Vec<(usize, usize)>,
}

/// True when a constructed ring system has stretched bonds, crowded atoms, or
/// atoms lying on bonds that refinement should remove. Geometry among pinned
/// atoms is intentional and never triggers refinement.
pub(crate) fn ring_system_needs_relaxation(
    molecule: &MoleculeGraph,
    system_atoms: &[usize],
    pinned: &[bool],
    coordinates: &[Vec2],
) -> bool {
    let model = RelaxationModel::new(molecule, system_atoms, pinned, coordinates);
    let position = |local_index: usize| coordinates[model.atoms[local_index]];
    let both_pinned = |first: usize, second: usize| !model.movable[first] && !model.movable[second];

    let has_distorted_bond = model.bonds.iter().any(|&(first, second)| {
        !both_pinned(first, second)
            && (position(first).distance_to(position(second)) - 1.0).abs() > BOND_LENGTH_TOLERANCE
    });
    let has_crowded_atoms = model.clash_pairs.iter().any(|&(first, second)| {
        !both_pinned(first, second)
            && position(first).distance_to(position(second)) < MIN_ATOM_SEPARATION
    });
    let has_pinched_angle = model.angle_pairs.iter().any(|&(first, second, _)| {
        !both_pinned(first, second)
            && position(first).distance_to(position(second)) < MIN_ATOM_SEPARATION
    });
    let has_atom_on_bond = (0..model.atoms.len()).any(|atom| {
        model.bonds.iter().any(|&(first, second)| {
            atom != first
                && atom != second
                && (model.movable[atom] || model.movable[first] || model.movable[second])
                && crate::geometry::point_to_segment_distance(
                    position(atom),
                    position(first),
                    position(second),
                ) < MIN_ATOM_TO_BOND_SEPARATION
        })
    });
    has_distorted_bond || has_crowded_atoms || has_pinched_angle || has_atom_on_bond
}

/// Refines the coordinates of one ring system in place.
pub(crate) fn relax_ring_system(
    molecule: &MoleculeGraph,
    system_atoms: &[usize],
    pinned: &[bool],
    coordinates: &mut [Vec2],
) {
    let model = RelaxationModel::new(molecule, system_atoms, pinned, coordinates);
    let mut positions: Vec<Vec2> = model.atoms.iter().map(|&atom| coordinates[atom]).collect();
    separate_coincident_atoms(&model, &mut positions);

    let (mut energy, mut gradient) = model.energy_and_gradient(&positions);
    let mut step_size = 0.05;
    for _ in 0..MAX_ITERATIONS {
        let largest_force = gradient
            .iter()
            .map(|force| force.length())
            .fold(0.0, f64::max);
        if largest_force < CONVERGED_FORCE {
            break;
        }

        let candidate = model.step_downhill(&positions, &gradient, step_size);
        let (candidate_energy, candidate_gradient) = model.energy_and_gradient(&candidate);
        if candidate_energy < energy {
            positions = candidate;
            energy = candidate_energy;
            gradient = candidate_gradient;
            step_size *= 1.2;
        } else {
            step_size *= 0.5;
            if step_size < 1e-10 {
                break;
            }
        }
    }

    for (local_index, &atom) in model.atoms.iter().enumerate() {
        coordinates[atom] = positions[local_index];
    }
}

/// Atoms built on the same spot give the springs no direction to push in.
/// Nudging the later atom outward from the system's center breaks the tie
/// the same way on every run.
fn separate_coincident_atoms(model: &RelaxationModel, positions: &mut [Vec2]) {
    let atom_count = positions.len();
    let center = positions
        .iter()
        .fold(Vec2::default(), |sum, &position| sum + position)
        * (1.0 / atom_count as f64);
    for first in 0..atom_count {
        for second in first + 1..atom_count {
            let bonded = model.bonds.contains(&(first, second));
            if bonded || positions[first].distance_to(positions[second]) > 1e-3 {
                continue;
            }
            let moved_atom = if model.movable[second] { second } else { first };
            if !model.movable[moved_atom] {
                continue;
            }
            let outward = positions[moved_atom] - center;
            let direction = if outward.length() > 1e-6 {
                outward.normalized()
            } else {
                Vec2::new(1.0, 0.0)
            };
            positions[moved_atom] = positions[moved_atom] + direction * 0.1;
        }
    }
}

impl RelaxationModel {
    fn new(
        molecule: &MoleculeGraph,
        system_atoms: &[usize],
        pinned: &[bool],
        coordinates: &[Vec2],
    ) -> Self {
        let atoms = system_atoms.to_vec();
        let mut local_index = vec![usize::MAX; molecule.n_atoms()];
        for (index, &atom) in atoms.iter().enumerate() {
            local_index[atom] = index;
        }
        let system_neighbors: Vec<Vec<usize>> = atoms
            .iter()
            .map(|&atom| {
                molecule.adj[atom]
                    .iter()
                    .map(|&(neighbor, _)| local_index[neighbor])
                    .filter(|&neighbor| neighbor != usize::MAX)
                    .collect()
            })
            .collect();

        let mut bonds = Vec::new();
        let mut angle_pairs = Vec::new();
        for (atom, neighbors) in system_neighbors.iter().enumerate() {
            for &neighbor in neighbors {
                if atom < neighbor {
                    bonds.push((atom, neighbor));
                }
            }
            for (pair_index, &first) in neighbors.iter().enumerate() {
                for &second in &neighbors[pair_index + 1..] {
                    let built_distance =
                        coordinates[atoms[first]].distance_to(coordinates[atoms[second]]);
                    angle_pairs.push((first, second, built_distance.max(MIN_ANGLE_PAIR_DISTANCE)));
                }
            }
        }

        let share_neighbor = |first: usize, second: usize| {
            system_neighbors[first]
                .iter()
                .any(|neighbor| system_neighbors[second].contains(neighbor))
        };
        let mut clash_pairs = Vec::new();
        for (first, first_neighbors) in system_neighbors.iter().enumerate() {
            for second in first + 1..atoms.len() {
                if !first_neighbors.contains(&second) && !share_neighbor(first, second) {
                    clash_pairs.push((first, second));
                }
            }
        }

        Self {
            movable: atoms.iter().map(|&atom| !pinned[atom]).collect(),
            atoms,
            bonds,
            angle_pairs,
            clash_pairs,
        }
    }

    /// Positions moved against the gradient by `step_size`, with each atom's
    /// move capped and pinned atoms left in place.
    fn step_downhill(&self, positions: &[Vec2], gradient: &[Vec2], step_size: f64) -> Vec<Vec2> {
        positions
            .iter()
            .zip(gradient)
            .zip(&self.movable)
            .map(|((&position, &force), &movable)| {
                if !movable {
                    return position;
                }
                let displacement = force * -step_size;
                let displacement_length = displacement.length();
                if displacement_length > MAX_STEP_DISPLACEMENT {
                    position + displacement * (MAX_STEP_DISPLACEMENT / displacement_length)
                } else {
                    position + displacement
                }
            })
            .collect()
    }

    /// Energy of the drawing and its gradient with respect to every atom
    /// position, computed in one pass over the springs.
    fn energy_and_gradient(&self, positions: &[Vec2]) -> (f64, Vec<Vec2>) {
        let mut energy = 0.0;
        let mut gradient = vec![Vec2::default(); positions.len()];
        let bond_springs = self
            .bonds
            .iter()
            .map(|&(first, second)| (first, second, 1.0, BOND_WEIGHT, false));
        let angle_springs = self
            .angle_pairs
            .iter()
            .map(|&(first, second, target)| (first, second, target, ANGLE_WEIGHT, false));
        let clash_springs = self
            .clash_pairs
            .iter()
            .map(|&(first, second)| (first, second, CLASH_DISTANCE, CLASH_WEIGHT, true));
        for (first, second, target, weight, repel_only) in
            bond_springs.chain(angle_springs).chain(clash_springs)
        {
            let offset = positions[first] - positions[second];
            let distance = offset.length();
            let deviation = distance - target;
            if repel_only && deviation > 0.0 {
                continue;
            }
            energy += weight * deviation * deviation;
            if distance > 1e-9 {
                let force = offset * (2.0 * weight * deviation / distance);
                gradient[first] = gradient[first] + force;
                gradient[second] = gradient[second] - force;
            }
        }
        (energy, gradient)
    }
}
