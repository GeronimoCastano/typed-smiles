//! Macrocycle outlines traced on a honeycomb lattice.
//!
//! A large ring drawn as a regular polygon becomes a wide circle whose nearly
//! straight angles hide its atoms. Chemists draw macrocycles with the 120°
//! zigzag of ordinary chains instead, which is exactly the outline of a
//! cluster of hexagons. This module traces such outlines, folds one corner
//! when the ring has an odd size, and chooses the rotation and direction that
//! put substituted atoms on outward corners and honor the ring's cis/trans
//! double bonds.

use std::collections::{HashMap, HashSet};

use crate::graph::{BondDirection, BondOrder, MoleculeGraph};
use crate::render::Vec2;

/// Rings of at least this size are drawn as lattice outlines.
pub(crate) const MACROCYCLE_MIN_ATOMS: usize = 12;

/// Smaller rings still use a lattice outline when one of their double bonds
/// is marked trans, which a regular polygon cannot show.
pub(crate) const TRANS_RING_MIN_ATOMS: usize = 9;

/// The smallest hexagon cluster with a ring-like outline is two fused
/// hexagons, whose outline has ten corners.
const SMALLEST_CLUSTER_OUTLINE: usize = 10;

/// Penalty for a cis/trans double bond drawn with the wrong geometry; large
/// enough to outweigh any number of crowded substituents.
const STEREO_VIOLATION_PENALTY: f64 = 1000.0;

/// Corner of the hexagon lattice in a doubled integer frame: x counts half
/// hexagon widths (√3/2 bond lengths) and y counts half bond lengths, so
/// every corner has exact integer coordinates.
type LatticeCorner = (i32, i32);

/// Hexagon of the lattice in axial coordinates.
type HexagonCell = (i32, i32);

const CELL_NEIGHBOR_OFFSETS: [(i32, i32); 6] = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)];

/// Whether a ring should be drawn as a lattice outline rather than a regular
/// polygon.
pub(crate) fn ring_uses_lattice_outline(molecule: &MoleculeGraph, ring: &[usize]) -> bool {
    ring.len() >= MACROCYCLE_MIN_ATOMS
        || (ring.len() >= TRANS_RING_MIN_ATOMS
            && ring_double_bond_requirements(molecule, ring)
                .iter()
                .any(|requirement| !requirement.same_side))
}

/// A ring double bond whose directional markers fix the relative side of
/// the two ring atoms flanking it.
#[derive(Clone, Copy)]
pub(crate) struct RingDoubleBondRequirement {
    /// Ring positions of the flanking atom, the two double-bonded atoms, and
    /// the other flanking atom, in ring order.
    pub flanking_before: usize,
    pub double_bond_start: usize,
    pub double_bond_end: usize,
    pub flanking_after: usize,
    /// True for cis flanking atoms, false for trans.
    pub same_side: bool,
}

/// Cis/trans requirements for every ring double bond whose both ends carry
/// directional markers.
pub(crate) fn ring_double_bond_requirements(
    molecule: &MoleculeGraph,
    ring: &[usize],
) -> Vec<RingDoubleBondRequirement> {
    let ring_size = ring.len();
    let mut requirements = Vec::new();
    for start in 0..ring_size {
        let end = (start + 1) % ring_size;
        let Some(bond_index) = molecule.bond_between(ring[start], ring[end]) else {
            continue;
        };
        if molecule.bonds[bond_index].order != BondOrder::Double {
            continue;
        }
        let before = (start + ring_size - 1) % ring_size;
        let after = (end + 1) % ring_size;
        let side_before =
            requested_side_of_ring_neighbor(molecule, ring[start], ring[end], ring[before]);
        let side_after =
            requested_side_of_ring_neighbor(molecule, ring[end], ring[start], ring[after]);
        if let (Some(side_before), Some(side_after)) = (side_before, side_after) {
            requirements.push(RingDoubleBondRequirement {
                flanking_before: before,
                double_bond_start: start,
                double_bond_end: end,
                flanking_after: after,
                same_side: side_before == side_after,
            });
        }
    }
    requirements
}

/// Side of the double-bond axis requested for `ring_neighbor` of `center` by
/// a directional marker on any single bond of `center`. A marker on an
/// exocyclic bond places the ring neighbor on the opposite side.
fn requested_side_of_ring_neighbor(
    molecule: &MoleculeGraph,
    center: usize,
    double_bond_partner: usize,
    ring_neighbor: usize,
) -> Option<i8> {
    molecule.adj[center]
        .iter()
        .find_map(|&(neighbor, bond_index)| {
            let bond = &molecule.bonds[bond_index];
            if neighbor == double_bond_partner || bond.direction == BondDirection::None {
                return None;
            }
            let marked_side = match bond.direction {
                BondDirection::Up => 1,
                BondDirection::Down => -1,
                BondDirection::None => 0,
            };
            let side = if bond.from == center {
                marked_side
            } else {
                -marked_side
            };
            Some(if neighbor == ring_neighbor {
                side
            } else {
                -side
            })
        })
}

/// True when the drawn flanking atoms of a ring double bond sit on the
/// requested sides of it.
pub(crate) fn requirement_is_drawn(
    requirement: &RingDoubleBondRequirement,
    position_of: impl Fn(usize) -> Vec2,
) -> bool {
    let start = position_of(requirement.double_bond_start);
    let axis = position_of(requirement.double_bond_end) - start;
    let side_before = axis.cross(position_of(requirement.flanking_before) - start);
    let side_after = axis.cross(position_of(requirement.flanking_after) - start);
    (side_before * side_after > 0.0) == requirement.same_side
}

/// Lattice-outline coordinates for a ring, indexed like `ring`.
/// `neighboring_rings` are the other rings of the same ring system.
pub(crate) fn macrocycle_ring_positions(
    molecule: &MoleculeGraph,
    ring: &[usize],
    neighboring_rings: &[&[usize]],
) -> Vec<Vec2> {
    let ring_size = ring.len();
    let crowding: Vec<usize> = ring
        .iter()
        .map(|&atom| {
            molecule.adj[atom]
                .iter()
                .filter(|(neighbor, _)| !ring.contains(neighbor))
                .count()
        })
        .collect();
    let inside_shared_path: Vec<bool> = (0..ring_size)
        .map(|ring_position| is_inside_shared_path(ring, ring_position, neighboring_rings))
        .collect();
    let requirements = ring_double_bond_requirements(molecule, ring);

    let mut best: Option<(f64, Vec<Vec2>)> = None;
    for outline in ring_outlines(ring_size) {
        for start in 0..ring_size {
            for direction in [1, ring_size - 1] {
                let positions: Vec<Vec2> = (0..ring_size)
                    .map(|ring_position| {
                        outline.corners[(start + direction * ring_position) % ring_size]
                    })
                    .collect();
                let corner_of =
                    |ring_position: usize| (start + direction * ring_position) % ring_size;
                let crowded_inward_corners: f64 = (0..ring_size)
                    .map(|ring_position| {
                        let corner = corner_of(ring_position);
                        let inward_penalty = if outline.points_outward[corner] {
                            0.0
                        } else {
                            1.0
                        };
                        let folded_penalty = if outline.folded_corner == Some(corner) {
                            0.5
                        } else {
                            0.0
                        };
                        let shared_path_penalty = if inside_shared_path[ring_position]
                            && outline.points_outward[corner]
                        {
                            1.0
                        } else {
                            0.0
                        };
                        crowding[ring_position] as f64 * (inward_penalty + folded_penalty)
                            + shared_path_penalty
                    })
                    .sum();
                let stereo_violations = requirements
                    .iter()
                    .filter(|requirement| {
                        !requirement_is_drawn(requirement, |ring_position| positions[ring_position])
                    })
                    .count() as f64;
                let penalty = crowded_inward_corners + STEREO_VIOLATION_PENALTY * stereo_violations;
                if best
                    .as_ref()
                    .is_none_or(|(best_penalty, _)| penalty < *best_penalty - 1e-9)
                {
                    best = Some((penalty, positions));
                }
            }
        }
    }
    best.map(|(_, positions)| positions)
        .unwrap_or_else(|| regular_polygon_positions(ring_size))
}

/// True when the atom and both of its ring neighbors also belong to one
/// other ring. That ring then shares a path through the atom and opens on
/// the far side, so the atom should point into the macrocycle.
fn is_inside_shared_path(
    ring: &[usize],
    ring_position: usize,
    neighboring_rings: &[&[usize]],
) -> bool {
    let ring_size = ring.len();
    let path = [
        ring[(ring_position + ring_size - 1) % ring_size],
        ring[ring_position],
        ring[(ring_position + 1) % ring_size],
    ];
    neighboring_rings
        .iter()
        .any(|other_ring| path.iter().all(|atom| other_ring.contains(atom)))
}

fn regular_polygon_positions(vertex_count: usize) -> Vec<Vec2> {
    let radius = 1.0 / (2.0 * (std::f64::consts::PI / vertex_count as f64).sin());
    (0..vertex_count)
        .map(|index| {
            Vec2::from_angle(2.0 * std::f64::consts::PI * index as f64 / vertex_count as f64)
                * radius
        })
        .collect()
}

/// A closed ring outline with its corners in counterclockwise order.
struct RingOutline {
    corners: Vec<Vec2>,
    /// Corners whose interior angle is 120°, leaving room outside for
    /// substituents.
    points_outward: Vec<bool>,
    /// The corner formed by folding two lattice corners into one, for odd
    /// ring sizes.
    folded_corner: Option<usize>,
}

/// Candidate outlines for a ring, roundest first.
fn ring_outlines(ring_size: usize) -> Vec<RingOutline> {
    let lattice_size = if ring_size.is_multiple_of(2) {
        ring_size.max(SMALLEST_CLUSTER_OUTLINE)
    } else {
        (ring_size + 1).max(SMALLEST_CLUSTER_OUTLINE)
    };
    let lattice_outlines: Vec<Vec<LatticeCorner>> = hexagon_clusters_with_outline(lattice_size)
        .iter()
        .filter_map(|cells| cluster_outline(cells))
        .collect();

    let mut outlines = Vec::new();
    for lattice_outline in lattice_outlines {
        let corners: Vec<Vec2> = lattice_outline
            .iter()
            .map(|&corner| lattice_point(corner))
            .collect();
        if corners.len() == ring_size {
            outlines.push(RingOutline {
                points_outward: corners_pointing_outward(&corners),
                corners,
                folded_corner: None,
            });
            continue;
        }
        if corners.len() != ring_size + 1 {
            continue;
        }
        let points_outward = corners_pointing_outward(&corners);
        for fold_start in 0..corners.len() {
            let fold_end = (fold_start + 1) % corners.len();
            if points_outward[fold_start] && !points_outward[fold_end] {
                outlines.push(fold_adjacent_corners(&corners, fold_start));
            }
        }
    }
    outlines
}

/// Replaces two neighboring corners by one corner pushed outward, turning a
/// lattice outline into one with a single odd-sized bend.
fn fold_adjacent_corners(corners: &[Vec2], fold_start: usize) -> RingOutline {
    let corner_count = corners.len();
    let fold_end = (fold_start + 1) % corner_count;
    let midpoint = (corners[fold_start] + corners[fold_end]) * 0.5;
    let before = corners[(fold_start + corner_count - 1) % corner_count];
    let after = corners[(fold_end + 1) % corner_count];
    // Outward lies to the right of the counterclockwise traversal.
    let outward_normal = (after - before).perpendicular().normalized() * -1.0;
    let folded_corner_position = midpoint + outward_normal * 0.45;

    let mut folded_corners = Vec::with_capacity(corner_count - 1);
    let mut folded_corner = 0;
    for (index, &corner) in corners.iter().enumerate() {
        if index == fold_end {
            continue;
        }
        if index == fold_start {
            folded_corner = folded_corners.len();
            folded_corners.push(folded_corner_position);
        } else {
            folded_corners.push(corner);
        }
    }
    RingOutline {
        points_outward: corners_pointing_outward(&folded_corners),
        corners: folded_corners,
        folded_corner: Some(folded_corner),
    }
}

fn corners_pointing_outward(corners: &[Vec2]) -> Vec<bool> {
    let corner_count = corners.len();
    (0..corner_count)
        .map(|index| {
            let previous = corners[(index + corner_count - 1) % corner_count];
            let next = corners[(index + 1) % corner_count];
            (corners[index] - previous).cross(next - corners[index]) > 0.0
        })
        .collect()
}

fn lattice_point(corner: LatticeCorner) -> Vec2 {
    Vec2::new(
        corner.0 as f64 * 3.0_f64.sqrt() / 2.0,
        corner.1 as f64 * 0.5,
    )
}

fn cell_corners(cell: HexagonCell) -> [LatticeCorner; 6] {
    let (column, row) = cell;
    let center = (2 * column + row, 3 * row);
    [(1, 1), (0, 2), (-1, 1), (-1, -1), (0, -2), (1, -1)]
        .map(|(offset_x, offset_y)| (center.0 + offset_x, center.1 + offset_y))
}

/// Outline of a hexagon cluster traced counterclockwise, or `None` when the
/// cluster encloses a hole or touches itself at a single corner.
fn cluster_outline(cells: &[HexagonCell]) -> Option<Vec<LatticeCorner>> {
    let mut cell_edges = HashSet::new();
    for &cell in cells {
        let corners = cell_corners(cell);
        for index in 0..6 {
            cell_edges.insert((corners[index], corners[(index + 1) % 6]));
        }
    }
    let mut next_corner: HashMap<LatticeCorner, LatticeCorner> = HashMap::new();
    for &(from, to) in &cell_edges {
        if !cell_edges.contains(&(to, from)) && next_corner.insert(from, to).is_some() {
            return None;
        }
    }

    let start = *next_corner.keys().min()?;
    let mut outline = vec![start];
    let mut current = next_corner[&start];
    while current != start {
        outline.push(current);
        current = *next_corner.get(&current)?;
        if outline.len() > next_corner.len() {
            return None;
        }
    }
    (outline.len() == next_corner.len()).then_some(outline)
}

/// Number of outline corners of a cluster without holes or pinches. Each
/// hexagon contributes six edges, and each pair of touching hexagons hides
/// two of them. Clusters whose count matches are traced afterwards to
/// confirm the outline is one simple loop.
fn outline_edge_count(cells: &[HexagonCell]) -> usize {
    let occupied: HashSet<HexagonCell> = cells.iter().copied().collect();
    let touching_pairs = cells
        .iter()
        .map(|&(column, row)| {
            CELL_NEIGHBOR_OFFSETS
                .iter()
                .filter(|&&(offset_column, offset_row)| {
                    occupied.contains(&(column + offset_column, row + offset_row))
                })
                .count()
        })
        .sum::<usize>()
        / 2;
    6 * cells.len() - 2 * touching_pairs
}

fn has_outline_of(cells: &[HexagonCell], corner_count: usize) -> bool {
    !cells.is_empty()
        && outline_edge_count(cells) == corner_count
        && cluster_outline(cells).is_some_and(|outline| outline.len() == corner_count)
}

/// Hexagon clusters whose outline has `corner_count` corners, roundest
/// first. Rounder clusters leave more room inside the ring and give
/// substituents more space outside it.
fn hexagon_clusters_with_outline(corner_count: usize) -> Vec<Vec<HexagonCell>> {
    let mut clusters: Vec<Vec<HexagonCell>> = Vec::new();
    if let Some(grown) = grown_cluster(corner_count) {
        clusters.push(grown);
    }
    for rows in 1..=corner_count / 4 {
        for row_length in 1..=corner_count / 2 {
            // A block of `rows` rows of `row_length` hexagons has an outline of
            // about 4 × (rows + row_length) − 2 corners; others cannot match.
            let estimated_corners = 4 * (rows + row_length) as i64 - 2;
            if (estimated_corners - corner_count as i64).abs() > 4 {
                continue;
            }
            for length_change in [0, 1, -1] {
                let cells = staggered_block(rows, row_length, length_change);
                if has_outline_of(&cells, corner_count) {
                    clusters.push(cells);
                }
            }
            let parallelogram: Vec<HexagonCell> = (0..rows as i32)
                .flat_map(|row| (0..row_length as i32).map(move |column| (column, row)))
                .collect();
            if has_outline_of(&parallelogram, corner_count) {
                clusters.push(parallelogram);
            }
        }
    }
    clusters.sort_by_key(|cluster| std::cmp::Reverse(cluster.len()));
    clusters.dedup();
    clusters
}

/// Rows of hexagons stacked in a brick pattern. Row lengths alternate by
/// `length_change` so the block can have straight or ragged sides.
fn staggered_block(rows: usize, row_length: usize, length_change: i32) -> Vec<HexagonCell> {
    let mut cells = Vec::new();
    for row in 0..rows as i32 {
        let length = if row % 2 == 1 {
            row_length as i32 + length_change
        } else {
            row_length as i32
        };
        if length <= 0 {
            return Vec::new();
        }
        let first_column = -(row / 2)
            - if row % 2 == 1 && length_change > 0 {
                1
            } else {
                0
            };
        cells.extend((0..length).map(|offset| (first_column + offset, row)));
    }
    cells
}

/// A compact cluster grown from two hexagons. Each step adds the hexagon
/// nearest the start that lengthens the outline by two corners, then fills
/// any notch that would not change the outline.
///
/// A hexagon touching the cluster along one unbroken run of `k` edges
/// changes the outline by `6 - 2k` corners and cannot enclose a hole, so
/// growth only inspects the run of occupied neighbors around each candidate.
fn grown_cluster(corner_count: usize) -> Option<Vec<HexagonCell>> {
    let mut cells: Vec<HexagonCell> = vec![(0, 0), (1, 0)];
    let mut occupied: HashSet<HexagonCell> = cells.iter().copied().collect();
    let mut current_length = SMALLEST_CLUSTER_OUTLINE;
    while current_length < corner_count {
        let mut candidates = free_neighbor_cells(&cells, &occupied);
        candidates.sort_by(|first, second| {
            cell_distance_from_origin(*first)
                .partial_cmp(&cell_distance_from_origin(*second))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(first.cmp(second))
        });
        let next_cell = candidates
            .into_iter()
            .find(|&candidate| contiguous_occupied_neighbors(candidate, &occupied) == Some(2))?;
        cells.push(next_cell);
        occupied.insert(next_cell);
        current_length += 2;
        fill_outline_notches(&mut cells, &mut occupied);
    }
    (current_length == corner_count).then_some(cells)
}

fn fill_outline_notches(cells: &mut Vec<HexagonCell>, occupied: &mut HashSet<HexagonCell>) {
    while let Some(notch) = free_neighbor_cells(cells, occupied)
        .into_iter()
        .find(|&candidate| contiguous_occupied_neighbors(candidate, occupied) == Some(3))
    {
        cells.push(notch);
        occupied.insert(notch);
    }
}

/// Unoccupied cells next to the cluster, in a fixed order.
fn free_neighbor_cells(cells: &[HexagonCell], occupied: &HashSet<HexagonCell>) -> Vec<HexagonCell> {
    let mut free: Vec<HexagonCell> = cells
        .iter()
        .flat_map(|&(column, row)| {
            CELL_NEIGHBOR_OFFSETS
                .iter()
                .map(move |&(offset_column, offset_row)| (column + offset_column, row + offset_row))
        })
        .filter(|cell| !occupied.contains(cell))
        .collect();
    free.sort_unstable();
    free.dedup();
    free
}

/// Number of occupied neighbors of `cell` when they form one unbroken run
/// around it, or `None` when they are split into separate runs.
fn contiguous_occupied_neighbors(
    cell: HexagonCell,
    occupied: &HashSet<HexagonCell>,
) -> Option<usize> {
    let (column, row) = cell;
    let neighbor_is_occupied: Vec<bool> = CELL_NEIGHBOR_OFFSETS
        .iter()
        .map(|&(offset_column, offset_row)| {
            occupied.contains(&(column + offset_column, row + offset_row))
        })
        .collect();
    let occupied_count = neighbor_is_occupied
        .iter()
        .filter(|&&is_occupied| is_occupied)
        .count();
    let run_starts = (0..6)
        .filter(|&index| neighbor_is_occupied[index] && !neighbor_is_occupied[(index + 5) % 6])
        .count();
    (occupied_count == 6 || run_starts == 1).then_some(occupied_count)
}

fn cell_distance_from_origin(cell: HexagonCell) -> f64 {
    let (column, row) = cell;
    let center = lattice_point((2 * column + row, 3 * row));
    let origin_pair_center = lattice_point((1, 0));
    center.distance_to(origin_pair_center)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layout_native;
    use crate::render::LayoutOutput;

    fn side_of(layout: &LayoutOutput, axis_start: usize, axis_end: usize, atom: usize) -> f64 {
        let start = layout.atoms[axis_start].pos;
        (layout.atoms[axis_end].pos - start).cross(layout.atoms[atom].pos - start)
    }

    #[test]
    fn large_ring_follows_the_hexagon_lattice() {
        let layout = layout_native("C1CCCCCCCCCCCCCCC1").unwrap();
        for center in 0..16 {
            let previous = layout.atoms[(center + 15) % 16].pos;
            let next = layout.atoms[(center + 1) % 16].pos;
            let flanking_distance = previous.distance_to(next);
            assert!(
                (flanking_distance - 3.0_f64.sqrt()).abs() < 1e-6,
                "angle at atom {center} is not 120° (flanking distance {flanking_distance:.3})"
            );
        }
    }

    #[test]
    fn marked_trans_ring_double_bond_is_drawn_trans() {
        // Atoms 5=6 with ring neighbors 4 and 7.
        let layout = layout_native("C1CCCC/C=C/CCCCC1").unwrap();
        assert!(side_of(&layout, 5, 6, 4) * side_of(&layout, 5, 6, 7) < 0.0);
    }

    #[test]
    fn marked_cis_ring_double_bond_is_drawn_cis() {
        let layout = layout_native("C1CCCC/C=C\\CCCCC1").unwrap();
        assert!(side_of(&layout, 5, 6, 4) * side_of(&layout, 5, 6, 7) > 0.0);
    }

    #[test]
    fn marker_on_ring_closure_digit_reads_from_the_closing_atom() {
        // Trans: the 8=7 double bond with ring neighbors 6 and 0.
        let layout = layout_native("C1CCCCCC/C=C/1").unwrap();
        assert!(side_of(&layout, 7, 8, 6) * side_of(&layout, 7, 8, 0) < 0.0);
        let layout = layout_native("C1CCCCCC/C=C\\1").unwrap();
        assert!(side_of(&layout, 7, 8, 6) * side_of(&layout, 7, 8, 0) > 0.0);
    }

    #[test]
    fn exocyclic_marker_sets_the_ring_geometry() {
        // The methyl is trans to atom 3, so ring atoms 11 and 3 are cis.
        let layout = layout_native("C/C1=C/CCCCCCCCCC1").unwrap();
        assert!(side_of(&layout, 1, 2, 12) * side_of(&layout, 1, 2, 3) > 0.0);
    }

    #[test]
    fn trans_double_bond_in_small_ring_is_reported_undepicted() {
        let layout = layout_native("C1CCC/C=C/CC1").unwrap();
        assert_eq!(layout.undepicted_stereo.len(), 1);
        let reason = &layout.undepicted_stereo[0].reason;
        assert!(reason.contains("atoms 4 and 5"), "{reason}");
        assert!(reason.contains("ring of 8 atoms"), "{reason}");
        assert!(reason.contains("at least 9 atoms"), "{reason}");
    }

    #[test]
    fn cis_double_bond_in_small_ring_has_no_size_advice() {
        let layout = layout_native("C1CCC/C=C\\CC1").unwrap();
        assert!(layout.undepicted_stereo.is_empty());
    }

    #[test]
    fn substituents_sit_on_outward_corners() {
        // 2-methylcyclopentadecanone: the methyl and carbonyl point out.
        let layout = layout_native("CC1CCCCCCCCCCCCCC1=O").unwrap();
        let ring: Vec<usize> = (1..=15).collect();
        let center = ring
            .iter()
            .fold(Vec2::default(), |sum, &atom| sum + layout.atoms[atom].pos)
            * (1.0 / ring.len() as f64);
        for (substituent, ring_atom) in [(0, 1), (16, 15)] {
            let outward = layout.atoms[ring_atom].pos - center;
            let bond = layout.atoms[substituent].pos - layout.atoms[ring_atom].pos;
            assert!(
                outward.normalized().cross(bond.normalized()).abs() < 0.9
                    && (outward.x * bond.x + outward.y * bond.y) > 0.0,
                "substituent {substituent} points into the ring"
            );
        }
    }

    #[test]
    fn odd_macrocycle_keeps_bonds_near_unit_length() {
        let layout = layout_native("C1CCCCCCCCCCCCCCCC1").unwrap();
        for bond in &layout.bonds {
            let length = layout.atoms[bond.from]
                .pos
                .distance_to(layout.atoms[bond.to].pos);
            assert!(
                (length - 1.0).abs() < 0.2,
                "bond {}-{} has length {length:.2}",
                bond.from,
                bond.to
            );
        }
    }

    #[test]
    fn lattice_outlines_exist_for_every_macrocycle_size() {
        for ring_size in MACROCYCLE_MIN_ATOMS..=40 {
            assert!(
                ring_outlines(ring_size)
                    .iter()
                    .any(|outline| outline.corners.len() == ring_size),
                "no outline for a {ring_size}-membered ring"
            );
        }
    }
}
