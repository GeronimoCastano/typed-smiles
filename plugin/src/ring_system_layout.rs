//! Local coordinates for one ring system.
//!
//! Rings that share atoms are drawn together as one unit, which the molecule
//! layout then rotates and moves into place. A system starts from its most
//! characteristic piece: a cage template, a lattice macrocycle, or the start
//! ring as a regular polygon. Every other ring grows from the atoms already
//! drawn. Its undrawn atoms form arcs between drawn ones, so a ring sharing
//! one bond becomes a regular polygon on that bond, while a ring sharing a
//! longer path, or touching the drawing at separate atoms, becomes a bridge.
//! A ring sharing a single atom opens away from it as a spiro ring. When these
//! pieces cannot all keep their ideal shapes, relaxation spreads the strain.

use std::f64::consts::PI;

use crate::geometry::{
    largest_angular_gap, normalize_angle, point_to_segment_distance, segments_cross,
};
use crate::graph::MoleculeGraph;
use crate::layout_relaxation::{relax_ring_system, ring_system_needs_relaxation};
use crate::macrocycles::{macrocycle_ring_positions, ring_uses_lattice_outline};
use crate::render::Vec2;
use crate::ring_templates::cage_template_positions;
use crate::rings::{ring_system_atoms, ring_system_rings};

/// Smallest turn at each atom along a drawn arc. A bridge spanning a wide gap
/// would otherwise run straight and hide its atoms.
const MIN_ARC_TURN_PER_BOND: f64 = PI / 6.0;
/// Clearance a newly drawn ring atom needs from atoms already drawn.
const ARC_ATOM_CLEARANCE: f64 = 0.5;
/// Clearance between a newly drawn atom or bond and the existing drawing.
const ARC_BOND_CLEARANCE: f64 = 0.3;
/// Bonds drawn more than this far from unit length make a drawing defective.
const STRETCHED_BOND_TOLERANCE: f64 = 0.2;

/// Coordinates for the atoms of one ring system in the system's own frame.
pub(crate) struct RingSystemLayout {
    /// Atoms of the system, sorted.
    pub atoms: Vec<usize>,
    /// Positions indexed by molecule atom; only system atoms are meaningful.
    pub positions: Vec<Vec2>,
}

impl RingSystemLayout {
    /// Angle pointing away from the system at `atom`: the middle of the
    /// widest gap between the atom's bonds inside the system.
    pub(crate) fn outward_angle(&self, molecule: &MoleculeGraph, atom: usize) -> f64 {
        let bond_angles: Vec<f64> = molecule.adj[atom]
            .iter()
            .filter(|(neighbor, _)| self.atoms.binary_search(neighbor).is_ok())
            .map(|&(neighbor, _)| (self.positions[neighbor] - self.positions[atom]).angle())
            .collect();
        match largest_angular_gap(&bond_angles) {
            Some((gap_start, gap_width)) => normalize_angle(gap_start + gap_width / 2.0),
            None => 0.0,
        }
    }
}

/// Lays out the ring system containing `start_ring`. The system starts
/// from a matching template or macrocycle, otherwise from `start_ring` drawn
/// as a regular polygon centered on the origin with its first atom straight
/// up. When that drawing still has crowded atoms or crossing bonds, plain
/// polygons and arcs are tried from the other rings of the system, in ring
/// order, until one is free of defects; otherwise the drawing with the
/// fewest defects is kept.
pub(crate) fn layout_ring_system(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    start_ring: usize,
) -> RingSystemLayout {
    let system_rings = ring_system_rings(rings, start_ring);
    let atoms = ring_system_atoms(rings, &system_rings);
    let system = SystemToDraw {
        molecule,
        rings,
        system_rings: &system_rings,
        atoms: &atoms,
    };

    let mut best = system.draw(start_ring, PresetShapes::Allowed);
    if !best.defects.is_clean() {
        for &alternative_start in &system_rings {
            if best.first_piece == FirstPiece::StartRing && alternative_start == start_ring {
                continue;
            }
            let candidate = system.draw(alternative_start, PresetShapes::Ignored);
            if candidate.defects < best.defects {
                best = candidate;
                if best.defects.is_clean() {
                    break;
                }
            }
        }
    }
    RingSystemLayout {
        atoms,
        positions: best.coordinates,
    }
}

/// The piece a ring system drawing grows from.
#[derive(Clone, Copy, PartialEq, Eq)]
enum FirstPiece {
    Template,
    Macrocycle,
    StartRing,
}

/// Whether cage templates and lattice macrocycle outlines may start a
/// drawing.
#[derive(Clone, Copy, PartialEq, Eq)]
enum PresetShapes {
    Allowed,
    Ignored,
}

/// One complete drawing of a ring system and how readable it is.
struct SystemDrawingResult {
    coordinates: Vec<Vec2>,
    defects: DrawingDefects,
    first_piece: FirstPiece,
}

struct SystemToDraw<'a> {
    molecule: &'a MoleculeGraph,
    rings: &'a [Vec<usize>],
    system_rings: &'a [usize],
    atoms: &'a [usize],
}

impl SystemToDraw<'_> {
    fn draw(&self, start_ring: usize, preset_shapes: PresetShapes) -> SystemDrawingResult {
        let atom_count = self.molecule.n_atoms();
        let mut drawing = SystemDrawing {
            molecule: self.molecule,
            coordinates: vec![Vec2::default(); atom_count],
            placed: vec![false; atom_count],
            pinned: vec![false; atom_count],
        };

        let first_piece = drawing.draw_first_piece(
            self.rings,
            self.system_rings,
            self.atoms,
            start_ring,
            preset_shapes,
        );
        drawing.attach_remaining_rings(self.rings, self.system_rings);
        if ring_system_needs_relaxation(
            self.molecule,
            self.atoms,
            &drawing.pinned,
            &drawing.coordinates,
        ) {
            relax_ring_system(
                self.molecule,
                self.atoms,
                &drawing.pinned,
                &mut drawing.coordinates,
            );
        }

        SystemDrawingResult {
            defects: DrawingDefects::measure(
                self.molecule,
                self.atoms,
                &drawing.pinned,
                &drawing.coordinates,
            ),
            coordinates: drawing.coordinates,
            first_piece,
        }
    }
}

/// Readability problems of a ring system drawing. Fields are ordered from
/// most to least disruptive, so comparing two values prefers the drawing
/// with fewer severe problems. Geometry among pinned template atoms is
/// intentional and not counted.
#[derive(PartialEq, PartialOrd)]
struct DrawingDefects {
    crowded_atom_pairs: usize,
    atoms_on_bonds: usize,
    crossing_bonds: usize,
    stretched_bonds: usize,
    bond_length_strain: f64,
}

impl DrawingDefects {
    fn measure(
        molecule: &MoleculeGraph,
        atoms: &[usize],
        pinned: &[bool],
        coordinates: &[Vec2],
    ) -> Self {
        let bonds: Vec<(usize, usize)> = molecule
            .bonds
            .iter()
            .filter(|bond| {
                atoms.binary_search(&bond.from).is_ok() && atoms.binary_search(&bond.to).is_ok()
            })
            .map(|bond| (bond.from, bond.to))
            .collect();
        let all_pinned = |members: &[usize]| members.iter().all(|&atom| pinned[atom]);

        let mut crossing_bonds = 0;
        for (index, &(first_from, first_to)) in bonds.iter().enumerate() {
            for &(second_from, second_to) in &bonds[index + 1..] {
                let members = [first_from, first_to, second_from, second_to];
                let shares_atom = first_from == second_from
                    || first_from == second_to
                    || first_to == second_from
                    || first_to == second_to;
                if !shares_atom
                    && !all_pinned(&members)
                    && segments_cross(
                        coordinates[first_from],
                        coordinates[first_to],
                        coordinates[second_from],
                        coordinates[second_to],
                    )
                {
                    crossing_bonds += 1;
                }
            }
        }

        let mut crowded_atom_pairs = 0;
        for (index, &first) in atoms.iter().enumerate() {
            for &second in &atoms[index + 1..] {
                if molecule.bond_between(first, second).is_none()
                    && !all_pinned(&[first, second])
                    && coordinates[first].distance_to(coordinates[second]) < ARC_ATOM_CLEARANCE
                {
                    crowded_atom_pairs += 1;
                }
            }
        }

        let atoms_on_bonds = atoms
            .iter()
            .map(|&atom| {
                bonds
                    .iter()
                    .filter(|&&(from, to)| {
                        atom != from
                            && atom != to
                            && !all_pinned(&[atom, from, to])
                            && point_to_segment_distance(
                                coordinates[atom],
                                coordinates[from],
                                coordinates[to],
                            ) < ARC_BOND_CLEARANCE
                    })
                    .count()
            })
            .sum();

        let bond_length_errors: Vec<f64> = bonds
            .iter()
            .filter(|&&(from, to)| !all_pinned(&[from, to]))
            .map(|&(from, to)| (coordinates[from].distance_to(coordinates[to]) - 1.0).abs())
            .collect();
        let stretched_bonds = bond_length_errors
            .iter()
            .filter(|&&error| error > STRETCHED_BOND_TOLERANCE)
            .count();
        let bond_length_strain = bond_length_errors.iter().sum();

        Self {
            crowded_atom_pairs,
            atoms_on_bonds,
            crossing_bonds,
            stretched_bonds,
            bond_length_strain,
        }
    }

    fn is_clean(&self) -> bool {
        self.crowded_atom_pairs == 0
            && self.atoms_on_bonds == 0
            && self.crossing_bonds == 0
            && self.stretched_bonds == 0
    }
}

/// A ring system drawing in progress. Only atoms of the system are ever
/// placed, so every decision depends on the system alone.
struct SystemDrawing<'a> {
    molecule: &'a MoleculeGraph,
    coordinates: Vec<Vec2>,
    placed: Vec<bool>,
    /// Template atoms whose perspective drawing relaxation must keep.
    pinned: Vec<bool>,
}

impl SystemDrawing<'_> {
    /// Draws the piece the rest of the system grows from.
    fn draw_first_piece(
        &mut self,
        rings: &[Vec<usize>],
        system_rings: &[usize],
        system_atoms: &[usize],
        start_ring: usize,
        preset_shapes: PresetShapes,
    ) -> FirstPiece {
        if preset_shapes == PresetShapes::Ignored {
            self.draw_regular_ring(&rings[start_ring], Vec2::default(), PI / 2.0);
            return FirstPiece::StartRing;
        }
        if let Some(template_positions) = cage_template_positions(self.molecule, system_atoms) {
            for (atom, position) in template_positions {
                self.place(atom, position);
                self.pinned[atom] = true;
            }
            return FirstPiece::Template;
        }

        let largest_macrocycle = system_rings
            .iter()
            .copied()
            .filter(|&ring_index| ring_uses_lattice_outline(self.molecule, &rings[ring_index]))
            .max_by_key(|&ring_index| (rings[ring_index].len(), std::cmp::Reverse(ring_index)));
        if let Some(ring_index) = largest_macrocycle {
            let ring = &rings[ring_index];
            let neighboring_rings: Vec<&[usize]> = system_rings
                .iter()
                .filter(|&&other_index| other_index != ring_index)
                .map(|&other_index| rings[other_index].as_slice())
                .collect();
            let positions = macrocycle_ring_positions(self.molecule, ring, &neighboring_rings);
            for (ring_position, &atom) in ring.iter().enumerate() {
                self.place(atom, positions[ring_position]);
            }
            return FirstPiece::Macrocycle;
        }

        self.draw_regular_ring(&rings[start_ring], Vec2::default(), PI / 2.0);
        FirstPiece::StartRing
    }

    fn place(&mut self, atom: usize, position: Vec2) {
        self.coordinates[atom] = position;
        self.placed[atom] = true;
    }

    /// Draws a ring as a regular polygon around `center`, with its first atom
    /// at `start_angle`.
    fn draw_regular_ring(&mut self, ring: &[usize], center: Vec2, start_angle: f64) {
        let angle_step = 2.0 * PI / ring.len() as f64;
        let radius = regular_polygon_radius(ring.len());
        for (ring_position, &atom) in ring.iter().enumerate() {
            if !self.placed[atom] {
                let angle = start_angle + angle_step * ring_position as f64;
                self.place(atom, center + Vec2::from_angle(angle) * radius);
            }
        }
    }

    fn attach_remaining_rings(&mut self, rings: &[Vec<usize>], system_rings: &[usize]) {
        loop {
            let mut progressed = false;
            for &ring_index in system_rings {
                let ring = &rings[ring_index];
                if self.is_drawn(ring) {
                    continue;
                }
                let segments = self.drawn_segments(ring);
                if segments.len() == 1 && segments[0] >= 2 {
                    self.draw_ring_gaps(ring);
                    progressed = true;
                }
            }
            if progressed {
                continue;
            }

            let bridging_ring = system_rings
                .iter()
                .copied()
                .filter(|&ring_index| {
                    !self.is_drawn(&rings[ring_index]) && self.drawn_count(&rings[ring_index]) >= 2
                })
                .max_by_key(|&ring_index| {
                    (
                        self.drawn_count(&rings[ring_index]),
                        std::cmp::Reverse(ring_index),
                    )
                });
            if let Some(ring_index) = bridging_ring {
                self.draw_ring_gaps(&rings[ring_index]);
                continue;
            }

            let spiro_ring = system_rings
                .iter()
                .copied()
                .find(|&ring_index| self.drawn_count(&rings[ring_index]) == 1);
            if let Some(ring_index) = spiro_ring {
                self.draw_spiro_ring(&rings[ring_index]);
                continue;
            }
            break;
        }
    }

    fn is_drawn(&self, ring: &[usize]) -> bool {
        ring.iter().all(|&atom| self.placed[atom])
    }

    fn drawn_count(&self, ring: &[usize]) -> usize {
        ring.iter().filter(|&&atom| self.placed[atom]).count()
    }

    /// Sizes of the runs of consecutive drawn atoms around a partly drawn
    /// ring.
    fn drawn_segments(&self, ring: &[usize]) -> Vec<usize> {
        let ring_size = ring.len();
        let Some(first_undrawn) = (0..ring_size).find(|&position| !self.placed[ring[position]])
        else {
            return vec![ring_size];
        };
        let mut segments = Vec::new();
        let mut run_length = 0;
        for step in 1..=ring_size {
            if self.placed[ring[(first_undrawn + step) % ring_size]] {
                run_length += 1;
            } else if run_length > 0 {
                segments.push(run_length);
                run_length = 0;
            }
        }
        segments
    }

    /// Draws every run of undrawn atoms in a ring as an arc between the drawn
    /// atoms at its ends.
    fn draw_ring_gaps(&mut self, ring: &[usize]) {
        let ring_size = ring.len();
        let Some(walk_start) = (0..ring_size).find(|&position| {
            self.placed[ring[position]] && !self.placed[ring[(position + 1) % ring_size]]
        }) else {
            return;
        };

        let mut gaps: Vec<(usize, Vec<usize>, usize)> = Vec::new();
        let mut gap_start = ring[walk_start];
        let mut gap_atoms = Vec::new();
        for step in 1..=ring_size {
            let atom = ring[(walk_start + step) % ring_size];
            if !self.placed[atom] {
                gap_atoms.push(atom);
                continue;
            }
            if !gap_atoms.is_empty() {
                gaps.push((gap_start, std::mem::take(&mut gap_atoms), atom));
            }
            gap_start = atom;
        }

        for (start_atom, arc_atoms, end_atom) in gaps {
            self.draw_arc(start_atom, &arc_atoms, end_atom);
        }
    }

    /// Places `arc_atoms` on a circular arc of equal bonds from `start_atom`
    /// to `end_atom`, on whichever side of them collides least with the
    /// drawing.
    fn draw_arc(&mut self, start_atom: usize, arc_atoms: &[usize], end_atom: usize) {
        let start = self.coordinates[start_atom];
        let end = self.coordinates[end_atom];
        let chord = end - start;
        let bond_count = arc_atoms.len() + 1;
        if chord.length() < 1e-6 {
            self.draw_closed_loop(start_atom, arc_atoms);
            return;
        }

        let shape = ArcShape::spanning(chord.length(), bond_count);
        let normal = chord.perpendicular().normalized();
        let candidates =
            [normal, normal * -1.0].map(|side| shape.points(start, end, side, arc_atoms.len()));
        let chosen = if self.arc_preference(start_atom, end_atom, &candidates[1])
            < self.arc_preference(start_atom, end_atom, &candidates[0])
        {
            &candidates[1]
        } else {
            &candidates[0]
        };
        for (&atom, &position) in arc_atoms.iter().zip(chosen) {
            self.place(atom, position);
        }
    }

    /// Ordering key for a candidate arc; smaller is better. Collisions count
    /// first; among equally clear arcs, the one farther from the drawing wins.
    fn arc_preference(&self, start_atom: usize, end_atom: usize, points: &[Vec2]) -> (usize, i64) {
        let drawn_atoms: Vec<usize> = (0..self.placed.len())
            .filter(|&atom| self.placed[atom])
            .collect();
        let arc_centroid = points
            .iter()
            .fold(Vec2::default(), |sum, &point| sum + point)
            * (1.0 / points.len() as f64);
        let distance_from_drawing: f64 = drawn_atoms
            .iter()
            .map(|&atom| arc_centroid.distance_to(self.coordinates[atom]))
            .sum();
        (
            self.arc_collisions(start_atom, end_atom, points),
            -(distance_from_drawing * 1e6).round() as i64,
        )
    }

    fn arc_collisions(&self, start_atom: usize, end_atom: usize, points: &[Vec2]) -> usize {
        let drawn_bonds: Vec<(usize, usize)> = self
            .molecule
            .bonds
            .iter()
            .filter(|bond| self.placed[bond.from] && self.placed[bond.to])
            .map(|bond| (bond.from, bond.to))
            .collect();
        let drawn_atoms: Vec<usize> = (0..self.placed.len())
            .filter(|&atom| self.placed[atom])
            .collect();

        let mut collisions = 0;
        for &point in points {
            let crowds_atom = drawn_atoms
                .iter()
                .any(|&atom| point.distance_to(self.coordinates[atom]) < ARC_ATOM_CLEARANCE);
            let lies_on_bond = drawn_bonds.iter().any(|&(from, to)| {
                point_to_segment_distance(point, self.coordinates[from], self.coordinates[to])
                    < ARC_BOND_CLEARANCE
            });
            collisions += usize::from(crowds_atom) + usize::from(lies_on_bond);
        }

        // Each new bond, with the drawn atom it starts or ends on, if any.
        let mut arc_path: Vec<(Vec2, Option<usize>)> =
            vec![(self.coordinates[start_atom], Some(start_atom))];
        arc_path.extend(points.iter().map(|&point| (point, None)));
        arc_path.push((self.coordinates[end_atom], Some(end_atom)));
        for segment in arc_path.windows(2) {
            let (segment_start, start_owner) = segment[0];
            let (segment_end, end_owner) = segment[1];
            let touches_segment =
                |atom: usize| Some(atom) == start_owner || Some(atom) == end_owner;
            for &(from, to) in &drawn_bonds {
                if !touches_segment(from)
                    && !touches_segment(to)
                    && segments_cross(
                        segment_start,
                        segment_end,
                        self.coordinates[from],
                        self.coordinates[to],
                    )
                {
                    collisions += 1;
                }
            }
            for &atom in &drawn_atoms {
                if !touches_segment(atom)
                    && point_to_segment_distance(self.coordinates[atom], segment_start, segment_end)
                        < ARC_BOND_CLEARANCE
                {
                    collisions += 1;
                }
            }
        }
        collisions
    }

    /// Draws a gap whose two ends coincide as a regular polygon through that
    /// point, opening away from the rest of the drawing.
    fn draw_closed_loop(&mut self, shared_atom: usize, loop_atoms: &[usize]) {
        let mut ring = vec![shared_atom];
        ring.extend_from_slice(loop_atoms);
        self.draw_ring_opening_away(&ring, 0);
    }

    /// Draws a ring that shares exactly one atom with the drawing as a regular
    /// polygon opening away from that atom's drawn neighbors.
    fn draw_spiro_ring(&mut self, ring: &[usize]) {
        if let Some(shared_position) = ring.iter().position(|&atom| self.placed[atom]) {
            self.draw_ring_opening_away(ring, shared_position);
        }
    }

    fn draw_ring_opening_away(&mut self, ring: &[usize], shared_position: usize) {
        let shared_atom = ring[shared_position];
        let shared_point = self.coordinates[shared_atom];
        let toward_drawn_neighbors = self.molecule.adj[shared_atom]
            .iter()
            .filter(|&&(neighbor, _)| self.placed[neighbor] && neighbor != shared_atom)
            .map(|&(neighbor, _)| (self.coordinates[neighbor] - shared_point).normalized())
            .fold(Vec2::default(), |sum, direction| sum + direction);
        let inward = if toward_drawn_neighbors.length() < 1e-6 {
            Vec2::new(0.0, -1.0)
        } else {
            toward_drawn_neighbors.normalized()
        };

        let radius = regular_polygon_radius(ring.len());
        let center = shared_point - inward * radius;
        let shared_angle = (shared_point - center).angle();
        let angle_step = 2.0 * PI / ring.len() as f64;
        for (ring_position, &atom) in ring.iter().enumerate() {
            if !self.placed[atom] {
                let offset = ring_position as f64 - shared_position as f64;
                self.place(
                    atom,
                    center + Vec2::from_angle(shared_angle + angle_step * offset) * radius,
                );
            }
        }
    }
}

/// A circular arc of equal bonds with a given total turn.
struct ArcShape {
    total_turn: f64,
    bond_length: f64,
    bond_count: usize,
}

impl ArcShape {
    /// The arc of `bond_count` unit bonds whose ends lie `span` apart. A span
    /// too wide for the minimum turn stretches the bonds instead of
    /// straightening the arc.
    fn spanning(span: f64, bond_count: usize) -> Self {
        let unit_bond_span = |total_turn: f64| {
            (total_turn / 2.0).sin() / (total_turn / (2.0 * bond_count as f64)).sin()
        };
        let min_turn = MIN_ARC_TURN_PER_BOND * bond_count as f64;
        let widest_unit_span = unit_bond_span(min_turn);
        if span >= widest_unit_span {
            return Self {
                total_turn: min_turn,
                bond_length: span / widest_unit_span,
                bond_count,
            };
        }

        // The span of an arc of unit bonds shrinks as its turn grows.
        let mut tighter = 2.0 * PI - 1e-9;
        let mut looser = min_turn;
        for _ in 0..80 {
            let middle = (tighter + looser) / 2.0;
            if unit_bond_span(middle) > span {
                looser = middle;
            } else {
                tighter = middle;
            }
        }
        Self {
            total_turn: (tighter + looser) / 2.0,
            bond_length: 1.0,
            bond_count,
        }
    }

    /// Positions of the `atom_count` interior arc atoms from `start` to `end`,
    /// bulging toward `side`.
    fn points(&self, start: Vec2, end: Vec2, side: Vec2, atom_count: usize) -> Vec<Vec2> {
        let turn_per_bond = self.total_turn / self.bond_count as f64;
        let radius = self.bond_length / (2.0 * (turn_per_bond / 2.0).sin());
        let midpoint = (start + end) * 0.5;
        let center = midpoint - side * (radius * (self.total_turn / 2.0).cos());
        let start_offset = start - center;
        let end_offset = end - center;
        let sweep = if start_offset
            .rotated(self.total_turn)
            .distance_to(end_offset)
            < start_offset
                .rotated(-self.total_turn)
                .distance_to(end_offset)
        {
            1.0
        } else {
            -1.0
        };
        (1..=atom_count)
            .map(|step| center + start_offset.rotated(sweep * turn_per_bond * step as f64))
            .collect()
    }
}

pub(crate) fn regular_polygon_radius(vertex_count: usize) -> f64 {
    1.0 / (2.0 * (PI / vertex_count as f64).sin())
}

#[cfg(test)]
mod tests {
    use crate::layout_native;
    use crate::render::LayoutOutput;

    fn bond_lengths(layout: &LayoutOutput) -> Vec<f64> {
        layout
            .bonds
            .iter()
            .filter(|bond| !bond.virtual_bond)
            .map(|bond| {
                layout.atoms[bond.from]
                    .pos
                    .distance_to(layout.atoms[bond.to].pos)
            })
            .collect()
    }

    fn closest_atoms(layout: &LayoutOutput) -> f64 {
        let mut closest = f64::INFINITY;
        for first in 0..layout.atoms.len() {
            for second in first + 1..layout.atoms.len() {
                if !layout.atoms[first].virtual_h && !layout.atoms[second].virtual_h {
                    closest = closest.min(
                        layout.atoms[first]
                            .pos
                            .distance_to(layout.atoms[second].pos),
                    );
                }
            }
        }
        closest
    }

    #[test]
    fn norbornane_bridge_is_drawn_inside_a_unit_bond_outline() {
        let layout = layout_native("C1CC2CCC1C2").unwrap();
        for length in bond_lengths(&layout) {
            assert!((length - 1.0).abs() < 0.05, "bond length {length:.3}");
        }
        assert!(closest_atoms(&layout) > 0.55);
    }

    #[test]
    fn peri_fused_rings_close_without_stretched_bonds() {
        for smiles in [
            "c1cc2ccc3ccc4ccc5ccc1c1c2c3c4c51",
            "CN1C[C@@H](C=C2[C@H]1Cc1c[nH]c3cccc2c13)C(=O)O",
        ] {
            let layout = layout_native(smiles).unwrap();
            for length in bond_lengths(&layout) {
                assert!(
                    (length - 1.0).abs() < 0.1,
                    "{smiles}: bond length {length:.3}"
                );
            }
        }
    }

    #[test]
    fn helicene_ends_do_not_overlap() {
        let layout = layout_native("c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c21").unwrap();
        assert!(
            closest_atoms(&layout) > 0.6,
            "closest atoms {:.3}",
            closest_atoms(&layout)
        );
    }

    #[test]
    fn complex_ring_layouts_are_deterministic() {
        for smiles in [
            "O=C1C[C@H]2OCC=C3CN4CC[C@]56[C@@H]4C[C@H]3[C@H]2[C@H]6N1c1ccccc15",
            "c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c21",
            "C1CCCCCCCCCCCCCCCC1",
        ] {
            let first = layout_native(smiles).unwrap();
            let second = layout_native(smiles).unwrap();
            for (first_atom, second_atom) in first.atoms.iter().zip(&second.atoms) {
                assert_eq!(
                    first_atom.pos.x.to_bits(),
                    second_atom.pos.x.to_bits(),
                    "{smiles}"
                );
                assert_eq!(
                    first_atom.pos.y.to_bits(),
                    second_atom.pos.y.to_bits(),
                    "{smiles}"
                );
            }
        }
    }

    #[test]
    fn ring_system_reached_by_a_chain_keeps_its_shape() {
        // The bridged system is placed from its anchor atom and must arrive
        // with the same bond lengths as when it starts the layout.
        let layout = layout_native("CCCC1CC2CCC1C2").unwrap();
        for length in bond_lengths(&layout) {
            assert!((length - 1.0).abs() < 0.05, "bond length {length:.3}");
        }
    }
}
