/// 2D coordinate generation for molecular graphs.
///
/// Algorithm (MVP):
///   1. Detect all rings using DFS cycle detection.
///   2. Identify ring systems (connected sets of rings sharing bonds/atoms).
///   3. Place each ring as a regular polygon. For fused rings, align the
///      shared edge first, then compute the polygon vertices.
///   4. Lay out acyclic chains via DFS, extending each bond at ±120° from
///      the incoming direction (standard organic chemistry depiction angle).
///   5. Translate the whole molecule so the centroid is at the origin.
use std::collections::HashSet;
use std::f64::consts::PI;

use crate::geometry::{angular_gaps, largest_angular_gap, normalize_angle, point_in_polygon};
use crate::graph::{AtomChirality, Bond, BondDirection, BondOrder, BondStereo, MoleculeGraph};
use crate::render::{
    AromaticRing, AtomOutput, BondOutput, LayoutOutput, UndepictedStereoOutput, Vec2,
};
use crate::ring_system_layout::layout_ring_system;
use crate::rings::{find_rings, ring_bond_set, ring_has_edge, rings_share_edge};
use crate::stereo::{depict_stereo, StereoDepiction};

/// Gap between the bounding boxes of dot-separated fragments, in bond lengths.
const FRAGMENT_GAP: f64 = 1.5;

pub fn compute_layout(molecule: &MoleculeGraph) -> Result<LayoutOutput, String> {
    if molecule.n_atoms() == 0 {
        return Ok(empty_layout_output());
    }

    let coordinates = layout_coordinates(molecule)?;
    let rings = find_rings(molecule);
    let ring_bonds = ring_bond_set(molecule, &rings);
    let stereo = depict_stereo(molecule, &coordinates, &rings, &ring_bonds);

    let mut atoms = build_atom_outputs(molecule, &coordinates, &stereo.hydrogen_stereo);
    let inner_directions = ring_inner_directions(molecule, &rings, &coordinates);
    let mut bonds = build_bond_outputs(molecule, &stereo, &inner_directions);
    append_virtual_hydrogen_outputs(
        molecule,
        &coordinates,
        &stereo.hydrogen_stereo,
        &mut atoms,
        &mut bonds,
    );

    let all_positions: Vec<Vec2> = atoms.iter().map(|atom| atom.pos).collect();
    let (bbox_width, bbox_height) = bounding_box(&all_positions);

    Ok(LayoutOutput {
        atoms,
        bonds,
        aromatic_rings: aromatic_ring_circles(molecule, &rings, &coordinates),
        undepicted_stereo: stereo
            .undepicted
            .into_iter()
            .map(|undepicted| UndepictedStereoOutput {
                atom: undepicted.atom,
                reason: undepicted.reason,
            })
            .collect(),
        bbox_width,
        bbox_height,
    })
}

fn empty_layout_output() -> LayoutOutput {
    LayoutOutput {
        atoms: Vec::new(),
        bonds: Vec::new(),
        aromatic_rings: Vec::new(),
        undepicted_stereo: Vec::new(),
        bbox_width: 0.0,
        bbox_height: 0.0,
    }
}

fn build_atom_outputs(
    molecule: &MoleculeGraph,
    coordinates: &[Vec2],
    hydrogen_stereo: &[Option<(BondStereo, Vec2)>],
) -> Vec<AtomOutput> {
    molecule
        .atoms
        .iter()
        .enumerate()
        .map(|(atom_index, atom)| {
            let lone_pairs = lone_pair_count(molecule, atom_index);
            AtomOutput {
                symbol: atom.symbol.clone(),
                pos: coordinates[atom_index],
                hcount: atom.hcount,
                implicit_h: implicit_h_count(molecule, atom_index),
                charge: atom.charge,
                isotope: atom.isotope.unwrap_or(0),
                lone_pairs,
                lone_pair_dirs: lone_pair_directions(
                    molecule,
                    atom_index,
                    coordinates,
                    lone_pairs as usize,
                ),
                abbrev: atom.abbrev.clone(),
                abbrev_style: atom.abbrev_style.clone(),
                abbrev_anchor: atom.abbrev_anchor,
                abbrev_anchor_len: atom.abbrev_anchor_len,
                abbrev_offset_x: atom.abbrev_offset_x,
                abbrev_offset_y: atom.abbrev_offset_y,
                atom_map: atom.atom_map,
                chirality: atom.chirality.as_str().to_string(),
                stereo_h: hydrogen_stereo[atom_index]
                    .map(|(stereo, _)| stereo.as_str().to_string())
                    .unwrap_or_else(|| "none".to_string()),
                stereo_h_dir: hydrogen_stereo[atom_index]
                    .map(|(_, direction)| direction)
                    .unwrap_or_default(),
                virtual_h: false,
            }
        })
        .collect()
}

fn build_bond_outputs(
    molecule: &MoleculeGraph,
    stereo: &StereoDepiction,
    inner_directions: &[(f64, f64)],
) -> Vec<BondOutput> {
    molecule
        .bonds
        .iter()
        .enumerate()
        .map(|(bond_index, bond)| BondOutput {
            from: bond.from,
            to: bond.to,
            order: bond.order.as_u8(),
            stereo: stereo.bond_stereo[bond_index].as_str().to_string(),
            stereo_tip: stereo.wedge_tips[bond_index],
            forced_stereo: bond.forced_stereo,
            direction: direction_as_str(bond.direction).to_string(),
            inner_x: inner_directions[bond_index].0,
            inner_y: inner_directions[bond_index].1,
            virtual_bond: false,
            aromatic: bond.aromatic,
        })
        .collect()
}

fn append_virtual_hydrogen_outputs(
    molecule: &MoleculeGraph,
    coordinates: &[Vec2],
    hydrogen_stereo: &[Option<(BondStereo, Vec2)>],
    atoms: &mut Vec<AtomOutput>,
    bonds: &mut Vec<BondOutput>,
) {
    for atom_index in 0..molecule.n_atoms() {
        let parent_position = coordinates[atom_index];
        let occupied_angles: Vec<f64> = molecule.adj[atom_index]
            .iter()
            .map(|&(neighbor, _)| {
                (coordinates[neighbor].y - parent_position.y)
                    .atan2(coordinates[neighbor].x - parent_position.x)
            })
            .collect();
        if should_add_virtual_hydrogen(molecule, atoms, hydrogen_stereo, atom_index) {
            let angle = hydrogen_label_angle(&occupied_angles);
            append_virtual_hydrogen(atoms, bonds, atom_index, angle);
        }
    }
}

fn should_add_virtual_hydrogen(
    molecule: &MoleculeGraph,
    atoms: &[AtomOutput],
    hydrogen_stereo: &[Option<(BondStereo, Vec2)>],
    atom_index: usize,
) -> bool {
    let atom = &molecule.atoms[atom_index];
    if atom.has_explicit_h {
        return atom.hcount > 0
            && hydrogen_stereo[atom_index].is_none()
            && !atom.chirality.is_tetrahedral();
    }

    let symbol = atoms[atom_index].symbol.as_str();
    let is_heteroatom = !matches!(symbol, "C" | "c" | "H" | "*");
    is_heteroatom && atoms[atom_index].implicit_h > 0
}

fn append_virtual_hydrogen(
    atoms: &mut Vec<AtomOutput>,
    bonds: &mut Vec<BondOutput>,
    parent_atom: usize,
    angle: f64,
) {
    let direction = Vec2::new(angle.cos(), angle.sin());
    let parent_position = atoms[parent_atom].pos;
    let hydrogen_index = atoms.len();
    atoms.push(AtomOutput {
        symbol: "H".to_string(),
        pos: Vec2::new(
            parent_position.x + direction.x * 0.35,
            parent_position.y + direction.y * 0.35,
        ),
        hcount: 0,
        implicit_h: 0,
        charge: 0,
        isotope: 0,
        lone_pairs: 0,
        lone_pair_dirs: Vec::new(),
        abbrev: String::new(),
        abbrev_style: String::new(),
        abbrev_anchor: 0,
        abbrev_anchor_len: 0,
        abbrev_offset_x: 0.0,
        abbrev_offset_y: 0.0,
        atom_map: 0,
        chirality: "none".to_string(),
        stereo_h: "none".to_string(),
        stereo_h_dir: Vec2::default(),
        virtual_h: true,
    });
    bonds.push(BondOutput {
        from: parent_atom,
        to: hydrogen_index,
        order: 1,
        stereo: "none".to_string(),
        stereo_tip: None,
        forced_stereo: false,
        direction: "none".to_string(),
        inner_x: 0.0,
        inner_y: 0.0,
        virtual_bond: true,
        aromatic: false,
    });
}

/// Circle parameters for every ring whose bonds were all aromatic in the
/// input, enabling the inscribed-circle depiction. The radius scales with the
/// ring's inradius so the circle stays clear of the bonds for any ring size.
fn aromatic_ring_circles(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    coordinates: &[Vec2],
) -> Vec<AromaticRing> {
    let mut circles = Vec::new();
    for ring in rings {
        let all_aromatic = (0..ring.len()).all(|ring_position| {
            let first_atom = ring[ring_position];
            let second_atom = ring[(ring_position + 1) % ring.len()];
            molecule
                .bond_between(first_atom, second_atom)
                .map(|index| molecule.bonds[index].aromatic)
                .unwrap_or(false)
        });
        if !all_aromatic {
            continue;
        }
        let ring_size = ring.len() as f64;
        let center_x = ring
            .iter()
            .map(|&atom_index| coordinates[atom_index].x)
            .sum::<f64>()
            / ring_size;
        let center_y = ring
            .iter()
            .map(|&atom_index| coordinates[atom_index].y)
            .sum::<f64>()
            / ring_size;
        let center = Vec2::new(center_x, center_y);
        let inradius = (0..ring.len())
            .map(|ring_position| {
                let first_position = coordinates[ring[ring_position]];
                let second_position = coordinates[ring[(ring_position + 1) % ring.len()]];
                let bond_midpoint = Vec2::new(
                    (first_position.x + second_position.x) / 2.0,
                    (first_position.y + second_position.y) / 2.0,
                );
                center.distance_to(bond_midpoint)
            })
            .fold(f64::INFINITY, f64::min);
        circles.push(AromaticRing {
            center,
            radius: inradius * 0.72,
        });
    }
    circles
}

/// Coordinates for every atom. Each connected fragment is laid out on its own;
/// dot-separated fragments are then arranged left to right in writing order,
/// vertically centered, with a fixed gap between their bounding boxes.
fn layout_coordinates(molecule: &MoleculeGraph) -> Result<Vec<Vec2>, String> {
    let components = connected_components(molecule);
    if components.len() <= 1 {
        return place_connected_molecule(molecule);
    }

    let mut coordinates = vec![Vec2::new(0.0, 0.0); molecule.n_atoms()];
    let mut cursor = 0.0;
    for (component_index, component) in components.iter().enumerate() {
        let component_molecule = component_subgraph(molecule, component);
        let component_coordinates = place_connected_molecule(&component_molecule)?;

        let min_x = component_coordinates
            .iter()
            .map(|position| position.x)
            .fold(f64::INFINITY, f64::min);
        let max_x = component_coordinates
            .iter()
            .map(|position| position.x)
            .fold(f64::NEG_INFINITY, f64::max);
        let min_y = component_coordinates
            .iter()
            .map(|position| position.y)
            .fold(f64::INFINITY, f64::min);
        let max_y = component_coordinates
            .iter()
            .map(|position| position.y)
            .fold(f64::NEG_INFINITY, f64::max);
        let center_y = (min_y + max_y) / 2.0;

        if component_index > 0 {
            cursor += FRAGMENT_GAP;
        }
        for (local_index, &global_index) in component.iter().enumerate() {
            coordinates[global_index] = Vec2::new(
                component_coordinates[local_index].x - min_x + cursor,
                component_coordinates[local_index].y - center_y,
            );
        }
        cursor += max_x - min_x;
    }

    center_coordinates(&mut coordinates);
    Ok(coordinates)
}

/// Lays out one connected molecule: rings first, then substituents, then
/// acyclic chains; centered and mirrored to the conventional handedness.
fn place_connected_molecule(molecule: &MoleculeGraph) -> Result<Vec<Vec2>, String> {
    let mut coordinates = vec![Vec2::new(0.0, 0.0); molecule.n_atoms()];
    let mut placed = vec![false; molecule.n_atoms()];

    // ── 1. Ring detection ──────────────────────────────────────────────────

    let rings = find_rings(molecule);

    // ── 2. Place ring systems first ───────────────────────────────────────

    if !rings.is_empty() {
        place_initial_ring_system(molecule, &rings, &mut coordinates, &mut placed);
    }

    // ── 3. Place ring substituents radially outward ──────────────────────

    place_substituents_for_placed_rings(molecule, &rings, &mut coordinates, &mut placed);

    // ── 4. Place remaining acyclic atoms (pure chain molecules) ──────────

    // For a ring-free molecule, grow the chain from the graph's center atom when
    // that center is a symmetric branch hub, so structures like a quaternary
    // carbon bearing four equal arms are laid out symmetrically instead of
    // lopsidedly from atom 0.
    let (root, symmetric_hub_root) = match placed.iter().position(|&is_placed| is_placed) {
        Some(atom_index) => (atom_index, false),
        None => acyclic_root(molecule),
    };
    if !placed[root] {
        coordinates[root] = Vec2::new(0.0, 0.0);
        placed[root] = true;
    }

    // A symmetric hub places its arms straddling the vertical axis so the figure
    // reads upright and mirror-symmetric; a plain chain starts at -30° so the
    // first bond is horizontal in the conventional zigzag.
    let initial_dir = if symmetric_hub_root {
        PI / 2.0 + PI / molecule.adj[root].len() as f64
    } else {
        -PI / 6.0
    };
    place_chain(
        molecule,
        root,
        initial_dir,
        symmetric_hub_root,
        &rings,
        &mut coordinates,
        &mut placed,
    );

    apply_curl_layout(molecule, &mut coordinates)?;
    apply_cis_trans_layout(molecule, &mut coordinates)?;
    separate_overlapping_branches(molecule, &rings, &mut coordinates);

    // ── 5. Center the molecule ────────────────────────────────────────────

    center_coordinates(&mut coordinates);

    // Mirror to match the layout handedness used by RDKit/PubChem, so wedges read
    // in the conventional orientation. Done before stereo assignment, so the
    // recomputed wedges still depict the correct enantiomer.
    for position in &mut coordinates {
        position.x = -position.x;
    }

    Ok(coordinates)
}

/// Connected components as sorted atom-index lists, ordered by first atom, so
/// fragments follow SMILES writing order.
fn connected_components(molecule: &MoleculeGraph) -> Vec<Vec<usize>> {
    let atom_count = molecule.n_atoms();
    let mut seen = vec![false; atom_count];
    let mut components = Vec::new();
    for start in 0..atom_count {
        if seen[start] {
            continue;
        }
        let mut component = vec![start];
        seen[start] = true;
        let mut stack = vec![start];
        while let Some(atom_index) = stack.pop() {
            for &(neighbor, _) in &molecule.adj[atom_index] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    component.push(neighbor);
                    stack.push(neighbor);
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components
}

/// Copy of one connected component with atom and bond indices renumbered to
/// 0..k, so the single-molecule placement can run on it unchanged.
fn component_subgraph(molecule: &MoleculeGraph, component: &[usize]) -> MoleculeGraph {
    let mut local_atom = vec![usize::MAX; molecule.n_atoms()];
    for (local_index, &global_index) in component.iter().enumerate() {
        local_atom[global_index] = local_index;
    }

    let mut local_bond_indices = vec![usize::MAX; molecule.bonds.len()];
    let mut bonds = Vec::new();
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        // Bonds never cross components, so checking one endpoint suffices.
        if local_atom[bond.from] != usize::MAX {
            local_bond_indices[bond_index] = bonds.len();
            let mut local_bond = bond.clone();
            local_bond.from = local_atom[bond.from];
            local_bond.to = local_atom[bond.to];
            bonds.push(local_bond);
        }
    }

    MoleculeGraph {
        atoms: component
            .iter()
            .map(|&global_index| molecule.atoms[global_index].clone())
            .collect(),
        bonds,
        adj: component
            .iter()
            .map(|&global_index| {
                molecule.adj[global_index]
                    .iter()
                    .map(|&(neighbor, bond_index)| {
                        (local_atom[neighbor], local_bond_indices[bond_index])
                    })
                    .collect()
            })
            .collect(),
        neighbor_bonds: component
            .iter()
            .map(|&global_index| {
                molecule.neighbor_bonds[global_index]
                    .iter()
                    .map(|&bond_index| local_bond_indices[bond_index])
                    .collect()
            })
            .collect(),
        has_preceding: component
            .iter()
            .map(|&global_index| molecule.has_preceding[global_index])
            .collect(),
        preceding_atom: component
            .iter()
            .map(|&global_index| {
                molecule.preceding_atom[global_index]
                    .map(|preceding_atom| local_atom[preceding_atom])
            })
            .collect(),
    }
}

/// Applies `!c` constraints after the automatic acyclic layout. For a written
/// path A-B-C!cD, every arm forward of C is reflected across the B-C axis when
/// needed so C-D repeats the A-B-C turn. Moving all forward arms together keeps
/// branch slots distinct at substituted centers.
fn apply_curl_layout(molecule: &MoleculeGraph, coordinates: &mut [Vec2]) -> Result<(), String> {
    let curl_bonds: Vec<usize> = molecule
        .bonds
        .iter()
        .enumerate()
        .filter_map(|(index, bond)| bond.curl.then_some(index))
        .collect();

    let mut curl_per_pivot = vec![0usize; molecule.n_atoms()];
    for &bond_index in &curl_bonds {
        curl_per_pivot[molecule.bonds[bond_index].from] += 1;
    }
    if let Some(pivot) = curl_per_pivot.iter().position(|&count| count > 1) {
        return Err(format!(
            "multiple !c bonds leave atom {pivot}; only one curl constraint is allowed per atom"
        ));
    }

    let rings = find_rings(molecule);
    let ring_bonds = ring_bond_set(molecule, &rings);

    for bond_index in curl_bonds {
        let bond = &molecule.bonds[bond_index];
        let pivot_atom = bond.from;
        let next_atom = bond.to;
        let preceding_atom = molecule.preceding_atom[pivot_atom].ok_or_else(|| {
            format!("!c on bond {pivot_atom}-{next_atom} needs two preceding chain bonds")
        })?;
        let first_atom = molecule.preceding_atom[preceding_atom].ok_or_else(|| {
            format!("!c on bond {pivot_atom}-{next_atom} needs two preceding chain bonds")
        })?;
        let incoming_bond = molecule
            .bond_between(preceding_atom, pivot_atom)
            .ok_or_else(|| format!("missing incoming bond {preceding_atom}-{pivot_atom} for !c"))?;
        if ring_bonds.contains(&bond_index) || ring_bonds.contains(&incoming_bond) {
            return Err("!c is not supported on a ring bond or directly after one".into());
        }

        let previous_horizontal = coordinates[preceding_atom].x - coordinates[first_atom].x;
        let previous_vertical = coordinates[preceding_atom].y - coordinates[first_atom].y;
        let incoming_horizontal = coordinates[pivot_atom].x - coordinates[preceding_atom].x;
        let incoming_vertical = coordinates[pivot_atom].y - coordinates[preceding_atom].y;
        let outgoing_horizontal = coordinates[next_atom].x - coordinates[pivot_atom].x;
        let outgoing_vertical = coordinates[next_atom].y - coordinates[pivot_atom].y;
        let previous_turn =
            previous_horizontal * incoming_vertical - previous_vertical * incoming_horizontal;
        let next_turn =
            incoming_horizontal * outgoing_vertical - incoming_vertical * outgoing_horizontal;
        if previous_turn.abs() < 1e-8 {
            return Err(format!(
                "!c on bond {pivot_atom}-{next_atom} has no preceding zigzag turn to repeat"
            ));
        }
        if next_turn.abs() < 1e-8 {
            return Err(format!(
                "!c on bond {pivot_atom}-{next_atom} cannot curl a linear continuation"
            ));
        }
        if previous_turn * next_turn > 0.0 {
            continue;
        }

        let forward_atoms = collect_subtree(molecule, pivot_atom, preceding_atom, preceding_atom);
        let reflection_axis = LineAxis {
            start: coordinates[preceding_atom],
            end: coordinates[pivot_atom],
        };
        for atom in forward_atoms {
            coordinates[atom] = reflect_point_across_line(coordinates[atom], reflection_axis);
        }
    }

    Ok(())
}

pub(crate) fn implicit_h_count(molecule: &MoleculeGraph, atom_index: usize) -> u8 {
    let atom = &molecule.atoms[atom_index];
    if atom.has_explicit_h {
        return 0;
    }

    let Some(normal_valences) = normal_valences(&atom.symbol) else {
        return 0;
    };

    let bond_order_sum: i16 = molecule.adj[atom_index]
        .iter()
        .map(|&(_, bond_index)| molecule.bonds[bond_index].order.as_u8() as i16)
        .sum();
    let required_valence = atom.charge as i16 + bond_order_sum;
    let valence = normal_valences
        .iter()
        .copied()
        .find(|&candidate| candidate >= required_valence)
        .unwrap_or(*normal_valences.last().unwrap());
    let remaining = valence - required_valence;
    remaining.max(0) as u8
}

/// Normal valences used when an unbracketed organic-subset atom needs implicit
/// hydrogens. The ordered alternatives let hypervalent P and S reach the next
/// valid valence instead of silently losing hydrogen count.
fn normal_valences(symbol: &str) -> Option<&'static [i16]> {
    match symbol {
        "C" | "Si" | "Sn" => Some(&[4]),
        "N" | "P" | "As" => Some(&[3, 5]),
        "O" => Some(&[2]),
        "S" | "Se" | "Te" => Some(&[2, 4, 6]),
        "B" => Some(&[3]),
        "F" | "Cl" | "Br" | "I" => Some(&[1]),
        _ => None,
    }
}

fn valence_electrons(symbol: &str) -> Option<i16> {
    match symbol {
        "B" => Some(3),
        "C" | "Si" | "Sn" => Some(4),
        "N" | "P" | "As" => Some(5),
        "O" | "S" | "Se" | "Te" => Some(6),
        "F" | "Cl" | "Br" | "I" => Some(7),
        _ => None,
    }
}

pub(crate) fn lone_pair_count(molecule: &MoleculeGraph, atom_index: usize) -> u8 {
    let atom = &molecule.atoms[atom_index];
    // Abbreviations hide their internal bonds, so lone pairs are never inferred
    // from the label text; they come only from an explicit `lp=N` modifier.
    if !atom.abbrev.is_empty() {
        return atom.abbrev_lone_pairs;
    }

    let Some(valence_electrons) = valence_electrons(&atom.symbol) else {
        return 0;
    };

    let bond_order_sum: i16 = molecule.adj[atom_index]
        .iter()
        .map(|&(_, bond_index)| molecule.bonds[bond_index].order.as_u8() as i16)
        .sum();
    let hydrogen_bonds = atom.hcount as i16 + implicit_h_count(molecule, atom_index) as i16;
    let nonbonding_electrons =
        valence_electrons - atom.charge as i16 - bond_order_sum - hydrogen_bonds;

    (nonbonding_electrons.max(0) / 2) as u8
}

fn lone_pair_directions(
    molecule: &MoleculeGraph,
    atom_index: usize,
    coordinates: &[Vec2],
    count: usize,
) -> Vec<Vec2> {
    if count == 0 {
        return Vec::new();
    }

    let occupied_angles: Vec<f64> = molecule.adj[atom_index]
        .iter()
        .filter_map(|&(neighbor, _)| {
            let horizontal_offset = coordinates[neighbor].x - coordinates[atom_index].x;
            let vertical_offset = coordinates[neighbor].y - coordinates[atom_index].y;
            (horizontal_offset.abs() > 1e-8 || vertical_offset.abs() > 1e-8)
                .then_some(vertical_offset.atan2(horizontal_offset))
        })
        .collect();

    let angles = if occupied_angles.is_empty() {
        spread_around_direction(PI / 2.0, count, PI)
    } else if occupied_angles.len() == 1 {
        spread_around_direction(
            normalize_angle(occupied_angles[0] + PI),
            count,
            PI * 2.0 / 3.0,
        )
    } else {
        let (best_start, best_gap) =
            largest_angular_gap(&occupied_angles).expect("occupied angles are nonempty");

        if count == 1 {
            vec![normalize_angle(best_start + best_gap / 2.0)]
        } else {
            let margin = (PI / 8.0).min(best_gap / 4.0);
            let usable_gap = (best_gap - 2.0 * margin).max(best_gap * 0.5);
            (0..count)
                .map(|direction_index| {
                    let position = (direction_index + 1) as f64 / (count + 1) as f64;
                    normalize_angle(best_start + margin + usable_gap * position)
                })
                .collect()
        }
    };

    angles
        .into_iter()
        .map(|angle| Vec2::new(angle.cos(), angle.sin()))
        .collect()
}

fn direction_as_str(direction: BondDirection) -> &'static str {
    match direction {
        BondDirection::None => "none",
        BondDirection::Up => "up",
        BondDirection::Down => "down",
    }
}

// ── OpenSMILES stereochemistry layout ────────────────────────────────────────

fn apply_cis_trans_layout(
    molecule: &MoleculeGraph,
    coordinates: &mut [Vec2],
) -> Result<(), String> {
    let mut handled_directional_bonds = HashSet::new();
    let rings = find_rings(molecule);

    // A directional bond between two double bonds describes both of them, so
    // each double bond reads every marked neighbor bond. Reorienting a shared
    // neighbor reflects across an axis through the shared atom, which keeps
    // the geometry of double bonds that were already placed.
    for double_bond in &molecule.bonds {
        if double_bond.order != BondOrder::Double {
            continue;
        }

        let left = directional_neighbors(molecule, double_bond.from, double_bond.to)?;
        let right = directional_neighbors(molecule, double_bond.to, double_bond.from)?;

        if left.is_empty() || right.is_empty() {
            continue;
        }

        // Inside a ring, the ring layout itself decides the side of each
        // neighbor; reflecting a neighbor would tear the ring apart. Stereo
        // depiction reports a ring the layout could not draw as written.
        if rings
            .iter()
            .any(|ring| ring_has_edge(ring, double_bond.from, double_bond.to))
        {
            handled_directional_bonds.extend(left.iter().map(|neighbor| neighbor.bond_index));
            handled_directional_bonds.extend(right.iter().map(|neighbor| neighbor.bond_index));
            continue;
        }

        let axis_from = coordinates[double_bond.from];
        let axis_to = coordinates[double_bond.to];
        let double_bond_axis = LineAxis {
            start: axis_from,
            end: axis_to,
        };
        for left_neighbor in left {
            handled_directional_bonds.insert(left_neighbor.bond_index);
            orient_subtree_to_side(
                molecule,
                coordinates,
                SubtreeBranch {
                    root: left_neighbor.neighbor,
                    parent: double_bond.from,
                    blocked_atom: double_bond.to,
                },
                double_bond_axis,
                left_neighbor.side,
            );
        }
        for right_neighbor in right {
            handled_directional_bonds.insert(right_neighbor.bond_index);
            orient_subtree_to_side(
                molecule,
                coordinates,
                SubtreeBranch {
                    root: right_neighbor.neighbor,
                    parent: double_bond.to,
                    blocked_atom: double_bond.from,
                },
                double_bond_axis,
                right_neighbor.side,
            );
        }
    }

    for (index, bond) in molecule.bonds.iter().enumerate() {
        if bond.direction == BondDirection::None || handled_directional_bonds.contains(&index) {
            continue;
        }
        if touches_double_bond(molecule, bond) {
            return Err("Directional / and \\ bonds must mark both ends of a double bond".into());
        }
        return Err("Directional / and \\ bonds are only supported around double bonds; use !w or !h for manual wedge drawing".into());
    }

    Ok(())
}

fn touches_double_bond(molecule: &MoleculeGraph, bond: &Bond) -> bool {
    [bond.from, bond.to].iter().any(|&atom| {
        molecule.adj[atom]
            .iter()
            .any(|&(_, bond_index)| molecule.bonds[bond_index].order == BondOrder::Double)
    })
}

#[derive(Clone, Copy)]
struct DirectionalNeighbor {
    bond_index: usize,
    neighbor: usize,
    side: i8,
}

#[derive(Clone, Copy)]
struct LineAxis {
    start: Vec2,
    end: Vec2,
}

#[derive(Clone, Copy)]
struct SubtreeBranch {
    root: usize,
    parent: usize,
    blocked_atom: usize,
}

fn directional_neighbors(
    molecule: &MoleculeGraph,
    center: usize,
    double_partner: usize,
) -> Result<Vec<DirectionalNeighbor>, String> {
    let mut found = Vec::new();
    for &(neighbor, bond_index) in &molecule.adj[center] {
        if neighbor == double_partner {
            continue;
        }
        let bond = &molecule.bonds[bond_index];
        if bond.direction == BondDirection::None {
            continue;
        }
        if bond.order != BondOrder::Single {
            return Err(
                "Directional / and \\ markers must be on single bonds adjacent to a double bond"
                    .into(),
            );
        }
        let raw = match bond.direction {
            BondDirection::Up => 1,
            BondDirection::Down => -1,
            BondDirection::None => 0,
        };
        let side = if bond.from == center { raw } else { -raw };
        let candidate = DirectionalNeighbor {
            bond_index,
            neighbor,
            side,
        };
        if found
            .iter()
            .any(|prev: &DirectionalNeighbor| prev.side == side)
        {
            return Err(
                "Conflicting or unsupported multiple directional bonds on one end of a double bond"
                    .into(),
            );
        }
        found.push(candidate);
    }
    Ok(found)
}

fn orient_subtree_to_side(
    molecule: &MoleculeGraph,
    coordinates: &mut [Vec2],
    branch: SubtreeBranch,
    axis: LineAxis,
    desired_side: i8,
) {
    let root_side = side_of_point(axis, coordinates[branch.root]);
    if root_side == 0 || root_side == desired_side {
        return;
    }

    let atoms = collect_subtree(molecule, branch.root, branch.parent, branch.blocked_atom);
    for atom in atoms {
        coordinates[atom] = reflect_point_across_line(coordinates[atom], axis);
    }
}

fn side_of_point(axis: LineAxis, point: Vec2) -> i8 {
    let cross = (axis.end.x - axis.start.x) * (point.y - axis.start.y)
        - (axis.end.y - axis.start.y) * (point.x - axis.start.x);
    if cross > 1e-8 {
        1
    } else if cross < -1e-8 {
        -1
    } else {
        0
    }
}

fn reflect_point_across_line(point: Vec2, axis: LineAxis) -> Vec2 {
    let horizontal_length = axis.end.x - axis.start.x;
    let vertical_length = axis.end.y - axis.start.y;
    let squared_length = horizontal_length * horizontal_length + vertical_length * vertical_length;
    if squared_length < 1e-10 {
        return point;
    }
    let projection_ratio = ((point.x - axis.start.x) * horizontal_length
        + (point.y - axis.start.y) * vertical_length)
        / squared_length;
    let projection = Vec2::new(
        axis.start.x + projection_ratio * horizontal_length,
        axis.start.y + projection_ratio * vertical_length,
    );
    Vec2::new(2.0 * projection.x - point.x, 2.0 * projection.y - point.y)
}

fn collect_subtree(
    molecule: &MoleculeGraph,
    root: usize,
    parent: usize,
    blocked: usize,
) -> Vec<usize> {
    let mut atoms = Vec::new();
    let mut seen = vec![false; molecule.n_atoms()];
    seen[parent] = true;
    seen[blocked] = true;

    let mut stack = vec![root];
    seen[root] = true;
    while let Some(atom) = stack.pop() {
        atoms.push(atom);
        for &(neighbor, _) in &molecule.adj[atom] {
            if !seen[neighbor] {
                seen[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }
    atoms
}

// ── Ring placement ────────────────────────────────────────────────────────────

/// Places the ring system containing the most central ring at the origin.
/// Other ring systems are anchored later, when a chain reaches them.
fn place_initial_ring_system(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    if rings.is_empty() {
        return;
    }

    let system = layout_ring_system(molecule, rings, initial_ring_index(rings));
    for &atom in &system.atoms {
        coordinates[atom] = system.positions[atom];
        placed[atom] = true;
    }
}

fn initial_ring_index(rings: &[Vec<usize>]) -> usize {
    let mut best_index = 0;
    let mut best_score = (0_usize, 0_usize);
    for (index, ring) in rings.iter().enumerate() {
        let fused_neighbors = rings
            .iter()
            .enumerate()
            .filter(|(other_index, other_ring)| {
                *other_index != index && rings_share_edge(ring, other_ring)
            })
            .count();
        let score = (fused_neighbors, ring.len());
        if score > best_score {
            best_index = index;
            best_score = score;
        }
    }
    best_index
}

/// Places the ring system containing `ring_index` so that it continues the
/// bond that reached `anchor` along `incoming_direction`: the system turns
/// until the side of `anchor` facing away from the rest of the system points
/// back along that bond.
fn place_ring_system_from_anchor(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    ring_index: usize,
    anchor: usize,
    incoming_direction: f64,
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    let system = layout_ring_system(molecule, rings, ring_index);
    let rotation = incoming_direction + PI - system.outward_angle(molecule, anchor);
    let anchor_position = coordinates[anchor];
    let local_anchor_position = system.positions[anchor];
    for &atom in &system.atoms {
        coordinates[atom] =
            anchor_position + (system.positions[atom] - local_anchor_position).rotated(rotation);
        placed[atom] = true;
    }
}

fn unfinished_ring_containing_atom(
    rings: &[Vec<usize>],
    atom: usize,
    placed: &[bool],
) -> Option<usize> {
    rings
        .iter()
        .position(|ring| ring.contains(&atom) && ring.iter().any(|&ring_atom| !placed[ring_atom]))
}

fn place_substituents_for_placed_rings(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    for ring in rings {
        if !ring.iter().all(|&atom_index| placed[atom_index]) {
            continue;
        }

        let ring_size = ring.len() as f64;
        let center_x = ring
            .iter()
            .map(|&atom_index| coordinates[atom_index].x)
            .sum::<f64>()
            / ring_size;
        let center_y = ring
            .iter()
            .map(|&atom_index| coordinates[atom_index].y)
            .sum::<f64>()
            / ring_size;
        for &atom in ring {
            let outward_direction =
                (coordinates[atom].y - center_y).atan2(coordinates[atom].x - center_x);
            place_ring_substituents(
                molecule,
                atom,
                outward_direction,
                rings,
                coordinates,
                placed,
            );
        }
    }
}

// ── Ring substituent placement ────────────────────────────────────────────────

/// Place all unplaced neighbors of a ring atom radially outward from the ring center,
/// then recursively extend any chains from those substituents.
fn place_ring_substituents(
    molecule: &MoleculeGraph,
    ring_atom: usize,
    outward_direction: f64,
    rings: &[Vec<usize>],
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    let unplaced_neighbors: Vec<usize> = molecule.adj[ring_atom]
        .iter()
        .filter(|&&(neighbor, _)| !placed[neighbor])
        .map(|&(neighbor, _)| neighbor)
        .collect();
    let directions = free_substituent_directions(
        molecule,
        ring_atom,
        unplaced_neighbors.len(),
        outward_direction,
        rings,
        coordinates,
        placed,
    );

    for (neighbor_index, &neighbor) in unplaced_neighbors.iter().enumerate() {
        let direction = directions[neighbor_index];
        let substituent_position = coordinates[ring_atom] + Vec2::from_angle(direction);
        if nearest_placed_distance(substituent_position, coordinates, placed) < CHAIN_CLEARANCE {
            mirror_blocking_branch(
                molecule,
                ring_atom,
                substituent_position,
                coordinates,
                placed,
            );
        }
        coordinates[neighbor] = substituent_position;
        placed[neighbor] = true;
        if let Some(ring_index) = unfinished_ring_containing_atom(rings, neighbor, placed) {
            place_ring_system_from_anchor(
                molecule,
                rings,
                ring_index,
                neighbor,
                direction,
                coordinates,
                placed,
            );
            place_substituents_for_placed_rings(molecule, rings, coordinates, placed);
        } else {
            place_chain(
                molecule,
                neighbor,
                direction,
                false,
                rings,
                coordinates,
                placed,
            );
        }
    }
}

fn free_substituent_directions(
    molecule: &MoleculeGraph,
    atom_index: usize,
    count: usize,
    fallback_direction: f64,
    rings: &[Vec<usize>],
    coordinates: &[Vec2],
    placed: &[bool],
) -> Vec<f64> {
    if count == 0 {
        return Vec::new();
    }

    let occupied_angles: Vec<f64> = molecule.adj[atom_index]
        .iter()
        .filter(|&&(neighbor, _)| placed[neighbor])
        .map(|&(neighbor, _)| {
            (coordinates[neighbor].y - coordinates[atom_index].y)
                .atan2(coordinates[neighbor].x - coordinates[atom_index].x)
        })
        .collect();

    if occupied_angles.len() < 2 {
        return spread_around_direction(fallback_direction, count, PI / 3.0);
    }

    let (best_start, best_gap) =
        widest_gap_outside_rings(atom_index, &occupied_angles, rings, coordinates, placed);

    let usable_gap = (best_gap - PI / 6.0).max(PI / 6.0);
    if count == 1 {
        vec![normalize_angle(best_start + best_gap / 2.0)]
    } else {
        (0..count)
            .map(|index| {
                let position = (index + 1) as f64 / (count + 1) as f64;
                normalize_angle(best_start + (best_gap - usable_gap) / 2.0 + usable_gap * position)
            })
            .collect()
    }
}

/// The widest angular gap between an atom's placed bonds. Among equally wide
/// gaps, such as the three 120° gaps at a fused-ring junction, it picks one
/// whose middle does not point into a ring containing the atom.
fn widest_gap_outside_rings(
    atom_index: usize,
    occupied_angles: &[f64],
    rings: &[Vec<usize>],
    coordinates: &[Vec2],
    placed: &[bool],
) -> (f64, f64) {
    let gaps = angular_gaps(occupied_angles);
    let widest = gaps.iter().map(|&(_, width)| width).fold(0.0, f64::max);
    let points_into_ring = |&(start, width): &(f64, f64)| {
        let probe = coordinates[atom_index] + Vec2::from_angle(start + width / 2.0) * 0.5;
        rings.iter().any(|ring| {
            ring.contains(&atom_index)
                && ring.iter().all(|&atom| placed[atom])
                && point_in_polygon(
                    probe,
                    &ring
                        .iter()
                        .map(|&atom| coordinates[atom])
                        .collect::<Vec<_>>(),
                )
        })
    };
    let widest_gaps: Vec<(f64, f64)> = gaps
        .into_iter()
        .filter(|&(_, width)| width > widest - 1e-6)
        .collect();
    widest_gaps
        .iter()
        .copied()
        .find(|gap| !points_into_ring(gap))
        .unwrap_or(widest_gaps[0])
}

fn spread_around_direction(center: f64, count: usize, spread: f64) -> Vec<f64> {
    if count == 1 {
        return vec![center];
    }
    (0..count)
        .map(|index| {
            normalize_angle(
                center - spread / 2.0 + spread * index as f64 / (count.saturating_sub(1)) as f64,
            )
        })
        .collect()
}

// ── Chain layout via DFS ──────────────────────────────────────────────────────

/// Picks the starting atom for a ring-free molecule and reports whether it is a
/// symmetric branch hub.
///
/// When the tree's center is a branch atom whose every arm carries an identical
/// subtree (e.g. a quaternary carbon bearing four equal chains), the layout
/// starts there and draws the arms with mirror symmetry. Other molecules keep
/// the conventional walk from atom 0, so their depiction is unchanged.
fn acyclic_root(molecule: &MoleculeGraph) -> (usize, bool) {
    let atom_count = molecule.n_atoms();
    if atom_count <= 1 || !is_connected(molecule) {
        return (0, false);
    }

    let mut degree: Vec<usize> = (0..atom_count)
        .map(|atom_index| molecule.adj[atom_index].len())
        .collect();
    let mut removed = vec![false; atom_count];
    let mut remaining = atom_count;
    let mut leaves: Vec<usize> = (0..atom_count)
        .filter(|&atom_index| degree[atom_index] <= 1)
        .collect();

    while remaining > 2 {
        let mut next = Vec::new();
        for &leaf in &leaves {
            if removed[leaf] {
                continue;
            }
            removed[leaf] = true;
            remaining -= 1;
            for &(neighbor, _) in &molecule.adj[leaf] {
                if !removed[neighbor] {
                    degree[neighbor] -= 1;
                    if degree[neighbor] == 1 {
                        next.push(neighbor);
                    }
                }
            }
        }
        if next.is_empty() {
            break;
        }
        leaves = next;
    }

    // One or two atoms survive: the tree center. Prefer the more-connected one,
    // breaking further ties toward the lower index for determinism.
    let center = (0..atom_count)
        .filter(|&atom_index| !removed[atom_index])
        .max_by(|&first_atom, &second_atom| {
            molecule.adj[first_atom]
                .len()
                .cmp(&molecule.adj[second_atom].len())
                .then(second_atom.cmp(&first_atom))
        })
        .unwrap_or(0);

    if is_symmetric_hub(molecule, center) {
        (center, true)
    } else {
        (0, false)
    }
}

/// True when `hub` has at least three arms and every arm leads to a subtree of
/// the same size. Such hubs are the ones whose evenly-spread arms would
/// otherwise be drawn as a rotational pinwheel.
fn is_symmetric_hub(molecule: &MoleculeGraph, hub: usize) -> bool {
    if molecule.adj[hub].len() < 3 {
        return false;
    }
    let sizes: Vec<usize> = molecule.adj[hub]
        .iter()
        .map(|&(neighbor, _)| arm_subtree_size(molecule, neighbor, hub))
        .collect();
    sizes.iter().all(|&size| size == sizes[0])
}

fn arm_subtree_size(molecule: &MoleculeGraph, start: usize, hub: usize) -> usize {
    let mut seen = vec![false; molecule.n_atoms()];
    seen[hub] = true;
    seen[start] = true;
    let mut stack = vec![start];
    let mut count = 0;
    while let Some(atom) = stack.pop() {
        count += 1;
        for &(neighbor, _) in &molecule.adj[atom] {
            if !seen[neighbor] {
                seen[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }
    count
}

fn is_connected(molecule: &MoleculeGraph) -> bool {
    let atom_count = molecule.n_atoms();
    if atom_count == 0 {
        return true;
    }
    let mut seen = vec![false; atom_count];
    let mut stack = vec![0];
    seen[0] = true;
    let mut count = 1;
    while let Some(atom) = stack.pop() {
        for &(neighbor, _) in &molecule.adj[atom] {
            if !seen[neighbor] {
                seen[neighbor] = true;
                count += 1;
                stack.push(neighbor);
            }
        }
    }
    count == atom_count
}

#[derive(Clone, Copy)]
struct PendingChainAtom {
    atom_index: usize,
    incoming_bond: Option<usize>,
    incoming_direction: f64,
    turn_sign: f64,
}

#[derive(Clone, Copy)]
struct UnplacedNeighbor {
    atom_index: usize,
    bond_index: usize,
    writing_order: usize,
}

struct ChainDirectionContext<'a> {
    molecule: &'a MoleculeGraph,
    atom_index: usize,
    incoming_direction: f64,
    turn_sign: f64,
    neighbor_count: usize,
    has_incoming_bond: bool,
    coordinates: &'a [Vec2],
    placed: &'a [bool],
}

fn place_chain(
    molecule: &MoleculeGraph,
    start: usize,
    incoming_angle: f64,
    symmetric_hub_root: bool,
    rings: &[Vec<usize>],
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    let mut pending_atoms = vec![PendingChainAtom {
        atom_index: start,
        incoming_bond: None,
        incoming_direction: incoming_angle,
        turn_sign: 1.0,
    }];

    while let Some(pending_atom) = pending_atoms.pop() {
        let atom_index = pending_atom.atom_index;
        let mut unplaced_neighbors: Vec<UnplacedNeighbor> = molecule.adj[atom_index]
            .iter()
            .enumerate()
            .filter(|(_, &(neighbor, _))| !placed[neighbor])
            .map(
                |(writing_order, &(neighbor, bond_index))| UnplacedNeighbor {
                    atom_index: neighbor,
                    bond_index,
                    writing_order,
                },
            )
            .collect();

        unplaced_neighbors.sort_by(|first, second| {
            let first_size = unplaced_subtree_size(molecule, first.atom_index, atom_index, placed);
            let second_size =
                unplaced_subtree_size(molecule, second.atom_index, atom_index, placed);
            second_size
                .cmp(&first_size)
                .then_with(|| first.writing_order.cmp(&second.writing_order))
        });

        let directions = square_planar_directions(
            molecule,
            atom_index,
            &unplaced_neighbors,
            coordinates,
            placed,
        )
        .unwrap_or_else(|| {
            chain_neighbor_directions(ChainDirectionContext {
                molecule,
                atom_index,
                incoming_direction: pending_atom.incoming_direction,
                turn_sign: pending_atom.turn_sign,
                neighbor_count: unplaced_neighbors.len(),
                has_incoming_bond: pending_atom.incoming_bond.is_some(),
                coordinates,
                placed,
            })
        });

        for (neighbor_index, neighbor) in unplaced_neighbors.iter().enumerate() {
            let selected_direction = if is_linear_atom(molecule, atom_index) {
                directions[neighbor_index]
            } else {
                resolve_chain_direction(
                    molecule,
                    atom_index,
                    directions[neighbor_index],
                    pending_atom.incoming_direction,
                    coordinates,
                    placed,
                )
            };

            coordinates[neighbor.atom_index] = Vec2::new(
                coordinates[atom_index].x + selected_direction.cos(),
                coordinates[atom_index].y + selected_direction.sin(),
            );
            placed[neighbor.atom_index] = true;

            if let Some(ring_index) =
                unfinished_ring_containing_atom(rings, neighbor.atom_index, placed)
            {
                place_ring_system_from_anchor(
                    molecule,
                    rings,
                    ring_index,
                    neighbor.atom_index,
                    selected_direction,
                    coordinates,
                    placed,
                );
                place_substituents_for_placed_rings(molecule, rings, coordinates, placed);
            } else {
                let next_turn_sign = if symmetric_hub_root && pending_atom.incoming_bond.is_none() {
                    if selected_direction.cos() >= 0.0 {
                        1.0
                    } else {
                        -1.0
                    }
                } else if normalize_angle(selected_direction - directions[neighbor_index]).abs()
                    > 1e-9
                {
                    let turn =
                        normalize_angle(selected_direction - pending_atom.incoming_direction);
                    if turn > 1e-9 {
                        -1.0
                    } else if turn < -1e-9 {
                        1.0
                    } else {
                        -pending_atom.turn_sign
                    }
                } else {
                    -pending_atom.turn_sign
                };
                pending_atoms.push(PendingChainAtom {
                    atom_index: neighbor.atom_index,
                    incoming_bond: Some(neighbor.bond_index),
                    incoming_direction: selected_direction,
                    turn_sign: next_turn_sign,
                });
            }
        }
    }
}

/// Minimum clearance between a newly placed chain atom and any atom already
/// placed elsewhere, in bond-length units.
const CHAIN_CLEARANCE: f64 = 0.55;

fn nearest_placed_distance(point: Vec2, coordinates: &[Vec2], placed: &[bool]) -> f64 {
    let mut minimum_distance = f64::INFINITY;
    for atom_index in 0..coordinates.len() {
        if placed[atom_index] {
            minimum_distance = minimum_distance.min(point.distance_to(coordinates[atom_index]));
        }
    }
    minimum_distance
}

/// Chooses the direction for the bond `u → next`. The natural (zigzag)
/// proposal is kept whenever its endpoint is clear of already-placed atoms.
/// On a collision, the preferred resolution is global: mirror the blocking
/// branch to the other side of its attachment bond, which frees the natural
/// spot and keeps both branches in textbook geometry. Otherwise the zigzag
/// turn flips to the opposite ideal slot if that one is free. Bond angles
/// are never bent to arbitrary values — every considered position keeps the
/// ideal angles a reference renderer would use — so if no ideal slot is
/// clear, the natural one is kept.
fn resolve_chain_direction(
    molecule: &MoleculeGraph,
    atom_index: usize,
    proposed_direction: f64,
    incoming_direction: f64,
    coordinates: &mut [Vec2],
    placed: &[bool],
) -> f64 {
    let origin = coordinates[atom_index];
    let endpoint =
        move |direction: f64| Vec2::new(origin.x + direction.cos(), origin.y + direction.sin());
    if nearest_placed_distance(endpoint(proposed_direction), coordinates, placed) >= CHAIN_CLEARANCE
    {
        return proposed_direction;
    }
    if mirror_blocking_branch(
        molecule,
        atom_index,
        endpoint(proposed_direction),
        coordinates,
        placed,
    ) {
        return proposed_direction;
    }

    let flipped_direction = normalize_angle(2.0 * incoming_direction - proposed_direction);
    if normalize_angle(flipped_direction - proposed_direction).abs() > 1e-9
        && nearest_placed_distance(endpoint(flipped_direction), coordinates, placed)
            >= CHAIN_CLEARANCE
    {
        let direction_is_taken = molecule.adj[atom_index].iter().any(|&(neighbor, _)| {
            placed[neighbor] && {
                let neighbor_angle = (coordinates[neighbor].y - coordinates[atom_index].y)
                    .atan2(coordinates[neighbor].x - coordinates[atom_index].x);
                normalize_angle(flipped_direction - neighbor_angle).abs() < PI / 6.0
            }
        });
        if !direction_is_taken {
            return flipped_direction;
        }
    }
    proposed_direction
}

/// Tries to clear the crowded spot `target` by reflecting the branch that
/// occupies it across its attachment bond (e.g. flipping an ortho substituent
/// to lean the other way). Candidate branches are subtrees hanging off a
/// single bond that contain every blocking atom, are fully placed, and do not
/// contain `u`. A candidate is accepted only if, after reflection, both the
/// target spot and every moved atom have full clearance; the smallest such
/// branch is flipped. Returns whether a reflection was applied.
fn mirror_blocking_branch(
    molecule: &MoleculeGraph,
    current_atom: usize,
    target: Vec2,
    coordinates: &mut [Vec2],
    placed: &[bool],
) -> bool {
    let atom_count = molecule.n_atoms();
    let blocking_atoms: Vec<usize> = (0..atom_count)
        .filter(|&atom_index| {
            placed[atom_index]
                && atom_index != current_atom
                && target.distance_to(coordinates[atom_index]) < CHAIN_CLEARANCE
        })
        .collect();
    if blocking_atoms.is_empty() {
        return false;
    }

    let mut best_reflection: Option<(Vec<usize>, usize, usize)> = None;
    for bond in &molecule.bonds {
        for (axis_start, axis_end) in [(bond.from, bond.to), (bond.to, bond.from)] {
            if !placed[axis_start] || !placed[axis_end] {
                continue;
            }
            let subtree_atoms = collect_subtree(molecule, axis_end, axis_start, axis_start);
            if subtree_atoms.contains(&current_atom) {
                continue;
            }
            if !subtree_atoms.iter().all(|&atom_index| placed[atom_index]) {
                continue;
            }
            if !blocking_atoms
                .iter()
                .all(|blocking_atom| subtree_atoms.contains(blocking_atom))
            {
                continue;
            }
            if best_reflection
                .as_ref()
                .is_some_and(|(previous_atoms, _, _)| previous_atoms.len() <= subtree_atoms.len())
            {
                continue;
            }

            let mut atom_in_subtree = vec![false; atom_count];
            for &atom_index in &subtree_atoms {
                atom_in_subtree[atom_index] = true;
            }
            let reflection_axis = LineAxis {
                start: coordinates[axis_start],
                end: coordinates[axis_end],
            };
            let reflect = |point: Vec2| reflect_point_across_line(point, reflection_axis);

            let target_clear = (0..atom_count)
                .filter(|&atom_index| placed[atom_index] && atom_index != current_atom)
                .all(|atom_index| {
                    let position = if atom_in_subtree[atom_index] {
                        reflect(coordinates[atom_index])
                    } else {
                        coordinates[atom_index]
                    };
                    target.distance_to(position) >= CHAIN_CLEARANCE
                });
            let branch_clear = target_clear
                && subtree_atoms.iter().all(|&subtree_atom| {
                    let reflected_position = reflect(coordinates[subtree_atom]);
                    (0..atom_count)
                        .filter(|&atom_index| {
                            placed[atom_index]
                                && !atom_in_subtree[atom_index]
                                && atom_index != axis_start
                        })
                        .all(|atom_index| {
                            reflected_position.distance_to(coordinates[atom_index])
                                >= CHAIN_CLEARANCE
                        })
                });
            if branch_clear {
                best_reflection = Some((subtree_atoms, axis_start, axis_end));
            }
        }
    }

    let Some((subtree_atoms, axis_start, axis_end)) = best_reflection else {
        return false;
    };
    let reflection_axis = LineAxis {
        start: coordinates[axis_start],
        end: coordinates[axis_end],
    };
    for &atom_index in &subtree_atoms {
        coordinates[atom_index] =
            reflect_point_across_line(coordinates[atom_index], reflection_axis);
    }
    true
}

/// Final collision pass. While atoms still overlap, ring-free branches that
/// contain an overlapping atom are reflected across the bond attaching them.
/// The single reflection that most reduces the total overlap is applied;
/// when none helps, pairs of reflections are tried, since two neighboring
/// branches sometimes have to turn together. Only when every ideal position
/// is taken is a branch turned slightly off its ideal angle. Branches
/// attached next to a marked cis/trans double bond or a `!c` curl keep their
/// side, since moving them would change the requested geometry.
fn separate_overlapping_branches(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    coordinates: &mut [Vec2],
) {
    let branches = flippable_branches(molecule, rings);
    for _ in 0..molecule.n_atoms() {
        let overlapping_atoms: HashSet<usize> = overlapping_atom_pairs(molecule, coordinates)
            .into_iter()
            .flat_map(|(first, second)| [first, second])
            .collect();
        if overlapping_atoms.is_empty() {
            return;
        }
        let candidates: Vec<&FlippableBranch> = branches
            .iter()
            .filter(|branch| {
                branch
                    .atoms
                    .iter()
                    .any(|atom| overlapping_atoms.contains(atom))
            })
            .collect();
        let current_overlap = total_overlap_depth(molecule, coordinates);

        let mut best = OverlapReduction::new(current_overlap);
        for branch in &candidates {
            best.consider(molecule, branch.reflected(coordinates));
        }
        if best.coordinates.is_none() {
            for first in &candidates {
                let after_first = first.reflected(coordinates);
                for second in &candidates {
                    if !std::ptr::eq(*first, *second) {
                        best.consider(molecule, second.reflected(&after_first));
                    }
                }
            }
        }
        if best.coordinates.is_none() {
            for branch in &candidates {
                for start in [coordinates.to_vec(), branch.reflected(coordinates)] {
                    for turn in BRANCH_ESCAPE_TURNS {
                        best.consider(molecule, branch.rotated(&start, turn));
                    }
                }
            }
        }

        let Some(reflected) = best.coordinates else {
            return;
        };
        coordinates.copy_from_slice(&reflected);
    }
}

/// The candidate drawing with the least total overlap seen so far, kept only
/// when it improves on the starting overlap.
struct OverlapReduction {
    overlap: f64,
    coordinates: Option<Vec<Vec2>>,
}

impl OverlapReduction {
    fn new(starting_overlap: f64) -> Self {
        Self {
            overlap: starting_overlap,
            coordinates: None,
        }
    }

    fn consider(&mut self, molecule: &MoleculeGraph, candidate: Vec<Vec2>) {
        let overlap = total_overlap_depth(molecule, &candidate);
        if overlap < self.overlap - 1e-9 {
            self.overlap = overlap;
            self.coordinates = Some(candidate);
        }
    }
}

/// Small turns of a branch about its attachment atom, tried when every
/// ideal position of the branch is already occupied.
const BRANCH_ESCAPE_TURNS: [f64; 4] = [PI / 12.0, -PI / 12.0, PI / 6.0, -PI / 6.0];

/// A ring-free branch hanging off one single bond.
struct FlippableBranch {
    attachment: usize,
    root: usize,
    atoms: Vec<usize>,
}

impl FlippableBranch {
    /// Coordinates with the branch reflected across its attachment bond.
    fn reflected(&self, coordinates: &[Vec2]) -> Vec<Vec2> {
        let axis = LineAxis {
            start: coordinates[self.attachment],
            end: coordinates[self.root],
        };
        let mut reflected = coordinates.to_vec();
        for &atom in &self.atoms {
            reflected[atom] = reflect_point_across_line(coordinates[atom], axis);
        }
        reflected
    }

    /// Coordinates with the branch turned by `angle` about its attachment
    /// atom.
    fn rotated(&self, coordinates: &[Vec2], angle: f64) -> Vec<Vec2> {
        let pivot = coordinates[self.attachment];
        let mut rotated = coordinates.to_vec();
        for &atom in &self.atoms {
            rotated[atom] = pivot + (coordinates[atom] - pivot).rotated(angle);
        }
        rotated
    }
}

fn flippable_branches(molecule: &MoleculeGraph, rings: &[Vec<usize>]) -> Vec<FlippableBranch> {
    let mut in_ring = vec![false; molecule.n_atoms()];
    for &atom in rings.iter().flatten() {
        in_ring[atom] = true;
    }
    let mut branches = Vec::new();
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        if bond.order != BondOrder::Single || is_geometry_constrained_bond(molecule, bond_index) {
            continue;
        }
        for (attachment, root) in [(bond.from, bond.to), (bond.to, bond.from)] {
            let atoms = collect_subtree(molecule, root, attachment, attachment);
            if !atoms.contains(&attachment) && atoms.iter().all(|&atom| !in_ring[atom]) {
                branches.push(FlippableBranch {
                    attachment,
                    root,
                    atoms,
                });
            }
        }
    }
    branches
}

/// Sum over overlapping atom pairs of how far they fall short of the chain
/// clearance.
fn total_overlap_depth(molecule: &MoleculeGraph, coordinates: &[Vec2]) -> f64 {
    overlapping_atom_pairs(molecule, coordinates)
        .iter()
        .map(|&(first, second)| {
            CHAIN_CLEARANCE - coordinates[first].distance_to(coordinates[second])
        })
        .sum()
}

fn overlapping_atom_pairs(molecule: &MoleculeGraph, coordinates: &[Vec2]) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for first in 0..coordinates.len() {
        for second in first + 1..coordinates.len() {
            if molecule.bond_between(first, second).is_none()
                && coordinates[first].distance_to(coordinates[second]) < CHAIN_CLEARANCE
            {
                pairs.push((first, second));
            }
        }
    }
    pairs
}

/// True for a bond whose neighborhood carries a drawing constraint: it is a
/// `!c` curl bond, or it touches an atom of a double bond that has
/// directional markers.
fn is_geometry_constrained_bond(molecule: &MoleculeGraph, bond_index: usize) -> bool {
    let bond = &molecule.bonds[bond_index];
    if bond.curl || bond.direction != BondDirection::None {
        return true;
    }
    [bond.from, bond.to].iter().any(|&atom| {
        molecule.adj[atom].iter().any(|&(_, neighbor_bond)| {
            let candidate = &molecule.bonds[neighbor_bond];
            candidate.curl
                || (candidate.order == BondOrder::Double
                    && [candidate.from, candidate.to].iter().any(|&double_atom| {
                        molecule.adj[double_atom].iter().any(|&(_, marked)| {
                            molecule.bonds[marked].direction != BondDirection::None
                        })
                    }))
        })
    })
}

/// Computes outgoing bond directions for the `count` unplaced neighbors of `u`.
///
/// Exact placement for a square-planar (`@SP`) center with four neighbors:
/// they sit at 90° steps around the atom, in the cyclic order given by the
/// shape class — the line traced through the neighbors in SMILES order reads
/// 'U' (@SP1), '4' (@SP2), or 'Z' (@SP3). An already-placed neighbor anchors
/// the rotation. Returns `None` (generic placement) for anything that is not
/// a clean four-coordinate @SP chain atom.
fn square_planar_directions(
    molecule: &MoleculeGraph,
    center: usize,
    unplaced: &[UnplacedNeighbor],
    coordinates: &[Vec2],
    placed: &[bool],
) -> Option<Vec<f64>> {
    let AtomChirality::SquarePlanar(class) = molecule.atoms[center].chirality else {
        return None;
    };
    let neighbor_bonds = &molecule.neighbor_bonds[center];
    if neighbor_bonds.len() != 4 {
        return None;
    }

    // Neighbor atoms in SMILES writing order.
    let writing: Vec<usize> = neighbor_bonds
        .iter()
        .map(|&b| {
            let bond = &molecule.bonds[b];
            if bond.from == center {
                bond.to
            } else {
                bond.from
            }
        })
        .collect();

    // corner_seq[j] = writing-order index occupying corner j; corners are 90°
    // apart. The three classes are the three ways to pair up trans neighbors.
    let corner_seq: [usize; 4] = match class {
        1 => [0, 1, 2, 3], // U: consecutive around the square
        2 => [0, 2, 1, 3], // 4: n1 trans n2
        _ => [0, 1, 3, 2], // Z: n1 trans n4
    };
    let mut corner_of = [0usize; 4];
    for (corner, &w) in corner_seq.iter().enumerate() {
        corner_of[w] = corner;
    }

    // An already-placed neighbor (the atom we were reached from, or a prior
    // fragment of the walk) fixes the square's rotation; a bare root center
    // defaults to an upright + cross.
    let base = writing
        .iter()
        .enumerate()
        .find(|&(_, &a)| placed[a])
        .map(|(w_idx, &a)| {
            let angle = (coordinates[a].y - coordinates[center].y)
                .atan2(coordinates[a].x - coordinates[center].x);
            angle - corner_of[w_idx] as f64 * (PI / 2.0)
        })
        .unwrap_or(PI);

    let directions = unplaced
        .iter()
        .map(|neighbor| {
            let w_idx = writing
                .iter()
                .position(|&atom| atom == neighbor.atom_index)?;
            Some(normalize_angle(base + corner_of[w_idx] as f64 * (PI / 2.0)))
        })
        .collect::<Option<Vec<f64>>>()?;
    Some(directions)
}

/// One or two substituents use the ±60° zigzag; three or more are spread across
/// the angular space left free by already-placed neighbors.
fn chain_neighbor_directions(context: ChainDirectionContext<'_>) -> Vec<f64> {
    if context.neighbor_count == 0 {
        return Vec::new();
    }

    if is_linear_atom(context.molecule, context.atom_index) && context.has_incoming_bond {
        return vec![context.incoming_direction; context.neighbor_count];
    }

    match context.neighbor_count {
        1 => vec![context.incoming_direction + context.turn_sign * (PI / 3.0)],
        2 => vec![
            context.incoming_direction + context.turn_sign * (PI / 3.0),
            context.incoming_direction - context.turn_sign * (PI / 3.0),
        ],
        _ => {
            let occupied: Vec<f64> = context.molecule.adj[context.atom_index]
                .iter()
                .filter(|&&(neighbor, _)| context.placed[neighbor])
                .map(|&(neighbor, _)| {
                    (context.coordinates[neighbor].y - context.coordinates[context.atom_index].y)
                        .atan2(
                            context.coordinates[neighbor].x
                                - context.coordinates[context.atom_index].x,
                        )
                })
                .collect();
            let slots = distribute_in_free_space(
                &occupied,
                context.neighbor_count,
                context.incoming_direction,
            );
            assign_slots_center_out(slots)
        }
    }
}

/// Reorders gap slots (in spatial order) so that the callers' size-sorted
/// neighbors are assigned outward from the middle of the free gap. The first
/// neighbor (largest subtree) takes the central slot, which continues the main
/// chain away from the already-placed atoms; smaller substituents fan out to the
/// sides. This keeps long branches from folding back over the rest of the
/// molecule at a crowded branch point.
fn assign_slots_center_out(slots: Vec<f64>) -> Vec<f64> {
    let slot_count = slots.len();
    let center = (slot_count as f64 - 1.0) / 2.0;
    let mut order: Vec<usize> = (0..slot_count).collect();
    order.sort_by(|&first_index, &second_index| {
        let first_distance = (first_index as f64 - center).abs();
        let second_distance = (second_index as f64 - center).abs();
        first_distance
            .partial_cmp(&second_distance)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(first_index.cmp(&second_index))
    });

    let mut assigned = vec![0.0; slot_count];
    for (neighbor_rank, &slot_index) in order.iter().enumerate() {
        assigned[neighbor_rank] = slots[slot_index];
    }
    assigned
}

/// Spreads `count` directions across the largest angular gap left by `occupied`
/// (the directions of already-placed neighbors). With nothing placed, the
/// directions are spaced evenly around the full circle starting at `fallback`.
fn distribute_in_free_space(
    occupied_angles: &[f64],
    count: usize,
    fallback_direction: f64,
) -> Vec<f64> {
    if count == 0 {
        return Vec::new();
    }

    if occupied_angles.is_empty() {
        return (0..count)
            .map(|index| {
                normalize_angle(fallback_direction + 2.0 * PI * index as f64 / count as f64)
            })
            .collect();
    }

    let (best_start, best_gap) =
        largest_angular_gap(occupied_angles).expect("occupied angles are nonempty");
    let segments = (count + 1) as f64;
    (1..=count)
        .map(|index| normalize_angle(best_start + best_gap * index as f64 / segments))
        .collect()
}

fn unplaced_subtree_size(
    molecule: &MoleculeGraph,
    root: usize,
    parent: usize,
    placed: &[bool],
) -> usize {
    let mut size = 0;
    let mut seen = vec![false; molecule.n_atoms()];
    seen[parent] = true;
    let mut stack = vec![root];
    seen[root] = true;

    while let Some(atom) = stack.pop() {
        if placed[atom] {
            continue;
        }
        size += 1;
        for &(neighbor, _) in &molecule.adj[atom] {
            if !seen[neighbor] {
                seen[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }

    size
}

fn is_linear_atom(molecule: &MoleculeGraph, atom_index: usize) -> bool {
    if molecule.adj[atom_index].len() != 2 {
        return false;
    }

    let mut double_count = 0;
    let mut has_triple = false;
    for &(_, bond_index) in &molecule.adj[atom_index] {
        match molecule.bonds[bond_index].order {
            BondOrder::Double => double_count += 1,
            BondOrder::Triple | BondOrder::Quadruple => has_triple = true,
            BondOrder::Single | BondOrder::Aromatic => {}
        }
    }

    has_triple || double_count == 2
}

// ── Virtual-H placement ───────────────────────────────────────────────────────

/// Returns the angle (radians) toward which the H-label group of a bracket atom
/// should be placed.  The label collapses all hydrogens into one glyph (e.g. "H₄"),
/// so a single direction is enough.  For atoms with no bonds the label sits to the
/// east (0°); for bonded atoms it points into the largest free angular gap.
fn hydrogen_label_angle(occupied_angles: &[f64]) -> f64 {
    if occupied_angles.is_empty() {
        return 0.0;
    }

    let (best_start, best_gap) =
        largest_angular_gap(occupied_angles).expect("occupied angles are nonempty");
    normalize_angle(best_start + best_gap / 2.0)
}

fn center_coordinates(coordinates: &mut [Vec2]) {
    if coordinates.is_empty() {
        return;
    }
    let center_x =
        coordinates.iter().map(|position| position.x).sum::<f64>() / coordinates.len() as f64;
    let center_y =
        coordinates.iter().map(|position| position.y).sum::<f64>() / coordinates.len() as f64;
    for position in coordinates {
        position.x -= center_x;
        position.y -= center_y;
    }
}

/// For each bond, returns the unit vector pointing from the bond midpoint
/// toward the centroid of the smallest ring containing that bond.
/// Returns (0.0, 0.0) for bonds not in any ring.
fn ring_inner_directions(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    coordinates: &[Vec2],
) -> Vec<(f64, f64)> {
    let mut directions = vec![(0.0_f64, 0.0_f64); molecule.bonds.len()];

    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        let best_ring = best_ring_for_inner_bond(molecule, rings, bond.from, bond.to);

        if let Some(ring) = best_ring {
            let ring_size = ring.len() as f64;
            let center_x = ring
                .iter()
                .map(|&atom_index| coordinates[atom_index].x)
                .sum::<f64>()
                / ring_size;
            let center_y = ring
                .iter()
                .map(|&atom_index| coordinates[atom_index].y)
                .sum::<f64>()
                / ring_size;

            let midpoint_x = (coordinates[bond.from].x + coordinates[bond.to].x) / 2.0;
            let midpoint_y = (coordinates[bond.from].y + coordinates[bond.to].y) / 2.0;

            let horizontal_direction = center_x - midpoint_x;
            let vertical_direction = center_y - midpoint_y;
            let direction_length = (horizontal_direction * horizontal_direction
                + vertical_direction * vertical_direction)
                .sqrt();
            if direction_length > 1e-6 {
                directions[bond_index] = (
                    horizontal_direction / direction_length,
                    vertical_direction / direction_length,
                );
            }
        }
    }

    directions
}

fn best_ring_for_inner_bond<'a>(
    molecule: &MoleculeGraph,
    rings: &'a [Vec<usize>],
    from: usize,
    to: usize,
) -> Option<&'a Vec<usize>> {
    rings
        .iter()
        .filter(|ring| ring_has_edge(ring, from, to))
        .max_by(|a, b| {
            ring_unsaturation_score(molecule, a)
                .cmp(&ring_unsaturation_score(molecule, b))
                .then_with(|| b.len().cmp(&a.len()))
        })
}

fn ring_unsaturation_score(molecule: &MoleculeGraph, ring: &[usize]) -> usize {
    (0..ring.len())
        .filter_map(|ring_index| {
            molecule.bond_between(ring[ring_index], ring[(ring_index + 1) % ring.len()])
        })
        .filter(|&bond_index| {
            matches!(
                molecule.bonds[bond_index].order,
                BondOrder::Double | BondOrder::Aromatic
            )
        })
        .count()
}

fn bounding_box(coordinates: &[Vec2]) -> (f64, f64) {
    if coordinates.is_empty() {
        return (0.0, 0.0);
    }
    let min_x = coordinates
        .iter()
        .map(|position| position.x)
        .fold(f64::INFINITY, f64::min);
    let max_x = coordinates
        .iter()
        .map(|position| position.x)
        .fold(f64::NEG_INFINITY, f64::max);
    let min_y = coordinates
        .iter()
        .map(|position| position.y)
        .fold(f64::INFINITY, f64::min);
    let max_y = coordinates
        .iter()
        .map(|position| position.y)
        .fold(f64::NEG_INFINITY, f64::max);
    (max_x - min_x + 1.0, max_y - min_y + 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring_atoms_with_substituent_inside(smiles: &str) -> Vec<usize> {
        let layout = crate::layout_native(smiles).unwrap();
        let molecule = crate::parse_molecule(smiles).unwrap();
        let rings = find_rings(&molecule);
        let in_ring: HashSet<usize> = rings.iter().flatten().copied().collect();
        (0..molecule.n_atoms())
            .filter(|atom| !in_ring.contains(atom))
            .filter(|&atom| {
                rings.iter().any(|ring| {
                    let corners: Vec<Vec2> = ring
                        .iter()
                        .map(|&ring_atom| layout.atoms[ring_atom].pos)
                        .collect();
                    crate::geometry::point_in_polygon(layout.atoms[atom].pos, &corners)
                })
            })
            .collect()
    }

    #[test]
    fn junction_methyls_point_away_from_the_fused_rings() {
        for smiles in [
            "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
            "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        ] {
            assert!(
                ring_atoms_with_substituent_inside(smiles).is_empty(),
                "{smiles}"
            );
        }
    }

    #[test]
    fn crowded_macrocycle_side_chains_are_separated() {
        let layout = crate::layout_native(
            "CCC1NC(=O)C(C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(C)NC(=O)C(CC(C)C)N(C)C(=O)C(C(C)C)NC(=O)C(CC(C)C)N(C)C(=O)CN(C)C1=O",
        )
        .unwrap();
        let heavy: Vec<Vec2> = layout
            .atoms
            .iter()
            .filter(|atom| !atom.virtual_h)
            .map(|atom| atom.pos)
            .collect();
        for first in 0..heavy.len() {
            for second in first + 1..heavy.len() {
                assert!(
                    heavy[first].distance_to(heavy[second]) > 0.5,
                    "atoms {first} and {second} overlap"
                );
            }
        }
    }

    #[test]
    fn steroid_ring_system_detects_four_rings() {
        let molecule =
            crate::parse_molecule("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
                .expect("steroid-like molecule should parse");
        let rings = find_rings(&molecule);
        assert_eq!(rings.len(), 4, "rings: {rings:?}");
    }

    #[test]
    fn ring_closure_after_branch_stays_on_branch_point() {
        let molecule = crate::parse_molecule("C1=CCCC(=O)1").expect("cyclopentenone should parse");
        let rings = find_rings(&molecule);
        assert!(rings.iter().any(|ring| ring.len() == 5), "rings: {rings:?}");

        let complex = "O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5";
        let molecule = crate::parse_molecule(complex).expect("complex fused system should parse");
        let rings = find_rings(&molecule);
        assert!(
            rings
                .iter()
                .any(|ring| ring.len() == 5 && ring.contains(&19)),
            "rings: {rings:?}"
        );
    }

    #[test]
    fn fused_double_bond_prefers_unsaturated_ring_side() {
        let molecule =
            crate::parse_molecule("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
                .expect("steroid-like molecule should parse");
        let rings = find_rings(&molecule);
        let mut checked_shared_double = false;

        for bond in molecule
            .bonds
            .iter()
            .filter(|bond| bond.order == BondOrder::Double)
        {
            let containing: Vec<&Vec<usize>> = rings
                .iter()
                .filter(|ring| ring_has_edge(ring, bond.from, bond.to))
                .collect();
            if containing.len() < 2 {
                continue;
            }

            let selected = best_ring_for_inner_bond(&molecule, &rings, bond.from, bond.to)
                .expect("ring missing");
            let selected_score = ring_unsaturation_score(&molecule, selected);
            let max_score = containing
                .iter()
                .map(|ring| ring_unsaturation_score(&molecule, ring))
                .max()
                .unwrap();
            assert_eq!(selected_score, max_score);
            checked_shared_double = true;
        }

        assert!(checked_shared_double, "expected a fused/shared double bond");
    }
}
