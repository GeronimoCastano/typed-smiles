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
use std::collections::{HashSet, VecDeque};
use std::f64::consts::PI;

use crate::graph::{AtomChirality, BondDirection, BondOrder, BondStereo, MoleculeGraph};
use crate::render::{AromaticRing, AtomOutput, BondOutput, LayoutOutput, Vec2};

/// Gap between the bounding boxes of dot-separated fragments, in bond lengths.
const FRAGMENT_GAP: f64 = 1.5;

pub fn compute_layout(molecule: &MoleculeGraph) -> Result<LayoutOutput, String> {
    if molecule.n_atoms() == 0 {
        return Ok(empty_layout_output());
    }

    let coordinates = layout_coordinates(molecule)?;
    let rings = find_rings(molecule);
    let ring_bonds = ring_bond_set(molecule, &rings);
    let mut rendered_stereo: Vec<BondStereo> =
        molecule.bonds.iter().map(|bond| bond.stereo).collect();
    let mut hydrogen_stereo = vec![None; molecule.n_atoms()];
    apply_tetrahedral_stereo(
        molecule,
        &coordinates,
        &ring_bonds,
        &mut rendered_stereo,
        &mut hydrogen_stereo,
    )?;

    let mut atoms = build_atom_outputs(molecule, &coordinates, &hydrogen_stereo);
    let inner_directions = ring_inner_directions(molecule, &rings, &coordinates);
    let mut bonds = build_bond_outputs(molecule, &rendered_stereo, &inner_directions);
    append_virtual_hydrogen_outputs(
        molecule,
        &coordinates,
        &hydrogen_stereo,
        &mut atoms,
        &mut bonds,
    );

    let all_positions: Vec<Vec2> = atoms.iter().map(|atom| atom.pos).collect();
    let (bbox_width, bbox_height) = bounding_box(&all_positions);

    Ok(LayoutOutput {
        atoms,
        bonds,
        aromatic_rings: aromatic_ring_circles(molecule, &rings, &coordinates),
        bbox_width,
        bbox_height,
    })
}

fn empty_layout_output() -> LayoutOutput {
    LayoutOutput {
        atoms: Vec::new(),
        bonds: Vec::new(),
        aromatic_rings: Vec::new(),
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
    rendered_stereo: &[BondStereo],
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
            stereo: rendered_stereo[bond_index].as_str().to_string(),
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
            bond_between(molecule, first_atom, second_atom)
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
        let incoming_bond = bond_between(molecule, preceding_atom, pivot_atom)
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

    let Some(valence) = standard_valence(&atom.symbol) else {
        return 0;
    };

    let bond_order_sum: i16 = molecule.adj[atom_index]
        .iter()
        .map(|&(_, bond_index)| molecule.bonds[bond_index].order.as_u8() as i16)
        .sum();
    let remaining = valence - atom.charge as i16 - bond_order_sum;
    remaining.max(0) as u8
}

fn standard_valence(symbol: &str) -> Option<i16> {
    match symbol {
        "C" | "Si" | "Sn" => Some(4),
        "N" | "P" | "As" => Some(3),
        "O" | "S" | "Se" | "Te" => Some(2),
        "B" => Some(3),
        "F" | "Cl" | "Br" | "I" => Some(1),
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

fn lone_pair_count(molecule: &MoleculeGraph, atom_index: usize) -> u8 {
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

    for double_bond in &molecule.bonds {
        if double_bond.order != BondOrder::Double {
            continue;
        }

        let left = directional_neighbors(
            molecule,
            double_bond.from,
            double_bond.to,
            &handled_directional_bonds,
        )?;
        let right = directional_neighbors(
            molecule,
            double_bond.to,
            double_bond.from,
            &handled_directional_bonds,
        )?;

        if left.is_empty() || right.is_empty() {
            if !left.is_empty() || !right.is_empty() {
                return Err(
                    "Directional / and \\ bonds must mark both ends of a double bond".into(),
                );
            }
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
        if bond.direction != BondDirection::None && !handled_directional_bonds.contains(&index) {
            return Err("Directional / and \\ bonds are only supported around double bonds; use !w or !h for manual wedge drawing".into());
        }
    }

    Ok(())
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
    handled: &HashSet<usize>,
) -> Result<Vec<DirectionalNeighbor>, String> {
    let mut found = Vec::new();
    for &(neighbor, bond_index) in &molecule.adj[center] {
        if neighbor == double_partner || handled.contains(&bond_index) {
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

/// Assigns wedge/hash bonds for tetrahedral stereocenters.
///
/// Chooses wedge vs. hash from the signed volume of the four neighbor directions
/// (in OpenSMILES order) so the depicted 3D structure reproduces the requested
/// chirality on our 2D layout. `@` requires a negative signed volume, `@@` a
/// positive one.
fn apply_tetrahedral_stereo(
    molecule: &MoleculeGraph,
    coordinates: &[Vec2],
    ring_bonds: &HashSet<usize>,
    rendered_stereo: &mut [BondStereo],
    stereo_h: &mut [Option<(BondStereo, Vec2)>],
) -> Result<(), String> {
    for (center, atom) in molecule.atoms.iter().enumerate() {
        let parity = match atom.chirality {
            // Square-planar centers are depicted exactly by the flat layout;
            // TB/OH/AL centers are accepted without stereo decoration.
            AtomChirality::None | AtomChirality::SquarePlanar(_) | AtomChirality::Undepicted => {
                continue
            }
            AtomChirality::Unsupported => {
                return Err("Unsupported chirality class".into());
            }
            AtomChirality::TetraAnti => -1.0, // @  ⇒ negative signed volume
            AtomChirality::TetraClockwise => 1.0, // @@ ⇒ positive signed volume
        };

        let neighbor_bonds = &molecule.neighbor_bonds[center];
        let hydrogen_count = (atom.hcount + implicit_h_count(molecule, center)) as usize;

        // Only clean tetrahedral centers (4 substituents, at most one of them H)
        // can be depicted unambiguously.
        if neighbor_bonds.len() + hydrogen_count != 4 || hydrogen_count > 1 {
            continue;
        }

        // Neighbors in OpenSMILES order, with the implicit/bracket hydrogen placed
        // after the "from" atom (or first if there is none).
        let mut neighbor_order: Vec<TetrahedralNeighbor> = neighbor_bonds
            .iter()
            .map(|&bond_index| {
                let bond = &molecule.bonds[bond_index];
                let atom_index = if bond.from == center {
                    bond.to
                } else {
                    bond.from
                };
                TetrahedralNeighbor::Bond {
                    bond_index,
                    atom_index,
                }
            })
            .collect();
        if hydrogen_count == 1 {
            let hydrogen_position = if molecule.has_preceding[center] { 1 } else { 0 };
            neighbor_order.insert(
                hydrogen_position.min(neighbor_order.len()),
                TetrahedralNeighbor::Hydrogen,
            );
        }

        let hydrogen_direction =
            stereochemical_hydrogen_direction(molecule, center, coordinates, None);

        let selected_bond =
            preferred_tetrahedral_bond(molecule, center, ring_bonds, rendered_stereo);
        let out_of_plane_index = if let Some(bond_index) = selected_bond {
            neighbor_order
                .iter()
                .position(|neighbor| {
                    matches!(
                        neighbor,
                        TetrahedralNeighbor::Bond {
                            bond_index: candidate,
                            ..
                        } if *candidate == bond_index
                    )
                })
                .expect("selected tetrahedral bond belongs to the center")
        } else if let Some(hydrogen_position) = neighbor_order
            .iter()
            .position(|neighbor| matches!(neighbor, TetrahedralNeighbor::Hydrogen))
        {
            hydrogen_position
        } else {
            continue;
        };

        let mut directions = [[0.0_f64; 3]; 4];
        for (neighbor_index, neighbor) in neighbor_order.iter().enumerate() {
            directions[neighbor_index] = if neighbor_index == out_of_plane_index {
                let direction = match neighbor {
                    TetrahedralNeighbor::Bond { atom_index, .. } => {
                        normalize_vector(vector_from(coordinates[center], coordinates[*atom_index]))
                    }
                    TetrahedralNeighbor::Hydrogen => hydrogen_direction,
                };
                [direction.x, direction.y, 1.0]
            } else {
                match neighbor {
                    TetrahedralNeighbor::Bond { atom_index, .. } => {
                        let direction = normalize_vector(vector_from(
                            coordinates[center],
                            coordinates[*atom_index],
                        ));
                        [direction.x, direction.y, 0.0]
                    }
                    TetrahedralNeighbor::Hydrogen => {
                        [hydrogen_direction.x, hydrogen_direction.y, -1.0]
                    }
                }
            };
        }

        let volume = signed_volume(&directions);
        if volume.abs() < 1e-9 {
            continue;
        }
        let stereo = if volume.signum() == parity {
            BondStereo::WedgeUp
        } else {
            BondStereo::WedgeDown
        };

        if let Some(bond_index) = selected_bond {
            rendered_stereo[bond_index] = stereo;
        } else {
            stereo_h[center] = Some((stereo, hydrogen_direction));
        }
    }
    Ok(())
}

/// One neighbor of a tetrahedral center in OpenSMILES ordering.
#[derive(Clone, Copy)]
enum TetrahedralNeighbor {
    Bond {
        bond_index: usize,
        atom_index: usize,
    },
    Hydrogen,
}

fn vector_from(origin: Vec2, destination: Vec2) -> Vec2 {
    Vec2::new(destination.x - origin.x, destination.y - origin.y)
}

fn normalize_vector(vector: Vec2) -> Vec2 {
    let length = (vector.x * vector.x + vector.y * vector.y).sqrt();
    if length > 1e-12 {
        Vec2::new(vector.x / length, vector.y / length)
    } else {
        vector
    }
}

/// Signed volume `(d1-d0)·((d2-d0)×(d3-d0))` of four 3D points.
pub(crate) fn signed_volume(points: &[[f64; 3]; 4]) -> f64 {
    let first_offset = [
        points[1][0] - points[0][0],
        points[1][1] - points[0][1],
        points[1][2] - points[0][2],
    ];
    let second_offset = [
        points[2][0] - points[0][0],
        points[2][1] - points[0][1],
        points[2][2] - points[0][2],
    ];
    let third_offset = [
        points[3][0] - points[0][0],
        points[3][1] - points[0][1],
        points[3][2] - points[0][2],
    ];
    first_offset[0] * (second_offset[1] * third_offset[2] - second_offset[2] * third_offset[1])
        - first_offset[1]
            * (second_offset[0] * third_offset[2] - second_offset[2] * third_offset[0])
        + first_offset[2]
            * (second_offset[0] * third_offset[1] - second_offset[1] * third_offset[0])
}

pub(crate) fn stereochemical_hydrogen_direction(
    molecule: &MoleculeGraph,
    atom_index: usize,
    coordinates: &[Vec2],
    preferred_direction: Option<f64>,
) -> Vec2 {
    let occupied_angles: Vec<f64> = molecule.adj[atom_index]
        .iter()
        .map(|&(neighbor, _)| {
            (coordinates[neighbor].y - coordinates[atom_index].y)
                .atan2(coordinates[neighbor].x - coordinates[atom_index].x)
        })
        .collect();

    if occupied_angles.is_empty() {
        return Vec2::new(0.0, -1.0);
    }

    if let Some(direction) = preferred_direction {
        if occupied_angles
            .iter()
            .all(|&angle| angle_delta(direction, angle).abs() > PI / 5.0)
        {
            return Vec2::new(direction.cos(), direction.sin());
        }
    }

    let (best_start, best_gap) =
        largest_angular_gap(&occupied_angles).expect("occupied angles are nonempty");
    let direction = normalize_angle(best_start + best_gap / 2.0);
    Vec2::new(direction.cos(), direction.sin())
}

fn preferred_tetrahedral_bond(
    molecule: &MoleculeGraph,
    atom_index: usize,
    ring_bonds: &HashSet<usize>,
    rendered_stereo: &[BondStereo],
) -> Option<usize> {
    let mut best: Option<(usize, i32)> = None;
    for &(neighbor, bond_index) in &molecule.adj[atom_index] {
        let bond = &molecule.bonds[bond_index];
        if bond.order != BondOrder::Single || rendered_stereo[bond_index] != BondStereo::None {
            continue;
        }

        let score = tetrahedral_bond_score(molecule, neighbor, bond_index, ring_bonds);
        if best
            .map(|(_, best_score)| score > best_score)
            .unwrap_or(true)
        {
            best = Some((bond_index, score));
        }
    }
    best.and_then(|(bond_index, score)| if score > 0 { Some(bond_index) } else { None })
}

fn tetrahedral_bond_score(
    molecule: &MoleculeGraph,
    neighbor: usize,
    bond_index: usize,
    ring_bonds: &HashSet<usize>,
) -> i32 {
    let neighbor_atom = &molecule.atoms[neighbor];
    let in_ring = ring_bonds.contains(&bond_index);
    let is_carbon = neighbor_atom.symbol == "C" || neighbor_atom.symbol == "c";
    let is_visible = !is_carbon || !neighbor_atom.abbrev.is_empty() || neighbor_atom.charge != 0;
    let is_terminal = molecule.adj[neighbor].len() == 1;

    let mut score = 0;
    if !in_ring {
        score += 100;
    }
    if is_visible {
        score += 50;
    }
    if is_terminal {
        score += 10;
    }
    if neighbor_atom.chirality != AtomChirality::None {
        score -= 20;
    }
    if in_ring {
        score -= 100;
    }
    score
}

fn ring_bond_set(molecule: &MoleculeGraph, rings: &[Vec<usize>]) -> HashSet<usize> {
    let mut ring_bonds = HashSet::new();
    for ring in rings {
        for (first_atom, second_atom) in ring_edges(ring) {
            if let Some(bond_index) = bond_between(molecule, first_atom, second_atom) {
                ring_bonds.insert(bond_index);
            }
        }
    }
    ring_bonds
}

// ── Ring detection ────────────────────────────────────────────────────────────

/// Returns cycles as lists of atom indices (the ring path).
fn find_rings(molecule: &MoleculeGraph) -> Vec<Vec<usize>> {
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
        if let Some(bond_index) = bond_between(molecule, first_atom, second_atom) {
            bits[bond_index / 64] |= 1_u64 << (bond_index % 64);
        }
    }
    bits
}

fn bond_between(molecule: &MoleculeGraph, first_atom: usize, second_atom: usize) -> Option<usize> {
    molecule.adj[first_atom]
        .iter()
        .find_map(|&(neighbor, bond_index)| {
            if neighbor == second_atom {
                Some(bond_index)
            } else {
                None
            }
        })
}

fn ring_edges(ring: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
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

// ── Ring placement ────────────────────────────────────────────────────────────

fn place_initial_ring_system(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    if rings.is_empty() {
        return;
    }

    // Place a central ring in the first ring system. Other standalone ring
    // systems in the same molecule are anchored later when the connecting
    // chain reaches them.
    let initial_ring = initial_ring_index(rings);
    place_regular_ring(
        &rings[initial_ring],
        Vec2::new(0.0, 0.0),
        90.0_f64.to_radians(),
        coordinates,
        placed,
    );

    let mut placed_rings: HashSet<usize> = HashSet::new();
    placed_rings.insert(initial_ring);
    place_connected_fused_rings(molecule, rings, &mut placed_rings, coordinates, placed);
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

fn rings_share_edge(first_ring: &[usize], second_ring: &[usize]) -> bool {
    shared_atoms(first_ring, second_ring)
        .windows(2)
        .any(|atom_pair| {
            ring_has_edge(first_ring, atom_pair[0], atom_pair[1])
                && ring_has_edge(second_ring, atom_pair[0], atom_pair[1])
        })
}

fn shared_atoms(first_ring: &[usize], second_ring: &[usize]) -> Vec<usize> {
    let mut atoms: Vec<usize> = first_ring
        .iter()
        .copied()
        .filter(|atom| second_ring.contains(atom))
        .collect();
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

fn ring_has_edge(ring: &[usize], first_atom: usize, second_atom: usize) -> bool {
    ring_edges(ring).any(|(edge_start, edge_end)| {
        (edge_start == first_atom && edge_end == second_atom)
            || (edge_start == second_atom && edge_end == first_atom)
    })
}

fn place_connected_fused_rings(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    placed_rings: &mut HashSet<usize>,
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    loop {
        let mut progressed = false;
        for (ring_index, ring) in rings.iter().enumerate() {
            if placed_rings.contains(&ring_index) {
                continue;
            }

            if placed_shared_edge(ring, placed).is_some() {
                place_fused_ring(ring, coordinates, placed);
                placed_rings.insert(ring_index);
                progressed = true;
            }
        }

        if progressed {
            continue;
        }

        // No edge-fused ring is ready: try a spiro ring, which joins the placed
        // structure at a single shared atom. Placing one may expose further
        // edge fusions, so re-enter the loop afterwards.
        let spiro = rings.iter().enumerate().find(|(index, ring)| {
            !placed_rings.contains(index)
                && ring
                    .iter()
                    .filter(|&&atom_index| placed[atom_index])
                    .count()
                    == 1
        });
        if let Some((ring_index, ring)) = spiro {
            place_spiro_ring(molecule, ring, coordinates, placed);
            placed_rings.insert(ring_index);
            continue;
        }

        break;
    }
}

/// Place a ring joined to the already-placed structure at a single shared
/// (spiro) atom, as a regular polygon opening away from that atom's placed
/// neighbors so the two rings do not overlap.
fn place_spiro_ring(
    molecule: &MoleculeGraph,
    ring: &[usize],
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    let ring_size = ring.len();
    let shared_position = (0..ring_size).find(|&ring_position| placed[ring[ring_position]]);
    let Some(shared_position) = shared_position else {
        place_regular_ring(ring, Vec2::new(0.0, 0.0), 0.0, coordinates, placed);
        return;
    };
    let spiro_atom = ring[shared_position];
    let spiro_position = coordinates[spiro_atom];

    let mut inward_horizontal = 0.0;
    let mut inward_vertical = 0.0;
    for bond in &molecule.bonds {
        let neighbor = if bond.from == spiro_atom {
            bond.to
        } else if bond.to == spiro_atom {
            bond.from
        } else {
            continue;
        };
        if placed[neighbor] {
            let horizontal_offset = coordinates[neighbor].x - spiro_position.x;
            let vertical_offset = coordinates[neighbor].y - spiro_position.y;
            let distance =
                (horizontal_offset * horizontal_offset + vertical_offset * vertical_offset).sqrt();
            if distance > 1e-6 {
                inward_horizontal += horizontal_offset / distance;
                inward_vertical += vertical_offset / distance;
            }
        }
    }
    let inward_length =
        (inward_horizontal * inward_horizontal + inward_vertical * inward_vertical).sqrt();
    let (inward_x, inward_y) = if inward_length < 1e-6 {
        (0.0, -1.0)
    } else {
        (
            inward_horizontal / inward_length,
            inward_vertical / inward_length,
        )
    };

    let radius = regular_polygon_radius(ring_size);
    let center = Vec2::new(
        spiro_position.x - inward_x * radius,
        spiro_position.y - inward_y * radius,
    );
    let start_angle = (spiro_position.y - center.y).atan2(spiro_position.x - center.x);
    let angle_step = 2.0 * PI / ring_size as f64;
    for (ring_position, &atom) in ring.iter().enumerate() {
        if !placed[atom] {
            let position_offset = ring_position as i64 - shared_position as i64;
            let angle = start_angle + angle_step * position_offset as f64;
            coordinates[atom] = Vec2::new(
                center.x + radius * angle.cos(),
                center.y + radius * angle.sin(),
            );
            placed[atom] = true;
        }
    }
}

/// Place a standalone ring as a regular n-gon centered at `center`.
/// `start_angle` is the angle (radians) of the first atom from the center.
fn place_regular_ring(
    ring: &[usize],
    center: Vec2,
    start_angle: f64,
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    let ring_size = ring.len() as f64;
    let radius = regular_polygon_radius(ring.len());

    for (ring_position, &atom) in ring.iter().enumerate() {
        if !placed[atom] {
            let angle = start_angle + (2.0 * PI / ring_size) * ring_position as f64;
            coordinates[atom] = Vec2::new(
                center.x + radius * angle.cos(),
                center.y + radius * angle.sin(),
            );
            placed[atom] = true;
        }
    }
}

fn regular_polygon_radius(vertex_count: usize) -> f64 {
    1.0 / (2.0 * (PI / vertex_count as f64).sin())
}

/// Place a fused ring where two atoms are already positioned.
fn place_fused_ring(ring: &[usize], coordinates: &mut [Vec2], placed: &mut [bool]) {
    let ring_size = ring.len();
    let Some((first_shared_position, second_shared_position)) = placed_shared_edge(ring, placed)
    else {
        place_regular_ring(ring, Vec2::new(0.0, 0.0), 0.0, coordinates, placed);
        return;
    };

    let first_shared_position_coordinates = coordinates[ring[first_shared_position]];
    let second_shared_position_coordinates = coordinates[ring[second_shared_position]];
    let midpoint = Vec2::new(
        (first_shared_position_coordinates.x + second_shared_position_coordinates.x) / 2.0,
        (first_shared_position_coordinates.y + second_shared_position_coordinates.y) / 2.0,
    );
    let bond_angle = (second_shared_position_coordinates.y - first_shared_position_coordinates.y)
        .atan2(second_shared_position_coordinates.x - first_shared_position_coordinates.x);

    let ring_size_as_float = ring_size as f64;
    let center_distance = 0.5 / (PI / ring_size_as_float).tan();

    let first_candidate_center = Vec2::new(
        midpoint.x + center_distance * (bond_angle + PI / 2.0).cos(),
        midpoint.y + center_distance * (bond_angle + PI / 2.0).sin(),
    );
    let second_candidate_center = Vec2::new(
        midpoint.x + center_distance * (bond_angle - PI / 2.0).cos(),
        midpoint.y + center_distance * (bond_angle - PI / 2.0).sin(),
    );
    let center = pick_farther_center(
        first_candidate_center,
        second_candidate_center,
        placed,
        coordinates,
    );

    let start_angle = (first_shared_position_coordinates.y - center.y)
        .atan2(first_shared_position_coordinates.x - center.x);
    let next_angle = (second_shared_position_coordinates.y - center.y)
        .atan2(second_shared_position_coordinates.x - center.x);
    let angle_step = 2.0 * PI / ring_size_as_float;
    let sweep = if angle_delta(start_angle + angle_step, next_angle).abs()
        <= angle_delta(start_angle - angle_step, next_angle).abs()
    {
        1.0
    } else {
        -1.0
    };

    let radius = regular_polygon_radius(ring_size);
    for (ring_position, &atom) in ring.iter().enumerate() {
        if !placed[atom] {
            let position_offset = ring_position as i64 - first_shared_position as i64;
            let angle = start_angle + sweep * angle_step * position_offset as f64;
            coordinates[atom] = Vec2::new(
                center.x + radius * angle.cos(),
                center.y + radius * angle.sin(),
            );
            placed[atom] = true;
        }
    }
}

fn angle_delta(first_angle: f64, second_angle: f64) -> f64 {
    let mut delta = first_angle - second_angle;
    while delta > PI {
        delta -= 2.0 * PI;
    }
    while delta < -PI {
        delta += 2.0 * PI;
    }
    delta
}

fn placed_shared_edge(ring: &[usize], placed: &[bool]) -> Option<(usize, usize)> {
    for first_position in 0..ring.len() {
        let second_position = (first_position + 1) % ring.len();
        if placed[ring[first_position]] && placed[ring[second_position]] {
            return Some((first_position, second_position));
        }
    }
    None
}

fn pick_farther_center(
    first_candidate: Vec2,
    second_candidate: Vec2,
    placed: &[bool],
    coordinates: &[Vec2],
) -> Vec2 {
    let mut first_total_distance = 0.0;
    let mut second_total_distance = 0.0;
    let mut placed_count = 0;
    for (atom_index, was_placed) in placed.iter().enumerate() {
        if *was_placed {
            first_total_distance += first_candidate.distance_to(coordinates[atom_index]);
            second_total_distance += second_candidate.distance_to(coordinates[atom_index]);
            placed_count += 1;
        }
    }
    if placed_count == 0 || first_total_distance >= second_total_distance {
        first_candidate
    } else {
        second_candidate
    }
}

fn place_ring_from_anchor(
    ring: &[usize],
    anchor: usize,
    center_direction: f64,
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    let Some(anchor_position) = ring.iter().position(|&atom_index| atom_index == anchor) else {
        return;
    };

    let ring_size = ring.len() as f64;
    let radius = regular_polygon_radius(ring.len());
    let center = Vec2::new(
        coordinates[anchor].x + radius * center_direction.cos(),
        coordinates[anchor].y + radius * center_direction.sin(),
    );
    let angle_step = 2.0 * PI / ring_size;
    let anchor_angle = (coordinates[anchor].y - center.y).atan2(coordinates[anchor].x - center.x);
    let start_angle = anchor_angle - angle_step * anchor_position as f64;

    place_regular_ring(ring, center, start_angle, coordinates, placed);
}

fn place_ring_system_from_anchor(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    ring_index: usize,
    anchor: usize,
    center_direction: f64,
    coordinates: &mut [Vec2],
    placed: &mut [bool],
) {
    place_ring_from_anchor(
        &rings[ring_index],
        anchor,
        center_direction,
        coordinates,
        placed,
    );

    let mut placed_rings: HashSet<usize> = HashSet::new();
    placed_rings.insert(ring_index);
    place_connected_fused_rings(molecule, rings, &mut placed_rings, coordinates, placed);
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
        coordinates,
        placed,
    );

    for (neighbor_index, &neighbor) in unplaced_neighbors.iter().enumerate() {
        let direction = directions[neighbor_index];
        coordinates[neighbor] = Vec2::new(
            coordinates[ring_atom].x + direction.cos(),
            coordinates[ring_atom].y + direction.sin(),
        );
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
        largest_angular_gap(&occupied_angles).expect("occupied angles are nonempty");

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

fn normalize_angle(angle: f64) -> f64 {
    let mut angle = angle;
    while angle > PI {
        angle -= 2.0 * PI;
    }
    while angle <= -PI {
        angle += 2.0 * PI;
    }
    angle
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

fn largest_angular_gap(angles: &[f64]) -> Option<(f64, f64)> {
    if angles.is_empty() {
        return None;
    }

    let mut sorted_angles = angles.to_vec();
    sorted_angles.sort_by(|first, second| {
        first
            .partial_cmp(second)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut best_start = sorted_angles[0];
    let mut best_gap = 0.0;
    for angle_index in 0..sorted_angles.len() {
        let start = sorted_angles[angle_index];
        let end = if angle_index + 1 < sorted_angles.len() {
            sorted_angles[angle_index + 1]
        } else {
            sorted_angles[0] + 2.0 * PI
        };
        let gap = end - start;
        if gap > best_gap {
            best_gap = gap;
            best_start = start;
        }
    }
    Some((best_start, best_gap))
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
            bond_between(
                molecule,
                ring[ring_index],
                ring[(ring_index + 1) % ring.len()],
            )
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
    use crate::graph::MoleculeGraph;

    #[test]
    fn steroid_ring_system_detects_four_rings() {
        let molecule = MoleculeGraph::from_smiles(
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
            Vec::new(),
            Vec::new(),
        )
        .expect("steroid-like molecule should parse");
        let rings = find_rings(&molecule);
        assert_eq!(rings.len(), 4, "rings: {rings:?}");
    }

    #[test]
    fn ring_closure_after_branch_stays_on_branch_point() {
        let pre = crate::preprocess_smiles("C1=CCCC(=O)1").expect("preprocess failed");
        let molecule = MoleculeGraph::from_smiles(
            &pre.smiles,
            pre.forced_direction_markers,
            pre.aromatic_atom_markers,
        )
        .expect("cyclopentenone should parse");
        let rings = find_rings(&molecule);
        assert!(rings.iter().any(|ring| ring.len() == 5), "rings: {rings:?}");

        let complex = "O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5";
        let pre = crate::preprocess_smiles(complex).expect("preprocess failed");
        let molecule = MoleculeGraph::from_smiles(
            &pre.smiles,
            pre.forced_direction_markers,
            pre.aromatic_atom_markers,
        )
        .expect("complex fused system should parse");
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
        let molecule = MoleculeGraph::from_smiles(
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
            Vec::new(),
            Vec::new(),
        )
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
