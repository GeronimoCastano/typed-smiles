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

pub fn compute_layout(mol: &MoleculeGraph) -> Result<LayoutOutput, String> {
    if mol.n_atoms() == 0 {
        return Ok(LayoutOutput {
            atoms: vec![],
            bonds: vec![],
            aromatic_rings: vec![],
            bbox_width: 0.0,
            bbox_height: 0.0,
        });
    }

    let coords = layout_coords(mol)?;

    // ── Build output ──────────────────────────────────────────────────────

    let rings = find_rings(mol);
    let ring_bonds = ring_bond_set(mol, &rings);
    let mut rendered_stereo: Vec<BondStereo> = mol.bonds.iter().map(|b| b.stereo).collect();
    let mut stereo_h: Vec<Option<(BondStereo, Vec2)>> = vec![None; mol.n_atoms()];
    apply_tetrahedral_stereo(
        mol,
        &coords,
        &ring_bonds,
        &mut rendered_stereo,
        &mut stereo_h,
    )?;

    let mut atoms: Vec<AtomOutput> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let lone_pairs = lone_pair_count(mol, i);
            AtomOutput {
                symbol: a.symbol.clone(),
                pos: coords[i],
                hcount: a.hcount,
                implicit_h: implicit_h_count(mol, i),
                charge: a.charge,
                isotope: a.isotope.unwrap_or(0),
                lone_pairs,
                lone_pair_dirs: lone_pair_directions(mol, i, &coords, lone_pairs as usize),
                abbrev: a.abbrev.clone(),
                abbrev_style: a.abbrev_style.clone(),
                abbrev_anchor: a.abbrev_anchor,
                abbrev_anchor_len: a.abbrev_anchor_len,
                chirality: a.chirality.as_str().to_string(),
                stereo_h: stereo_h[i]
                    .map(|(stereo, _)| stereo.as_str().to_string())
                    .unwrap_or_else(|| "none".to_string()),
                stereo_h_dir: stereo_h[i].map(|(_, dir)| dir).unwrap_or_default(),
                virtual_h: false,
            }
        })
        .collect();

    let inner_dirs = ring_inner_directions(mol, &rings, &coords);

    let mut bonds: Vec<BondOutput> = mol
        .bonds
        .iter()
        .enumerate()
        .map(|(i, b)| BondOutput {
            from: b.from,
            to: b.to,
            order: b.order.as_u8(),
            stereo: rendered_stereo[i].as_str().to_string(),
            forced_stereo: b.forced_stereo,
            direction: direction_as_str(b.direction).to_string(),
            inner_x: inner_dirs[i].0,
            inner_y: inner_dirs[i].1,
            virtual_bond: false,
            aromatic: b.aromatic,
        })
        .collect();

    // Append virtual H atoms so that bracket-notation H groups and implicit
    // heteroatom H labels are addressable by index.  Each parent atom gets at
    // most one virtual H regardless of hcount, because all hydrogens collapse
    // into a single rendered glyph ("H₄", "H₂", etc.).
    let push_virtual_h = |atoms: &mut Vec<AtomOutput>,
                          bonds: &mut Vec<BondOutput>,
                          parent_idx: usize,
                          dir: Vec2| {
        let parent_pos = atoms[parent_idx].pos;
        let h_idx = atoms.len();
        atoms.push(AtomOutput {
            symbol: "H".to_string(),
            // ~0.35 bond-lengths: roughly where the H glyph sits within an
            // inline label at default font size.
            pos: Vec2::new(parent_pos.x + dir.x * 0.35, parent_pos.y + dir.y * 0.35),
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
            chirality: "none".to_string(),
            stereo_h: "none".to_string(),
            stereo_h_dir: Vec2::default(),
            virtual_h: true,
        });
        bonds.push(BondOutput {
            from: parent_idx,
            to: h_idx,
            order: 1,
            stereo: "none".to_string(),
            forced_stereo: false,
            direction: "none".to_string(),
            inner_x: 0.0,
            inner_y: 0.0,
            virtual_bond: true,
            aromatic: false,
        });
    };

    for i in 0..mol.n_atoms() {
        let mol_atom = &mol.atoms[i];
        let parent_pos = coords[i];
        let existing_adj: Vec<f64> = mol.adj[i]
            .iter()
            .map(|&(j, _)| (coords[j].y - parent_pos.y).atan2(coords[j].x - parent_pos.x))
            .collect();

        if mol_atom.has_explicit_h && mol_atom.hcount > 0 {
            // Bracket atom with explicit H (e.g. [OH-], [NH4+]).
            // Skip atoms whose H is drawn as a stereo label (stereo_h set) or
            // absorbed into a tetrahedral stereo bond assignment (which causes
            // preferred_tetrahedral_bond to carry the stereo on a neighbor
            // bond instead of drawing an H label). Non-tetrahedral chirality
            // classes never consume the H, so their label still renders.
            if stereo_h[i].is_some() || mol_atom.chirality.is_tetrahedral() {
                continue;
            }
            let angle = h_label_angle(&existing_adj);
            let dir = Vec2::new(angle.cos(), angle.sin());
            push_virtual_h(&mut atoms, &mut bonds, i, dir);
        } else if !mol_atom.has_explicit_h {
            // Non-bracket heteroatom with implicit H (e.g. plain N, O, S).
            // Only emit when show-hetero-h would render a label (non-carbon,
            // non-wildcard, has implicit H).
            let sym = &atoms[i].symbol;
            let is_hetero = sym != "C" && sym != "c" && sym != "H" && sym != "*";
            let ih = atoms[i].implicit_h;
            if !is_hetero || ih == 0 {
                continue;
            }
            let angle = h_label_angle(&existing_adj);
            let dir = Vec2::new(angle.cos(), angle.sin());
            push_virtual_h(&mut atoms, &mut bonds, i, dir);
        }
    }

    let all_positions: Vec<Vec2> = atoms.iter().map(|a| a.pos).collect();
    let (bbox_w, bbox_h) = bounding_box(&all_positions);

    Ok(LayoutOutput {
        atoms,
        bonds,
        aromatic_rings: aromatic_ring_circles(mol, &rings, &coords),
        bbox_width: bbox_w,
        bbox_height: bbox_h,
    })
}

/// Circle parameters for every ring whose bonds were all aromatic in the
/// input, enabling the inscribed-circle depiction. The radius scales with the
/// ring's inradius so the circle stays clear of the bonds for any ring size.
fn aromatic_ring_circles(
    mol: &MoleculeGraph,
    rings: &[Vec<usize>],
    coords: &[Vec2],
) -> Vec<AromaticRing> {
    let mut circles = Vec::new();
    for ring in rings {
        let all_aromatic = (0..ring.len()).all(|k| {
            let a = ring[k];
            let b = ring[(k + 1) % ring.len()];
            bond_between(mol, a, b)
                .map(|idx| mol.bonds[idx].aromatic)
                .unwrap_or(false)
        });
        if !all_aromatic {
            continue;
        }
        let n = ring.len() as f64;
        let cx = ring.iter().map(|&i| coords[i].x).sum::<f64>() / n;
        let cy = ring.iter().map(|&i| coords[i].y).sum::<f64>() / n;
        let center = Vec2::new(cx, cy);
        // Inradius: distance from the centroid to the nearest bond midpoint.
        let inradius = (0..ring.len())
            .map(|k| {
                let a = coords[ring[k]];
                let b = coords[ring[(k + 1) % ring.len()]];
                center.dist(Vec2::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0))
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
fn layout_coords(mol: &MoleculeGraph) -> Result<Vec<Vec2>, String> {
    let components = connected_components(mol);
    if components.len() <= 1 {
        return place_molecule(mol);
    }

    let mut coords = vec![Vec2::new(0.0, 0.0); mol.n_atoms()];
    let mut cursor = 0.0;
    for (k, comp) in components.iter().enumerate() {
        let sub = fragment_subgraph(mol, comp);
        let sub_coords = place_molecule(&sub)?;

        let min_x = sub_coords.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let max_x = sub_coords.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let min_y = sub_coords.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let max_y = sub_coords.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        let mid_y = (min_y + max_y) / 2.0;

        if k > 0 {
            cursor += FRAGMENT_GAP;
        }
        for (local, &global) in comp.iter().enumerate() {
            coords[global] = Vec2::new(
                sub_coords[local].x - min_x + cursor,
                sub_coords[local].y - mid_y,
            );
        }
        cursor += max_x - min_x;
    }

    center_coords(&mut coords);
    Ok(coords)
}

/// Lays out one connected molecule: rings first, then substituents, then
/// acyclic chains; centered and mirrored to the conventional handedness.
fn place_molecule(mol: &MoleculeGraph) -> Result<Vec<Vec2>, String> {
    let mut coords = vec![Vec2::new(0.0, 0.0); mol.n_atoms()];
    let mut placed = vec![false; mol.n_atoms()];

    // ── 1. Ring detection ──────────────────────────────────────────────────

    let rings = find_rings(mol);

    // ── 2. Place ring systems first ───────────────────────────────────────

    if !rings.is_empty() {
        place_initial_ring_system(mol, &rings, &mut coords, &mut placed);
    }

    // ── 3. Place ring substituents radially outward ──────────────────────

    place_substituents_for_placed_rings(mol, &rings, &mut coords, &mut placed);

    // ── 4. Place remaining acyclic atoms (pure chain molecules) ──────────

    // For a ring-free molecule, grow the chain from the graph's center atom when
    // that center is a symmetric branch hub, so structures like a quaternary
    // carbon bearing four equal arms are laid out symmetrically instead of
    // lopsidedly from atom 0.
    let (root, symmetric_hub_root) = match placed.iter().position(|&p| p) {
        Some(p) => (p, false),
        None => acyclic_root(mol),
    };
    if !placed[root] {
        coords[root] = Vec2::new(0.0, 0.0);
        placed[root] = true;
    }

    // A symmetric hub places its arms straddling the vertical axis so the figure
    // reads upright and mirror-symmetric; a plain chain starts at -30° so the
    // first bond is horizontal in the conventional zigzag.
    let initial_dir = if symmetric_hub_root {
        PI / 2.0 + PI / mol.adj[root].len() as f64
    } else {
        -PI / 6.0
    };
    place_chain(
        mol,
        root,
        initial_dir,
        symmetric_hub_root,
        &rings,
        &mut coords,
        &mut placed,
    );

    apply_cis_trans_layout(mol, &mut coords)?;

    // ── 5. Center the molecule ────────────────────────────────────────────

    center_coords(&mut coords);

    // Mirror to match the layout handedness used by RDKit/PubChem, so wedges read
    // in the conventional orientation. Done before stereo assignment, so the
    // recomputed wedges still depict the correct enantiomer.
    for c in coords.iter_mut() {
        c.x = -c.x;
    }

    Ok(coords)
}

/// Connected components as sorted atom-index lists, ordered by first atom, so
/// fragments follow SMILES writing order.
fn connected_components(mol: &MoleculeGraph) -> Vec<Vec<usize>> {
    let n = mol.n_atoms();
    let mut seen = vec![false; n];
    let mut components = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut component = vec![start];
        seen[start] = true;
        let mut stack = vec![start];
        while let Some(u) = stack.pop() {
            for &(v, _) in &mol.adj[u] {
                if !seen[v] {
                    seen[v] = true;
                    component.push(v);
                    stack.push(v);
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
fn fragment_subgraph(mol: &MoleculeGraph, component: &[usize]) -> MoleculeGraph {
    let mut local_atom = vec![usize::MAX; mol.n_atoms()];
    for (local, &global) in component.iter().enumerate() {
        local_atom[global] = local;
    }

    let mut local_bond = vec![usize::MAX; mol.bonds.len()];
    let mut bonds = Vec::new();
    for (b_idx, bond) in mol.bonds.iter().enumerate() {
        // Bonds never cross components, so checking one endpoint suffices.
        if local_atom[bond.from] != usize::MAX {
            local_bond[b_idx] = bonds.len();
            let mut b = bond.clone();
            b.from = local_atom[bond.from];
            b.to = local_atom[bond.to];
            bonds.push(b);
        }
    }

    MoleculeGraph {
        atoms: component.iter().map(|&g| mol.atoms[g].clone()).collect(),
        bonds,
        adj: component
            .iter()
            .map(|&g| {
                mol.adj[g]
                    .iter()
                    .map(|&(v, b)| (local_atom[v], local_bond[b]))
                    .collect()
            })
            .collect(),
        neighbor_bonds: component
            .iter()
            .map(|&g| mol.neighbor_bonds[g].iter().map(|&b| local_bond[b]).collect())
            .collect(),
        has_preceding: component.iter().map(|&g| mol.has_preceding[g]).collect(),
    }
}

pub(crate) fn implicit_h_count(mol: &MoleculeGraph, atom_idx: usize) -> u8 {
    let atom = &mol.atoms[atom_idx];
    if atom.has_explicit_h {
        return 0;
    }

    let Some(valence) = standard_valence(&atom.symbol) else {
        return 0;
    };

    let bond_order_sum: i16 = mol.adj[atom_idx]
        .iter()
        .map(|&(_, bond_idx)| mol.bonds[bond_idx].order.as_u8() as i16)
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

fn lone_pair_count(mol: &MoleculeGraph, atom_idx: usize) -> u8 {
    let atom = &mol.atoms[atom_idx];
    if atom.abbrev != "" {
        return 0;
    }

    let Some(valence_electrons) = valence_electrons(&atom.symbol) else {
        return 0;
    };

    let bond_order_sum: i16 = mol.adj[atom_idx]
        .iter()
        .map(|&(_, bond_idx)| mol.bonds[bond_idx].order.as_u8() as i16)
        .sum();
    let hydrogen_bonds = atom.hcount as i16 + implicit_h_count(mol, atom_idx) as i16;
    let nonbonding_electrons =
        valence_electrons - atom.charge as i16 - bond_order_sum - hydrogen_bonds;

    (nonbonding_electrons.max(0) / 2) as u8
}

fn lone_pair_directions(
    mol: &MoleculeGraph,
    atom_idx: usize,
    coords: &[Vec2],
    count: usize,
) -> Vec<Vec2> {
    if count == 0 {
        return Vec::new();
    }

    let mut occupied: Vec<f64> = mol.adj[atom_idx]
        .iter()
        .filter_map(|&(neighbor, _)| {
            let dx = coords[neighbor].x - coords[atom_idx].x;
            let dy = coords[neighbor].y - coords[atom_idx].y;
            (dx.abs() > 1e-8 || dy.abs() > 1e-8).then_some(dy.atan2(dx))
        })
        .collect();

    let angles = if occupied.is_empty() {
        spread_around_direction(PI / 2.0, count, PI)
    } else if occupied.len() == 1 {
        spread_around_direction(normalize_angle(occupied[0] + PI), count, PI * 2.0 / 3.0)
    } else {
        occupied.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut best_start = occupied[0];
        let mut best_gap = 0.0;
        for i in 0..occupied.len() {
            let start = occupied[i];
            let end = if i + 1 < occupied.len() {
                occupied[i + 1]
            } else {
                occupied[0] + 2.0 * PI
            };
            let gap = end - start;
            if gap > best_gap {
                best_gap = gap;
                best_start = start;
            }
        }

        if count == 1 {
            vec![normalize_angle(best_start + best_gap / 2.0)]
        } else {
            let margin = (PI / 8.0).min(best_gap / 4.0);
            let usable_gap = (best_gap - 2.0 * margin).max(best_gap * 0.5);
            (0..count)
                .map(|i| {
                    let t = (i + 1) as f64 / (count + 1) as f64;
                    normalize_angle(best_start + margin + usable_gap * t)
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

fn apply_cis_trans_layout(mol: &MoleculeGraph, coords: &mut Vec<Vec2>) -> Result<(), String> {
    let mut handled_directional_bonds = HashSet::new();

    for double_bond in &mol.bonds {
        if double_bond.order != BondOrder::Double {
            continue;
        }

        let left = directional_neighbors(
            mol,
            double_bond.from,
            double_bond.to,
            &handled_directional_bonds,
        )?;
        let right = directional_neighbors(
            mol,
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

        let axis_from = coords[double_bond.from];
        let axis_to = coords[double_bond.to];
        for left in left {
            handled_directional_bonds.insert(left.bond_idx);
            orient_subtree_to_side(
                mol,
                coords,
                left.neighbor,
                double_bond.from,
                double_bond.to,
                axis_from,
                axis_to,
                left.side,
            );
        }
        for right in right {
            handled_directional_bonds.insert(right.bond_idx);
            orient_subtree_to_side(
                mol,
                coords,
                right.neighbor,
                double_bond.to,
                double_bond.from,
                axis_from,
                axis_to,
                right.side,
            );
        }
    }

    for (idx, bond) in mol.bonds.iter().enumerate() {
        if bond.direction != BondDirection::None && !handled_directional_bonds.contains(&idx) {
            return Err("Directional / and \\ bonds are only supported around double bonds; use !w or !h for manual wedge drawing".into());
        }
    }

    Ok(())
}

#[derive(Clone, Copy)]
struct DirectionalNeighbor {
    bond_idx: usize,
    neighbor: usize,
    side: i8,
}

fn directional_neighbors(
    mol: &MoleculeGraph,
    center: usize,
    double_partner: usize,
    handled: &HashSet<usize>,
) -> Result<Vec<DirectionalNeighbor>, String> {
    let mut found = Vec::new();
    for &(neighbor, bond_idx) in &mol.adj[center] {
        if neighbor == double_partner || handled.contains(&bond_idx) {
            continue;
        }
        let bond = &mol.bonds[bond_idx];
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
            bond_idx,
            neighbor,
            side,
        };
        if found.iter().any(|prev: &DirectionalNeighbor| prev.side == side) {
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
    mol: &MoleculeGraph,
    coords: &mut Vec<Vec2>,
    root: usize,
    parent: usize,
    blocked: usize,
    axis_from: Vec2,
    axis_to: Vec2,
    desired_side: i8,
) {
    let root_side = side_of_point(axis_from, axis_to, coords[root]);
    if root_side == 0 || root_side == desired_side {
        return;
    }

    let atoms = collect_subtree(mol, root, parent, blocked);
    for atom in atoms {
        coords[atom] = reflect_point_across_line(coords[atom], axis_from, axis_to);
    }
}

fn side_of_point(a: Vec2, b: Vec2, p: Vec2) -> i8 {
    let cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    if cross > 1e-8 {
        1
    } else if cross < -1e-8 {
        -1
    } else {
        0
    }
}

fn reflect_point_across_line(p: Vec2, a: Vec2, b: Vec2) -> Vec2 {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-10 {
        return p;
    }
    let t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
    let proj = Vec2::new(a.x + t * dx, a.y + t * dy);
    Vec2::new(2.0 * proj.x - p.x, 2.0 * proj.y - p.y)
}

fn collect_subtree(mol: &MoleculeGraph, root: usize, parent: usize, blocked: usize) -> Vec<usize> {
    let mut atoms = Vec::new();
    let mut seen = vec![false; mol.n_atoms()];
    seen[parent] = true;
    seen[blocked] = true;

    let mut stack = vec![root];
    seen[root] = true;
    while let Some(atom) = stack.pop() {
        atoms.push(atom);
        for &(neighbor, _) in &mol.adj[atom] {
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
    mol: &MoleculeGraph,
    coords: &[Vec2],
    ring_bonds: &HashSet<usize>,
    rendered_stereo: &mut Vec<BondStereo>,
    stereo_h: &mut Vec<Option<(BondStereo, Vec2)>>,
) -> Result<(), String> {
    for (center, atom) in mol.atoms.iter().enumerate() {
        let parity = match atom.chirality {
            // Square-planar centers are depicted exactly by the flat layout;
            // TB/OH/AL centers are accepted without stereo decoration.
            AtomChirality::None
            | AtomChirality::SquarePlanar(_)
            | AtomChirality::Undepicted => continue,
            AtomChirality::Unsupported => {
                return Err("Unsupported chirality class".into());
            }
            AtomChirality::TetraAnti => -1.0,      // @  ⇒ negative signed volume
            AtomChirality::TetraClockwise => 1.0,  // @@ ⇒ positive signed volume
        };

        let neighbor_bonds = &mol.neighbor_bonds[center];
        let n_h = (atom.hcount + implicit_h_count(mol, center)) as usize;

        // Only clean tetrahedral centers (4 substituents, at most one of them H)
        // can be depicted unambiguously.
        if neighbor_bonds.len() + n_h != 4 || n_h > 1 {
            continue;
        }

        // Neighbors in OpenSMILES order, with the implicit/bracket hydrogen placed
        // after the "from" atom (or first if there is none).
        let mut order: Vec<Slot> = neighbor_bonds
            .iter()
            .map(|&bond| {
                let b = &mol.bonds[bond];
                let other = if b.from == center { b.to } else { b.from };
                Slot::Bond { bond, other }
            })
            .collect();
        if n_h == 1 {
            let h_pos = if mol.has_preceding[center] { 1 } else { 0 };
            order.insert(h_pos.min(order.len()), Slot::Hydrogen);
        }

        // Direction the implicit hydrogen points (largest angular gap).
        let h_dir = stereo_h_direction(mol, center, coords, None);

        // Draw one neighbor out of plane: prefer an exocyclic substituent, else the H.
        let selected_bond = preferred_tetrahedral_bond(mol, center, ring_bonds, rendered_stereo);
        let out_idx = if let Some(bond_idx) = selected_bond {
            order
                .iter()
                .position(|s| matches!(s, Slot::Bond { bond, .. } if *bond == bond_idx))
                .unwrap()
        } else if let Some(pos) = order.iter().position(|s| matches!(s, Slot::Hydrogen)) {
            pos
        } else {
            continue; // No exocyclic substituent and no hydrogen: cannot depict.
        };

        // Trial geometry with the out-of-plane neighbor toward the viewer (z = +1);
        // an undrawn implicit hydrogen sits on the far side (z = -1).
        let mut dirs = [[0.0_f64; 3]; 4];
        for (i, slot) in order.iter().enumerate() {
            dirs[i] = if i == out_idx {
                let d = match slot {
                    Slot::Bond { other, .. } => unit(sub(coords[*other], coords[center])),
                    Slot::Hydrogen => h_dir,
                };
                [d.x, d.y, 1.0]
            } else {
                match slot {
                    Slot::Bond { other, .. } => {
                        let d = unit(sub(coords[*other], coords[center]));
                        [d.x, d.y, 0.0]
                    }
                    Slot::Hydrogen => [h_dir.x, h_dir.y, -1.0],
                }
            };
        }

        let volume = signed_volume(&dirs);
        if volume.abs() < 1e-9 {
            continue; // Degenerate 2D geometry; leave undecorated.
        }
        // z = +1 reproduces the requested parity ⇒ the bond points at the viewer.
        let stereo = if volume.signum() == parity {
            BondStereo::WedgeUp
        } else {
            BondStereo::WedgeDown
        };

        if let Some(bond_idx) = selected_bond {
            rendered_stereo[bond_idx] = stereo;
        } else {
            stereo_h[center] = Some((stereo, h_dir));
        }
    }
    Ok(())
}

/// One neighbor of a tetrahedral center in OpenSMILES ordering.
#[derive(Clone, Copy)]
enum Slot {
    Bond { bond: usize, other: usize },
    Hydrogen,
}

fn sub(a: Vec2, b: Vec2) -> Vec2 {
    Vec2::new(a.x - b.x, a.y - b.y)
}

fn unit(v: Vec2) -> Vec2 {
    let len = (v.x * v.x + v.y * v.y).sqrt();
    if len > 1e-12 {
        Vec2::new(v.x / len, v.y / len)
    } else {
        v
    }
}

/// Signed volume `(d1-d0)·((d2-d0)×(d3-d0))` of four 3D points.
pub(crate) fn signed_volume(d: &[[f64; 3]; 4]) -> f64 {
    let a = [d[1][0] - d[0][0], d[1][1] - d[0][1], d[1][2] - d[0][2]];
    let b = [d[2][0] - d[0][0], d[2][1] - d[0][1], d[2][2] - d[0][2]];
    let c = [d[3][0] - d[0][0], d[3][1] - d[0][1], d[3][2] - d[0][2]];
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

pub(crate) fn stereo_h_direction(
    mol: &MoleculeGraph,
    atom_idx: usize,
    coords: &[Vec2],
    preferred_dir: Option<f64>,
) -> Vec2 {
    let mut occupied: Vec<f64> = mol.adj[atom_idx]
        .iter()
        .map(|&(neighbor, _)| {
            (coords[neighbor].y - coords[atom_idx].y).atan2(coords[neighbor].x - coords[atom_idx].x)
        })
        .collect();

    if occupied.is_empty() {
        return Vec2::new(0.0, -1.0);
    }

    if let Some(dir) = preferred_dir {
        if occupied
            .iter()
            .all(|&angle| angle_delta(dir, angle).abs() > PI / 5.0)
        {
            return Vec2::new(dir.cos(), dir.sin());
        }
    }

    occupied.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut best_start = occupied[0];
    let mut best_gap = 0.0;
    for i in 0..occupied.len() {
        let start = occupied[i];
        let end = if i + 1 < occupied.len() {
            occupied[i + 1]
        } else {
            occupied[0] + 2.0 * PI
        };
        let gap = end - start;
        if gap > best_gap {
            best_gap = gap;
            best_start = start;
        }
    }

    let dir = normalize_angle(best_start + best_gap / 2.0);
    Vec2::new(dir.cos(), dir.sin())
}

fn preferred_tetrahedral_bond(
    mol: &MoleculeGraph,
    atom_idx: usize,
    ring_bonds: &HashSet<usize>,
    rendered_stereo: &[BondStereo],
) -> Option<usize> {
    let mut best: Option<(usize, i32)> = None;
    for &(neighbor, bond_idx) in &mol.adj[atom_idx] {
        let bond = &mol.bonds[bond_idx];
        if bond.order != BondOrder::Single || rendered_stereo[bond_idx] != BondStereo::None {
            continue;
        }

        let score = tetrahedral_bond_score(mol, neighbor, bond_idx, ring_bonds);
        if best.map(|(_, best_score)| score > best_score).unwrap_or(true) {
            best = Some((bond_idx, score));
        }
    }
    best.and_then(|(bond_idx, score)| if score > 0 { Some(bond_idx) } else { None })
}

fn tetrahedral_bond_score(
    mol: &MoleculeGraph,
    neighbor: usize,
    bond_idx: usize,
    ring_bonds: &HashSet<usize>,
) -> i32 {
    let neighbor_atom = &mol.atoms[neighbor];
    let in_ring = ring_bonds.contains(&bond_idx);
    let is_carbon = neighbor_atom.symbol == "C" || neighbor_atom.symbol == "c";
    let is_visible = !is_carbon || neighbor_atom.abbrev != "" || neighbor_atom.charge != 0;
    let is_terminal = mol.adj[neighbor].len() == 1;

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

fn ring_bond_set(mol: &MoleculeGraph, rings: &[Vec<usize>]) -> HashSet<usize> {
    let mut ring_bonds = HashSet::new();
    for ring in rings {
        for i in 0..ring.len() {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            if let Some(bond_idx) = bond_between(mol, a, b) {
                ring_bonds.insert(bond_idx);
            }
        }
    }
    ring_bonds
}

// ── Ring detection ────────────────────────────────────────────────────────────

/// Returns cycles as lists of atom indices (the ring path).
fn find_rings(mol: &MoleculeGraph) -> Vec<Vec<usize>> {
    let target_count = cycle_rank(mol);
    if target_count == 0 {
        return Vec::new();
    }

    let bit_words = (mol.bonds.len() + 63) / 64;
    let mut seen_cycles: HashSet<Vec<u64>> = HashSet::new();
    let mut candidates: Vec<(Vec<usize>, Vec<u64>)> = Vec::new();

    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        if let Some(ring) = shortest_path_excluding_bond(mol, bond.from, bond.to, bond_idx) {
            if ring.len() < 3 {
                continue;
            }
            let bits = ring_bond_bits(mol, &ring, bit_words);
            if seen_cycles.insert(bits.clone()) {
                candidates.push((ring, bits));
            }
        }
    }

    candidates.sort_by(|(ring_a, bits_a), (ring_b, bits_b)| {
        ring_a
            .len()
            .cmp(&ring_b.len())
            .then_with(|| bits_a.cmp(bits_b))
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

fn cycle_rank(mol: &MoleculeGraph) -> usize {
    if mol.n_atoms() == 0 {
        return 0;
    }

    let mut seen = vec![false; mol.n_atoms()];
    let mut components = 0;
    for start in 0..mol.n_atoms() {
        if seen[start] {
            continue;
        }
        components += 1;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(atom) = stack.pop() {
            for &(neighbor, _) in &mol.adj[atom] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
    }

    mol.bonds.len() + components - mol.n_atoms()
}

/// BFS from `from` to `to`, intentionally skipping one bond so the path plus
/// that skipped bond is a ring candidate.
fn shortest_path_excluding_bond(
    mol: &MoleculeGraph,
    from: usize,
    to: usize,
    excluded_bond: usize,
) -> Option<Vec<usize>> {
    let n = mol.n_atoms();
    let mut bfs_parent: Vec<Option<usize>> = vec![None; n];

    bfs_parent[from] = Some(from); // root sentinel
    let mut queue = VecDeque::new();
    queue.push_back(from);

    while let Some(u) = queue.pop_front() {
        for &(v, bond_idx) in &mol.adj[u] {
            if bond_idx == excluded_bond {
                continue;
            }
            if bfs_parent[v].is_some() {
                continue;
            }
            bfs_parent[v] = Some(u);

            if v == to {
                let mut path = Vec::new();
                let mut cur = to;
                loop {
                    path.push(cur);
                    if cur == from {
                        break;
                    }
                    cur = match bfs_parent[cur] {
                        Some(p) => p,
                        None => return None,
                    };
                }
                path.reverse();
                return Some(path);
            }
            queue.push_back(v);
        }
    }

    None
}

fn ring_bond_bits(mol: &MoleculeGraph, ring: &[usize], bit_words: usize) -> Vec<u64> {
    let mut bits = vec![0_u64; bit_words];
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if let Some(bond_idx) = bond_between(mol, a, b) {
            bits[bond_idx / 64] |= 1_u64 << (bond_idx % 64);
        }
    }
    bits
}

fn bond_between(mol: &MoleculeGraph, a: usize, b: usize) -> Option<usize> {
    mol.adj[a]
        .iter()
        .find_map(|&(neighbor, bond_idx)| if neighbor == b { Some(bond_idx) } else { None })
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
    mol: &MoleculeGraph,
    rings: &[Vec<usize>],
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
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
        coords,
        placed,
    );

    let mut placed_rings: HashSet<usize> = HashSet::new();
    placed_rings.insert(initial_ring);
    place_connected_fused_rings(mol, rings, &mut placed_rings, coords, placed);
}

fn initial_ring_index(rings: &[Vec<usize>]) -> usize {
    let mut best_idx = 0;
    let mut best_score = (0_usize, 0_usize);
    for (idx, ring) in rings.iter().enumerate() {
        let fused_neighbors = rings
            .iter()
            .enumerate()
            .filter(|(other_idx, other)| *other_idx != idx && rings_share_edge(ring, other))
            .count();
        let score = (fused_neighbors, ring.len());
        if score > best_score {
            best_idx = idx;
            best_score = score;
        }
    }
    best_idx
}

fn rings_share_edge(a: &[usize], b: &[usize]) -> bool {
    shared_atoms(a, b)
        .windows(2)
        .any(|pair| ring_has_edge(a, pair[0], pair[1]) && ring_has_edge(b, pair[0], pair[1]))
}

fn shared_atoms(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut atoms: Vec<usize> = a.iter().copied().filter(|atom| b.contains(atom)).collect();
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

fn ring_has_edge(ring: &[usize], a: usize, b: usize) -> bool {
    for i in 0..ring.len() {
        let x = ring[i];
        let y = ring[(i + 1) % ring.len()];
        if (x == a && y == b) || (x == b && y == a) {
            return true;
        }
    }
    false
}

fn place_connected_fused_rings(
    mol: &MoleculeGraph,
    rings: &[Vec<usize>],
    placed_rings: &mut HashSet<usize>,
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    loop {
        let mut progressed = false;
        for (ring_idx, ring) in rings.iter().enumerate() {
            if placed_rings.contains(&ring_idx) {
                continue;
            }

            if placed_shared_edge(ring, placed).is_some() {
                place_fused_ring(mol, ring, coords, placed);
                placed_rings.insert(ring_idx);
                progressed = true;
            }
        }

        if progressed {
            continue;
        }

        // No edge-fused ring is ready: try a spiro ring, which joins the placed
        // structure at a single shared atom. Placing one may expose further
        // edge fusions, so re-enter the loop afterwards.
        let spiro = rings.iter().enumerate().find(|(idx, ring)| {
            !placed_rings.contains(idx)
                && ring.iter().filter(|&&a| placed[a]).count() == 1
        });
        if let Some((ring_idx, ring)) = spiro {
            place_spiro_ring(mol, ring, coords, placed);
            placed_rings.insert(ring_idx);
            continue;
        }

        break;
    }
}

/// Place a ring joined to the already-placed structure at a single shared
/// (spiro) atom, as a regular polygon opening away from that atom's placed
/// neighbors so the two rings do not overlap.
fn place_spiro_ring(
    mol: &MoleculeGraph,
    ring: &[usize],
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    let n = ring.len();
    let Some(si) = (0..n).find(|&i| placed[ring[i]]) else {
        place_regular_ring(ring, Vec2::new(0.0, 0.0), 0.0, coords, placed);
        return;
    };
    let s = ring[si];
    let p = coords[s];

    // Mean unit vector from the spiro atom toward its already-placed neighbors
    // points into the existing ring system; the new ring is centered opposite.
    let mut dx = 0.0;
    let mut dy = 0.0;
    for bond in &mol.bonds {
        let other = if bond.from == s {
            bond.to
        } else if bond.to == s {
            bond.from
        } else {
            continue;
        };
        if placed[other] {
            let vx = coords[other].x - p.x;
            let vy = coords[other].y - p.y;
            let len = (vx * vx + vy * vy).sqrt();
            if len > 1e-6 {
                dx += vx / len;
                dy += vy / len;
            }
        }
    }
    let into_len = (dx * dx + dy * dy).sqrt();
    let (ix, iy) = if into_len < 1e-6 {
        (0.0, -1.0)
    } else {
        (dx / into_len, dy / into_len)
    };

    let r = 1.0 / (2.0 * (PI / n as f64).sin());
    let center = Vec2::new(p.x - ix * r, p.y - iy * r);
    let start = (p.y - center.y).atan2(p.x - center.x);
    let step = 2.0 * PI / n as f64;
    for (i, &atom) in ring.iter().enumerate() {
        if !placed[atom] {
            let offset = i as i64 - si as i64;
            let angle = start + step * offset as f64;
            coords[atom] = Vec2::new(center.x + r * angle.cos(), center.y + r * angle.sin());
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
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    let n = ring.len() as f64;
    // For a regular n-gon with unit bond length, the circumradius is:
    //   R = 1 / (2 * sin(π/n))
    let r = 1.0 / (2.0 * (PI / n).sin());

    for (i, &atom) in ring.iter().enumerate() {
        if !placed[atom] {
            let angle = start_angle + (2.0 * PI / n) * i as f64;
            coords[atom] = Vec2::new(center.x + r * angle.cos(), center.y + r * angle.sin());
            placed[atom] = true;
        }
    }
}

/// Place a fused ring where two atoms are already positioned.
fn place_fused_ring(
    _mol: &MoleculeGraph,
    ring: &[usize],
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    let n = ring.len();

    // Find the first pair of consecutive already-placed atoms (the shared edge)
    let Some((ei, ej)) = placed_shared_edge(ring, placed) else {
        // Fallback: place as standalone (shouldn't happen)
        place_regular_ring(ring, Vec2::new(0.0, 0.0), 0.0, coords, placed);
        return;
    };

    // The shared bond goes from ring[ei] to ring[ej].
    // Place remaining atoms of the n-gon on the opposite side from any already-placed atoms.
    let p1 = coords[ring[ei]];
    let p2 = coords[ring[ej]];

    // Midpoint and perpendicular direction
    let mid = Vec2::new((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
    let bond_angle = (p2.y - p1.y).atan2(p2.x - p1.x);

    // For a regular n-gon, the center is at distance d from the midpoint of each edge:
    //   d = (1/2) * cot(π/n)  [in units of bond length]
    let fn_ = n as f64;
    let center_dist = 0.5 / (PI / fn_).tan();

    // Determine which side to place the center (away from existing atoms)
    let perp1 = Vec2::new(
        mid.x + center_dist * (bond_angle + PI / 2.0).cos(),
        mid.y + center_dist * (bond_angle + PI / 2.0).sin(),
    );
    let perp2 = Vec2::new(
        mid.x + center_dist * (bond_angle - PI / 2.0).cos(),
        mid.y + center_dist * (bond_angle - PI / 2.0).sin(),
    );

    // Pick the side farther from existing placed atoms (avoid overlap)
    let center = pick_farther_center(perp1, perp2, placed, coords);

    // The angle from center to ring[ei]
    let start_angle = (p1.y - center.y).atan2(p1.x - center.x);
    let next_angle = (p2.y - center.y).atan2(p2.x - center.x);
    let step = 2.0 * PI / fn_;
    let sweep = if angle_delta(start_angle + step, next_angle).abs()
        <= angle_delta(start_angle - step, next_angle).abs()
    {
        1.0
    } else {
        -1.0
    };

    let r = 1.0 / (2.0 * (PI / fn_).sin());
    for (i, &atom) in ring.iter().enumerate() {
        if !placed[atom] {
            // Compute offset from ei position
            let offset = i as i64 - ei as i64;
            let angle = start_angle + sweep * step * offset as f64;
            coords[atom] = Vec2::new(center.x + r * angle.cos(), center.y + r * angle.sin());
            placed[atom] = true;
        }
    }
}

fn angle_delta(a: f64, b: f64) -> f64 {
    let mut delta = a - b;
    while delta > PI {
        delta -= 2.0 * PI;
    }
    while delta < -PI {
        delta += 2.0 * PI;
    }
    delta
}

fn placed_shared_edge(ring: &[usize], placed: &[bool]) -> Option<(usize, usize)> {
    for i in 0..ring.len() {
        let j = (i + 1) % ring.len();
        if placed[ring[i]] && placed[ring[j]] {
            return Some((i, j));
        }
    }
    None
}

fn pick_farther_center(
    c1: Vec2,
    c2: Vec2,
    placed: &[bool],
    coords: &[Vec2],
) -> Vec2 {
    let mut sum1 = 0.0_f64;
    let mut sum2 = 0.0_f64;
    let mut count = 0;
    for (atom_idx, was_placed) in placed.iter().enumerate() {
        if *was_placed {
            sum1 += c1.dist(coords[atom_idx]);
            sum2 += c2.dist(coords[atom_idx]);
            count += 1;
        }
    }
    if count == 0 || sum1 >= sum2 {
        c1
    } else {
        c2
    }
}

fn place_ring_from_anchor(
    ring: &[usize],
    anchor: usize,
    center_dir: f64,
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    let Some(anchor_pos) = ring.iter().position(|&a| a == anchor) else {
        return;
    };

    let n = ring.len() as f64;
    let r = 1.0 / (2.0 * (PI / n).sin());
    let center = Vec2::new(
        coords[anchor].x + r * center_dir.cos(),
        coords[anchor].y + r * center_dir.sin(),
    );
    let step = 2.0 * PI / n;
    let anchor_angle = (coords[anchor].y - center.y).atan2(coords[anchor].x - center.x);
    let start_angle = anchor_angle - step * anchor_pos as f64;

    place_regular_ring(ring, center, start_angle, coords, placed);
}

fn place_ring_system_from_anchor(
    mol: &MoleculeGraph,
    rings: &[Vec<usize>],
    ring_idx: usize,
    anchor: usize,
    center_dir: f64,
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    place_ring_from_anchor(&rings[ring_idx], anchor, center_dir, coords, placed);

    let mut placed_rings: HashSet<usize> = HashSet::new();
    placed_rings.insert(ring_idx);
    place_connected_fused_rings(mol, rings, &mut placed_rings, coords, placed);
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
    mol: &MoleculeGraph,
    rings: &[Vec<usize>],
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    for ring in rings {
        if !ring.iter().all(|&a| placed[a]) {
            continue;
        }

        let n = ring.len() as f64;
        let cx: f64 = ring.iter().map(|&a| coords[a].x).sum::<f64>() / n;
        let cy: f64 = ring.iter().map(|&a| coords[a].y).sum::<f64>() / n;
        for &atom in ring {
            let outward = (coords[atom].y - cy).atan2(coords[atom].x - cx);
            place_ring_substituents(mol, atom, outward, rings, coords, placed);
        }
    }
}

// ── Ring substituent placement ────────────────────────────────────────────────

/// Place all unplaced neighbors of a ring atom radially outward from the ring center,
/// then recursively extend any chains from those substituents.
fn place_ring_substituents(
    mol: &MoleculeGraph,
    ring_atom: usize,
    outward_dir: f64,
    rings: &[Vec<usize>],
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    let unplaced: Vec<usize> = mol.adj[ring_atom]
        .iter()
        .filter(|&&(v, _)| !placed[v])
        .map(|&(v, _)| v)
        .collect();
    let dirs = free_substituent_directions(
        mol,
        ring_atom,
        unplaced.len(),
        outward_dir,
        coords,
        placed,
    );

    for (i, &v) in unplaced.iter().enumerate() {
        let dir = dirs[i];
        coords[v] = Vec2::new(
            coords[ring_atom].x + dir.cos(),
            coords[ring_atom].y + dir.sin(),
        );
        placed[v] = true;
        if let Some(ring_idx) = unfinished_ring_containing_atom(rings, v, placed) {
            place_ring_system_from_anchor(mol, rings, ring_idx, v, dir, coords, placed);
            place_substituents_for_placed_rings(mol, rings, coords, placed);
        } else {
            place_chain(mol, v, dir, false, rings, coords, placed);
        }
    }
}

fn free_substituent_directions(
    mol: &MoleculeGraph,
    atom_idx: usize,
    count: usize,
    fallback_dir: f64,
    coords: &[Vec2],
    placed: &[bool],
) -> Vec<f64> {
    if count == 0 {
        return Vec::new();
    }

    let mut occupied: Vec<f64> = mol.adj[atom_idx]
        .iter()
        .filter(|&&(neighbor, _)| placed[neighbor])
        .map(|&(neighbor, _)| {
            (coords[neighbor].y - coords[atom_idx].y)
                .atan2(coords[neighbor].x - coords[atom_idx].x)
        })
        .collect();

    if occupied.len() < 2 {
        return spread_around_direction(fallback_dir, count, PI / 3.0);
    }

    occupied.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut best_start = occupied[0];
    let mut best_gap = 0.0;
    for i in 0..occupied.len() {
        let start = occupied[i];
        let end = if i + 1 < occupied.len() {
            occupied[i + 1]
        } else {
            occupied[0] + 2.0 * PI
        };
        let gap = end - start;
        if gap > best_gap {
            best_gap = gap;
            best_start = start;
        }
    }

    let usable_gap = (best_gap - PI / 6.0).max(PI / 6.0);
    if count == 1 {
        vec![normalize_angle(best_start + best_gap / 2.0)]
    } else {
        (0..count)
            .map(|i| {
                let t = (i + 1) as f64 / (count + 1) as f64;
                normalize_angle(best_start + (best_gap - usable_gap) / 2.0 + usable_gap * t)
            })
            .collect()
    }
}

fn spread_around_direction(center: f64, count: usize, spread: f64) -> Vec<f64> {
    if count == 1 {
        return vec![center];
    }
    (0..count)
        .map(|i| {
            normalize_angle(
                center - spread / 2.0 + spread * i as f64 / (count.saturating_sub(1)) as f64,
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
fn acyclic_root(mol: &MoleculeGraph) -> (usize, bool) {
    let n = mol.n_atoms();
    if n <= 1 || !is_connected(mol) {
        return (0, false);
    }

    let mut degree: Vec<usize> = (0..n).map(|i| mol.adj[i].len()).collect();
    let mut removed = vec![false; n];
    let mut remaining = n;
    let mut leaves: Vec<usize> = (0..n).filter(|&i| degree[i] <= 1).collect();

    while remaining > 2 {
        let mut next = Vec::new();
        for &leaf in &leaves {
            if removed[leaf] {
                continue;
            }
            removed[leaf] = true;
            remaining -= 1;
            for &(neighbor, _) in &mol.adj[leaf] {
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
    let center = (0..n)
        .filter(|&i| !removed[i])
        .max_by(|&a, &b| mol.adj[a].len().cmp(&mol.adj[b].len()).then(b.cmp(&a)))
        .unwrap_or(0);

    if is_symmetric_hub(mol, center) {
        (center, true)
    } else {
        (0, false)
    }
}

/// True when `hub` has at least three arms and every arm leads to a subtree of
/// the same size. Such hubs are the ones whose evenly-spread arms would
/// otherwise be drawn as a rotational pinwheel.
fn is_symmetric_hub(mol: &MoleculeGraph, hub: usize) -> bool {
    if mol.adj[hub].len() < 3 {
        return false;
    }
    let sizes: Vec<usize> = mol.adj[hub]
        .iter()
        .map(|&(neighbor, _)| arm_subtree_size(mol, neighbor, hub))
        .collect();
    sizes.iter().all(|&s| s == sizes[0])
}

fn arm_subtree_size(mol: &MoleculeGraph, start: usize, hub: usize) -> usize {
    let mut seen = vec![false; mol.n_atoms()];
    seen[hub] = true;
    seen[start] = true;
    let mut stack = vec![start];
    let mut count = 0;
    while let Some(atom) = stack.pop() {
        count += 1;
        for &(neighbor, _) in &mol.adj[atom] {
            if !seen[neighbor] {
                seen[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }
    count
}

fn is_connected(mol: &MoleculeGraph) -> bool {
    let n = mol.n_atoms();
    if n == 0 {
        return true;
    }
    let mut seen = vec![false; n];
    let mut stack = vec![0];
    seen[0] = true;
    let mut count = 1;
    while let Some(atom) = stack.pop() {
        for &(neighbor, _) in &mol.adj[atom] {
            if !seen[neighbor] {
                seen[neighbor] = true;
                count += 1;
                stack.push(neighbor);
            }
        }
    }
    count == n
}

/// DFS traversal that extends unplaced atoms at ±120° from the incoming direction.
fn place_chain(
    mol: &MoleculeGraph,
    start: usize,
    incoming_angle: f64,
    symmetric_hub_root: bool,
    rings: &[Vec<usize>],
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    // Iterative DFS to avoid stack overflow on long chains
    // Stack entries: (atom_idx, incoming_bond_idx, incoming_direction_radians, turn_sign)
    // turn_sign alternates to produce a zigzag: +1 = turn left, -1 = turn right
    let mut stack: Vec<(usize, Option<usize>, f64, f64)> = Vec::new();
    stack.push((start, None, incoming_angle, 1.0));

    while let Some((u, incoming_bond_idx, dir, sign)) = stack.pop() {
        let mut unplaced_neighbors: Vec<(usize, usize, usize)> = mol.adj[u]
            .iter()
            .enumerate()
            .filter(|(_, &(v, _))| !placed[v])
            .map(|(order, &(v, bond_idx))| (v, bond_idx, order))
            .collect();

        unplaced_neighbors.sort_by(|a, b| {
            let size_a = unplaced_subtree_size(mol, a.0, u, placed);
            let size_b = unplaced_subtree_size(mol, b.0, u, placed);
            size_b.cmp(&size_a).then_with(|| a.2.cmp(&b.2))
        });

        let dirs = square_planar_directions(mol, u, &unplaced_neighbors, coords, placed)
            .unwrap_or_else(|| {
                chain_neighbor_directions(
                    mol,
                    u,
                    dir,
                    sign,
                    unplaced_neighbors.len(),
                    incoming_bond_idx.is_some(),
                    coords,
                    placed,
                )
            });

        for (i, &(v, bond_idx, _)) in unplaced_neighbors.iter().enumerate() {
            // Independent branches can converge on the same lattice point (e.g.
            // ortho ring substituents growing toward each other). If the
            // proposed spot is already crowded, deflect toward a clear one.
            let new_dir = if is_linear_atom(mol, u) {
                dirs[i]
            } else {
                resolve_chain_direction(mol, u, dirs[i], dir, coords, placed)
            };

            coords[v] = Vec2::new(coords[u].x + new_dir.cos(), coords[u].y + new_dir.sin());
            placed[v] = true;

            if let Some(ring_idx) = unfinished_ring_containing_atom(rings, v, placed) {
                place_ring_system_from_anchor(mol, rings, ring_idx, v, new_dir, coords, placed);
                place_substituents_for_placed_rings(mol, rings, coords, placed);
            } else {
                // Arms leaving a symmetric hub bend away from the vertical axis
                // (sign keyed to the horizontal lean), making the equal arms
                // mirror images instead of a rotational pinwheel. Every other
                // step alternates the turn sign for the usual zigzag.
                let next_sign = if symmetric_hub_root && incoming_bond_idx.is_none() {
                    if new_dir.cos() >= 0.0 {
                        1.0
                    } else {
                        -1.0
                    }
                } else if (normalize_angle(new_dir - dirs[i])).abs() > 1e-9 {
                    // A deflected bond changed which way this step turned, so
                    // the child's zigzag continues from the actual turn taken.
                    let turn = normalize_angle(new_dir - dir);
                    if turn > 1e-9 {
                        -1.0
                    } else if turn < -1e-9 {
                        1.0
                    } else {
                        -sign
                    }
                } else {
                    -sign
                };
                stack.push((v, Some(bond_idx), new_dir, next_sign));
            }
        }
    }
}

/// Minimum clearance between a newly placed chain atom and any atom already
/// placed elsewhere, in bond-length units.
const CHAIN_CLEARANCE: f64 = 0.55;

/// Distance from `p` to the nearest already-placed atom.
fn nearest_placed_distance(p: Vec2, coords: &[Vec2], placed: &[bool]) -> f64 {
    let mut min = f64::INFINITY;
    for i in 0..coords.len() {
        if placed[i] {
            min = min.min(p.dist(coords[i]));
        }
    }
    min
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
    mol: &MoleculeGraph,
    u: usize,
    proposed: f64,
    incoming: f64,
    coords: &mut Vec<Vec2>,
    placed: &[bool],
) -> f64 {
    let origin = coords[u];
    let end_of = move |d: f64| Vec2::new(origin.x + d.cos(), origin.y + d.sin());
    if nearest_placed_distance(end_of(proposed), coords, placed) >= CHAIN_CLEARANCE {
        return proposed;
    }
    if mirror_blocking_branch(mol, u, end_of(proposed), coords, placed) {
        return proposed;
    }

    // Flip the zigzag turn: the proposal mirrored across the incoming
    // direction is the other ideal slot at this atom.
    let flipped = normalize_angle(2.0 * incoming - proposed);
    if normalize_angle(flipped - proposed).abs() > 1e-9
        && nearest_placed_distance(end_of(flipped), coords, placed) >= CHAIN_CLEARANCE
    {
        let taken = mol.adj[u].iter().any(|&(v2, _)| {
            placed[v2] && {
                let na = (coords[v2].y - coords[u].y).atan2(coords[v2].x - coords[u].x);
                normalize_angle(flipped - na).abs() < PI / 6.0
            }
        });
        if !taken {
            return flipped;
        }
    }
    proposed
}

/// Tries to clear the crowded spot `target` by reflecting the branch that
/// occupies it across its attachment bond (e.g. flipping an ortho substituent
/// to lean the other way). Candidate branches are subtrees hanging off a
/// single bond that contain every blocking atom, are fully placed, and do not
/// contain `u`. A candidate is accepted only if, after reflection, both the
/// target spot and every moved atom have full clearance; the smallest such
/// branch is flipped. Returns whether a reflection was applied.
fn mirror_blocking_branch(
    mol: &MoleculeGraph,
    u: usize,
    target: Vec2,
    coords: &mut Vec<Vec2>,
    placed: &[bool],
) -> bool {
    let n = mol.n_atoms();
    let blockers: Vec<usize> = (0..n)
        .filter(|&i| placed[i] && i != u && target.dist(coords[i]) < CHAIN_CLEARANCE)
        .collect();
    if blockers.is_empty() {
        return false;
    }

    let mut best: Option<(Vec<usize>, usize, usize)> = None;
    for bond in &mol.bonds {
        for (a, b) in [(bond.from, bond.to), (bond.to, bond.from)] {
            if !placed[a] || !placed[b] {
                continue;
            }
            let subtree = collect_subtree(mol, b, a, a);
            if subtree.contains(&u) {
                continue;
            }
            if !subtree.iter().all(|&t| placed[t]) {
                continue;
            }
            if !blockers.iter().all(|bl| subtree.contains(bl)) {
                continue;
            }
            if best.as_ref().is_some_and(|(prev, _, _)| prev.len() <= subtree.len()) {
                continue;
            }

            let mut in_subtree = vec![false; n];
            for &t in &subtree {
                in_subtree[t] = true;
            }
            let reflect = |p: Vec2| reflect_point_across_line(p, coords[a], coords[b]);

            // The target spot must actually clear once the branch flips.
            let target_clear = (0..n).filter(|&i| placed[i] && i != u).all(|i| {
                let p = if in_subtree[i] { reflect(coords[i]) } else { coords[i] };
                target.dist(p) >= CHAIN_CLEARANCE
            });
            // The flipped branch must land in open space: full clearance from
            // every atom outside it (internal distances are preserved).
            let branch_clear = target_clear
                && subtree.iter().all(|&t| {
                    let tp = reflect(coords[t]);
                    (0..n)
                        .filter(|&i| placed[i] && !in_subtree[i] && i != a)
                        .all(|i| tp.dist(coords[i]) >= CHAIN_CLEARANCE)
                });
            if branch_clear {
                best = Some((subtree, a, b));
            }
        }
    }

    let Some((subtree, a, b)) = best else {
        return false;
    };
    let (axis_a, axis_b) = (coords[a], coords[b]);
    for &t in &subtree {
        coords[t] = reflect_point_across_line(coords[t], axis_a, axis_b);
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
    mol: &MoleculeGraph,
    center: usize,
    unplaced: &[(usize, usize, usize)],
    coords: &[Vec2],
    placed: &[bool],
) -> Option<Vec<f64>> {
    let AtomChirality::SquarePlanar(class) = mol.atoms[center].chirality else {
        return None;
    };
    let neighbor_bonds = &mol.neighbor_bonds[center];
    if neighbor_bonds.len() != 4 {
        return None;
    }

    // Neighbor atoms in SMILES writing order.
    let writing: Vec<usize> = neighbor_bonds
        .iter()
        .map(|&b| {
            let bond = &mol.bonds[b];
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
            let angle = (coords[a].y - coords[center].y).atan2(coords[a].x - coords[center].x);
            angle - corner_of[w_idx] as f64 * (PI / 2.0)
        })
        .unwrap_or(PI);

    let dirs = unplaced
        .iter()
        .map(|&(atom, _, _)| {
            let w_idx = writing.iter().position(|&a| a == atom)?;
            Some(normalize_angle(base + corner_of[w_idx] as f64 * (PI / 2.0)))
        })
        .collect::<Option<Vec<f64>>>()?;
    Some(dirs)
}

/// One or two substituents use the ±60° zigzag; three or more are spread across
/// the angular space left free by already-placed neighbors.
fn chain_neighbor_directions(
    mol: &MoleculeGraph,
    u: usize,
    dir: f64,
    sign: f64,
    count: usize,
    has_incoming: bool,
    coords: &[Vec2],
    placed: &[bool],
) -> Vec<f64> {
    if count == 0 {
        return Vec::new();
    }

    // A linear atom (cumulene/alkyne center) keeps the incoming direction so the
    // chain stays straight through it.
    if is_linear_atom(mol, u) && has_incoming {
        return vec![dir; count];
    }

    match count {
        1 => vec![dir + sign * (PI / 3.0)],
        2 => vec![dir + sign * (PI / 3.0), dir - sign * (PI / 3.0)],
        _ => {
            let occupied: Vec<f64> = mol.adj[u]
                .iter()
                .filter(|&&(neighbor, _)| placed[neighbor])
                .map(|&(neighbor, _)| {
                    (coords[neighbor].y - coords[u].y).atan2(coords[neighbor].x - coords[u].x)
                })
                .collect();
            let slots = distribute_in_free_space(&occupied, count, dir);
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
    let n = slots.len();
    let center = (n as f64 - 1.0) / 2.0;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let da = (a as f64 - center).abs();
        let db = (b as f64 - center).abs();
        da.partial_cmp(&db)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });

    let mut assigned = vec![0.0; n];
    for (neighbor_rank, &slot_idx) in order.iter().enumerate() {
        assigned[neighbor_rank] = slots[slot_idx];
    }
    assigned
}

/// Spreads `count` directions across the largest angular gap left by `occupied`
/// (the directions of already-placed neighbors). With nothing placed, the
/// directions are spaced evenly around the full circle starting at `fallback`.
fn distribute_in_free_space(occupied: &[f64], count: usize, fallback: f64) -> Vec<f64> {
    if count == 0 {
        return Vec::new();
    }

    if occupied.is_empty() {
        return (0..count)
            .map(|i| normalize_angle(fallback + 2.0 * PI * i as f64 / count as f64))
            .collect();
    }

    let mut occ = occupied.to_vec();
    occ.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut best_start = occ[0];
    let mut best_gap = 0.0;
    for i in 0..occ.len() {
        let start = occ[i];
        let end = if i + 1 < occ.len() {
            occ[i + 1]
        } else {
            occ[0] + 2.0 * PI
        };
        let gap = end - start;
        if gap > best_gap {
            best_gap = gap;
            best_start = start;
        }
    }

    // Split the gap into `count + 1` segments so bonds keep a margin from the edges.
    let segments = (count + 1) as f64;
    (1..=count)
        .map(|i| normalize_angle(best_start + best_gap * i as f64 / segments))
        .collect()
}

fn unplaced_subtree_size(
    mol: &MoleculeGraph,
    root: usize,
    parent: usize,
    placed: &[bool],
) -> usize {
    let mut size = 0;
    let mut seen = vec![false; mol.n_atoms()];
    seen[parent] = true;
    let mut stack = vec![root];
    seen[root] = true;

    while let Some(atom) = stack.pop() {
        if placed[atom] {
            continue;
        }
        size += 1;
        for &(neighbor, _) in &mol.adj[atom] {
            if !seen[neighbor] {
                seen[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }

    size
}

fn is_linear_atom(mol: &MoleculeGraph, atom_idx: usize) -> bool {
    if mol.adj[atom_idx].len() != 2 {
        return false;
    }

    let mut double_count = 0;
    let mut has_triple = false;
    for &(_, bond_idx) in &mol.adj[atom_idx] {
        match mol.bonds[bond_idx].order {
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
fn h_label_angle(existing_adj: &[f64]) -> f64 {
    if existing_adj.is_empty() {
        return 0.0;
    }

    let mut occ = existing_adj.to_vec();
    occ.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut best_start = occ[0];
    let mut best_gap = 0.0_f64;
    for i in 0..occ.len() {
        let start = occ[i];
        let end = if i + 1 < occ.len() {
            occ[i + 1]
        } else {
            occ[0] + 2.0 * PI
        };
        let gap = end - start;
        if gap > best_gap {
            best_gap = gap;
            best_start = start;
        }
    }

    normalize_angle(best_start + best_gap / 2.0)
}

// ── Utilities ─────────────────────────────────────────────────────────────────

fn center_coords(coords: &mut Vec<Vec2>) {
    if coords.is_empty() {
        return;
    }
    let cx = coords.iter().map(|p| p.x).sum::<f64>() / coords.len() as f64;
    let cy = coords.iter().map(|p| p.y).sum::<f64>() / coords.len() as f64;
    for p in coords.iter_mut() {
        p.x -= cx;
        p.y -= cy;
    }
}

/// For each bond, returns the unit vector pointing from the bond midpoint
/// toward the centroid of the smallest ring containing that bond.
/// Returns (0.0, 0.0) for bonds not in any ring.
fn ring_inner_directions(
    mol: &MoleculeGraph,
    rings: &[Vec<usize>],
    coords: &[Vec2],
) -> Vec<(f64, f64)> {
    let mut dirs = vec![(0.0_f64, 0.0_f64); mol.bonds.len()];

    for (b_idx, bond) in mol.bonds.iter().enumerate() {
        let best_ring = best_ring_for_inner_bond(mol, rings, bond.from, bond.to);

        if let Some(ring) = best_ring {
            let n = ring.len() as f64;
            let cx = ring.iter().map(|&a| coords[a].x).sum::<f64>() / n;
            let cy = ring.iter().map(|&a| coords[a].y).sum::<f64>() / n;

            let mx = (coords[bond.from].x + coords[bond.to].x) / 2.0;
            let my = (coords[bond.from].y + coords[bond.to].y) / 2.0;

            let vx = cx - mx;
            let vy = cy - my;
            let vlen = (vx * vx + vy * vy).sqrt();
            if vlen > 1e-6 {
                dirs[b_idx] = (vx / vlen, vy / vlen);
            }
        }
    }

    dirs
}

fn best_ring_for_inner_bond<'a>(
    mol: &MoleculeGraph,
    rings: &'a [Vec<usize>],
    from: usize,
    to: usize,
) -> Option<&'a Vec<usize>> {
    rings
        .iter()
        .filter(|ring| ring_has_edge(ring, from, to))
        .max_by(|a, b| {
            ring_unsaturation_score(mol, a)
                .cmp(&ring_unsaturation_score(mol, b))
                .then_with(|| b.len().cmp(&a.len()))
        })
}

fn ring_unsaturation_score(mol: &MoleculeGraph, ring: &[usize]) -> usize {
    (0..ring.len())
        .filter_map(|i| bond_between(mol, ring[i], ring[(i + 1) % ring.len()]))
        .filter(|&bond_idx| {
            matches!(
                mol.bonds[bond_idx].order,
                BondOrder::Double | BondOrder::Aromatic
            )
        })
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::MoleculeGraph;

    #[test]
    fn steroid_ring_system_detects_four_rings() {
        let mol = MoleculeGraph::from_smiles(
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
            Vec::new(),
            Vec::new(),
        )
        .expect("steroid-like molecule should parse");
        let rings = find_rings(&mol);
        assert_eq!(rings.len(), 4, "rings: {rings:?}");
    }

    #[test]
    fn ring_closure_after_branch_stays_on_branch_point() {
        let pre = crate::preprocess_smiles("C1=CCCC(=O)1").expect("preprocess failed");
        let mol = MoleculeGraph::from_smiles(
            &pre.smiles,
            pre.forced_direction_markers,
            pre.aromatic_atom_markers,
        )
        .expect("cyclopentenone should parse");
        let rings = find_rings(&mol);
        assert!(rings.iter().any(|ring| ring.len() == 5), "rings: {rings:?}");

        let complex = "O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5";
        let pre = crate::preprocess_smiles(complex).expect("preprocess failed");
        let mol = MoleculeGraph::from_smiles(
            &pre.smiles,
            pre.forced_direction_markers,
            pre.aromatic_atom_markers,
        )
        .expect("complex fused system should parse");
        let rings = find_rings(&mol);
        assert!(
            rings
                .iter()
                .any(|ring| ring.len() == 5 && ring.contains(&19)),
            "rings: {rings:?}"
        );
    }

    #[test]
    fn fused_double_bond_prefers_unsaturated_ring_side() {
        let mol = MoleculeGraph::from_smiles(
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
            Vec::new(),
            Vec::new(),
        )
        .expect("steroid-like molecule should parse");
        let rings = find_rings(&mol);
        let mut checked_shared_double = false;

        for bond in mol.bonds.iter().filter(|bond| bond.order == BondOrder::Double) {
            let containing: Vec<&Vec<usize>> = rings
                .iter()
                .filter(|ring| ring_has_edge(ring, bond.from, bond.to))
                .collect();
            if containing.len() < 2 {
                continue;
            }

            let selected =
                best_ring_for_inner_bond(&mol, &rings, bond.from, bond.to).expect("ring missing");
            let selected_score = ring_unsaturation_score(&mol, selected);
            let max_score = containing
                .iter()
                .map(|ring| ring_unsaturation_score(&mol, ring))
                .max()
                .unwrap();
            assert_eq!(selected_score, max_score);
            checked_shared_double = true;
        }

        assert!(checked_shared_double, "expected a fused/shared double bond");
    }
}

fn bounding_box(coords: &[Vec2]) -> (f64, f64) {
    if coords.is_empty() {
        return (0.0, 0.0);
    }
    let min_x = coords.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let max_x = coords.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let min_y = coords.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let max_y = coords.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    (max_x - min_x + 1.0, max_y - min_y + 1.0)
}
