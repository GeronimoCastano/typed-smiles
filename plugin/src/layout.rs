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

use crate::graph::MoleculeGraph;
use crate::render::{AtomOutput, BondOutput, LayoutOutput, Vec2};

pub fn compute_layout(mol: &MoleculeGraph) -> Result<LayoutOutput, String> {
    if mol.n_atoms() == 0 {
        return Ok(LayoutOutput {
            atoms: vec![],
            bonds: vec![],
            bbox_width: 0.0,
            bbox_height: 0.0,
        });
    }

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

    let root = placed.iter().position(|&p| p).unwrap_or(0);
    if !placed[root] {
        coords[root] = Vec2::new(0.0, 0.0);
        placed[root] = true;
    }

    // Start at -30° so the first bond is horizontal in zigzag depiction
    let initial_dir = -PI / 6.0;
    place_chain(mol, root, initial_dir, &rings, &mut coords, &mut placed);

    // ── 5. Center the molecule ────────────────────────────────────────────

    center_coords(&mut coords);

    // ── 5. Build output ───────────────────────────────────────────────────

    let atoms: Vec<AtomOutput> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| AtomOutput {
            symbol: a.symbol.clone(),
            pos: coords[i],
            hcount: a.hcount,
            implicit_h: implicit_h_count(mol, i),
            charge: a.charge,
            abbrev: a.abbrev.clone(),
        })
        .collect();

    let inner_dirs = ring_inner_directions(mol, &rings, &coords);

    let bonds: Vec<BondOutput> = mol
        .bonds
        .iter()
        .enumerate()
        .map(|(i, b)| BondOutput {
            from: b.from,
            to: b.to,
            order: b.order.as_u8(),
            stereo: b.stereo.as_str().to_string(),
            inner_x: inner_dirs[i].0,
            inner_y: inner_dirs[i].1,
        })
        .collect();

    let (bbox_w, bbox_h) = bounding_box(&coords);

    Ok(LayoutOutput {
        atoms,
        bonds,
        bbox_width: bbox_w,
        bbox_height: bbox_h,
    })
}

fn implicit_h_count(mol: &MoleculeGraph, atom_idx: usize) -> u8 {
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

// ── Ring detection (simple DFS cycle finder) ─────────────────────────────────

/// Returns cycles as lists of atom indices (the ring path).
fn find_rings(mol: &MoleculeGraph) -> Vec<Vec<usize>> {
    let n = mol.n_atoms();
    let mut visited = vec![false; n];
    let mut rings: Vec<Vec<usize>> = Vec::new();

    for start in 0..n {
        if !visited[start] {
            dfs_rings(mol, start, usize::MAX, &mut visited, &mut rings);
        }
    }

    // Deduplicate and keep only SSSR candidates (shortest rings)
    rings.sort_by_key(|r| r.len());
    rings.dedup();
    rings
}

fn dfs_rings(
    mol: &MoleculeGraph,
    u: usize,
    prev: usize,
    visited: &mut Vec<bool>,
    rings: &mut Vec<Vec<usize>>,
) {
    visited[u] = true;
    for &(v, _) in &mol.adj[u] {
        if v == prev {
            continue;
        }
        if visited[v] {
            // Back edge u→v. Find the SMALLEST ring through this edge via BFS
            // on the full graph so that fused bicyclics are detected correctly.
            // (Walking the DFS parent tree instead would give the outer perimeter
            // of fused systems, not each individual ring.)
            let ring = bfs_shortest_ring(mol, u, v);
            if ring.len() >= 3 {
                rings.push(ring);
            }
        } else {
            dfs_rings(mol, v, u, visited, rings);
        }
    }
}

/// BFS from `back_to` to `back_from`, intentionally skipping the direct
/// `back_to → back_from` edge so that the path goes through the ring body.
fn bfs_shortest_ring(mol: &MoleculeGraph, back_from: usize, back_to: usize) -> Vec<usize> {
    let n = mol.n_atoms();
    let mut bfs_parent: Vec<Option<usize>> = vec![None; n];

    bfs_parent[back_to] = Some(back_to); // root sentinel
    let mut queue = VecDeque::new();
    queue.push_back(back_to);

    while let Some(u) = queue.pop_front() {
        for &(v, _) in &mol.adj[u] {
            // Skip the direct back edge so BFS must find the ring path
            if u == back_to && v == back_from {
                continue;
            }
            if bfs_parent[v].is_some() {
                continue;
            }
            bfs_parent[v] = Some(u);

            if v == back_from {
                // Reconstruct ring atoms from back_from back to back_to
                let mut ring = Vec::new();
                let mut cur = back_from;
                loop {
                    ring.push(cur);
                    if cur == back_to {
                        break;
                    }
                    cur = match bfs_parent[cur] {
                        Some(p) => p,
                        None => break,
                    };
                }
                return ring;
            }
            queue.push_back(v);
        }
    }

    vec![] // unreachable for a real ring
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

    // Place the first ring system only. Other standalone ring systems in the
    // same molecule are anchored later when the connecting chain reaches them.
    place_regular_ring(
        &rings[0],
        Vec2::new(0.0, 0.0),
        90.0_f64.to_radians(),
        coords,
        placed,
    );

    let mut placed_rings: HashSet<usize> = HashSet::new();
    placed_rings.insert(0);
    place_connected_fused_rings(mol, rings, &mut placed_rings, coords, placed);
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

            let placed_count = ring.iter().filter(|&&a| placed[a]).count();
            if placed_count >= 2 {
                place_fused_ring(mol, ring, coords, placed);
                placed_rings.insert(ring_idx);
                progressed = true;
            }
        }

        if !progressed {
            break;
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
    let mut shared_edge: Option<(usize, usize)> = None; // (ring_pos_i, ring_pos_j)
    for i in 0..n {
        let j = (i + 1) % n;
        if placed[ring[i]] && placed[ring[j]] {
            shared_edge = Some((i, j));
            break;
        }
    }

    let Some((ei, ej)) = shared_edge else {
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
    let center = pick_farther_center(perp1, perp2, ring, placed, coords);

    // The angle from center to ring[ei]
    let start_angle = (p1.y - center.y).atan2(p1.x - center.x);

    let r = 1.0 / (2.0 * (PI / fn_).sin());
    for (i, &atom) in ring.iter().enumerate() {
        if !placed[atom] {
            // Compute offset from ei position
            let offset = i as i64 - ei as i64;
            let angle = start_angle + (2.0 * PI / fn_) * offset as f64;
            coords[atom] = Vec2::new(center.x + r * angle.cos(), center.y + r * angle.sin());
            placed[atom] = true;
        }
    }
}

fn pick_farther_center(
    c1: Vec2,
    c2: Vec2,
    ring: &[usize],
    placed: &[bool],
    coords: &[Vec2],
) -> Vec2 {
    let mut sum1 = 0.0_f64;
    let mut sum2 = 0.0_f64;
    let mut count = 0;
    for &a in ring {
        if placed[a] {
            sum1 += c1.dist(coords[a]);
            sum2 += c2.dist(coords[a]);
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
    rings.iter().position(|ring| {
        ring.contains(&atom) && ring.iter().any(|&ring_atom| !placed[ring_atom])
    })
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

    for (i, &v) in unplaced.iter().enumerate() {
        let turn = if unplaced.len() == 1 {
            0.0
        } else {
            let spread = PI / 3.0;
            -spread / 2.0 + spread * i as f64 / (unplaced.len() as f64 - 1.0)
        };
        let dir = outward_dir + turn;
        coords[v] = Vec2::new(
            coords[ring_atom].x + dir.cos(),
            coords[ring_atom].y + dir.sin(),
        );
        placed[v] = true;
        if let Some(ring_idx) = unfinished_ring_containing_atom(rings, v, placed) {
            place_ring_system_from_anchor(mol, rings, ring_idx, v, dir, coords, placed);
            place_substituents_for_placed_rings(mol, rings, coords, placed);
        } else {
            place_chain(mol, v, dir, rings, coords, placed);
        }
    }
}

// ── Chain layout via DFS ──────────────────────────────────────────────────────

/// DFS traversal that extends unplaced atoms at ±120° from the incoming direction.
fn place_chain(
    mol: &MoleculeGraph,
    start: usize,
    incoming_angle: f64,
    rings: &[Vec<usize>],
    coords: &mut Vec<Vec2>,
    placed: &mut Vec<bool>,
) {
    // Iterative DFS to avoid stack overflow on long chains
    // Stack entries: (atom_idx, incoming_direction_radians, turn_sign)
    // turn_sign alternates to produce a zigzag: +1 = turn left, -1 = turn right
    let mut stack: Vec<(usize, f64, f64)> = Vec::new();
    stack.push((start, incoming_angle, 1.0));

    while let Some((u, dir, sign)) = stack.pop() {
        let unplaced_neighbors: Vec<usize> = mol.adj[u]
            .iter()
            .filter(|&&(v, _)| !placed[v])
            .map(|&(v, _)| v)
            .collect();

        for (i, &v) in unplaced_neighbors.iter().enumerate() {
            // Alternate ±60° to produce a standard chemistry zigzag (120° bond angles)
            let turn = if i == 0 {
                sign * (PI / 3.0)
            } else {
                -sign * (PI / 3.0)
            };
            let new_dir = dir + turn;

            coords[v] = Vec2::new(coords[u].x + new_dir.cos(), coords[u].y + new_dir.sin());
            placed[v] = true;

            if let Some(ring_idx) = unfinished_ring_containing_atom(rings, v, placed) {
                place_ring_system_from_anchor(mol, rings, ring_idx, v, new_dir, coords, placed);
                place_substituents_for_placed_rings(mol, rings, coords, placed);
            } else {
                // Alternate turn sign for zigzag
                stack.push((v, new_dir, -sign));
            }
        }
    }
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
        // Find the smallest ring containing both endpoints of this bond
        let best_ring = rings
            .iter()
            .filter(|r| r.contains(&bond.from) && r.contains(&bond.to))
            .min_by_key(|r| r.len());

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
