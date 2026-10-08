/// Decides how the stereochemistry written in the SMILES appears in the 2D
/// drawing, and reports every configuration the drawing cannot show.
///
/// Coordinates are fixed before this module runs. Tetrahedral centers receive
/// a wedge or hash whose narrow tip sits on the stereocenter (IUPAC 2006,
/// ST-1.1), so a reader interprets each wedge at exactly one atom. Square-planar
/// centers and cis/trans double bonds are shown by geometry alone, so their
/// drawn geometry is checked against the written configuration. Allene,
/// trigonal-bipyramidal, and octahedral configurations are reported as
/// undepicted rather than drawn misleadingly flat.
use std::collections::HashSet;
use std::f64::consts::PI;

use crate::geometry::{angle_delta, largest_angular_gap, normalize_angle};
use crate::graph::{AtomChirality, BondDirection, BondOrder, BondStereo, MoleculeGraph};
use crate::layout::{implicit_h_count, lone_pair_count};
use crate::macrocycles::{ring_double_bond_requirements, TRANS_RING_MIN_ATOMS};
use crate::render::Vec2;
use crate::rings::ring_has_edge;

/// A written stereo configuration that the drawing does not show.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct UndepictedStereo {
    pub atom: usize,
    pub reason: String,
}

pub(crate) struct StereoDepiction {
    /// Per bond: the wedge, hash, wavy, or dashed style to draw.
    pub bond_stereo: Vec<BondStereo>,
    /// Per bond: the atom at the narrow end of a wedge or hash.
    pub wedge_tips: Vec<Option<usize>>,
    /// Per atom: a wedged or hashed hydrogen and its drawing direction.
    pub hydrogen_stereo: Vec<Option<(BondStereo, Vec2)>>,
    pub undepicted: Vec<UndepictedStereo>,
}

pub(crate) fn depict_stereo(
    molecule: &MoleculeGraph,
    coordinates: &[Vec2],
    rings: &[Vec<usize>],
    ring_bonds: &HashSet<usize>,
) -> StereoDepiction {
    let mut depiction = StereoDepiction {
        bond_stereo: molecule.bonds.iter().map(|bond| bond.stereo).collect(),
        // Drawing extensions are written after the atom at their narrow end.
        wedge_tips: molecule
            .bonds
            .iter()
            .map(|bond| bond.stereo.is_wedge().then_some(bond.from))
            .collect(),
        hydrogen_stereo: vec![None; molecule.n_atoms()],
        undepicted: Vec::new(),
    };

    for (center, atom) in molecule.atoms.iter().enumerate() {
        let outcome = match atom.chirality {
            AtomChirality::None => Ok(()),
            AtomChirality::TetraAnti | AtomChirality::TetraClockwise => {
                depict_tetrahedral_center(molecule, center, coordinates, ring_bonds, &mut depiction)
            }
            AtomChirality::SquarePlanar(class) => {
                check_square_planar_geometry(molecule, center, class, coordinates)
            }
            AtomChirality::Allenal(_) => {
                Err("allene (extended tetrahedral) configurations are not drawn".to_string())
            }
            AtomChirality::TrigonalBipyramidal(_) => {
                Err("trigonal-bipyramidal configurations are not drawn".to_string())
            }
            AtomChirality::Octahedral(_) => {
                Err("octahedral configurations are not drawn".to_string())
            }
        };
        if let Err(reason) = outcome {
            depiction.undepicted.push(UndepictedStereo {
                atom: center,
                reason: format!("{}: {reason}", describe_stereo_atom(molecule, center)),
            });
        }
    }

    report_undrawn_double_bond_configurations(
        molecule,
        coordinates,
        rings,
        &mut depiction.undepicted,
    );
    depiction
}

fn describe_stereo_atom(molecule: &MoleculeGraph, atom_index: usize) -> String {
    let atom = &molecule.atoms[atom_index];
    format!(
        "`{}` written {} at character {} (atom {atom_index})",
        atom.written_symbol(),
        atom.chirality.notation(),
        atom.source_position
    )
}

// ── Tetrahedral centers ──────────────────────────────────────────────────────

/// One neighbor of a tetrahedral center in OpenSMILES ordering.
#[derive(Clone, Copy, PartialEq)]
enum TetrahedralNeighbor {
    Bond {
        bond_index: usize,
        atom_index: usize,
    },
    /// A counted (implicit or bracket) hydrogen, which can carry the wedge.
    Hydrogen,
    /// The lone pair of a three-coordinate center such as a sulfoxide sulfur.
    LonePair,
}

/// Assigns a wedge or hash that makes the drawing reproduce the written
/// tetrahedral configuration.
///
/// The four neighbor directions (in OpenSMILES order) form a signed volume:
/// `@` requires a negative one, `@@` a positive one. The wedged neighbor points
/// toward the viewer and an undrawn hydrogen or lone pair away from it.
fn depict_tetrahedral_center(
    molecule: &MoleculeGraph,
    center: usize,
    coordinates: &[Vec2],
    ring_bonds: &HashSet<usize>,
    depiction: &mut StereoDepiction,
) -> Result<(), String> {
    let implicit_neighbor = implicit_tetrahedral_neighbor(molecule, center)?;
    let neighbor_order = tetrahedral_neighbor_order(molecule, center, implicit_neighbor);
    let implicit_direction = stereochemical_hydrogen_direction(molecule, center, coordinates, None);
    let parity = if molecule.atoms[center].chirality == AtomChirality::TetraAnti {
        -1.0
    } else {
        1.0
    };

    let candidates = wedge_candidates(
        molecule,
        center,
        ring_bonds,
        &depiction.bond_stereo,
        implicit_neighbor,
    );
    for candidate in candidates {
        let directions = tetrahedral_directions(
            &neighbor_order,
            candidate,
            center,
            coordinates,
            implicit_direction,
        );
        let volume = signed_volume(&directions);
        if volume.abs() < 1e-9 {
            continue;
        }
        let stereo = if volume.signum() == parity {
            BondStereo::WedgeUp
        } else {
            BondStereo::WedgeDown
        };
        match candidate {
            TetrahedralNeighbor::Bond { bond_index, .. } => {
                depiction.bond_stereo[bond_index] = stereo;
                depiction.wedge_tips[bond_index] = Some(center);
            }
            TetrahedralNeighbor::Hydrogen => {
                depiction.hydrogen_stereo[center] = Some((stereo, implicit_direction));
            }
            TetrahedralNeighbor::LonePair => unreachable!("a lone pair is never wedged"),
        }
        return Ok(());
    }
    Err("the layout leaves no bond that can show this configuration".to_string())
}

/// The fourth tetrahedral neighbor that is not a written bond, if any, or the
/// reason the atom cannot be a tetrahedral stereocenter.
fn implicit_tetrahedral_neighbor(
    molecule: &MoleculeGraph,
    center: usize,
) -> Result<Option<TetrahedralNeighbor>, String> {
    let atom = &molecule.atoms[center];
    let bond_count = molecule.neighbor_bonds[center].len();
    let hydrogen_count = usize::from(atom.hcount + implicit_h_count(molecule, center));

    if hydrogen_count >= 2 {
        return Err(format!(
            "it carries {hydrogen_count} hydrogens, so it is not a stereocenter"
        ));
    }
    if bond_count + hydrogen_count == 4 {
        return Ok((hydrogen_count == 1).then_some(TetrahedralNeighbor::Hydrogen));
    }
    if bond_count == 3 && hydrogen_count == 0 && lone_pair_count(molecule, center) >= 1 {
        return Ok(Some(TetrahedralNeighbor::LonePair));
    }
    if is_cumulated_center(molecule, center) {
        return Err(
            "on the central atom of an allene it describes an extended tetrahedral \
             configuration, which is not drawn"
                .to_string(),
        );
    }
    Err(format!(
        "a tetrahedral stereocenter needs four neighbors, or three and a lone pair, but it has \
         {}",
        bond_count + hydrogen_count
    ))
}

fn is_cumulated_center(molecule: &MoleculeGraph, center: usize) -> bool {
    let bonds = &molecule.adj[center];
    bonds.len() == 2
        && bonds
            .iter()
            .all(|&(_, bond_index)| molecule.bonds[bond_index].order == BondOrder::Double)
}

/// Neighbors in OpenSMILES order. A counted hydrogen follows the preceding
/// atom, or comes first when the center starts the SMILES (OpenSMILES 3.9).
/// The specification leaves the lone pair's place open; it counts as the last
/// neighbor, as in RDKit, which agrees with the hydrogen rule whenever the
/// center has a preceding atom.
fn tetrahedral_neighbor_order(
    molecule: &MoleculeGraph,
    center: usize,
    implicit_neighbor: Option<TetrahedralNeighbor>,
) -> Vec<TetrahedralNeighbor> {
    let mut neighbor_order: Vec<TetrahedralNeighbor> = molecule.neighbor_bonds[center]
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
    match implicit_neighbor {
        Some(TetrahedralNeighbor::LonePair) => neighbor_order.push(TetrahedralNeighbor::LonePair),
        Some(hydrogen) => {
            let hydrogen_position = usize::from(molecule.has_preceding[center]);
            neighbor_order.insert(hydrogen_position.min(neighbor_order.len()), hydrogen);
        }
        None => {}
    }
    neighbor_order
}

/// Wedge candidates, most readable first: preferred substituent bonds, then a
/// counted hydrogen, then the remaining single bonds, including ring bonds for
/// spiro and ring-fusion centers that have nothing else to wedge.
fn wedge_candidates(
    molecule: &MoleculeGraph,
    center: usize,
    ring_bonds: &HashSet<usize>,
    bond_stereo: &[BondStereo],
    implicit_neighbor: Option<TetrahedralNeighbor>,
) -> Vec<TetrahedralNeighbor> {
    let mut scored_bonds: Vec<(i32, TetrahedralNeighbor)> = molecule.adj[center]
        .iter()
        .filter(|&&(_, bond_index)| {
            molecule.bonds[bond_index].order == BondOrder::Single
                && bond_stereo[bond_index] == BondStereo::None
        })
        .map(|&(neighbor, bond_index)| {
            (
                tetrahedral_bond_score(molecule, neighbor, bond_index, ring_bonds),
                TetrahedralNeighbor::Bond {
                    bond_index,
                    atom_index: neighbor,
                },
            )
        })
        .collect();
    // A stable sort keeps writing order among equally readable bonds.
    scored_bonds.sort_by_key(|&(score, _)| std::cmp::Reverse(score));

    let (preferred, remaining): (Vec<_>, Vec<_>) =
        scored_bonds.into_iter().partition(|&(score, _)| score > 0);
    let mut candidates: Vec<TetrahedralNeighbor> = preferred
        .into_iter()
        .map(|(_, neighbor)| neighbor)
        .collect();
    if implicit_neighbor == Some(TetrahedralNeighbor::Hydrogen) {
        candidates.push(TetrahedralNeighbor::Hydrogen);
    }
    candidates.extend(remaining.into_iter().map(|(_, neighbor)| neighbor));
    candidates
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

/// Pseudo-3D directions for the signed volume: drawn bonds lie in the page,
/// the wedged neighbor points toward the viewer, and an undrawn hydrogen or
/// lone pair points away.
fn tetrahedral_directions(
    neighbor_order: &[TetrahedralNeighbor],
    wedged_neighbor: TetrahedralNeighbor,
    center: usize,
    coordinates: &[Vec2],
    implicit_direction: Vec2,
) -> [[f64; 3]; 4] {
    let mut directions = [[0.0_f64; 3]; 4];
    for (slot, neighbor) in neighbor_order.iter().enumerate() {
        let direction = match neighbor {
            TetrahedralNeighbor::Bond { atom_index, .. } => {
                normalize_vector(vector_from(coordinates[center], coordinates[*atom_index]))
            }
            TetrahedralNeighbor::Hydrogen | TetrahedralNeighbor::LonePair => implicit_direction,
        };
        let depth = if *neighbor == wedged_neighbor {
            1.0
        } else if matches!(neighbor, TetrahedralNeighbor::Bond { .. }) {
            0.0
        } else {
            -1.0
        };
        directions[slot] = [direction.x, direction.y, depth];
    }
    directions
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

// ── Square-planar centers ────────────────────────────────────────────────────

/// Neighbors sharpest to `trans` within this angle still read as opposite
/// corners of a square-planar center.
const TRANS_ANGLE_TOLERANCE: f64 = PI / 6.0;

/// Checks that the drawn square shows the written shape class: the neighbors
/// written first and third are trans for 'U' (@SP1), first and second for '4'
/// (@SP2), and first and fourth for 'Z' (@SP3).
fn check_square_planar_geometry(
    molecule: &MoleculeGraph,
    center: usize,
    class: u8,
    coordinates: &[Vec2],
) -> Result<(), String> {
    let neighbor_bonds = &molecule.neighbor_bonds[center];
    if neighbor_bonds.len() != 4 {
        return Err(format!(
            "a square-planar center needs four written neighbors, but it has {}",
            neighbor_bonds.len()
        ));
    }
    let angles: Vec<f64> = neighbor_bonds
        .iter()
        .map(|&bond_index| {
            let bond = &molecule.bonds[bond_index];
            let neighbor = if bond.from == center {
                bond.to
            } else {
                bond.from
            };
            (coordinates[neighbor].y - coordinates[center].y)
                .atan2(coordinates[neighbor].x - coordinates[center].x)
        })
        .collect();
    let trans_pairs: [(usize, usize); 2] = match class {
        1 => [(0, 2), (1, 3)],
        2 => [(0, 1), (2, 3)],
        _ => [(0, 3), (1, 2)],
    };
    let pairs_are_opposite = trans_pairs.iter().all(|&(first, second)| {
        angle_delta(angles[first], angles[second]).abs() >= PI - TRANS_ANGLE_TOLERANCE
    });
    if pairs_are_opposite {
        return Ok(());
    }
    Err(
        "the layout cannot place its four neighbors at the corners of a square with the \
         written trans pairs"
            .to_string(),
    )
}

// ── Double bonds ─────────────────────────────────────────────────────────────

/// Reports cis/trans double bonds whose drawn substituents do not lie on the
/// sides their directional `/` and `\` bonds require, for example when a small
/// ring forces the opposite geometry.
fn report_undrawn_double_bond_configurations(
    molecule: &MoleculeGraph,
    coordinates: &[Vec2],
    rings: &[Vec<usize>],
    undepicted: &mut Vec<UndepictedStereo>,
) {
    for double_bond in molecule
        .bonds
        .iter()
        .filter(|bond| bond.order == BondOrder::Double)
    {
        let left = marked_substituent_sides(molecule, double_bond.from, double_bond.to);
        let right = marked_substituent_sides(molecule, double_bond.to, double_bond.from);
        if left.is_empty() || right.is_empty() {
            continue;
        }
        // Cis/trans depends only on whether substituents share a side, so a
        // drawing mirrored as a whole still shows the written configuration.
        let axis_start = coordinates[double_bond.from];
        let axis_end = coordinates[double_bond.to];
        let side_agreements: Vec<i8> = left
            .iter()
            .chain(right.iter())
            .map(|&(neighbor, side)| {
                side_of_axis(axis_start, axis_end, coordinates[neighbor]) * side
            })
            .collect();
        let shows_written_configuration = side_agreements.iter().all(|&agreement| agreement == 1)
            || side_agreements.iter().all(|&agreement| agreement == -1);
        if shows_written_configuration {
            continue;
        }
        let from_atom = &molecule.atoms[double_bond.from];
        let to_atom = &molecule.atoms[double_bond.to];
        undepicted.push(UndepictedStereo {
            atom: double_bond.from,
            reason: format!(
                "the double bond between `{}` at character {} and `{}` at character {} (atoms \
                 {} and {}): the layout cannot place its `/` and `\\` substituents on the \
                 written sides{}",
                from_atom.written_symbol(),
                from_atom.source_position,
                to_atom.written_symbol(),
                to_atom.source_position,
                double_bond.from,
                double_bond.to,
                small_ring_trans_advice(molecule, rings, double_bond.from, double_bond.to)
            ),
        });
    }
}

/// Explains why a trans double bond written inside a small ring cannot be
/// drawn, or returns an empty string for any other double bond.
fn small_ring_trans_advice(
    molecule: &MoleculeGraph,
    rings: &[Vec<usize>],
    first_atom: usize,
    second_atom: usize,
) -> String {
    let Some(smallest_ring) = rings
        .iter()
        .filter(|ring| ring_has_edge(ring, first_atom, second_atom))
        .min_by_key(|ring| ring.len())
    else {
        return String::new();
    };
    let requests_trans = ring_double_bond_requirements(molecule, smallest_ring)
        .iter()
        .any(|requirement| {
            let bond_atoms = [
                smallest_ring[requirement.double_bond_start],
                smallest_ring[requirement.double_bond_end],
            ];
            bond_atoms.contains(&first_atom)
                && bond_atoms.contains(&second_atom)
                && !requirement.same_side
        });
    if !requests_trans || smallest_ring.len() >= TRANS_RING_MIN_ATOMS {
        return String::new();
    }
    format!(
        "; it lies in a ring of {} atoms, and a flat drawing can show a trans ring double bond \
         only in a ring of at least {TRANS_RING_MIN_ATOMS} atoms",
        smallest_ring.len()
    )
}

/// Directional substituents of one double-bond end with the side of the bond
/// axis each must occupy (+1 or -1), read in the bond's writing orientation.
/// Only the relative sides of substituents carry meaning.
fn marked_substituent_sides(
    molecule: &MoleculeGraph,
    end_atom: usize,
    partner_atom: usize,
) -> Vec<(usize, i8)> {
    molecule.adj[end_atom]
        .iter()
        .filter(|&&(neighbor, _)| neighbor != partner_atom)
        .filter_map(|&(neighbor, bond_index)| {
            let bond = &molecule.bonds[bond_index];
            let written_side = match bond.direction {
                BondDirection::None => return None,
                BondDirection::Up => 1,
                BondDirection::Down => -1,
            };
            let side = if bond.from == end_atom {
                written_side
            } else {
                -written_side
            };
            Some((neighbor, side))
        })
        .collect()
}

fn side_of_axis(axis_start: Vec2, axis_end: Vec2, point: Vec2) -> i8 {
    let cross = (axis_end.x - axis_start.x) * (point.y - axis_start.y)
        - (axis_end.y - axis_start.y) * (point.x - axis_start.x);
    if cross > 1e-9 {
        1
    } else if cross < -1e-9 {
        -1
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use crate::layout_native;
    use crate::render::LayoutOutput;

    fn layout(smiles: &str) -> LayoutOutput {
        layout_native(smiles).unwrap_or_else(|error| panic!("{smiles}: {error}"))
    }

    fn undepicted_reasons(smiles: &str) -> Vec<String> {
        layout(smiles)
            .undepicted_stereo
            .into_iter()
            .map(|undepicted| undepicted.reason)
            .collect()
    }

    /// Every wedge or hash drawn for a stereocenter has its narrow tip on it.
    fn assert_wedges_start_at_stereocenters(smiles: &str) {
        let layout_output = layout(smiles);
        for bond in &layout_output.bonds {
            if bond.stereo != "wedge_up" && bond.stereo != "wedge_down" {
                continue;
            }
            let tip = bond.stereo_tip.expect("every wedge names its tip");
            assert!(tip == bond.from || tip == bond.to, "{smiles}");
            if !bond.forced_stereo {
                assert_ne!(layout_output.atoms[tip].chirality, "none", "{smiles}");
            }
        }
    }

    #[test]
    fn written_stereo_in_repository_molecules_is_drawn() {
        let library = include_str!("../../src/molecule/library.typ");
        let library_smiles = library
            .lines()
            .filter_map(|line| line.split_once(": \"").map(|(_, rest)| rest))
            .filter_map(|rest| rest.strip_suffix("\","))
            .map(|smiles| smiles.replace("\\\\", "\\"));
        for smiles in library_smiles {
            let layout_output = layout(&smiles);
            assert!(
                layout_output.undepicted_stereo.is_empty(),
                "{smiles}: {:?}",
                layout_output.undepicted_stereo
            );
            assert_wedges_start_at_stereocenters(&smiles);
        }
    }

    #[test]
    fn heteroatom_stereocenters_put_the_wedge_tip_on_the_heteroatom() {
        for smiles in [
            "C[N@+](CC)(CCC)Cc1ccccc1",
            "C[P@](=O)(CC)c1ccccc1",
            "OC(=O)[C@@H]1CCCN1",
            "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O",
        ] {
            assert_wedges_start_at_stereocenters(smiles);
            assert!(undepicted_reasons(smiles).is_empty(), "{smiles}");
        }
        let ammonium = layout("C[N@+](CC)(CCC)Cc1ccccc1");
        let wedge = ammonium
            .bonds
            .iter()
            .find(|bond| bond.stereo != "none")
            .expect("the ammonium center is wedged");
        assert_eq!(wedge.stereo_tip, Some(1));
    }

    #[test]
    fn lone_pair_stereocenters_are_drawn() {
        for smiles in [
            "C[S@](=O)c1ccccc1",
            "COc1ccc2[nH]c(nc2c1)[S@](=O)Cc1ncc(C)c(OC)c1C",
            "C[P@](CC)c1ccccc1",
            "C[N@](CC)CCC",
            "C[C@-](F)Cl",
        ] {
            let layout_output = layout(smiles);
            assert!(layout_output.undepicted_stereo.is_empty(), "{smiles}");
            assert!(
                layout_output.bonds.iter().any(|bond| bond.stereo != "none"),
                "{smiles} has no wedge"
            );
        }
    }

    #[test]
    fn inverting_a_lone_pair_center_flips_its_wedge() {
        let first = layout("C[S@](=O)c1ccccc1");
        let second = layout("C[S@@](=O)c1ccccc1");
        let wedge = |layout_output: &LayoutOutput| {
            layout_output
                .bonds
                .iter()
                .find(|bond| bond.stereo != "none")
                .map(|bond| (bond.from, bond.to, bond.stereo.clone()))
        };
        let (first_from, first_to, first_stereo) = wedge(&first).unwrap();
        let (second_from, second_to, second_stereo) = wedge(&second).unwrap();
        assert_eq!((first_from, first_to), (second_from, second_to));
        assert_ne!(first_stereo, second_stereo);
    }

    #[test]
    fn spiro_and_ring_fusion_centers_wedge_a_ring_bond() {
        for smiles in ["CC[C@H](O1)CC[C@@]12CCCO2", "C[C@@]12CCCC[C@@]1(C)CCCC2"] {
            assert!(undepicted_reasons(smiles).is_empty(), "{smiles}");
            assert_wedges_start_at_stereocenters(smiles);
        }
    }

    #[test]
    fn non_stereocenters_with_tetrahedral_marks_are_reported() {
        let reasons = undepicted_reasons("[C@H2](F)Cl");
        assert!(reasons[0].contains("carries 2 hydrogens"));
        assert!(reasons[0].contains("character 1"));
        let reasons = undepicted_reasons("C[C@](F)Cl");
        assert!(reasons[0].contains("needs four neighbors, or three and a lone pair"));
        let reasons = undepicted_reasons("NC(Br)=[C@]=C(O)C");
        assert!(reasons[0].contains("central atom of an allene"));
    }

    #[test]
    fn square_planar_geometry_is_checked() {
        for class in 1..=3 {
            let smiles = format!("N[Pt@SP{class}](N)(Cl)Cl");
            assert!(undepicted_reasons(&smiles).is_empty(), "{smiles}");
        }
        assert!(undepicted_reasons("N1CC[NH2][Pt@SP1]1(Cl)Br").is_empty());
        let reasons = undepicted_reasons("N1CC[NH2][Pt@SP2]1(Cl)Br");
        assert!(reasons[0].contains("corners of a square"));
        let reasons = undepicted_reasons("[Pt@SP1](Cl)(Cl)N");
        assert!(reasons[0].contains("needs four written neighbors"));
    }

    #[test]
    fn cis_trans_bonds_the_layout_cannot_honor_are_reported() {
        for smiles in [
            "C/C=C/C",
            "C/C=C\\C",
            "F/C=C1.Cl/1",
            "C1CCC/C=C\\CC1",
            "C/C=C/C=C/C",
        ] {
            assert!(undepicted_reasons(smiles).is_empty(), "{smiles}");
        }
        let reasons = undepicted_reasons("C1CCC/C=C/CC1");
        assert!(reasons[0].contains("double bond between `C` at character 7"));
    }
}
