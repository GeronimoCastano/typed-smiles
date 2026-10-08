//! Rigid scaffold alignment of independently laid-out molecules.
//!
//! Each molecule keeps the coordinates chosen by the layout engine. A rotation,
//! optionally combined with a reflection, then orients every molecule so its
//! corresponding atoms best overlay the reference molecule's atoms. The core
//! geometry itself is not changed, so two cores only coincide exactly when the
//! layout engine drew them identically.
use crate::graph::{Bond, BondOrder, MoleculeGraph};
use crate::layout::compute_layout;
use crate::render::Vec2;
use crate::substructure::find_embeddings;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Fits whose root-mean-square deviations differ by less than this many bond
/// lengths are treated as equally good, so floating-point noise never decides
/// between symmetric correspondences.
const EQUAL_FIT_TOLERANCE: f64 = 1e-6;

/// Distinct scaffold occurrences listed in an ambiguity diagnostic.
const LISTED_OCCURRENCE_LIMIT: usize = 6;

#[derive(Debug, Deserialize)]
pub struct AlignmentRequest {
    pub molecules: Vec<String>,
    /// SMARTS pattern shared by every molecule, when correspondences are found
    /// by matching rather than given explicitly.
    #[serde(default)]
    pub scaffold: Option<String>,
    /// Per-molecule explicit correspondences; `None` entries use the scaffold.
    #[serde(default)]
    pub atoms: Vec<Option<Vec<usize>>>,
    #[serde(default)]
    pub reference: usize,
    #[serde(default = "reflection_allowed_by_default")]
    pub allow_reflection: bool,
}

fn reflection_allowed_by_default() -> bool {
    true
}

/// A 2x2 orthogonal matrix acting on layout coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct OrthogonalTransform {
    pub xx: f64,
    pub xy: f64,
    pub yx: f64,
    pub yy: f64,
}

impl OrthogonalTransform {
    const IDENTITY: Self = Self {
        xx: 1.0,
        xy: 0.0,
        yx: 0.0,
        yy: 1.0,
    };

    /// The rotation whose cosine and sine are proportional to `dot` and `cross`.
    fn rotation_from_components(dot: f64, cross: f64) -> Self {
        let length = dot.hypot(cross);
        if length == 0.0 {
            return Self::IDENTITY;
        }
        let (cosine, sine) = (dot / length, cross / length);
        Self {
            xx: cosine,
            xy: -sine,
            yx: sine,
            yy: cosine,
        }
    }

    /// This transform applied after mirroring the x axis.
    fn after_x_mirror(self) -> Self {
        Self {
            xx: -self.xx,
            xy: self.xy,
            yx: -self.yx,
            yy: self.yy,
        }
    }

    fn apply(self, point: Vec2) -> Vec2 {
        Vec2::new(
            self.xx * point.x + self.xy * point.y,
            self.yx * point.x + self.yy * point.y,
        )
    }

    fn is_reflection(self) -> bool {
        self.xx * self.yy - self.xy * self.yx < 0.0
    }
}

/// How one molecule is oriented onto the reference.
#[derive(Debug, Clone, Serialize)]
pub struct MoleculeAlignment {
    /// Molecule atoms in correspondence order: entry `k` overlays entry `k`
    /// of the reference molecule's atoms.
    pub atoms: Vec<usize>,
    /// Maps this molecule's layout coordinates into the reference frame.
    pub transform: OrthogonalTransform,
    pub reflected: bool,
    /// Largest distance, in bond lengths, between a corresponding atom and its
    /// reference atom after alignment.
    pub deviation: f64,
}

struct PreparedMolecule {
    graph: MoleculeGraph,
    positions: Vec<Vec2>,
}

/// Correspondences a molecule may use: one explicit list, or every embedding
/// of the scaffold's single occurrence.
struct CorrespondenceCandidates {
    embeddings: Vec<Vec<usize>>,
}

struct CandidateFit {
    atoms: Vec<usize>,
    attachment_agreement: usize,
    transform: OrthogonalTransform,
    root_mean_square_deviation: f64,
    deviation: f64,
}

pub fn align_molecules(request: &AlignmentRequest) -> Result<Vec<MoleculeAlignment>, String> {
    validate_request_shape(request)?;
    let molecules = request
        .molecules
        .iter()
        .enumerate()
        .map(|(molecule_index, smiles)| prepare_molecule(molecule_index, smiles))
        .collect::<Result<Vec<_>, String>>()?;
    let candidates = molecules
        .iter()
        .enumerate()
        .map(|(molecule_index, molecule)| {
            correspondence_candidates(request, molecule_index, molecule)
        })
        .collect::<Result<Vec<_>, String>>()?;

    let reference_atoms = candidates[request.reference].embeddings[0].clone();
    if reference_atoms.len() < 2 {
        return Err(format!(
            "alignment needs at least two corresponding atoms to fix an orientation, got {}; use a scaffold or atom list with two or more atoms",
            reference_atoms.len()
        ));
    }
    let reference = &molecules[request.reference];
    molecules
        .iter()
        .zip(&candidates)
        .enumerate()
        .map(|(molecule_index, (molecule, molecule_candidates))| {
            if molecule_index == request.reference {
                return Ok(MoleculeAlignment {
                    atoms: reference_atoms.clone(),
                    transform: OrthogonalTransform::IDENTITY,
                    reflected: false,
                    deviation: 0.0,
                });
            }
            let best_fit = best_candidate_fit(
                reference,
                &reference_atoms,
                molecule,
                molecule_candidates,
                request.allow_reflection,
            );
            Ok(MoleculeAlignment {
                atoms: best_fit.atoms,
                transform: best_fit.transform,
                reflected: best_fit.transform.is_reflection(),
                deviation: best_fit.deviation,
            })
        })
        .collect()
}

fn validate_request_shape(request: &AlignmentRequest) -> Result<(), String> {
    if request.molecules.len() < 2 {
        return Err(format!(
            "alignment needs at least two molecules, got {}; pass the molecules to compare",
            request.molecules.len()
        ));
    }
    if request.reference >= request.molecules.len() {
        return Err(format!(
            "reference molecule {} does not exist; use an index from 0 to {}",
            request.reference,
            request.molecules.len() - 1
        ));
    }
    if !request.atoms.is_empty() && request.atoms.len() != request.molecules.len() {
        return Err(format!(
            "atoms has {} entries for {} molecules; give one entry (auto or an index list) per molecule",
            request.atoms.len(),
            request.molecules.len()
        ));
    }
    if request.scaffold.is_none() && explicit_atoms(request, 0).is_none() {
        return Err(
            "alignment needs a scaffold pattern or an explicit atom list for every molecule".into(),
        );
    }
    Ok(())
}

fn explicit_atoms(request: &AlignmentRequest, molecule_index: usize) -> Option<&Vec<usize>> {
    request.atoms.get(molecule_index).and_then(Option::as_ref)
}

fn prepare_molecule(molecule_index: usize, smiles: &str) -> Result<PreparedMolecule, String> {
    let graph = crate::parse_molecule(smiles)
        .map_err(|error| format!("molecule {molecule_index}: {error}"))?;
    let layout =
        compute_layout(&graph).map_err(|error| format!("molecule {molecule_index}: {error}"))?;
    let positions = layout.atoms[..graph.n_atoms()]
        .iter()
        .map(|atom| atom.pos)
        .collect();
    Ok(PreparedMolecule { graph, positions })
}

/// Without a scaffold, an atom list is an exact atom-by-atom correspondence.
/// With a scaffold, an atom list only selects which occurrence to use, and the
/// correspondence within it is chosen like any symmetric scaffold match.
fn correspondence_candidates(
    request: &AlignmentRequest,
    molecule_index: usize,
    molecule: &PreparedMolecule,
) -> Result<CorrespondenceCandidates, String> {
    let selected_atoms = explicit_atoms(request, molecule_index);
    if let Some(atoms) = selected_atoms {
        validate_explicit_atoms(request, molecule_index, molecule, atoms)?;
    }
    let Some(scaffold) = &request.scaffold else {
        let Some(atoms) = selected_atoms else {
            return Err(format!(
                "molecule {molecule_index} has no atom list and no scaffold is given; pass a scaffold or an atom list for every molecule"
            ));
        };
        return Ok(CorrespondenceCandidates {
            embeddings: vec![atoms.clone()],
        });
    };
    let embeddings = find_embeddings(&molecule.graph, scaffold)
        .map_err(|error| format!("molecule {molecule_index}: {error}"))?;
    if embeddings.is_empty() {
        return Err(format!(
            "scaffold {scaffold:?} does not occur in molecule {molecule_index} ({:?}); use a scaffold shared by every molecule",
            request.molecules[molecule_index]
        ));
    }
    let mut embeddings = match selected_atoms {
        Some(atoms) => selected_occurrence_embeddings(scaffold, molecule_index, embeddings, atoms)?,
        None => single_symmetry_class_embeddings(scaffold, molecule_index, molecule, embeddings)?,
    };
    embeddings.sort_unstable();
    Ok(CorrespondenceCandidates { embeddings })
}

fn selected_occurrence_embeddings(
    scaffold: &str,
    molecule_index: usize,
    embeddings: Vec<Vec<usize>>,
    atoms: &[usize],
) -> Result<Vec<Vec<usize>>, String> {
    let selected_occurrence = sorted(atoms);
    let occurrence_embeddings: Vec<Vec<usize>> = embeddings
        .into_iter()
        .filter(|embedding| sorted(embedding) == selected_occurrence)
        .collect();
    if occurrence_embeddings.is_empty() {
        return Err(format!(
            "the atom list {} for molecule {molecule_index} is not an occurrence of scaffold {scaffold:?}; list the atoms of one place where the scaffold occurs, or omit the scaffold for an exact correspondence",
            index_list(atoms)
        ));
    }
    Ok(occurrence_embeddings)
}

fn single_symmetry_class_embeddings(
    scaffold: &str,
    molecule_index: usize,
    molecule: &PreparedMolecule,
    embeddings: Vec<Vec<usize>>,
) -> Result<Vec<Vec<usize>>, String> {
    let mut occurrence_groups = symmetry_distinct_occurrences(molecule, embeddings);
    if occurrence_groups.len() > 1 {
        return Err(ambiguous_occurrence_message(
            scaffold,
            molecule_index,
            &occurrence_groups,
        ));
    }
    Ok(occurrence_groups.remove(0))
}

fn validate_explicit_atoms(
    request: &AlignmentRequest,
    molecule_index: usize,
    molecule: &PreparedMolecule,
    atoms: &[usize],
) -> Result<(), String> {
    let atom_count = molecule.graph.n_atoms();
    if let Some(&missing_atom) = atoms.iter().find(|&&atom| atom >= atom_count) {
        return Err(format!(
            "atom {missing_atom} does not exist in molecule {molecule_index}; use indices from 0 to {}",
            atom_count - 1
        ));
    }
    if sorted(atoms).windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(format!(
            "the atom list {} for molecule {molecule_index} repeats an atom; list each corresponding atom once",
            index_list(atoms)
        ));
    }
    if let Some(first_atoms) = request.atoms.iter().flatten().next() {
        if first_atoms.len() != atoms.len() && request.scaffold.is_none() {
            return Err(format!(
                "the atom list for molecule {molecule_index} has {} atoms but another list has {}; give every molecule the same number of corresponding atoms",
                atoms.len(),
                first_atoms.len()
            ));
        }
    }
    Ok(())
}

/// Groups scaffold embeddings into occurrences that molecular symmetry makes
/// interchangeable, such as the two methyl groups of an isopropyl group.
/// Each group holds the embeddings of all its occurrences, ordered by its
/// lowest embedding.
fn symmetry_distinct_occurrences(
    molecule: &PreparedMolecule,
    embeddings: Vec<Vec<usize>>,
) -> Vec<Vec<Vec<usize>>> {
    let classes = symmetry_classes(&molecule.graph);
    let mut occurrences: BTreeMap<Vec<usize>, Vec<Vec<usize>>> = BTreeMap::new();
    for embedding in embeddings {
        occurrences
            .entry(sorted(&embedding))
            .or_default()
            .push(embedding);
    }
    let mut groups: BTreeMap<Vec<usize>, Vec<Vec<usize>>> = BTreeMap::new();
    for occurrence_embeddings in occurrences.into_values() {
        let symmetry_key = occurrence_embeddings
            .iter()
            .map(|embedding| {
                embedding
                    .iter()
                    .map(|&atom| classes[atom])
                    .collect::<Vec<_>>()
            })
            .min()
            .expect("every occurrence has an embedding");
        groups
            .entry(symmetry_key)
            .or_default()
            .extend(occurrence_embeddings);
    }
    let mut groups: Vec<Vec<Vec<usize>>> = groups.into_values().collect();
    for group in &mut groups {
        group.sort_unstable();
    }
    groups.sort_unstable_by(|first, second| first[0].cmp(&second[0]));
    groups
}

/// Atom classes shared exactly by atoms that are indistinguishable from their
/// element, charge, hydrogens, and connectivity. Initial atom invariants are
/// refined by neighbor classes until no class splits any further.
fn symmetry_classes(graph: &MoleculeGraph) -> Vec<usize> {
    let initial_invariants: Vec<_> = (0..graph.n_atoms())
        .map(|atom| {
            let properties = &graph.atoms[atom];
            (
                properties.symbol.as_str(),
                properties.abbrev.as_str(),
                properties.aromatic,
                properties.charge,
                properties.isotope,
                crate::hydrogen_count(graph, atom),
                graph.adj[atom].len(),
            )
        })
        .collect();
    let mut classes = dense_ranks(&initial_invariants);
    loop {
        let neighborhood_signatures: Vec<(usize, Vec<(u8, usize)>)> = (0..graph.n_atoms())
            .map(|atom| {
                let mut neighbors: Vec<(u8, usize)> = graph.adj[atom]
                    .iter()
                    .map(|&(neighbor, bond)| (bond_class(&graph.bonds[bond]), classes[neighbor]))
                    .collect();
                neighbors.sort_unstable();
                (classes[atom], neighbors)
            })
            .collect();
        let refined_classes = dense_ranks(&neighborhood_signatures);
        if class_count(&refined_classes) == class_count(&classes) {
            return refined_classes;
        }
        classes = refined_classes;
    }
}

fn bond_class(bond: &Bond) -> u8 {
    if bond.aromatic {
        return 0;
    }
    match bond.order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 0,
    }
}

/// Replaces each value by its rank among the distinct values.
fn dense_ranks<T: Ord + Clone>(values: &[T]) -> Vec<usize> {
    let mut distinct_values = values.to_vec();
    distinct_values.sort_unstable();
    distinct_values.dedup();
    values
        .iter()
        .map(|value| {
            distinct_values
                .binary_search(value)
                .expect("every value is ranked")
        })
        .collect()
}

fn class_count(classes: &[usize]) -> usize {
    classes
        .iter()
        .max()
        .map_or(0, |highest_class| highest_class + 1)
}

fn ambiguous_occurrence_message(
    scaffold: &str,
    molecule_index: usize,
    occurrence_groups: &[Vec<Vec<usize>>],
) -> String {
    let listed: Vec<String> = occurrence_groups
        .iter()
        .take(LISTED_OCCURRENCE_LIMIT)
        .map(|group| index_list(&sorted(&group[0])))
        .collect();
    let remainder = if occurrence_groups.len() > LISTED_OCCURRENCE_LIMIT {
        ", …"
    } else {
        ""
    };
    format!(
        "scaffold {scaffold:?} occurs at {} chemically distinct places in molecule {molecule_index}, for example at atoms {}{remainder}; choose one by passing its atoms as that molecule's atoms entry",
        occurrence_groups.len(),
        listed.join(", ")
    )
}

/// Formats atom indices as a Typst array so diagnostics can be copied into
/// an `atoms:` argument.
fn index_list(atoms: &[usize]) -> String {
    let indices: Vec<String> = atoms.iter().map(usize::to_string).collect();
    format!("({})", indices.join(", "))
}

fn sorted(atoms: &[usize]) -> Vec<usize> {
    let mut sorted_atoms = atoms.to_vec();
    sorted_atoms.sort_unstable();
    sorted_atoms
}

/// Chooses the correspondence and orientation for one molecule. Agreement of
/// substituent attachment points decides first, so a symmetric scaffold keeps
/// its substituents on the same side as the reference; the closer overlay
/// decides next, then a rotation is preferred over a reflection, and finally
/// the lowest atom list keeps the choice deterministic.
fn best_candidate_fit(
    reference: &PreparedMolecule,
    reference_atoms: &[usize],
    molecule: &PreparedMolecule,
    candidates: &CorrespondenceCandidates,
    allow_reflection: bool,
) -> CandidateFit {
    let mut fits: Vec<CandidateFit> = Vec::new();
    for atoms in &candidates.embeddings {
        let attachment_agreement =
            attachment_agreement(reference, reference_atoms, molecule, atoms);
        let reference_points: Vec<Vec2> = reference_atoms
            .iter()
            .map(|&atom| reference.positions[atom])
            .collect();
        let molecule_points: Vec<Vec2> =
            atoms.iter().map(|&atom| molecule.positions[atom]).collect();
        let mut transforms = vec![best_rotation(&molecule_points, &reference_points)];
        if allow_reflection {
            transforms.push(best_reflected_rotation(&molecule_points, &reference_points));
        }
        for transform in transforms {
            let deviations = overlay_deviations(transform, &molecule_points, &reference_points);
            fits.push(CandidateFit {
                atoms: atoms.clone(),
                attachment_agreement,
                transform,
                root_mean_square_deviation: root_mean_square(&deviations),
                deviation: deviations.iter().copied().fold(0.0, f64::max),
            });
        }
    }
    fits.into_iter()
        .min_by(|first, second| {
            second
                .attachment_agreement
                .cmp(&first.attachment_agreement)
                .then_with(|| fit_quality_order(first, second))
                .then_with(|| {
                    first
                        .transform
                        .is_reflection()
                        .cmp(&second.transform.is_reflection())
                })
                .then_with(|| first.atoms.cmp(&second.atoms))
        })
        .expect("every molecule has at least one correspondence candidate")
}

fn fit_quality_order(first: &CandidateFit, second: &CandidateFit) -> std::cmp::Ordering {
    let difference = first.root_mean_square_deviation - second.root_mean_square_deviation;
    if difference.abs() < EQUAL_FIT_TOLERANCE {
        std::cmp::Ordering::Equal
    } else {
        difference.total_cmp(&0.0)
    }
}

/// Scores how well substituent attachment points line up: one point for each
/// corresponding scaffold atom substituted in both molecules, and one more when
/// the attached atoms are the same elements or labels.
fn attachment_agreement(
    reference: &PreparedMolecule,
    reference_atoms: &[usize],
    molecule: &PreparedMolecule,
    atoms: &[usize],
) -> usize {
    reference_atoms
        .iter()
        .zip(atoms)
        .map(|(&reference_atom, &atom)| {
            let reference_substituents =
                substituent_labels(reference, reference_atoms, reference_atom);
            let molecule_substituents = substituent_labels(molecule, atoms, atom);
            if reference_substituents.is_empty() || molecule_substituents.is_empty() {
                0
            } else if reference_substituents == molecule_substituents {
                2
            } else {
                1
            }
        })
        .sum()
}

/// Sorted labels of the non-hydrogen neighbors outside the corresponded atoms.
fn substituent_labels(
    molecule: &PreparedMolecule,
    corresponded_atoms: &[usize],
    atom: usize,
) -> Vec<String> {
    let graph = &molecule.graph;
    let mut labels: Vec<String> = graph.adj[atom]
        .iter()
        .map(|&(neighbor, _)| neighbor)
        .filter(|neighbor| !corresponded_atoms.contains(neighbor))
        .map(|neighbor| &graph.atoms[neighbor])
        .filter(|neighbor_atom| neighbor_atom.symbol != "H")
        .map(|neighbor_atom| {
            if neighbor_atom.abbrev.is_empty() {
                neighbor_atom.symbol.clone()
            } else {
                neighbor_atom.abbrev.clone()
            }
        })
        .collect();
    labels.sort();
    labels
}

fn centroid(points: &[Vec2]) -> Vec2 {
    let count = points.len() as f64;
    Vec2::new(
        points.iter().map(|point| point.x).sum::<f64>() / count,
        points.iter().map(|point| point.y).sum::<f64>() / count,
    )
}

fn centered(points: &[Vec2]) -> Vec<Vec2> {
    let center = centroid(points);
    points
        .iter()
        .map(|point| Vec2::new(point.x - center.x, point.y - center.y))
        .collect()
}

/// The least-squares rotation carrying `moving` onto `fixed` after both are
/// centered (the two-dimensional Kabsch solution).
fn best_rotation(moving: &[Vec2], fixed: &[Vec2]) -> OrthogonalTransform {
    let moving = centered(moving);
    let fixed = centered(fixed);
    let (cross_sum, dot_sum) = moving.iter().zip(&fixed).fold(
        (0.0, 0.0),
        |(cross_sum, dot_sum), (moving_point, fixed_point)| {
            (
                cross_sum + moving_point.x * fixed_point.y - moving_point.y * fixed_point.x,
                dot_sum + moving_point.x * fixed_point.x + moving_point.y * fixed_point.y,
            )
        },
    );
    OrthogonalTransform::rotation_from_components(dot_sum, cross_sum)
}

fn best_reflected_rotation(moving: &[Vec2], fixed: &[Vec2]) -> OrthogonalTransform {
    let mirrored: Vec<Vec2> = moving
        .iter()
        .map(|point| Vec2::new(-point.x, point.y))
        .collect();
    best_rotation(&mirrored, fixed).after_x_mirror()
}

/// Distances between transformed moving points and fixed points once both
/// centroids coincide.
fn overlay_deviations(transform: OrthogonalTransform, moving: &[Vec2], fixed: &[Vec2]) -> Vec<f64> {
    let moving = centered(moving);
    let fixed = centered(fixed);
    moving
        .iter()
        .zip(&fixed)
        .map(|(&moving_point, &fixed_point)| transform.apply(moving_point).distance_to(fixed_point))
        .collect()
}

fn root_mean_square(values: &[f64]) -> f64 {
    (values.iter().map(|value| value * value).sum::<f64>() / values.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request(molecules: &[&str], scaffold: Option<&str>) -> AlignmentRequest {
        AlignmentRequest {
            molecules: molecules.iter().map(|smiles| smiles.to_string()).collect(),
            scaffold: scaffold.map(str::to_string),
            atoms: Vec::new(),
            reference: 0,
            allow_reflection: true,
        }
    }

    fn aligned_positions(
        request: &AlignmentRequest,
        alignment: &MoleculeAlignment,
        molecule_index: usize,
    ) -> Vec<Vec2> {
        let molecule =
            prepare_molecule(molecule_index, &request.molecules[molecule_index]).unwrap();
        let points: Vec<Vec2> = alignment
            .atoms
            .iter()
            .map(|&atom| molecule.positions[atom])
            .collect();
        centered(&points)
            .into_iter()
            .map(|point| alignment.transform.apply(point))
            .collect()
    }

    #[test]
    fn reference_keeps_its_orientation() {
        let alignments = align_molecules(&request(
            &["CC(=O)c1ccccc1", "CC(O)c1ccccc1"],
            Some("c1ccccc1"),
        ))
        .unwrap();
        assert_eq!(alignments[0].transform, OrthogonalTransform::IDENTITY);
        assert!(!alignments[0].reflected);
        assert_eq!(alignments[0].deviation, 0.0);
    }

    #[test]
    fn differently_written_molecules_overlay_their_scaffold() {
        let alignment_request = request(
            &["c1ccccc1C(=O)O", "O=C(O)c1ccccc1"],
            Some("c1ccccc1C(=O)O"),
        );
        let alignments = align_molecules(&alignment_request).unwrap();
        assert!(
            alignments[1].deviation < 1e-6,
            "deviation {}",
            alignments[1].deviation
        );
        let reference = aligned_positions(&alignment_request, &alignments[0], 0);
        let molecule = aligned_positions(&alignment_request, &alignments[1], 1);
        for (reference_point, molecule_point) in reference.iter().zip(&molecule) {
            assert!(reference_point.distance_to(*molecule_point) < 1e-6);
        }
    }

    #[test]
    fn symmetric_scaffold_keeps_substituents_on_the_reference_side() {
        let alignment_request = request(&["Clc1ccccc1OC", "COc1ccccc1Br"], Some("c1ccccc1"));
        let alignments = align_molecules(&alignment_request).unwrap();
        let reference = prepare_molecule(0, &alignment_request.molecules[0]).unwrap();
        let molecule = prepare_molecule(1, &alignment_request.molecules[1]).unwrap();
        let methoxy_ring_atom = |prepared: &PreparedMolecule, atoms: &[usize]| {
            atoms
                .iter()
                .position(|&atom| substituent_labels(prepared, atoms, atom) == ["O"])
                .unwrap()
        };
        assert_eq!(
            methoxy_ring_atom(&reference, &alignments[0].atoms),
            methoxy_ring_atom(&molecule, &alignments[1].atoms)
        );
        assert!(alignments[1].deviation < 1e-6);
    }

    #[test]
    fn rotation_is_preferred_when_a_reflection_fits_equally() {
        let alignments =
            align_molecules(&request(&["Oc1ccccc1", "Nc1ccccc1"], Some("c1ccccc1"))).unwrap();
        assert!(!alignments[1].reflected);
    }

    #[test]
    fn mirror_image_correspondence_needs_a_reflection() {
        let mut alignment_request = request(&["CC(O)N", "CC(O)N"], None);
        alignment_request.atoms = vec![Some(vec![0, 1, 2, 3]), Some(vec![0, 1, 3, 2])];
        let reflected = align_molecules(&alignment_request).unwrap();
        assert!(reflected[1].reflected);
        assert!(
            reflected[1].deviation < 1e-6,
            "deviation {}",
            reflected[1].deviation
        );

        alignment_request.allow_reflection = false;
        let rotated = align_molecules(&alignment_request).unwrap();
        assert!(!rotated[1].reflected);
        assert!(
            rotated[1].deviation > 0.1,
            "deviation {}",
            rotated[1].deviation
        );
    }

    #[test]
    fn explicit_atoms_choose_between_distinct_scaffold_occurrences() {
        let mut alignment_request =
            request(&["c1ccccc1", "c1ccccc1-c1ccc(Cl)cc1"], Some("c1ccccc1"));
        let error = align_molecules(&alignment_request).unwrap_err();
        assert!(
            error.contains("occurs at 2 chemically distinct places in molecule 1"),
            "{error}"
        );
        assert!(
            error.contains("at atoms (0, 1, 2, 3, 4, 5), (6, 7, 8, 9, 11, 12);"),
            "{error}"
        );
        alignment_request.atoms = vec![None, Some(vec![12, 11, 9, 8, 7, 6])];
        let alignments = align_molecules(&alignment_request).unwrap();
        assert_eq!(sorted(&alignments[1].atoms), vec![6, 7, 8, 9, 11, 12]);
    }

    #[test]
    fn selected_occurrence_still_matches_substituents() {
        let mut alignment_request =
            request(&["Clc1ccccc1", "c1ccccc1-c1ccc(Cl)cc1"], Some("c1ccccc1"));
        alignment_request.atoms = vec![None, Some(vec![6, 7, 8, 9, 11, 12])];
        let alignments = align_molecules(&alignment_request).unwrap();
        let chlorinated_position = |atoms: &[usize], chlorinated_atom: usize| {
            atoms
                .iter()
                .position(|&atom| atom == chlorinated_atom)
                .unwrap()
        };
        assert_eq!(
            chlorinated_position(&alignments[0].atoms, 1),
            chlorinated_position(&alignments[1].atoms, 9)
        );
    }

    #[test]
    fn symmetry_equivalent_occurrences_are_chosen_automatically() {
        let alignments = align_molecules(&request(
            &["CC(C)Cc1ccccc1", "CC(C)(O)Cc1ccccc1"],
            Some("CCCc1ccccc1"),
        ))
        .unwrap();
        assert_eq!(alignments[0].atoms[1..], [1, 3, 4, 5, 6, 7, 8, 9]);
        let alignments = align_molecules(&request(
            &["c1ccccc1", "c1ccccc1-c1ccccc1"],
            Some("c1ccccc1"),
        ))
        .unwrap();
        assert!(alignments[1].deviation < 1e-6);
    }

    #[test]
    fn explicit_atoms_align_without_a_scaffold() {
        let mut alignment_request = request(&["CCO", "OCC"], None);
        alignment_request.atoms = vec![Some(vec![0, 1, 2]), Some(vec![2, 1, 0])];
        let alignments = align_molecules(&alignment_request).unwrap();
        assert!(alignments[1].deviation < 1e-6);
    }

    #[test]
    fn invalid_requests_are_diagnostics() {
        let error_for =
            |alignment_request: AlignmentRequest| align_molecules(&alignment_request).unwrap_err();

        assert!(error_for(request(&["CCO"], Some("CC"))).contains("at least two molecules"));
        assert!(error_for(request(&["CCO", "CCN"], None)).contains("needs a scaffold pattern"));
        assert!(error_for(request(&["CCO", "CCN"], Some("c1ccccc1")))
            .contains("does not occur in molecule 0"));
        assert!(error_for(request(&["CCO", "CCN"], Some("[Q]")))
            .contains("molecule 0: typed-smiles: invalid SMARTS"));
        assert!(
            error_for(request(&["CCO", "C1CC"], Some("CC"))).contains("molecule 1: invalid SMILES")
        );
        assert!(error_for(request(&["CO", "CN"], Some("C")))
            .contains("at least two corresponding atoms"));

        let mut out_of_range_reference = request(&["CCO", "CCN"], Some("CC"));
        out_of_range_reference.reference = 2;
        assert!(error_for(out_of_range_reference).contains("reference molecule 2 does not exist"));

        let mut wrong_entry_count = request(&["CC", "CC"], Some("CC"));
        wrong_entry_count.atoms = vec![None];
        assert!(error_for(wrong_entry_count).contains("atoms has 1 entries for 2 molecules"));

        let mut missing_atom = request(&["CC", "CC"], Some("CC"));
        missing_atom.atoms = vec![None, Some(vec![0, 5])];
        assert!(error_for(missing_atom).contains("atom 5 does not exist in molecule 1"));

        let mut repeated_atom = request(&["CC", "CC"], None);
        repeated_atom.atoms = vec![Some(vec![0, 1]), Some(vec![1, 1])];
        assert!(error_for(repeated_atom).contains("repeats an atom"));

        let mut unequal_lists = request(&["CCC", "CCC"], None);
        unequal_lists.atoms = vec![Some(vec![0, 1, 2]), Some(vec![0, 1])];
        assert!(error_for(unequal_lists).contains("same number of corresponding atoms"));

        let mut non_matching_list = request(&["CCO", "CCO"], Some("CO"));
        non_matching_list.atoms = vec![None, Some(vec![0, 1])];
        assert!(error_for(non_matching_list).contains("is not an occurrence of scaffold"));

        let mut partially_explicit = request(&["CC", "CC"], None);
        partially_explicit.atoms = vec![Some(vec![0, 1]), None];
        assert!(
            error_for(partially_explicit).contains("molecule 1 has no atom list and no scaffold")
        );
    }
}
