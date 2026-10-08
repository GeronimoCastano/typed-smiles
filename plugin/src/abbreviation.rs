//! Automatic abbreviation of terminal functional groups.
//!
//! Catalogue groups are recognized on the parsed molecular graph, never on the
//! SMILES text, so every atom order that encodes the same structure behaves the
//! same way. Contraction is a display transformation: the molecule used for
//! formulas and masses is never modified, the layout runs on a copy in which
//! each group is a single label atom, and the result is mapped back to the
//! original atom indices so atom references keep their meaning.

use crate::graph::{Atom, AtomChirality, BondDirection, BondOrder, BondStereo, MoleculeGraph};
use crate::layout::{compute_layout, implicit_h_count};
use crate::render::{AbbreviationGroupOutput, AbbreviationLabelOutput, AtomOutput, LayoutOutput};

/// One atom of a catalogue group, described outward from the attachment atom.
struct GroupAtom {
    element: &'static str,
    charge: i8,
    hydrogens: u8,
    substituents: &'static [GroupBond],
}

struct GroupBond {
    order: BondOrder,
    atom: GroupAtom,
}

pub(crate) struct GroupDefinition {
    name: &'static str,
    /// Label read with the rest of the molecule on its left, in `{label}`
    /// notation; `>` marks the glyph drawn over the attachment atom.
    label: &'static str,
    /// Label read with the rest of the molecule on its right.
    reversed_label: &'static str,
    /// Accepted structures rooted at the attachment atom. Each group attaches
    /// to the rest of the molecule through one single bond.
    forms: &'static [GroupAtom],
}

const fn terminal_atom(element: &'static str, charge: i8, hydrogens: u8) -> GroupAtom {
    GroupAtom {
        element,
        charge,
        hydrogens,
        substituents: &[],
    }
}

const fn bonded(order: BondOrder, atom: GroupAtom) -> GroupBond {
    GroupBond { order, atom }
}

const METHYL: GroupAtom = terminal_atom("C", 0, 3);
const FLUORINE: GroupAtom = terminal_atom("F", 0, 0);
const CARBONYL_OXYGEN: GroupAtom = terminal_atom("O", 0, 0);

/// The catalogue in priority order. When candidates overlap, the group listed
/// first wins, so the selection never depends on the order of the request.
/// Alkoxy groups precede acetyl so esters keep their carbonyl drawn.
const CATALOGUE: &[GroupDefinition] = &[
    GroupDefinition {
        name: "tBu",
        label: "tBu",
        reversed_label: "tBu",
        forms: &[GroupAtom {
            element: "C",
            charge: 0,
            hydrogens: 0,
            substituents: &[
                bonded(BondOrder::Single, METHYL),
                bonded(BondOrder::Single, METHYL),
                bonded(BondOrder::Single, METHYL),
            ],
        }],
    },
    GroupDefinition {
        name: "CF3",
        label: ">CF_3",
        reversed_label: "F_3>C",
        forms: &[GroupAtom {
            element: "C",
            charge: 0,
            hydrogens: 0,
            substituents: &[
                bonded(BondOrder::Single, FLUORINE),
                bonded(BondOrder::Single, FLUORINE),
                bonded(BondOrder::Single, FLUORINE),
            ],
        }],
    },
    GroupDefinition {
        name: "NO2",
        label: ">NO_2",
        reversed_label: "O_2>N",
        forms: &[
            GroupAtom {
                element: "N",
                charge: 1,
                hydrogens: 0,
                substituents: &[
                    bonded(BondOrder::Double, CARBONYL_OXYGEN),
                    bonded(BondOrder::Single, terminal_atom("O", -1, 0)),
                ],
            },
            GroupAtom {
                element: "N",
                charge: 0,
                hydrogens: 0,
                substituents: &[
                    bonded(BondOrder::Double, CARBONYL_OXYGEN),
                    bonded(BondOrder::Double, CARBONYL_OXYGEN),
                ],
            },
        ],
    },
    GroupDefinition {
        name: "CN",
        label: ">CN",
        reversed_label: "N>C",
        forms: &[GroupAtom {
            element: "C",
            charge: 0,
            hydrogens: 0,
            substituents: &[bonded(BondOrder::Triple, terminal_atom("N", 0, 0))],
        }],
    },
    GroupDefinition {
        name: "OEt",
        label: ">OEt",
        reversed_label: "Et>O",
        forms: &[GroupAtom {
            element: "O",
            charge: 0,
            hydrogens: 0,
            substituents: &[bonded(
                BondOrder::Single,
                GroupAtom {
                    element: "C",
                    charge: 0,
                    hydrogens: 2,
                    substituents: &[bonded(BondOrder::Single, METHYL)],
                },
            )],
        }],
    },
    GroupDefinition {
        name: "OMe",
        label: ">OMe",
        reversed_label: "Me>O",
        forms: &[GroupAtom {
            element: "O",
            charge: 0,
            hydrogens: 0,
            substituents: &[bonded(BondOrder::Single, METHYL)],
        }],
    },
    GroupDefinition {
        name: "Ac",
        label: "Ac",
        reversed_label: "Ac",
        forms: &[GroupAtom {
            element: "C",
            charge: 0,
            hydrogens: 0,
            substituents: &[
                bonded(BondOrder::Double, CARBONYL_OXYGEN),
                bonded(BondOrder::Single, METHYL),
            ],
        }],
    },
];

/// Catalogue names in priority order.
pub(crate) fn catalogue_names() -> Vec<&'static str> {
    CATALOGUE.iter().map(|definition| definition.name).collect()
}

/// Resolves a request — `all`, or comma-separated catalogue names — into
/// definitions in catalogue priority order.
pub(crate) fn requested_definitions(
    request: &str,
) -> Result<Vec<&'static GroupDefinition>, String> {
    if request.trim() == "all" {
        return Ok(CATALOGUE.iter().collect());
    }

    let mut requested_names: Vec<&str> = Vec::new();
    for raw_name in request.split(',') {
        let name = raw_name.trim();
        if !CATALOGUE.iter().any(|definition| definition.name == name) {
            return Err(format!(
                "unknown automatic abbreviation {name:?}; available abbreviations are {}",
                catalogue_names().join(", ")
            ));
        }
        if requested_names.contains(&name) {
            return Err(format!(
                "automatic abbreviation {name:?} is requested more than once"
            ));
        }
        requested_names.push(name);
    }

    Ok(CATALOGUE
        .iter()
        .filter(|definition| requested_names.contains(&definition.name))
        .collect())
}

/// A catalogue group found in the molecule, in original atom indices.
struct GroupMatch {
    definition: &'static GroupDefinition,
    attachment_atom: usize,
    external_atom: usize,
    /// Every group atom, including the attachment atom, ascending.
    atoms: Vec<usize>,
    /// Bond indices inside the group.
    bonds: Vec<usize>,
}

impl GroupMatch {
    fn hidden_atoms(&self) -> impl Iterator<Item = usize> + '_ {
        self.atoms
            .iter()
            .copied()
            .filter(|&atom_index| atom_index != self.attachment_atom)
    }

    fn conflicts_with(&self, other: &GroupMatch) -> bool {
        let shares_atoms = self
            .atoms
            .iter()
            .any(|atom_index| other.atoms.contains(atom_index));
        // Two labels bonded directly to each other would hide the whole bond
        // between them, so a group may not attach to another group.
        let attaches_to_other =
            other.atoms.contains(&self.external_atom) || self.atoms.contains(&other.external_atom);
        shares_atoms || attaches_to_other
    }
}

/// Groups to contract, chosen greedily in catalogue priority order and then by
/// ascending attachment and external atom indices.
fn select_groups(
    molecule: &MoleculeGraph,
    definitions: &[&'static GroupDefinition],
) -> Vec<GroupMatch> {
    let mut selected_groups: Vec<GroupMatch> = Vec::new();
    for definition in definitions {
        for candidate in find_group_matches(molecule, definition) {
            let conflicts = selected_groups
                .iter()
                .any(|selected_group| candidate.conflicts_with(selected_group));
            if !conflicts {
                selected_groups.push(candidate);
            }
        }
    }
    selected_groups
}

fn find_group_matches(
    molecule: &MoleculeGraph,
    definition: &'static GroupDefinition,
) -> Vec<GroupMatch> {
    let mut matches = Vec::new();
    for attachment_atom in 0..molecule.n_atoms() {
        let mut attachment_bonds = molecule.adj[attachment_atom].clone();
        attachment_bonds.sort_unstable();
        for (external_atom, attachment_bond) in attachment_bonds {
            if !is_attachment_bond(molecule, attachment_bond) {
                continue;
            }
            let group_match = definition
                .forms
                .iter()
                .find_map(|form| match_group_atom(molecule, form, attachment_atom, external_atom));
            if let Some(mut group_match) = group_match {
                group_match.atoms.sort_unstable();
                group_match.atoms.dedup();
                // A tree pattern must cover distinct atoms; anything else is a
                // ring the catalogue does not describe.
                if group_match.atoms.len() != group_match.bonds.len() + 1 {
                    continue;
                }
                matches.push(GroupMatch {
                    definition,
                    attachment_atom,
                    external_atom,
                    atoms: group_match.atoms,
                    bonds: group_match.bonds,
                });
            }
        }
    }
    matches
}

/// Atoms and bonds covered by a matched subtree.
struct SubtreeMatch {
    atoms: Vec<usize>,
    bonds: Vec<usize>,
}

/// Matches `pattern` at `atom_index`, whose only neighbor outside the pattern
/// is `parent_atom`.
fn match_group_atom(
    molecule: &MoleculeGraph,
    pattern: &GroupAtom,
    atom_index: usize,
    parent_atom: usize,
) -> Option<SubtreeMatch> {
    if !atom_matches(molecule, pattern, atom_index) {
        return None;
    }
    if molecule.adj[atom_index].len() != pattern.substituents.len() + 1 {
        return None;
    }
    let children: Vec<(usize, usize)> = molecule.adj[atom_index]
        .iter()
        .copied()
        .filter(|&(neighbor, _)| neighbor != parent_atom)
        .collect();
    let mut subtree = SubtreeMatch {
        atoms: vec![atom_index],
        bonds: Vec::new(),
    };
    let mut used_children = vec![false; children.len()];
    match_substituents(
        molecule,
        pattern.substituents,
        atom_index,
        &children,
        &mut used_children,
        &mut subtree,
    )
    .then_some(subtree)
}

/// Assigns each substituent pattern to a distinct child, backtracking when an
/// early assignment leaves a later pattern unmatched.
fn match_substituents(
    molecule: &MoleculeGraph,
    substituents: &[GroupBond],
    atom_index: usize,
    children: &[(usize, usize)],
    used_children: &mut [bool],
    subtree: &mut SubtreeMatch,
) -> bool {
    let Some((substituent, remaining_substituents)) = substituents.split_first() else {
        return true;
    };
    for (child_position, &(child_atom, child_bond)) in children.iter().enumerate() {
        if used_children[child_position]
            || !internal_bond_matches(molecule, substituent, child_bond)
        {
            continue;
        }
        let Some(child_subtree) =
            match_group_atom(molecule, &substituent.atom, child_atom, atom_index)
        else {
            continue;
        };

        let atom_count = subtree.atoms.len();
        let bond_count = subtree.bonds.len();
        used_children[child_position] = true;
        subtree.atoms.extend(child_subtree.atoms);
        subtree.bonds.push(child_bond);
        subtree.bonds.extend(child_subtree.bonds);
        if match_substituents(
            molecule,
            remaining_substituents,
            atom_index,
            children,
            used_children,
            subtree,
        ) {
            return true;
        }
        used_children[child_position] = false;
        subtree.atoms.truncate(atom_count);
        subtree.bonds.truncate(bond_count);
    }
    false
}

/// A group atom must carry exactly the meaning its label expresses. Isotopes,
/// stereo marks, unexpected charges, aromaticity, bracket hydrogens, and
/// manual labels all keep the atom expanded.
fn atom_matches(molecule: &MoleculeGraph, pattern: &GroupAtom, atom_index: usize) -> bool {
    let atom = &molecule.atoms[atom_index];
    let has_bracket_hydrogens = atom.has_explicit_h && atom.hcount > 0;
    let hydrogen_count = atom.hcount + implicit_h_count(molecule, atom_index);
    atom.symbol == pattern.element
        && atom.abbrev.is_empty()
        && !atom.aromatic
        && atom.isotope.is_none()
        && atom.charge == pattern.charge
        && atom.chirality == AtomChirality::None
        && !has_bracket_hydrogens
        && hydrogen_count == pattern.hydrogens
}

fn internal_bond_matches(
    molecule: &MoleculeGraph,
    substituent: &GroupBond,
    bond_index: usize,
) -> bool {
    let bond = &molecule.bonds[bond_index];
    bond.order == substituent.order && is_unmarked(molecule, bond_index)
}

/// The bond joining a group to the rest of the molecule stays drawn, so it may
/// carry stereo or direction marks; it only needs to be an ordinary single bond.
fn is_attachment_bond(molecule: &MoleculeGraph, bond_index: usize) -> bool {
    let bond = &molecule.bonds[bond_index];
    bond.order == BondOrder::Single && !bond.aromatic
}

fn is_unmarked(molecule: &MoleculeGraph, bond_index: usize) -> bool {
    let bond = &molecule.bonds[bond_index];
    !bond.aromatic
        && bond.stereo == BondStereo::None
        && bond.direction == BondDirection::None
        && !bond.forced_stereo
        && !bond.curl
}

/// Lays out `molecule` with every selected catalogue group drawn as one label.
/// The returned layout uses the original atom indices; atoms hidden inside a
/// label are marked `contracted` and share the label's position.
pub(crate) fn layout_with_abbreviations(
    molecule: &MoleculeGraph,
    definitions: &[&'static GroupDefinition],
) -> Result<LayoutOutput, String> {
    let groups = select_groups(molecule, definitions);
    let kept_atoms = kept_atom_indices(molecule, &groups);
    let display_molecule = contracted_molecule(molecule, &groups, &kept_atoms)?;
    let display_layout = compute_layout(&display_molecule)?;
    expand_to_original_atoms(molecule, &groups, &kept_atoms, display_layout)
}

fn kept_atom_indices(molecule: &MoleculeGraph, groups: &[GroupMatch]) -> Vec<usize> {
    let mut is_hidden = vec![false; molecule.n_atoms()];
    for group in groups {
        for hidden_atom in group.hidden_atoms() {
            is_hidden[hidden_atom] = true;
        }
    }
    (0..molecule.n_atoms())
        .filter(|&atom_index| !is_hidden[atom_index])
        .collect()
}

/// The display graph: hidden group atoms removed and each attachment atom
/// replaced by a label atom, exactly as if the group had been written as a
/// manual `{label}`.
fn contracted_molecule(
    molecule: &MoleculeGraph,
    groups: &[GroupMatch],
    kept_atoms: &[usize],
) -> Result<MoleculeGraph, String> {
    let mut display_molecule = molecule.induced_subgraph(kept_atoms);
    for group in groups {
        let display_index = kept_atoms
            .binary_search(&group.attachment_atom)
            .map_err(|_| {
                format!(
                    "abbreviation {} lost its attachment atom",
                    group.definition.name
                )
            })?;
        let attachment_symbol = &molecule.atoms[group.attachment_atom].symbol;
        display_molecule.atoms[display_index] = label_atom(group.definition, attachment_symbol)?;
    }
    Ok(display_molecule)
}

fn label_atom(definition: &GroupDefinition, attachment_symbol: &str) -> Result<Atom, String> {
    let label = parsed_label(definition.label)?;
    // Heteroatom-anchored labels take the attachment element's color, like the
    // atom label they replace.
    let style = if attachment_symbol == "C" {
        String::new()
    } else {
        attachment_symbol.to_string()
    };
    Ok(Atom {
        symbol: "*".to_string(),
        aromatic: false,
        hcount: 0,
        has_explicit_h: true,
        isotope: None,
        charge: 0,
        chirality: AtomChirality::None,
        abbrev: label.text,
        abbrev_style: style,
        abbrev_anchor: label.anchor,
        abbrev_anchor_len: label.anchor_len,
        abbrev_lone_pairs: 0,
        abbrev_offset_x: 0.0,
        abbrev_offset_y: 0.0,
    })
}

fn parsed_label(raw_label: &str) -> Result<AbbreviationLabelOutput, String> {
    let label = crate::parse_abbreviation_label(raw_label)?;
    Ok(AbbreviationLabelOutput {
        text: label.text,
        anchor: label.anchor,
        anchor_len: label.anchor_len,
    })
}

/// Maps a display-graph layout back to original atom indices. Virtual
/// hydrogens keep following the real atoms.
fn expand_to_original_atoms(
    molecule: &MoleculeGraph,
    groups: &[GroupMatch],
    kept_atoms: &[usize],
    display_layout: LayoutOutput,
) -> Result<LayoutOutput, String> {
    let original_count = molecule.n_atoms();
    let kept_count = kept_atoms.len();
    let original_index = |display_index: usize| {
        if display_index < kept_count {
            kept_atoms[display_index]
        } else {
            original_count + (display_index - kept_count)
        }
    };

    let mut display_index_of = vec![None; original_count];
    for (display_index, &atom_index) in kept_atoms.iter().enumerate() {
        display_index_of[atom_index] = Some(display_index);
    }
    let mut group_of = vec![None; original_count];
    for (group_index, group) in groups.iter().enumerate() {
        for &atom_index in &group.atoms {
            group_of[atom_index] = Some(group_index);
        }
    }

    let mut atoms = Vec::with_capacity(display_layout.atoms.len() + original_count - kept_count);
    for atom_index in 0..original_count {
        let atom_output = match display_index_of[atom_index] {
            Some(display_index) => AtomOutput {
                abbreviation_group: group_of[atom_index],
                ..display_layout.atoms[display_index].clone()
            },
            None => {
                let group_index = group_of[atom_index].ok_or_else(|| {
                    format!("atom {atom_index} vanished outside any abbreviation")
                })?;
                let attachment_atom = groups[group_index].attachment_atom;
                let label_display_index = display_index_of[attachment_atom].ok_or_else(|| {
                    format!("abbreviation attachment atom {attachment_atom} was not kept")
                })?;
                let label_position = display_layout.atoms[label_display_index].pos;
                contracted_atom_output(molecule, atom_index, group_index, label_position)
            }
        };
        atoms.push(atom_output);
    }
    atoms.extend(display_layout.atoms[kept_count..].iter().cloned());

    let bonds = display_layout
        .bonds
        .into_iter()
        .map(|bond_output| {
            let from = original_index(bond_output.from);
            let to = original_index(bond_output.to);
            crate::render::BondOutput {
                from,
                to,
                ..bond_output
            }
        })
        .collect();

    let abbreviation_groups = groups
        .iter()
        .map(|group| group_output(molecule, group))
        .collect::<Result<Vec<_>, String>>()?;

    Ok(LayoutOutput {
        atoms,
        bonds,
        abbreviation_groups,
        ..display_layout
    })
}

fn contracted_atom_output(
    molecule: &MoleculeGraph,
    atom_index: usize,
    group_index: usize,
    label_position: crate::render::Vec2,
) -> AtomOutput {
    let atom = &molecule.atoms[atom_index];
    AtomOutput {
        symbol: atom.symbol.clone(),
        pos: label_position,
        hcount: atom.hcount,
        implicit_h: implicit_h_count(molecule, atom_index),
        charge: atom.charge,
        isotope: atom.isotope.unwrap_or(0),
        lone_pairs: 0,
        lone_pair_dirs: Vec::new(),
        abbrev: String::new(),
        abbrev_style: String::new(),
        abbrev_anchor: 0,
        abbrev_anchor_len: 0,
        abbrev_offset_x: 0.0,
        abbrev_offset_y: 0.0,
        chirality: atom.chirality.as_str().to_string(),
        stereo_h: "none".to_string(),
        stereo_h_dir: crate::render::Vec2::default(),
        virtual_h: false,
        contracted: true,
        abbreviation_group: Some(group_index),
    }
}

fn group_output(
    molecule: &MoleculeGraph,
    group: &GroupMatch,
) -> Result<AbbreviationGroupOutput, String> {
    Ok(AbbreviationGroupOutput {
        name: group.definition.name.to_string(),
        attachment_atom: group.attachment_atom,
        external_atom: group.external_atom,
        atoms: group.atoms.clone(),
        bonds: group
            .bonds
            .iter()
            .map(|&bond_index| {
                let bond = &molecule.bonds[bond_index];
                [bond.from, bond.to]
            })
            .collect(),
        label: parsed_label(group.definition.label)?,
        reversed_label: parsed_label(group.definition.reversed_label)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{layout_abbreviated_native, layout_native, parse_molecule};

    fn abbreviated(smiles: &str, request: &str) -> LayoutOutput {
        layout_abbreviated_native(smiles, request)
            .unwrap_or_else(|error| panic!("{smiles} with {request}: {error}"))
    }

    fn group_names(layout_output: &LayoutOutput) -> Vec<&str> {
        layout_output
            .abbreviation_groups
            .iter()
            .map(|group| group.name.as_str())
            .collect()
    }

    fn only_group(layout_output: &LayoutOutput) -> &AbbreviationGroupOutput {
        assert_eq!(
            layout_output.abbreviation_groups.len(),
            1,
            "{:?}",
            group_names(layout_output)
        );
        &layout_output.abbreviation_groups[0]
    }

    #[test]
    fn every_catalogue_group_is_recognized_on_a_benzene_ring() {
        let cases = [
            ("CC(C)(C)c1ccccc1", "tBu", 1, 4),
            ("FC(F)(F)c1ccccc1", "CF3", 1, 4),
            ("[O-][N+](=O)c1ccccc1", "NO2", 1, 3),
            ("O=N(=O)c1ccccc1", "NO2", 1, 3),
            ("N#Cc1ccccc1", "CN", 1, 2),
            ("CCOc1ccccc1", "OEt", 2, 3),
            ("COc1ccccc1", "OMe", 1, 2),
            ("CC(=O)Nc1ccccc1", "Ac", 1, 3),
        ];
        for (smiles, name, attachment_atom, group_size) in cases {
            let layout_output = abbreviated(smiles, name);
            let group = only_group(&layout_output);
            assert_eq!(group.name, name, "{smiles}");
            assert_eq!(group.attachment_atom, attachment_atom, "{smiles}");
            assert_eq!(group.atoms.len(), group_size, "{smiles}");
            assert_eq!(group.bonds.len(), group_size - 1, "{smiles}");
        }
    }

    #[test]
    fn equivalent_atom_orders_contract_the_same_group() {
        for smiles in [
            "FC(F)(F)c1ccccc1",
            "c1ccccc1C(F)(F)F",
            "c1ccc(cc1)C(F)(F)F",
            "C(F)(c1ccccc1)(F)F",
        ] {
            let layout_output = abbreviated(smiles, "CF3");
            let group = only_group(&layout_output);
            let molecule = parse_molecule(smiles).unwrap();
            let hidden_symbols: Vec<&str> = group
                .atoms
                .iter()
                .filter(|&&atom_index| atom_index != group.attachment_atom)
                .map(|&atom_index| molecule.atoms[atom_index].symbol.as_str())
                .collect();
            assert_eq!(hidden_symbols, ["F", "F", "F"], "{smiles}");
            assert_eq!(molecule.atoms[group.attachment_atom].symbol, "C");
            assert!(molecule.atoms[group.external_atom].aromatic, "{smiles}");
        }
    }

    #[test]
    fn contraction_keeps_original_atom_indices_and_virtual_hydrogens() {
        let smiles = "COc1ccc(cc1)[NH3+]";
        let full_layout = layout_native(smiles).unwrap();
        let layout_output = abbreviated(smiles, "OMe");
        assert_eq!(layout_output.atoms.len(), full_layout.atoms.len());

        for (atom_index, (abbreviated_atom, full_atom)) in layout_output
            .atoms
            .iter()
            .zip(&full_layout.atoms)
            .enumerate()
        {
            assert_eq!(
                abbreviated_atom.virtual_h, full_atom.virtual_h,
                "atom {atom_index}"
            );
            if abbreviated_atom.abbreviation_group.is_none() {
                assert_eq!(
                    abbreviated_atom.symbol, full_atom.symbol,
                    "atom {atom_index}"
                );
            }
        }

        let label_position = layout_output.atoms[1].pos;
        assert!(layout_output.atoms[0].contracted);
        assert_eq!(layout_output.atoms[0].abbreviation_group, Some(0));
        assert!(layout_output.atoms[0].pos.dist(label_position) < 1e-12);
        assert!(!layout_output.atoms[1].contracted);
        assert_eq!(layout_output.atoms[1].abbreviation_group, Some(0));
        assert_eq!(layout_output.atoms[1].abbrev, "OMe");
        assert_eq!(layout_output.atoms[1].abbrev_style, "O");
    }

    #[test]
    fn hidden_bonds_leave_the_drawn_bond_list() {
        let layout_output = abbreviated("FC(F)(F)c1ccc(cc1)OC", "CF3,OMe");
        let contracted_atoms: Vec<usize> = (0..layout_output.atoms.len())
            .filter(|&atom_index| layout_output.atoms[atom_index].contracted)
            .collect();
        assert_eq!(contracted_atoms, [0, 2, 3, 11]);
        for bond_output in &layout_output.bonds {
            assert!(!contracted_atoms.contains(&bond_output.from));
            assert!(!contracted_atoms.contains(&bond_output.to));
        }
        let attachment_bond_drawn = layout_output.bonds.iter().any(|bond_output| {
            (bond_output.from, bond_output.to) == (1, 4)
                || (bond_output.from, bond_output.to) == (4, 1)
        });
        assert!(attachment_bond_drawn);
    }

    #[test]
    fn groups_whose_meaning_a_label_would_hide_stay_expanded() {
        let expanded_cases = [
            ("FC(F)c1ccccc1", "CF3"),
            ("[13CH3]Oc1ccccc1", "OMe"),
            ("[CH3]Oc1ccccc1", "OMe"),
            ("C[O+](C)c1ccccc1", "OMe"),
            ("CO", "OMe"),
            ("C!wOc1ccccc1", "OMe"),
            ("C1CCOC1", "OMe"),
            ("[O-][N+]([O-])c1ccccc1", "NO2"),
            ("[N+]#[C-]", "CN"),
            ("{OMe}c1ccccc1", "OMe"),
            ("C1=CC=CC=C1", "all"),
            ("CC(C)(C)[C@H](F)Cl", "OMe"),
        ];
        for (smiles, request) in expanded_cases {
            let layout_output = abbreviated(smiles, request);
            assert!(
                layout_output.abbreviation_groups.is_empty(),
                "{smiles}: {:?}",
                group_names(&layout_output)
            );
            assert!(
                layout_output.atoms.iter().all(|atom| !atom.contracted),
                "{smiles}"
            );
        }
    }

    #[test]
    fn marked_attachment_bonds_keep_their_marks() {
        let layout_output = abbreviated("F/C=C/C(F)(F)F", "CF3");
        assert_eq!(group_names(&layout_output), ["CF3"]);
        let directional_bond = layout_output
            .bonds
            .iter()
            .find(|bond_output| (bond_output.from, bond_output.to) == (2, 3))
            .expect("attachment bond stays drawn");
        assert_eq!(directional_bond.direction, "up");
    }

    #[test]
    fn overlapping_candidates_follow_catalogue_priority() {
        let ether = abbreviated("CCOC", "OMe,OEt");
        assert_eq!(group_names(&ether), ["OEt"]);
        assert_eq!(group_names(&abbreviated("CCOC", "OEt,OMe")), ["OEt"]);

        let methyl_acetate = abbreviated("COC(C)=O", "Ac,OMe");
        assert_eq!(group_names(&methyl_acetate), ["OMe"]);

        let aryl_acetate = abbreviated("CC(=O)Oc1ccccc1", "Ac,OMe");
        assert_eq!(group_names(&aryl_acetate), ["Ac"]);
    }

    #[test]
    fn groups_never_attach_to_another_group() {
        let hexafluoroethane = abbreviated("FC(F)(F)C(F)(F)F", "CF3");
        let group = only_group(&hexafluoroethane);
        assert_eq!(group.attachment_atom, 1);
        assert_eq!(group.external_atom, 4);
        assert!(!hexafluoroethane.atoms[4].contracted);

        let neopentane = abbreviated("CC(C)(C)C", "tBu");
        let group = only_group(&neopentane);
        assert_eq!(group.atoms, [1, 2, 3, 4]);
        assert_eq!(group.external_atom, 0);
    }

    #[test]
    fn several_groups_contract_together() {
        let layout_output = abbreviated("COc1ccc(cc1)C(F)(F)F", "all");
        assert_eq!(group_names(&layout_output), ["CF3", "OMe"]);
        let mut covered_atoms: Vec<usize> = layout_output
            .abbreviation_groups
            .iter()
            .flat_map(|group| group.atoms.iter().copied())
            .collect();
        covered_atoms.sort_unstable();
        assert_eq!(covered_atoms, [0, 1, 8, 9, 10, 11]);
    }

    #[test]
    fn requests_name_catalogue_groups_once() {
        assert!(requested_definitions("all").is_ok());
        assert_eq!(
            requested_definitions("OMe, CF3")
                .unwrap()
                .iter()
                .map(|definition| definition.name)
                .collect::<Vec<_>>(),
            ["CF3", "OMe"]
        );
        let unknown = requested_definitions("Boc").err().unwrap();
        assert!(
            unknown.contains("unknown automatic abbreviation \"Boc\""),
            "{unknown}"
        );
        assert!(
            unknown.contains("tBu, CF3, NO2, CN, OEt, OMe, Ac"),
            "{unknown}"
        );
        let repeated = requested_definitions("OMe,OMe").err().unwrap();
        assert!(repeated.contains("more than once"), "{repeated}");
        assert!(requested_definitions("").is_err());
    }

    #[test]
    fn catalogue_labels_anchor_the_attachment_element() {
        let attachment_elements = ["C", "C", "N", "C", "O", "O", "C"];
        for (definition, element) in CATALOGUE.iter().zip(attachment_elements) {
            for raw_label in [definition.label, definition.reversed_label] {
                let label = parsed_label(raw_label).unwrap();
                if label.anchor_len == 0 {
                    continue;
                }
                let glyph: String = label
                    .text
                    .chars()
                    .skip(label.anchor)
                    .take(label.anchor_len)
                    .collect();
                assert_eq!(glyph, element, "{}", definition.name);
            }
        }
    }

    /// Signed volume of the depicted tetrahedron around `center`, taking
    /// neighbors in SMILES order and lifting wedged bonds out of the plane.
    fn depicted_handedness(layout_output: &LayoutOutput, smiles: &str, center: usize) -> f64 {
        let molecule = parse_molecule(smiles).unwrap();
        let center_position = layout_output.atoms[center].pos;
        let mut directions = [[0.0; 3]; 4];
        for (slot, &bond_index) in molecule.neighbor_bonds[center].iter().enumerate() {
            let bond = &molecule.bonds[bond_index];
            let neighbor = if bond.from == center {
                bond.to
            } else {
                bond.from
            };
            let neighbor_position = layout_output.atoms[neighbor].pos;
            let distance = neighbor_position.dist(center_position);
            let depth = layout_output
                .bonds
                .iter()
                .find(|bond_output| {
                    (bond_output.from, bond_output.to) == (center, neighbor)
                        || (bond_output.from, bond_output.to) == (neighbor, center)
                })
                .map(|bond_output| match bond_output.stereo.as_str() {
                    "wedge_up" => 1.0,
                    "wedge_down" => -1.0,
                    _ => 0.0,
                })
                .unwrap();
            directions[slot] = [
                (neighbor_position.x - center_position.x) / distance,
                (neighbor_position.y - center_position.y) / distance,
                depth,
            ];
        }
        crate::layout::signed_volume(&directions)
    }

    #[test]
    fn stereocenters_next_to_a_label_keep_their_configuration() {
        for smiles in [
            "F[C@](Cl)(Br)C(F)(F)F",
            "F[C@@](Cl)(Br)C(F)(F)F",
            "Cl[C@](F)(OC)C#N",
        ] {
            let full_handedness = depicted_handedness(&layout_native(smiles).unwrap(), smiles, 1);
            let abbreviated_handedness =
                depicted_handedness(&abbreviated(smiles, "all"), smiles, 1);
            assert!(full_handedness.abs() > 1e-6, "{smiles}");
            assert_eq!(
                full_handedness.signum(),
                abbreviated_handedness.signum(),
                "{smiles}"
            );
        }
    }
}
