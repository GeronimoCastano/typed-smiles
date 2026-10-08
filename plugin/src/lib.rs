mod abbreviation;
mod alignment;
#[cfg(test)]
mod conformance_tests;
mod error;
mod geometry;
mod graph;
mod kekulize;
mod label;
mod layout;
#[cfg(test)]
mod layout_quality;
mod layout_relaxation;
mod macrocycles;
mod parser;
#[cfg(test)]
mod performance_tests;
mod render;
mod ring_system_layout;
mod ring_templates;
mod rings;
mod stereo;
mod substructure;

pub use alignment::{AlignmentRequest, MoleculeAlignment};
pub use render::LayoutOutput;
pub use substructure::SubstructureMatch;

use graph::MoleculeGraph;
use layout::compute_layout;
use ptable::Element;
use std::collections::BTreeMap;

fn parse_molecule(smiles: &str) -> Result<MoleculeGraph, String> {
    parser::parse_smiles(smiles).map_err(|error| format!("invalid SMILES {smiles:?}: {error}"))
}

/// Looks up an element by its exact symbol, such as "Cl". The table search in
/// `ptable` assumes well-formed symbols, so arbitrary user text is compared
/// against every element instead.
pub(crate) fn element_from_symbol(symbol: &str) -> Option<Element> {
    (1..=118)
        .filter_map(Element::from_atomic_number)
        .find(|element| element.get_symbol() == symbol)
}

fn atomic_mass(atom: &graph::Atom) -> Result<f64, String> {
    if !atom.abbrev.is_empty() {
        return Err(format!(
            "cannot compute molecular weight: abbreviation {{{}}} has no defined composition",
            atom.abbrev
        ));
    }
    if atom.symbol == "*" {
        return Err("cannot compute molecular weight: wildcard atom `*` has no mass".to_string());
    }
    if let Some(isotope) = atom.isotope {
        return Err(format!(
            "cannot compute molecular weight: isotope [{isotope}{}] needs a nuclide mass, \
             not a standard atomic weight",
            atom.symbol
        ));
    }

    element_from_symbol(&atom.symbol)
        .map(|element| element.get_atomic_mass() as f64)
        .ok_or_else(|| {
            format!(
                "cannot compute molecular weight: unknown element {}",
                atom.symbol
            )
        })
}

fn hydrogen_count(molecule: &MoleculeGraph, atom_index: usize) -> u8 {
    molecule.atoms[atom_index].hcount + layout::implicit_h_count(molecule, atom_index)
}

fn compute_molecular_weight(molecule: &MoleculeGraph) -> Result<f64, String> {
    let hydrogen_mass = Element::Hydrogen.get_atomic_mass() as f64;
    let mut molecular_weight = 0.0;

    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        molecular_weight += atomic_mass(atom)?;
        molecular_weight += f64::from(hydrogen_count(molecule, atom_index)) * hydrogen_mass;
    }

    Ok(molecular_weight)
}

fn compute_molecular_formula(molecule: &MoleculeGraph) -> Result<String, String> {
    let mut element_counts = BTreeMap::<String, u32>::new();
    let mut formal_charge = 0i16;

    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        if !atom.abbrev.is_empty() {
            return Err(format!(
                "cannot compute molecular formula: abbreviation {{{}}} has no defined composition",
                atom.abbrev
            ));
        }
        if atom.symbol == "*" {
            return Err(
                "cannot compute molecular formula: wildcard atom `*` has no element".to_string(),
            );
        }
        if let Some(isotope) = atom.isotope {
            return Err(format!(
                "cannot compute molecular formula: isotope [{isotope}{}] needs isotope-aware formula output",
                atom.symbol
            ));
        }
        if element_from_symbol(&atom.symbol).is_none() {
            return Err(format!(
                "cannot compute molecular formula: unknown element {}",
                atom.symbol
            ));
        }

        *element_counts.entry(atom.symbol.clone()).or_default() += 1;
        let hydrogen_count =
            u32::from(atom.hcount) + u32::from(layout::implicit_h_count(molecule, atom_index));
        if hydrogen_count > 0 {
            *element_counts.entry("H".to_string()).or_default() += hydrogen_count;
        }
        formal_charge += i16::from(atom.charge);
    }

    let mut ordered_symbols = Vec::new();
    if element_counts.contains_key("C") {
        ordered_symbols.push("C".to_string());
        if element_counts.contains_key("H") {
            ordered_symbols.push("H".to_string());
        }
        ordered_symbols.extend(
            element_counts
                .keys()
                .filter(|symbol| symbol.as_str() != "C" && symbol.as_str() != "H")
                .cloned(),
        );
    } else {
        ordered_symbols.extend(element_counts.keys().cloned());
    }

    let mut formula = String::new();
    for symbol in ordered_symbols {
        formula.push_str(&symbol);
        let count = element_counts
            .get(&symbol)
            .expect("formula symbol came from the element-count map");
        if *count > 1 {
            formula.push_str(&count.to_string());
        }
    }

    if formal_charge != 0 {
        formula.push('^');
        if formal_charge.abs() != 1 {
            formula.push_str(&formal_charge.abs().to_string());
        }
        formula.push(if formal_charge > 0 { '+' } else { '-' });
    }

    Ok(formula)
}

// ── WASM / Typst plugin entrypoint ──────────────────────────────────────────

#[cfg(target_arch = "wasm32")]
mod wasm_entrypoint {
    use super::*;
    use wasm_minimal_protocol::*;

    initiate_protocol!();

    /// Called from Typst as `smiles-plugin.layout(bytes(smiles-str))`.
    /// Returns JSON-encoded `LayoutOutput`.
    /// Accepts extended SMILES with `{label}` abbreviation syntax.
    #[wasm_func]
    pub fn layout(smiles: &[u8]) -> Result<Vec<u8>, String> {
        let smiles =
            core::str::from_utf8(smiles).map_err(|error| format!("UTF-8 error: {error}"))?;
        let molecule = parse_molecule(smiles)?;
        let layout = compute_layout(&molecule)?;
        serde_json::to_vec(&layout).map_err(|error| format!("JSON error: {error}"))
    }

    /// Called from Typst as `smiles-plugin.layout_abbreviated(bytes(smiles-str),
    /// bytes(request))`, where the request is `all` or comma-separated
    /// catalogue names. Returns JSON-encoded `LayoutOutput` in original atom
    /// indices, with matched terminal groups drawn as labels.
    #[wasm_func]
    pub fn layout_abbreviated(smiles: &[u8], request: &[u8]) -> Result<Vec<u8>, String> {
        let smiles =
            core::str::from_utf8(smiles).map_err(|error| format!("UTF-8 error: {error}"))?;
        let request =
            core::str::from_utf8(request).map_err(|error| format!("UTF-8 error: {error}"))?;
        let definitions = abbreviation::requested_definitions(request)?;
        let molecule = parse_molecule(smiles)?;
        let layout = abbreviation::layout_with_abbreviations(&molecule, &definitions)?;
        serde_json::to_vec(&layout).map_err(|error| format!("JSON error: {error}"))
    }

    /// Catalogue names accepted by `layout_abbreviated`, as a JSON array in
    /// priority order.
    #[wasm_func]
    pub fn abbreviation_names() -> Result<Vec<u8>, String> {
        serde_json::to_vec(&abbreviation::catalogue_names())
            .map_err(|error| format!("JSON error: {error}"))
    }

    #[wasm_func]
    pub fn substructure_matches(smiles: &[u8], pattern: &[u8]) -> Result<Vec<u8>, String> {
        let smiles =
            core::str::from_utf8(smiles).map_err(|error| format!("UTF-8 error: {error}"))?;
        let pattern =
            core::str::from_utf8(pattern).map_err(|error| format!("UTF-8 error: {error}"))?;
        let molecule = parse_molecule(smiles)?;
        let matches = substructure::find_matches(&molecule, pattern)?;
        serde_json::to_vec(&matches).map_err(|error| format!("JSON error: {error}"))
    }

    /// Called from Typst with a JSON-encoded `AlignmentRequest`.
    /// Returns a JSON array of `MoleculeAlignment`, one per molecule.
    #[wasm_func]
    pub fn align_molecules(request: &[u8]) -> Result<Vec<u8>, String> {
        let request: AlignmentRequest = serde_json::from_slice(request)
            .map_err(|error| format!("alignment request error: {error}"))?;
        let alignments = alignment::align_molecules(&request)
            .map_err(|error| format!("typed-smiles: align-molecules: {error}"))?;
        serde_json::to_vec(&alignments).map_err(|error| format!("JSON error: {error}"))
    }

    /// Called from Typst as `smiles-plugin.mol_weight(bytes(smiles-str))`.
    /// Returns the molecular weight in g/mol as a JSON number.
    #[wasm_func]
    pub fn mol_weight(smiles: &[u8]) -> Result<Vec<u8>, String> {
        let smiles =
            core::str::from_utf8(smiles).map_err(|error| format!("UTF-8 error: {error}"))?;
        let molecule = parse_molecule(smiles)?;
        let molecular_weight = compute_molecular_weight(&molecule)?;
        serde_json::to_vec(&molecular_weight).map_err(|error| format!("JSON error: {error}"))
    }

    /// Called from Typst as `smiles-plugin.mol_formula(bytes(smiles-str))`.
    /// Returns a Hill-ordered molecular formula string with a net-charge suffix.
    #[wasm_func]
    pub fn mol_formula(smiles: &[u8]) -> Result<Vec<u8>, String> {
        let smiles =
            core::str::from_utf8(smiles).map_err(|error| format!("UTF-8 error: {error}"))?;
        let molecule = parse_molecule(smiles)?;
        let formula = compute_molecular_formula(&molecule)?;
        serde_json::to_vec(&formula).map_err(|error| format!("JSON error: {error}"))
    }
}

// ── Native entrypoint for tests / CLI ───────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
pub fn layout_native(smiles: &str) -> Result<LayoutOutput, String> {
    let molecule = parse_molecule(smiles)?;
    compute_layout(&molecule)
}

#[cfg(not(target_arch = "wasm32"))]
pub fn align_molecules_native(
    request: &AlignmentRequest,
) -> Result<Vec<MoleculeAlignment>, String> {
    alignment::align_molecules(request)
}

#[cfg(not(target_arch = "wasm32"))]
pub fn layout_abbreviated_native(smiles: &str, request: &str) -> Result<LayoutOutput, String> {
    let definitions = abbreviation::requested_definitions(request)?;
    let molecule = parse_molecule(smiles)?;
    abbreviation::layout_with_abbreviations(&molecule, &definitions)
}

#[cfg(not(target_arch = "wasm32"))]
pub fn mol_weight_native(smiles: &str) -> Result<f64, String> {
    let molecule = parse_molecule(smiles)?;
    compute_molecular_weight(&molecule)
}

#[cfg(not(target_arch = "wasm32"))]
pub fn mol_formula_native(smiles: &str) -> Result<String, String> {
    let molecule = parse_molecule(smiles)?;
    compute_molecular_formula(&molecule)
}

#[cfg(not(target_arch = "wasm32"))]
pub fn substructure_matches_native(
    smiles: &str,
    pattern: &str,
) -> Result<Vec<SubstructureMatch>, String> {
    substructure::find_matches(&parse_molecule(smiles)?, pattern)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Smallest distance between any two distinct non-virtual atoms.
    fn min_atom_distance(layout_output: &LayoutOutput) -> f64 {
        let mut min = f64::INFINITY;
        for i in 0..layout_output.atoms.len() {
            for j in (i + 1)..layout_output.atoms.len() {
                if layout_output.atoms[i].virtual_h || layout_output.atoms[j].virtual_h {
                    continue;
                }
                min = min.min(layout_output.atoms[i].pos.dist(layout_output.atoms[j].pos));
            }
        }
        min
    }

    #[test]
    fn ortho_ring_substituents_do_not_collide() {
        // Aspirin: the acetyl C=O and the carboxyl OH grow from ortho ring
        // positions toward each other and must not land on the same point.
        let layout_output = layout_native("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap();
        assert!(
            min_atom_distance(&layout_output) > 0.5,
            "atoms overlap: min distance {}",
            min_atom_distance(&layout_output)
        );

        // The carboxyl carbon (atom 10: ring C9, =O 11, OH 12) must keep its
        // textbook trigonal geometry: the conflict is resolved by flipping
        // the acetyl branch aside, not by squeezing the carboxyl bonds into a
        // narrow fan or a straight line.
        let carboxyl_center = layout_output.atoms[10].pos;
        let angles: Vec<f64> = [9, 11, 12]
            .iter()
            .map(|&neighbor_index| {
                let neighbor_position = layout_output.atoms[neighbor_index].pos;
                (neighbor_position.y - carboxyl_center.y)
                    .atan2(neighbor_position.x - carboxyl_center.x)
            })
            .collect();
        for first_angle_index in 0..angles.len() {
            for second_angle_index in (first_angle_index + 1)..angles.len() {
                let mut delta = (angles[first_angle_index] - angles[second_angle_index]).abs();
                if delta > std::f64::consts::PI {
                    delta = 2.0 * std::f64::consts::PI - delta;
                }
                let deg = delta.to_degrees();
                assert!(
                    (deg - 120.0).abs() < 1.0,
                    "carboxyl bond angle {deg:.1} degrees is not the ideal 120"
                );
            }
        }
    }

    #[test]
    fn crowded_branch_chains_do_not_collide() {
        // Two ortho substituent chains long enough to sweep past each other.
        let layout_output = layout_native("CCCC1=CC=CC=C1CCC").unwrap();
        assert!(min_atom_distance(&layout_output) > 0.5);
    }

    #[test]
    fn spiro_rings_share_one_atom_without_overlap() {
        // 1,6-dioxaspiro[4.4]nonane: two five-membered rings joined at a single
        // spiro atom. Each ring must close as a regular pentagon rather than
        // unravel into a chain, so no atoms collide.
        let layout_output = layout_native("CC[C@H](O1)CC[C@@]12CCCO2").unwrap();
        assert!(
            min_atom_distance(&layout_output) > 0.5,
            "spiro atoms overlap: min distance {}",
            min_atom_distance(&layout_output)
        );
        // The two ring oxygens (atoms 3 and 10) both neighbor the spiro carbon
        // (atom 6) and must sit one bond length away from it.
        let spiro = layout_output.atoms[6].pos;
        for oxygen_index in [3usize, 10] {
            let distance = layout_output.atoms[oxygen_index].pos.dist(spiro);
            assert!(
                (distance - 1.0).abs() < 0.05,
                "ring O-spiro bond length {} for atom {}",
                distance,
                oxygen_index
            );
        }
    }

    #[test]
    fn benzene_kekule() {
        let layout_output = layout_native("C1=CC=CC=C1").expect("benzene layout failed");
        assert_eq!(layout_output.atoms.len(), 6);
        assert_eq!(layout_output.bonds.len(), 6);
    }

    #[test]
    fn ethanol() {
        // CCO: 3 real atoms + 1 virtual H for terminal O (implicit OH)
        let layout_output = layout_native("CCO").expect("ethanol layout failed");
        assert_eq!(layout_output.atoms.len(), 4);
        assert!(layout_output.atoms[3].virtual_h);
    }

    #[test]
    fn implicit_h_counts_ethanol() {
        // implicit_h counts for real atoms only (virtual H has implicit_h = 0)
        let layout_output = layout_native("CCO").expect("ethanol layout failed");
        let counts: Vec<u8> = layout_output
            .atoms
            .iter()
            .map(|atom| atom.implicit_h)
            .collect();
        assert_eq!(counts, vec![3, 2, 1, 0]);
    }

    #[test]
    fn explicit_h_suppresses_implicit_h() {
        let layout_output = layout_native("[NH4+]").expect("ammonium layout failed");
        assert_eq!(layout_output.atoms[0].hcount, 4);
        assert_eq!(layout_output.atoms[0].implicit_h, 0);
    }

    #[test]
    fn bracket_hydrogen_atoms_fold_into_neighbor() {
        // Methane written with four explicit [H] atoms collapses into one
        // carbon carrying the folded hydrogen count; the only remaining H is the
        // virtual label placeholder, never a drawn atom.
        let layout_output = layout_native("C([H])([H])([H])[H]").expect("methane layout failed");
        assert_eq!(layout_output.atoms[0].symbol, "C");
        assert_eq!(layout_output.atoms[0].hcount, 4);
        assert_eq!(layout_output.atoms[0].implicit_h, 0);
        let drawn_hydrogen_count = layout_output
            .atoms
            .iter()
            .filter(|atom| atom.symbol == "H" && !atom.virtual_h)
            .count();
        assert_eq!(drawn_hydrogen_count, 0);
    }

    #[test]
    fn folded_hydrogens_drop_from_neighbor_ordering() {
        // The dichloromethyl carbon keeps only its three heavy neighbors after
        // the explicit hydrogen is folded away.
        let layout_output = layout_native("ClC([H])Cl").expect("dichloromethane layout failed");
        let heavy_atom_count = layout_output
            .atoms
            .iter()
            .filter(|atom| atom.symbol != "H")
            .count();
        assert_eq!(heavy_atom_count, 3);
        assert_eq!(
            layout_output
                .bonds
                .iter()
                .filter(|bond| !bond.virtual_bond)
                .count(),
            2
        );
    }

    #[test]
    fn isotopic_and_charged_hydrogens_are_kept() {
        // Deuterium is a real, drawn atom; only the plain [H] atoms fold away.
        let layout_output =
            layout_native("[2H]C([H])([H])[H]").expect("deuteromethane layout failed");
        let drawn_hydrogen_count = layout_output
            .atoms
            .iter()
            .filter(|atom| atom.symbol == "H" && !atom.virtual_h)
            .count();
        assert_eq!(drawn_hydrogen_count, 1);
    }

    #[test]
    fn double_bond() {
        let layout_output = layout_native("C=C").expect("ethylene layout failed");
        assert_eq!(layout_output.bonds[0].order, 2);
    }

    #[test]
    fn cumulated_double_bonds_are_linear() {
        let layout_output = layout_native("O=C=O").expect("carbon dioxide layout failed");
        assert!(atoms_are_collinear(&layout_output, 0, 1, 2));
    }

    #[test]
    fn single_triple_chain_is_linear_at_middle_atom() {
        let layout_output = layout_native("CC#N").expect("nitrile layout failed");
        assert!(atoms_are_collinear(&layout_output, 0, 1, 2));
    }

    #[test]
    fn saturated_two_neighbor_atom_keeps_zigzag() {
        let layout_output = layout_native("CCC").expect("propane layout failed");
        assert!(!atoms_are_collinear(&layout_output, 0, 1, 2));
    }

    #[test]
    fn cyclohexane() {
        let layout_output = layout_native("C1CCCCC1").expect("cyclohexane layout failed");
        assert_eq!(layout_output.atoms.len(), 6);
        assert_eq!(layout_output.bonds.len(), 6);
    }

    #[test]
    fn separated_ring_systems_do_not_overlap() {
        let layout_output = layout_native("CC(N)C(=O)OCCC1=CC=CC=C1NCC1=CC=CC=C1")
            .expect("two-ring molecule layout failed");
        assert!(max_bond_length(&layout_output) < 1.2);
    }

    #[test]
    fn symmetric_hub_is_mirror_symmetric() {
        // Tetraethylmethane's four equal arms must be drawn as a mirror-symmetric
        // figure (the "tweezers" depiction), not a rotational pinwheel. The layout
        // is built about the vertical axis, so reflecting every atom across that
        // axis must reproduce the same set of positions.
        let layout_output =
            layout_native("CCC(CC)(CC)CC").expect("tetraethylmethane layout failed");
        let cx = layout_output
            .atoms
            .iter()
            .map(|atom| atom.pos.x)
            .sum::<f64>()
            / layout_output.atoms.len() as f64;
        for atom in &layout_output.atoms {
            let mirrored_x = 2.0 * cx - atom.pos.x;
            let found = layout_output.atoms.iter().any(|b| {
                (b.pos.x - mirrored_x).abs() < 1e-6 && (b.pos.y - atom.pos.y).abs() < 1e-6
            });
            assert!(
                found,
                "no mirror partner for atom at ({:.3}, {:.3})",
                atom.pos.x, atom.pos.y
            );
        }
    }

    #[test]
    fn branch_point_routes_largest_subtree_straight_ahead() {
        // The central acetal carbon has four substituents; the long ester chain
        // must continue away from the rest of the molecule instead of folding
        // back over the other ester group.
        let layout_output =
            layout_native("BrCC(=O)OC(C)(O)OC(=O)C").expect("acetal diester layout failed");
        assert!(
            min_nonbonded_distance(&layout_output) > 0.5,
            "branches overlap: min non-bonded distance {}",
            min_nonbonded_distance(&layout_output)
        );
    }

    #[test]
    fn steroid_ring_system_has_regular_bond_lengths() {
        let layout_output =
            layout_native("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
                .expect("steroid-like molecule layout failed");
        assert!(max_bond_length(&layout_output) < 1.25);
    }

    #[test]
    fn isobutane() {
        let layout_output = layout_native("CC(C)C").expect("isobutane layout failed");
        assert_eq!(layout_output.atoms.len(), 4);
        assert_eq!(layout_output.bonds.len(), 3);
    }

    #[test]
    fn naphthalene_kekule() {
        let layout_output =
            layout_native("C1=CC2=CC=CC=C2C=C1").expect("naphthalene layout failed");
        assert_eq!(layout_output.atoms.len(), 10);
        assert_eq!(layout_output.bonds.len(), 11);
    }

    // ── Dot-disconnected structures ──────────────────────────────────────────

    #[test]
    fn dot_creates_no_bond() {
        let layout_output = layout_native("CCO.CCO").expect("two ethanols failed");
        assert_eq!(
            layout_output.atoms.iter().filter(|a| !a.virtual_h).count(),
            6
        );
        // Two C-C-O fragments: 4 real bonds, none crossing the dot.
        let real: Vec<_> = layout_output
            .bonds
            .iter()
            .filter(|b| !b.virtual_bond)
            .collect();
        assert_eq!(real.len(), 4);
        assert!(real.iter().all(|b| (b.from < 3) == (b.to < 3)));
    }

    #[test]
    fn salt_fragments_are_disconnected_and_ordered() {
        // Sodium acetate: no bond may touch the sodium ion, and fragments keep
        // SMILES writing order left to right with a visible gap.
        let layout_output = layout_native("CC(=O)[O-].[Na+]").expect("sodium acetate failed");
        let sodium_index = 4;
        assert_eq!(layout_output.atoms[sodium_index].symbol, "Na");
        assert_eq!(layout_output.atoms[sodium_index].charge, 1);
        assert!(layout_output
            .bonds
            .iter()
            .all(|bond| bond.from != sodium_index && bond.to != sodium_index));

        let max_acetate_x = (0..4)
            .map(|atom_index| layout_output.atoms[atom_index].pos.x)
            .fold(f64::MIN, f64::max);
        assert!(
            layout_output.atoms[sodium_index].pos.x >= max_acetate_x + 1.0,
            "Na+ should sit clearly right of the acetate fragment"
        );
    }

    #[test]
    fn bare_ion_pair() {
        let layout_output = layout_native("[Na+].[Cl-]").expect("sodium chloride failed");
        assert_eq!(layout_output.atoms.len(), 2);
        assert!(layout_output.bonds.is_empty());
        assert!(layout_output.atoms[1].pos.x > layout_output.atoms[0].pos.x);
    }

    #[test]
    fn three_fragments() {
        let layout_output = layout_native("O.O.O").expect("three waters failed");
        let heavy: Vec<_> = layout_output
            .atoms
            .iter()
            .filter(|a| !a.virtual_h)
            .collect();
        assert_eq!(heavy.len(), 3);
        assert!(layout_output.bonds.iter().all(|bond| bond.virtual_bond));
        assert!(heavy[0].pos.x < heavy[1].pos.x && heavy[1].pos.x < heavy[2].pos.x);
    }

    #[test]
    fn ring_closure_across_dot_joins_fragments() {
        // OpenSMILES: "C1.C1" is ethane — the ring-closure bond still forms
        // even though a dot separates the digits.
        let layout_output = layout_native("C1.C1").expect("dot ring closure failed");
        assert_eq!(layout_output.atoms.len(), 2);
        assert_eq!(layout_output.bonds.len(), 1);
        assert_eq!(layout_output.atoms[0].implicit_h, 3);
    }

    #[test]
    fn aromatic_fragments_kekulize_independently() {
        let layout_output = layout_native("c1ccccc1.c1ccccc1").expect("two benzenes failed");
        assert_eq!(layout_output.atoms.len(), 12);
        assert_eq!(order_counts(&layout_output), (6, 6));
        assert!(layout_output
            .bonds
            .iter()
            .all(|b| (b.from < 6) == (b.to < 6)));
    }

    // ── Aromatic (lowercase) SMILES and kekulization ─────────────────────────

    fn order_counts(layout_output: &LayoutOutput) -> (usize, usize) {
        let real = layout_output.bonds.iter().filter(|b| !b.virtual_bond);
        (
            real.clone().filter(|bond| bond.order == 1).count(),
            real.filter(|bond| bond.order == 2).count(),
        )
    }

    #[test]
    fn benzene_aromatic() {
        let layout_output = layout_native("c1ccccc1").expect("aromatic benzene failed");
        assert_eq!(layout_output.atoms.len(), 6);
        assert_eq!(order_counts(&layout_output), (3, 3));
        assert!(layout_output.atoms.iter().all(|atom| atom.implicit_h == 1));
    }

    #[test]
    fn aromatic_atom_indices_match_writing_order() {
        // Kekulization must not reorder atoms: index N is the Nth atom token,
        // so show-indices / highlight / arrow references keep working.
        let layout_output = layout_native("Cc1ccncc1").expect("4-methylpyridine failed");
        let symbols: Vec<&str> = layout_output
            .atoms
            .iter()
            .filter(|a| !a.virtual_h)
            .map(|atom| atom.symbol.as_str())
            .collect();
        assert_eq!(symbols, vec!["C", "C", "C", "C", "N", "C", "C"]);
    }

    #[test]
    fn pyridine_aromatic() {
        let layout_output = layout_native("c1ccncc1").expect("pyridine failed");
        assert_eq!(layout_output.atoms[3].symbol, "N");
        assert_eq!(layout_output.atoms[3].implicit_h, 0);
        assert_eq!(order_counts(&layout_output), (3, 3));
    }

    #[test]
    fn pyrrole_aromatic() {
        let layout_output = layout_native("c1cc[nH]c1").expect("pyrrole failed");
        let nitrogen = layout_output
            .atoms
            .iter()
            .find(|atom| atom.symbol == "N")
            .unwrap();
        assert_eq!(nitrogen.hcount, 1);
        assert_eq!(order_counts(&layout_output), (3, 2));
    }

    #[test]
    fn furan_and_thiophene_aromatic() {
        for smiles in ["c1occc1", "c1sccc1"] {
            let layout_output = layout_native(smiles).expect("5-ring heteroaromatic failed");
            assert_eq!(
                order_counts(&layout_output),
                (3, 2),
                "wrong kekulization for {smiles}"
            );
            let hetero = &layout_output.atoms[1];
            assert_eq!(hetero.implicit_h, 0);
        }
    }

    #[test]
    fn imidazole_aromatic() {
        let layout_output = layout_native("c1cnc[nH]1").expect("imidazole failed");
        assert_eq!(order_counts(&layout_output), (3, 2));
    }

    #[test]
    fn n_methylpyrrole_aromatic() {
        // A three-connected aromatic n carries no H and no double bond.
        let layout_output = layout_native("Cn1cccc1").expect("N-methylpyrrole failed");
        assert_eq!(layout_output.atoms[1].symbol, "N");
        assert_eq!(layout_output.atoms[1].implicit_h, 0);
        assert_eq!(order_counts(&layout_output), (4, 2));
    }

    #[test]
    fn naphthalene_aromatic() {
        let layout_output = layout_native("c1ccc2ccccc2c1").expect("aromatic naphthalene failed");
        assert_eq!(layout_output.atoms.len(), 10);
        assert_eq!(order_counts(&layout_output), (6, 5));
    }

    #[test]
    fn indane_mixed_aromatic_aliphatic() {
        // Spec example: aromatic ring fused to an aliphatic ring.
        let layout_output = layout_native("c1ccc2CCCc2c1").expect("indane failed");
        assert_eq!(layout_output.atoms.len(), 9);
        assert_eq!(order_counts(&layout_output), (7, 3));
    }

    #[test]
    fn biphenyl_explicit_and_implicit_single_link() {
        for smiles in ["c1ccccc1-c1ccccc1", "c1ccccc1c1ccccc1"] {
            let layout_output = layout_native(smiles).expect("biphenyl failed");
            assert_eq!(layout_output.atoms.len(), 12);
            assert_eq!(
                order_counts(&layout_output),
                (7, 6),
                "wrong kekulization for {smiles}"
            );
            // The inter-ring bond stays single.
            let link = layout_output
                .bonds
                .iter()
                .find(|b| (b.from < 6) != (b.to < 6))
                .unwrap();
            assert_eq!(link.order, 1);
        }
    }

    #[test]
    fn pyridinone_exocyclic_double_bond() {
        // The carbonyl carbon already has its double bond outside the ring and
        // must not receive another one during kekulization.
        let layout_output = layout_native("O=c1cccc[nH]1").expect("2-pyridinone failed");
        assert_eq!(order_counts(&layout_output), (4, 3));
    }

    #[test]
    fn explicit_aromatic_bond_symbol() {
        let layout_output = layout_native("c1:c:c:c:c:c1").expect("explicit ':' benzene failed");
        assert_eq!(order_counts(&layout_output), (3, 3));
    }

    #[test]
    fn charged_aromatic_rings() {
        // Pyridinium: the protonated nitrogen still takes part in a double bond.
        let pyridinium = layout_native("c1cc[nH+]cc1").expect("pyridinium failed");
        assert_eq!(order_counts(&pyridinium), (3, 3));

        // Pyrylium: positively charged oxygen participates in a double bond.
        let pyrylium = layout_native("c1cc[o+]cc1").expect("pyrylium failed");
        assert_eq!(order_counts(&pyrylium), (3, 3));
    }

    #[test]
    fn wildcard_in_aromatic_ring() {
        let layout_output = layout_native("c1cc*cc1").expect("aromatic ring with wildcard failed");
        let (_, doubles) = order_counts(&layout_output);
        assert!(
            doubles >= 2,
            "expected an alternating pattern, got {doubles} double bonds"
        );
    }

    #[test]
    fn azulene_aromatic() {
        // Non-alternant 5-7 fused system: every carbon needs a double bond.
        let layout_output = layout_native("c1ccc2cccc2cc1").expect("azulene failed");
        assert_eq!(layout_output.atoms.len(), 10);
        assert_eq!(order_counts(&layout_output), (6, 5));
    }

    #[test]
    fn caffeine_aromatic() {
        let layout_output = layout_native("Cn1cnc2c1c(=O)n(C)c(=O)n2C").expect("caffeine failed");
        // Purine core: the imidazole C=N plus the C4=C5 bridge double bond,
        // and the two exocyclic carbonyls.
        assert_eq!(
            layout_output.atoms.iter().filter(|a| !a.virtual_h).count(),
            14
        );
        let (_, doubles) = order_counts(&layout_output);
        assert_eq!(doubles, 4);
    }

    #[test]
    fn aromatic_bond_symbol_requires_aromatic_atoms() {
        let err = layout_native("C:C").expect_err("':' between aliphatic atoms should fail");
        assert!(err.contains("aromatic"));
    }

    #[test]
    fn aromatic_atom_outside_ring_errors() {
        let err = layout_native("Cc").expect_err("acyclic aromatic atom should fail");
        assert!(err.contains("ring"));
    }

    #[test]
    fn unkekulizable_ring_errors() {
        // Pyrrole written without its hydrogen has five atoms all demanding a
        // double bond; no perfect matching exists.
        let err = layout_native("c1ccnc1").expect_err("H-less pyrrole should fail");
        assert!(err.contains("cannot assign alternating double bonds"));
        assert!(err.contains("starting with `c` at character 1"));
        assert!(err.contains("c1cc[nH]c1"));

        let err = layout_native("c1ccccc1.c1ccnc1").expect_err("second system should fail");
        assert!(err.contains("starting with `c` at character 10"));
    }

    // ── Standards-first stereochemistry and drawing extensions ───────────────

    #[test]
    fn forced_wedge_up_bond() {
        let layout_output = layout_native("C!wN").expect("forced wedge up failed");
        assert_eq!(layout_output.bonds[0].stereo, "wedge_up");
        assert_eq!(layout_output.bonds[0].direction, "none");
        assert!(layout_output.bonds[0].forced_stereo);
    }

    #[test]
    fn forced_wedge_down_bond() {
        let layout_output = layout_native("C!hN").expect("forced wedge down failed");
        assert_eq!(layout_output.bonds[0].stereo, "wedge_down");
        assert_eq!(layout_output.bonds[0].direction, "none");
        assert!(layout_output.bonds[0].forced_stereo);
    }

    #[test]
    fn forced_wavy_bond() {
        let layout_output = layout_native("C!sN").expect("forced wavy bond failed");
        assert_eq!(layout_output.bonds[0].stereo, "wavy");
        assert_eq!(layout_output.bonds[0].direction, "none");
        assert!(!layout_output.bonds[0].forced_stereo);
    }

    #[test]
    fn forced_dashed_bond() {
        let layout_output = layout_native("C!dN").expect("forced dashed bond failed");
        assert_eq!(layout_output.bonds[0].stereo, "dashed");
        assert_eq!(layout_output.bonds[0].direction, "none");
        assert!(!layout_output.bonds[0].forced_stereo);
    }

    #[test]
    fn forced_wavy_and_dashed_in_chain() {
        let layout_output = layout_native("CC!sO!dN").expect("wavy/dashed chain failed");
        assert_eq!(layout_output.bonds[0].stereo, "none");
        assert_eq!(layout_output.bonds[1].stereo, "wavy");
        assert_eq!(layout_output.bonds[2].stereo, "dashed");
    }

    #[test]
    fn forced_wavy_does_not_disturb_real_directional_bonds() {
        // A genuine trans alkene next to a forced wavy bond: the wavy marker
        // must not consume or shift the cis/trans direction tokens.
        let layout_output = layout_native("F/C=C/C!sN").expect("mixed directional/wavy failed");
        let wavy = layout_output
            .bonds
            .iter()
            .filter(|bond| bond.stereo == "wavy")
            .count();
        assert_eq!(wavy, 1);
        let directional = layout_output
            .bonds
            .iter()
            .filter(|bond| bond.direction != "none")
            .count();
        assert_eq!(directional, 2);
    }

    #[test]
    fn slash_without_double_bond_is_invalid() {
        let err = layout_native("C/N").expect_err("isolated slash should fail");
        assert!(err.contains("Directional"));
    }

    #[test]
    fn forced_wedge_in_chain() {
        let layout_output = layout_native("CC!wN").expect("forced wedge in chain failed");
        assert_eq!(layout_output.bonds[0].stereo, "none");
        assert_eq!(layout_output.bonds[1].stereo, "wedge_up");
        assert!(layout_output.bonds[1].forced_stereo);
    }

    #[test]
    fn alkene_directional_bonds_do_not_render_as_wedges() {
        let layout_output = layout_native("F/C=C/F").expect("E-alkene failed");
        let double = layout_output
            .bonds
            .iter()
            .find(|bond| bond.order == 2)
            .unwrap();
        assert_eq!(double.order, 2);
        assert!(layout_output.bonds.iter().all(|bond| bond.stereo == "none"));
        assert_eq!(
            layout_output
                .bonds
                .iter()
                .filter(|bond| bond.direction != "none")
                .count(),
            2
        );
    }

    #[test]
    fn trans_alkene_substituents_are_opposite() {
        let layout_output = layout_native("F/C=C/F").expect("trans alkene failed");
        assert_eq!(alkene_substituent_side_product(&layout_output), -1);
    }

    #[test]
    fn cis_alkene_substituents_are_same_side() {
        let layout_output = layout_native("F/C=C\\F").expect("cis alkene failed");
        assert_eq!(alkene_substituent_side_product(&layout_output), 1);
    }

    #[test]
    fn branch_direction_matches_opensmiles_examples() {
        let trans = layout_native("C(\\F)=C/F").expect("branch trans failed");
        let cis = layout_native("C(/F)=C/F").expect("branch cis failed");
        assert_eq!(alkene_substituent_side_product(&trans), -1);
        assert_eq!(alkene_substituent_side_product(&cis), 1);
    }

    #[test]
    fn pyrethroid_like_smiles_with_multiple_alkene_markers_parses() {
        let smiles = "CC1=C(C(=O)C[C@@H]1OC(=O)[C@@H]2[C@H](C2(C)C)/C=C(\\C)/C(=O)OC)C/C=C\\C=C";
        let layout_output = layout_native(smiles).expect("pyrethroid-like molecule failed");
        assert_eq!(
            layout_output
                .atoms
                .iter()
                .filter(|atom| !atom.virtual_h)
                .count(),
            27
        );
        assert_eq!(
            layout_output
                .bonds
                .iter()
                .filter(|bond| bond.direction != "none")
                .count(),
            5
        );
    }

    #[test]
    fn conjugated_diene_shares_directional_bond_between_double_bonds() {
        let trans_trans = layout_native("C/C=C/C=C/C").expect("trans,trans-diene failed");
        assert_eq!(
            double_bond_substituent_side_product(&trans_trans, 0, (1, 2), 3),
            -1
        );
        assert_eq!(
            double_bond_substituent_side_product(&trans_trans, 2, (3, 4), 5),
            -1
        );

        let trans_cis = layout_native("C/C=C/C=C\\C").expect("trans,cis-diene failed");
        assert_eq!(
            double_bond_substituent_side_product(&trans_cis, 0, (1, 2), 3),
            -1
        );
        assert_eq!(
            double_bond_substituent_side_product(&trans_cis, 2, (3, 4), 5),
            1
        );
    }

    #[test]
    fn directional_bond_next_to_carbonyl_marks_only_the_alkene() {
        // Atoms: O0 C1 O2 C3 C4 C5 O6 O7; the alkene is C3=C4.
        let fumaric = layout_native("OC(=O)/C=C/C(=O)O").expect("fumaric acid failed");
        assert_eq!(
            double_bond_substituent_side_product(&fumaric, 1, (3, 4), 5),
            -1
        );

        let maleic = layout_native("OC(=O)/C=C\\C(=O)O").expect("maleic acid failed");
        assert_eq!(
            double_bond_substituent_side_product(&maleic, 1, (3, 4), 5),
            1
        );

        layout_native("CC(=O)/C=C/C").expect("trans enone failed");
        layout_native("C=C/C=C/C").expect("terminal diene failed");
    }

    #[test]
    fn one_sided_directional_bond_errors() {
        let err = layout_native("F/C=CF").expect_err("one-sided marker should fail");
        assert!(err.contains("must mark both ends"));
    }

    #[test]
    fn directional_bond_away_from_double_bonds_errors() {
        let err = layout_native("F/CC").expect_err("stray marker should fail");
        assert!(err.contains("only supported around double bonds"));
    }

    #[test]
    fn conflicting_cis_trans_markers_error() {
        let err = layout_native("C/C(\\F)=C/F").expect_err("conflicting markers should fail");
        assert!(err.contains("Conflicting"));
    }

    #[test]
    fn tetrahedral_chirality_adds_rendered_stereo() {
        let layout_output = layout_native("N[C@@H](C)C(=O)O").expect("chiral alanine failed");
        assert!(layout_output
            .atoms
            .iter()
            .any(|atom| atom.chirality == "tetra_clockwise"));
        // Exactly one wedge for the single stereocenter; direction checked below.
        let wedge_bonds = layout_output
            .bonds
            .iter()
            .filter(|bond| bond.stereo != "none")
            .count();
        let wedge_h = layout_output
            .atoms
            .iter()
            .filter(|atom| atom.stereo_h != "none")
            .count();
        assert_eq!(wedge_bonds + wedge_h, 1);
        assert!(chirality_matches_smiles("N[C@@H](C)C(=O)O"));
    }

    #[test]
    fn inverting_chirality_flips_the_wedge() {
        // Same skeleton/layout, opposite chirality token ⇒ the wedge on the chosen
        // bond must flip. Guards against regressing to a fixed @→up / @@→down map.
        let clockwise_layout = layout_native("N[C@@H](C)C(=O)O").expect("R failed");
        let anticlockwise_layout = layout_native("N[C@H](C)C(=O)O").expect("S failed");
        let stereo = |layout_output: &LayoutOutput| -> String {
            layout_output
                .bonds
                .iter()
                .find(|bond| bond.stereo != "none")
                .map(|bond| bond.stereo.clone())
                .or_else(|| {
                    layout_output
                        .atoms
                        .iter()
                        .find(|atom| atom.stereo_h != "none")
                        .map(|atom| atom.stereo_h.clone())
                })
                .unwrap_or_else(|| "none".to_string())
        };
        assert_ne!(stereo(&clockwise_layout), "none");
        assert_ne!(stereo(&anticlockwise_layout), "none");
        assert_ne!(stereo(&clockwise_layout), stereo(&anticlockwise_layout));
        assert!(chirality_matches_smiles("N[C@@H](C)C(=O)O"));
        assert!(chirality_matches_smiles("N[C@H](C)C(=O)O"));
    }

    #[test]
    fn tetrahedral_stereo_prefers_oh_substituent() {
        let layout_output = layout_native("CC[C@@H](O)CC/C=C/CO").expect("chiral alcohol failed");
        assert!(chiral_oxygen_bond_has_stereo(&layout_output));
    }

    #[test]
    fn chiral_alcohol_keeps_long_chain_as_continuation() {
        let layout_output = layout_native("CC[C@@H](O)CC/C=C/CO").expect("chiral alcohol failed");
        // The long chain leaves the stereocenter (atom 2) as a straight
        // continuation; the short ethyl tail (atom 0) sits on the opposite side.
        // Direction-agnostic so it survives the conventional horizontal mirror.
        let continuation_direction = layout_output.atoms[4].pos.x - layout_output.atoms[2].pos.x;
        assert!(continuation_direction.abs() > 1e-6);
        assert!(
            (layout_output.atoms[9].pos.x - layout_output.atoms[4].pos.x) * continuation_direction
                > 0.0,
            "terminal alcohol should remain on the zig-zag continuation"
        );
        assert!(
            (layout_output.atoms[0].pos.x - layout_output.atoms[2].pos.x) * continuation_direction
                < 0.0,
            "short ethyl branch should sit opposite the long continuation"
        );
        assert!(layout_output
            .atoms
            .iter()
            .all(|atom| atom.stereo_h == "none"));
    }

    #[test]
    fn steroid_stereo_prefers_exocyclic_oh_over_ring_bond() {
        let layout_output =
            layout_native("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
                .expect("steroid-like molecule layout failed");
        assert!(chiral_oxygen_bond_has_stereo(&layout_output));
    }

    #[test]
    fn steroid_ring_chiral_hydrogens_are_rendered() {
        let smiles = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O";
        let layout_output = layout_native(smiles).expect("steroid-like molecule layout failed");

        // The three ring-fusion stereocenters (no exocyclic substituent) render
        // an explicit wedge/hash hydrogen. The 17-OH carbon and the quaternary
        // C13 instead wedge their exocyclic substituent, so they get no H label.
        let stereo_h: Vec<&str> = layout_output
            .atoms
            .iter()
            .filter_map(|atom| (atom.stereo_h != "none").then_some(atom.stereo_h.as_str()))
            .collect();
        assert_eq!(stereo_h.len(), 3);

        // Adjacent ring-fusion stereocenters must point to opposite faces.
        assert_ne!(layout_output.atoms[5].stereo_h, "none");
        assert_ne!(layout_output.atoms[6].stereo_h, "none");
        assert_ne!(
            layout_output.atoms[5].stereo_h,
            layout_output.atoms[6].stereo_h
        );

        // The 17-OH bond is wedged (not the hydrogen).
        assert!(chiral_oxygen_bond_has_stereo(&layout_output));

        // Every stereocenter is depicted with the geometrically correct handedness.
        assert!(chirality_matches_smiles(smiles));
    }

    // ── Extended chirality classes and quadruple bonds ───────────────────────

    /// Cosine of the angle neighbor-a → center → neighbor-b.
    fn bond_angle_cos(
        layout_output: &LayoutOutput,
        center: usize,
        first_neighbor: usize,
        second_neighbor: usize,
    ) -> f64 {
        let center_position = layout_output.atoms[center].pos;
        let (first_offset_x, first_offset_y) = (
            layout_output.atoms[first_neighbor].pos.x - center_position.x,
            layout_output.atoms[first_neighbor].pos.y - center_position.y,
        );
        let (second_offset_x, second_offset_y) = (
            layout_output.atoms[second_neighbor].pos.x - center_position.x,
            layout_output.atoms[second_neighbor].pos.y - center_position.y,
        );
        (first_offset_x * second_offset_x + first_offset_y * second_offset_y)
            / ((first_offset_x * first_offset_x + first_offset_y * first_offset_y).sqrt()
                * (second_offset_x * second_offset_x + second_offset_y * second_offset_y).sqrt())
    }

    #[test]
    fn square_planar_cis_and_trans() {
        // Neighbors in writing order: N(0), N(2), Cl(3), Cl(4) around Pt(1).
        // @SP1 (U shape) puts the two Cl on adjacent corners: cis, 90° apart.
        let cis = layout_native("N[Pt@SP1](N)(Cl)Cl").expect("cisplatin failed");
        assert!(bond_angle_cos(&cis, 1, 3, 4).abs() < 1e-6);
        // @SP2 (4 shape) pairs Cl trans to Cl: 180° apart.
        let trans = layout_native("N[Pt@SP2](N)(Cl)Cl").expect("transplatin failed");
        assert!((bond_angle_cos(&trans, 1, 3, 4) + 1.0).abs() < 1e-6);
        // @SP3 (Z shape) pairs N1 trans to Cl4, so the Cl pair is cis again.
        let third_configuration = layout_native("N[Pt@SP3](N)(Cl)Cl").expect("SP3 failed");
        assert!(bond_angle_cos(&third_configuration, 1, 3, 4).abs() < 1e-6);
        assert!((bond_angle_cos(&third_configuration, 1, 0, 4) + 1.0).abs() < 1e-6);
    }

    #[test]
    fn square_planar_all_neighbors_at_right_angles() {
        let layout_output = layout_native("C[Fe@SP1](Cl)(Br)I").expect("SP iron failed");
        assert_eq!(layout_output.atoms[1].chirality, "square_planar");
        for (a, b) in [(0, 2), (2, 3), (3, 4)] {
            assert!(bond_angle_cos(&layout_output, 1, a, b).abs() < 1e-6);
        }
    }

    #[test]
    fn extended_stereo_classes_are_reported_as_undepicted() {
        for (smiles, chirality) in [
            ("S[As@TB1](F)(Cl)(Br)N", "trigonal_bipyramidal"),
            ("C[Co@OH1](F)(Cl)(Br)(I)N", "octahedral"),
            ("C[Co@OH28](F)(Cl)(Br)(I)N", "octahedral"),
            ("NC(Br)=[C@AL1]=C(O)C", "allenal"),
        ] {
            let layout_output = layout_native(smiles).expect("extended chirality should parse");
            let center = layout_output
                .atoms
                .iter()
                .position(|atom| atom.chirality == chirality)
                .unwrap_or_else(|| panic!("missing {chirality} center for {smiles}"));
            assert_eq!(layout_output.undepicted_stereo.len(), 1, "{smiles}");
            assert_eq!(layout_output.undepicted_stereo[0].atom, center);
            assert!(layout_output.undepicted_stereo[0]
                .reason
                .contains("not drawn"));
            assert!(layout_output.bonds.iter().all(|bond| bond.stereo == "none"));
            assert!(layout_output
                .atoms
                .iter()
                .all(|atom| atom.stereo_h == "none"));
        }
    }

    #[test]
    fn quadruple_bond_order() {
        // The classic metal-metal quadruple bond, e.g. in [Re2Cl8]2-.
        let layout_output = layout_native("[Re]$[Re]").expect("quadruple bond failed");
        assert_eq!(layout_output.bonds[0].order, 4);
    }

    #[test]
    fn quadruple_bond_counts_toward_valence() {
        let layout_output = layout_native("C$C").expect("C$C failed");
        assert_eq!(layout_output.bonds[0].order, 4);
        assert!(layout_output.atoms.iter().all(|atom| atom.implicit_h == 0));
    }

    // ── Abbreviation tests ────────────────────────────────────────────────────

    #[test]
    fn abbreviation_label_becomes_wildcard_atom() {
        let molecule = parse_molecule("C({PPh3})=O").expect("abbreviation failed");
        assert_eq!(molecule.atoms.len(), 3);
        assert_eq!(molecule.atoms[1].symbol, "*");
        assert!(molecule.atoms[1].has_explicit_h);
        assert_eq!(molecule.atoms[1].abbrev, "PPh3");
        assert_eq!(molecule.atoms[1].abbrev_style, "");
        assert_eq!(molecule.atoms[1].abbrev_anchor_len, 0);
    }

    #[test]
    fn abbreviation_labels_keep_writing_order() {
        let molecule = parse_molecule("{OEt}C(=O){NHR}").expect("abbreviations failed");
        assert_eq!(molecule.atoms[0].abbrev, "OEt");
        assert_eq!(molecule.atoms[3].abbrev, "NHR");
    }

    #[test]
    fn abbreviation_styles_reach_their_atoms() {
        let molecule = parse_molecule("{PPh3|P}C({LG|red})=O").expect("styled labels failed");
        assert_eq!(molecule.atoms[0].abbrev, "PPh3");
        assert_eq!(molecule.atoms[0].abbrev_style, "P");
        assert_eq!(molecule.atoms[2].abbrev, "LG");
        assert_eq!(molecule.atoms[2].abbrev_style, "red");
    }

    #[test]
    fn abbreviation_anchor_positions() {
        let molecule = parse_molecule("{>CAT}C({C>AT}){CA>T}").expect("anchors failed");
        for (atom_index, anchor) in [(0, 0), (2, 1), (3, 2)] {
            assert_eq!(molecule.atoms[atom_index].abbrev, "CAT");
            assert_eq!(molecule.atoms[atom_index].abbrev_anchor, anchor);
            assert_eq!(molecule.atoms[atom_index].abbrev_anchor_len, 1);
        }
    }

    #[test]
    fn arrow_marker_outside_abbreviation_is_rejected() {
        assert!(parse_molecule("C>C").is_err());
        assert!(parse_molecule("{CAT>}C").is_err());
    }

    #[test]
    fn malformed_delimiters_report_actionable_errors() {
        let unclosed_label = layout_native("C{OH").expect_err("unclosed label should fail");
        assert!(unclosed_label.contains("unclosed custom label"));
        assert!(unclosed_label.contains("add `}`"));

        let unclosed_bracket = layout_native("C[NH2").expect_err("unclosed bracket should fail");
        assert!(unclosed_bracket.contains("unclosed bracket atom"));
        assert!(unclosed_bracket.contains("add `]`"));

        let unmatched_label_end =
            layout_native("CO}").expect_err("unmatched label end should fail");
        assert!(unmatched_label_end.contains("unmatched `}`"));
    }

    #[test]
    fn empty_smiles_reports_actionable_error() {
        let error = layout_native("   ").expect_err("empty SMILES should fail");
        assert!(error.contains("expression is empty"));
        assert!(error.contains("at least one atom"));
    }

    #[test]
    fn unclosed_ring_reports_ring_number_and_correction() {
        let error = layout_native("C1CC").expect_err("unclosed ring should fail");
        assert!(error.contains("ring closure 1 opened at character 2"));
        assert!(error.contains("never closed"));
        assert!(error.contains("repeat the ring number"));
    }

    #[test]
    fn unknown_abbreviation_style_reports_supported_forms() {
        let error = layout_native("C{LG|chartreuse}").expect_err("unknown style should fail");
        assert!(error.contains("unknown abbreviation style `chartreuse`"));
        assert!(error.contains("element symbol"));
        assert!(error.contains("#RRGGBB"));
    }

    #[test]
    fn forced_wedge_extensions_mark_their_bonds() {
        let molecule = parse_molecule("C!wN!hO").expect("forced wedges failed");
        assert_eq!(molecule.bonds[0].stereo, graph::BondStereo::WedgeUp);
        assert_eq!(molecule.bonds[1].stereo, graph::BondStereo::WedgeDown);
        assert!(molecule.bonds.iter().all(|bond| bond.forced_stereo));
        assert!(molecule
            .bonds
            .iter()
            .all(|bond| bond.direction == graph::BondDirection::None));
    }

    #[test]
    fn curl_extensions_combine_with_styles_and_orders() {
        let molecule = parse_molecule("CCC!cC").expect("plain curl failed");
        assert!(molecule.bonds[2].curl);
        assert_eq!(molecule.bonds[2].stereo, graph::BondStereo::None);
        assert_eq!(molecule.bonds[2].order, graph::BondOrder::Single);

        let molecule = parse_molecule("CCC!c!wC").expect("curl wedge failed");
        assert!(molecule.bonds[2].curl);
        assert_eq!(molecule.bonds[2].stereo, graph::BondStereo::WedgeUp);

        let molecule = parse_molecule("CCC!c=C").expect("curl double failed");
        assert!(molecule.bonds[2].curl);
        assert_eq!(molecule.bonds[2].order, graph::BondOrder::Double);
    }

    #[test]
    fn ring_closures_after_branches_keep_writing_order() {
        let molecule = parse_molecule("C1=CCCC(=O)1").expect("cyclopentenone failed");
        let carbonyl_carbon = 4;
        let written_neighbors: Vec<usize> = molecule.neighbor_bonds[carbonyl_carbon]
            .iter()
            .map(|&bond_index| {
                let bond = &molecule.bonds[bond_index];
                bond.from + bond.to - carbonyl_carbon
            })
            .collect();
        assert_eq!(written_neighbors, vec![3, 5, 0]);

        let molecule = parse_molecule("C(=O)(O)1N.C1").expect("trailing ring bond failed");
        assert_eq!(molecule.bonds.len(), 4);
        assert!(molecule.adj[0].iter().any(|&(neighbor, _)| neighbor == 4));
    }

    #[test]
    fn organic_atoms_record_aromaticity() {
        let molecule = parse_molecule("CCO").expect("ethanol failed");
        assert!(molecule.atoms.iter().all(|atom| !atom.aromatic));

        let molecule = parse_molecule("Clc1ccccc1").expect("chlorobenzene failed");
        let aromatic: Vec<bool> = molecule.atoms.iter().map(|atom| atom.aromatic).collect();
        assert_eq!(aromatic, vec![false, true, true, true, true, true, true]);
        assert_eq!(molecule.atoms[0].symbol, "Cl");
    }

    #[test]
    fn bracket_symbols_read_whole_elements() {
        // [nH] is an aromatic nitrogen, and the `c` in [Sc] belongs to scandium.
        let molecule = parse_molecule("c1cc[nH]c1[Sc]").expect("bracket atoms failed");
        assert!(molecule.atoms[3].aromatic);
        assert_eq!(molecule.atoms[3].hcount, 1);
        assert_eq!(molecule.atoms[5].symbol, "Sc");
        assert!(!molecule.atoms[5].aromatic);
    }

    #[test]
    fn abbrev_assigned_in_layout() {
        let layout_output = layout_native("{PPh3}C=O").expect("abbrev layout failed");
        assert_eq!(layout_output.atoms.len(), 3);
        assert_eq!(layout_output.atoms[0].abbrev, "PPh3");
        assert_eq!(layout_output.atoms[1].abbrev, "");
    }

    #[test]
    fn abbrev_multiple_in_layout() {
        // [*]C(=O)[*] → atoms: 0=[*](OEt), 1=C, 2=O (branch), 3=[*](NHR)
        let layout_output = layout_native("{OEt}C(=O){NHR}").expect("multi abbrev layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "OEt");
        assert_eq!(layout_output.atoms[3].abbrev, "NHR");
    }

    #[test]
    fn abbrev_style_assigned_in_layout() {
        let layout_output = layout_native("{PPh3|P}C=O").expect("styled abbrev layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "PPh3");
        assert_eq!(layout_output.atoms[0].abbrev_style, "P");
    }

    #[test]
    fn abbrev_anchor_assigned_in_layout() {
        let layout_output = layout_native("{>PPh3}C=O").expect("anchored abbrev layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "PPh3");
        assert_eq!(layout_output.atoms[0].abbrev_anchor, 0);
        assert_eq!(layout_output.atoms[0].abbrev_anchor_len, 1);
    }

    // ── Abbreviation lone-pair modifier (`lp=N`) ─────────────────────────────

    #[test]
    fn abbrev_lp_parses_all_counts() {
        for n in 1u8..=4 {
            let raw = format!("Cl|Cl|lp={n}");
            let label = crate::label::parse_abbreviation_label(&raw).expect("lp parse failed");
            assert_eq!(label.text, "Cl");
            assert_eq!(label.style, "Cl");
            assert_eq!(label.lone_pairs, Some(n));
        }
    }

    #[test]
    fn abbrev_lp_without_style() {
        let label =
            crate::label::parse_abbreviation_label("OR|lp=2").expect("lp without style failed");
        assert_eq!(label.text, "OR");
        assert_eq!(label.style, "");
        assert_eq!(label.lone_pairs, Some(2));
    }

    #[test]
    fn abbrev_lp_with_anchor_marker() {
        let label =
            crate::label::parse_abbreviation_label(">PPh_3|P|lp=1").expect("anchored lp failed");
        assert_eq!(label.text, "PPh_3");
        assert_eq!(label.anchor, 0);
        assert_eq!(label.anchor_len, 1);
        assert_eq!(label.style, "P");
        assert_eq!(label.lone_pairs, Some(1));
    }

    #[test]
    fn abbrev_lp_with_anchor_no_style() {
        let label = crate::label::parse_abbreviation_label(">OR|lp=2")
            .expect("anchored lp no style failed");
        assert_eq!(label.text, "OR");
        assert_eq!(label.anchor_len, 1);
        assert_eq!(label.style, "");
        assert_eq!(label.lone_pairs, Some(2));
    }

    #[test]
    fn abbrev_existing_syntax_preserved() {
        assert_eq!(
            crate::label::parse_abbreviation_label("PPh3")
                .unwrap()
                .lone_pairs,
            None
        );
        assert_eq!(
            crate::label::parse_abbreviation_label("PPh3|P")
                .unwrap()
                .style,
            "P"
        );
        assert_eq!(
            crate::label::parse_abbreviation_label("PPh3|P")
                .unwrap()
                .lone_pairs,
            None
        );
        let anchored = crate::label::parse_abbreviation_label(">PPh3").unwrap();
        assert_eq!(anchored.anchor_len, 1);
        assert_eq!(anchored.lone_pairs, None);
        let anchored_styled = crate::label::parse_abbreviation_label(">PPh3|red").unwrap();
        assert_eq!(anchored_styled.style, "red");
        assert_eq!(anchored_styled.lone_pairs, None);
    }

    #[test]
    fn abbrev_lp_rejects_out_of_range_and_malformed() {
        assert!(crate::label::parse_abbreviation_label("Cl|Cl|lp=0").is_err());
        assert!(crate::label::parse_abbreviation_label("Cl|Cl|lp=5").is_err());
        assert!(crate::label::parse_abbreviation_label("Cl|Cl|lp=-1").is_err());
        assert!(crate::label::parse_abbreviation_label("Cl|Cl|lp=2.5").is_err());
        assert!(crate::label::parse_abbreviation_label("Cl|Cl|lp=two").is_err());
        assert!(crate::label::parse_abbreviation_label("Cl|Cl|lp=").is_err());
        // Duplicate lp modifiers.
        assert!(crate::label::parse_abbreviation_label("Cl|Cl|lp=1|lp=2").is_err());
        // Unknown named modifier.
        assert!(crate::label::parse_abbreviation_label("Cl|Cl|foo=1").is_err());
        // Style after a named modifier.
        assert!(crate::label::parse_abbreviation_label("Cl|lp=1|Cl").is_err());
    }

    #[test]
    fn abbrev_lp_count_reaches_layout_output() {
        let layout_output = layout_native("{>Cl|Cl|lp=3}C").expect("lp layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "Cl");
        assert_eq!(layout_output.atoms[0].lone_pairs, 3);
        // A fallback direction record is emitted per declared pair so lp()
        // references stay resolvable.
        assert_eq!(layout_output.atoms[0].lone_pair_dirs.len(), 3);
    }

    // ── Abbreviation rendering offset (`offset=(x,y)`) ─────────────────────

    #[test]
    fn abbrev_offset_parses_with_style_and_lone_pairs() {
        let label = crate::label::parse_abbreviation_label("H|grey|lp=1|offset=(0.1, -0.2)")
            .expect("offset parse failed");
        assert_eq!(label.text, "H");
        assert_eq!(label.style, "grey");
        assert_eq!(label.lone_pairs, Some(1));
        assert_eq!(label.offset, Some((0.1, -0.2)));

        let reversed = crate::label::parse_abbreviation_label("H|offset=(-.25, .4)|lp=2")
            .expect("reordered modifiers failed");
        assert_eq!(reversed.style, "");
        assert_eq!(reversed.lone_pairs, Some(2));
        assert_eq!(reversed.offset, Some((-0.25, 0.4)));
    }

    #[test]
    fn abbrev_offset_rejects_malformed_values() {
        assert!(crate::label::parse_abbreviation_label("H|offset=0.1,0.2").is_err());
        assert!(crate::label::parse_abbreviation_label("H|offset=(0.1)").is_err());
        assert!(crate::label::parse_abbreviation_label("H|offset=(0.1,0.2,0.3)").is_err());
        assert!(crate::label::parse_abbreviation_label("H|offset=(left,0.2)").is_err());
        assert!(crate::label::parse_abbreviation_label("H|offset=(NaN,0.2)").is_err());
        assert!(
            crate::label::parse_abbreviation_label("H|offset=(0.1,0.2)|offset=(0.3,0.4)").is_err()
        );
    }

    #[test]
    fn abbrev_offset_reaches_output_without_changing_layout_coordinates() {
        let baseline = layout_native("C{H}").expect("baseline layout failed");
        let displaced = layout_native("C{H|offset=(0.3,-0.2)}").expect("offset layout failed");

        assert_eq!(baseline.atoms[0].pos.x, displaced.atoms[0].pos.x);
        assert_eq!(baseline.atoms[0].pos.y, displaced.atoms[0].pos.y);
        assert_eq!(baseline.atoms[1].pos.x, displaced.atoms[1].pos.x);
        assert_eq!(baseline.atoms[1].pos.y, displaced.atoms[1].pos.y);
        assert_eq!(baseline.bbox_width, displaced.bbox_width);
        assert_eq!(baseline.bbox_height, displaced.bbox_height);
        assert_eq!(displaced.atoms[1].abbrev_offset_x, 0.3);
        assert_eq!(displaced.atoms[1].abbrev_offset_y, -0.2);
    }

    #[test]
    fn abbrev_without_lp_has_no_pairs() {
        let layout_output = layout_native("{PPh3}C=O").expect("no-lp layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "PPh3");
        assert_eq!(layout_output.atoms[0].lone_pairs, 0);
        assert_eq!(layout_output.atoms[0].lone_pair_dirs.len(), 0);
    }

    #[test]
    fn bracket_atom_is_not_literal_label() {
        let bracket = layout_native("[N]").expect("bracket nitrogen failed");
        let label = layout_native("{N}").expect("literal label failed");
        assert_eq!(bracket.atoms[0].symbol, "N");
        assert_eq!(bracket.atoms[0].abbrev, "");
        assert_eq!(label.atoms[0].symbol, "*");
        assert_eq!(label.atoms[0].abbrev, "N");
    }

    // ── Virtual H atoms for bracket-notation hydrogens ───────────────────────

    #[test]
    fn bracket_h_group_is_one_addressable_index() {
        // [OH-]: one O heavy atom + one virtual H group → 2 atoms total.
        let layout_output = layout_native("[OH-]").expect("[OH-] layout failed");
        assert_eq!(layout_output.atoms.len(), 2);
        assert_eq!(layout_output.atoms[0].symbol, "O");
        assert!(layout_output.atoms[1].virtual_h);
        assert_eq!(layout_output.atoms[1].symbol, "H");
        assert_eq!(layout_output.bonds.len(), 1);
        assert!(layout_output.bonds[0].virtual_bond);
    }

    #[test]
    fn ammonium_h_group_is_one_index() {
        // [NH4+]: despite hcount=4 the H-label is one glyph → exactly 1 virtual H.
        let layout_output = layout_native("[NH4+]").expect("[NH4+] layout failed");
        assert_eq!(layout_output.atoms.len(), 2); // 1 N + 1 virtual H group
        assert!(layout_output.atoms[1].virtual_h);
        assert_eq!(layout_output.bonds.len(), 1);
        assert!(layout_output.bonds[0].virtual_bond);
    }

    #[test]
    fn virtual_h_has_valid_position() {
        let layout_output = layout_native("[OH-]").expect("[OH-] layout failed");
        let oxygen_position = layout_output.atoms[0].pos;
        let hydrogen_position = layout_output.atoms[1].pos;
        let distance = ((hydrogen_position.x - oxygen_position.x).powi(2)
            + (hydrogen_position.y - oxygen_position.y).powi(2))
        .sqrt();
        assert!(
            (distance - 0.35).abs() < 1e-6,
            "H should be 0.35 bond lengths from O, got {distance}"
        );
    }

    // ── Implicit H for expanded valence table ────────────────────────────────

    #[test]
    fn implicit_h_boron() {
        // B (organic subset, valence 3): B bonded to one C → 2 implicit H
        let layout_output = layout_native("BC").expect("boron layout failed");
        assert_eq!(layout_output.atoms[0].implicit_h, 2);
    }

    #[test]
    fn lone_pair_counts_common_organic_atoms() {
        let alcohol = layout_native("CCO").expect("alcohol layout failed");
        assert_eq!(alcohol.atoms[2].lone_pairs, 2);
        assert_eq!(alcohol.atoms[2].lone_pair_dirs.len(), 2);

        let amine = layout_native("CCN").expect("amine layout failed");
        assert_eq!(amine.atoms[2].lone_pairs, 1);

        let chloride = layout_native("CCl").expect("chloride layout failed");
        assert_eq!(chloride.atoms[1].lone_pairs, 3);
    }

    #[test]
    fn lone_pair_counts_respect_charge_and_explicit_hydrogens() {
        let ammonium = layout_native("[NH4+]").expect("ammonium layout failed");
        assert_eq!(ammonium.atoms[0].lone_pairs, 0);

        let formate = layout_native("[O-]C=O").expect("formate layout failed");
        assert_eq!(formate.atoms[0].lone_pairs, 3);
        assert_eq!(formate.atoms[2].lone_pairs, 2);
    }

    fn max_bond_length(layout_output: &LayoutOutput) -> f64 {
        layout_output
            .bonds
            .iter()
            .map(|bond| {
                let from = layout_output.atoms[bond.from].pos;
                let to = layout_output.atoms[bond.to].pos;
                from.dist(to)
            })
            .fold(0.0, f64::max)
    }

    /// Smallest distance between any two atoms that are not bonded to each other.
    /// A value well below the ~1.0 bond length signals two parts of the molecule
    /// being laid out on top of one another.
    fn min_nonbonded_distance(layout_output: &LayoutOutput) -> f64 {
        let bonded: std::collections::HashSet<(usize, usize)> = layout_output
            .bonds
            .iter()
            .map(|b| (b.from.min(b.to), b.from.max(b.to)))
            .collect();
        let mut min = f64::INFINITY;
        for i in 0..layout_output.atoms.len() {
            for j in (i + 1)..layout_output.atoms.len() {
                if bonded.contains(&(i, j)) {
                    continue;
                }
                min = min.min(layout_output.atoms[i].pos.dist(layout_output.atoms[j].pos));
            }
        }
        min
    }

    fn atoms_are_collinear(
        layout_output: &LayoutOutput,
        first_atom: usize,
        middle_atom: usize,
        last_atom: usize,
    ) -> bool {
        let first_position = layout_output.atoms[first_atom].pos;
        let middle_position = layout_output.atoms[middle_atom].pos;
        let last_position = layout_output.atoms[last_atom].pos;
        let cross = (middle_position.x - first_position.x) * (last_position.y - middle_position.y)
            - (middle_position.y - first_position.y) * (last_position.x - middle_position.x);
        cross.abs() < 1e-8
    }

    /// End-to-end check: reconstructs the depicted 3D geometry from the rendered
    /// output and verifies every stereocenter's signed volume matches `@`/`@@`.
    fn chirality_matches_smiles(smiles: &str) -> bool {
        use crate::layout::implicit_h_count;
        use crate::stereo::{signed_volume, stereochemical_hydrogen_direction};

        let molecule = parse_molecule(smiles).expect("graph build failed");
        let layout_output = compute_layout(&molecule).expect("layout failed");
        let coordinates: Vec<crate::render::Vec2> =
            layout_output.atoms.iter().map(|atom| atom.pos).collect();

        for center in 0..molecule.atoms.len() {
            let parity = match layout_output.atoms[center].chirality.as_str() {
                "tetra_anti" => -1.0_f64,
                "tetra_clockwise" => 1.0_f64,
                _ => continue,
            };
            let neighbor_bonds = &molecule.neighbor_bonds[center];
            let hydrogen_count =
                (molecule.atoms[center].hcount + implicit_h_count(&molecule, center)) as usize;
            if neighbor_bonds.len() + hydrogen_count != 4 || hydrogen_count > 1 {
                continue;
            }

            // Neighbor order with the implicit hydrogen inserted (same rule as the
            // renderer): after the "from" atom, or first if there is none.
            #[derive(Clone, Copy)]
            enum TetrahedralNeighbor {
                Bond(usize, usize),
                Hydrogen,
            }
            let mut neighbor_order: Vec<TetrahedralNeighbor> = neighbor_bonds
                .iter()
                .map(|&bond_index| {
                    let bond = &molecule.bonds[bond_index];
                    let other = if bond.from == center {
                        bond.to
                    } else {
                        bond.from
                    };
                    TetrahedralNeighbor::Bond(bond_index, other)
                })
                .collect();
            if hydrogen_count == 1 {
                let insertion_position = if molecule.has_preceding[center] { 1 } else { 0 };
                neighbor_order.insert(
                    insertion_position.min(neighbor_order.len()),
                    TetrahedralNeighbor::Hydrogen,
                );
            }

            // Which bond (if any) carries the rendered wedge, and its z sign.
            let wedge_bond = neighbor_bonds
                .iter()
                .find(|&&bond_index| layout_output.bonds[bond_index].stereo != "none")
                .copied();
            let hydrogen_direction =
                stereochemical_hydrogen_direction(&molecule, center, &coordinates, None);

            let normalize = |offset_x: f64, offset_y: f64| {
                let length = (offset_x * offset_x + offset_y * offset_y).sqrt();
                if length > 1e-12 {
                    (offset_x / length, offset_y / length)
                } else {
                    (offset_x, offset_y)
                }
            };
            let stereo_depth = |stereo: &str| if stereo == "wedge_up" { 1.0 } else { -1.0 };

            let mut directions = [[0.0_f64; 3]; 4];
            for (neighbor_index, neighbor) in neighbor_order.iter().enumerate() {
                directions[neighbor_index] = match neighbor {
                    TetrahedralNeighbor::Bond(bond_index, other) => {
                        let (direction_x, direction_y) = normalize(
                            coordinates[*other].x - coordinates[center].x,
                            coordinates[*other].y - coordinates[center].y,
                        );
                        let depth = if Some(*bond_index) == wedge_bond {
                            stereo_depth(&layout_output.bonds[*bond_index].stereo)
                        } else {
                            0.0
                        };
                        [direction_x, direction_y, depth]
                    }
                    TetrahedralNeighbor::Hydrogen => {
                        // If the H itself is wedged, use that; otherwise it sits on
                        // the face opposite the wedged heavy substituent.
                        let depth = if layout_output.atoms[center].stereo_h != "none" {
                            stereo_depth(&layout_output.atoms[center].stereo_h)
                        } else if let Some(bond_index) = wedge_bond {
                            -stereo_depth(&layout_output.bonds[bond_index].stereo)
                        } else {
                            0.0
                        };
                        let direction = if layout_output.atoms[center].stereo_h != "none" {
                            layout_output.atoms[center].stereo_h_dir
                        } else {
                            hydrogen_direction
                        };
                        [direction.x, direction.y, depth]
                    }
                };
            }

            let volume = signed_volume(&directions);
            if volume.abs() < 1e-9 || volume.signum() != parity {
                return false;
            }
        }
        true
    }

    fn chiral_oxygen_bond_has_stereo(layout_output: &LayoutOutput) -> bool {
        layout_output.bonds.iter().any(|bond| {
            let from = &layout_output.atoms[bond.from];
            let to = &layout_output.atoms[bond.to];
            bond.stereo != "none"
                && ((from.chirality != "none" && to.symbol == "O")
                    || (to.chirality != "none" && from.symbol == "O"))
        })
    }

    fn alkene_substituent_side_product(layout_output: &LayoutOutput) -> i8 {
        let double = layout_output
            .bonds
            .iter()
            .find(|bond| bond.order == 2)
            .unwrap();
        let first_alkene_atom = double.from;
        let second_alkene_atom = double.to;
        let left = layout_output
            .bonds
            .iter()
            .find(|bond| {
                bond.order == 1 && (bond.from == first_alkene_atom || bond.to == first_alkene_atom)
            })
            .unwrap();
        let right = layout_output
            .bonds
            .iter()
            .find(|bond| {
                bond.order == 1
                    && (bond.from == second_alkene_atom || bond.to == second_alkene_atom)
            })
            .unwrap();
        let left_neighbor = if left.from == first_alkene_atom {
            left.to
        } else {
            left.from
        };
        let right_neighbor = if right.from == second_alkene_atom {
            right.to
        } else {
            right.from
        };
        side(
            layout_output.atoms[first_alkene_atom].pos,
            layout_output.atoms[second_alkene_atom].pos,
            layout_output.atoms[left_neighbor].pos,
        ) * side(
            layout_output.atoms[first_alkene_atom].pos,
            layout_output.atoms[second_alkene_atom].pos,
            layout_output.atoms[right_neighbor].pos,
        )
    }

    /// Returns 1 when the two substituents sit on the same side of the
    /// double bond (cis) and -1 when they sit on opposite sides (trans).
    fn double_bond_substituent_side_product(
        layout_output: &LayoutOutput,
        first_substituent: usize,
        double_bond_atoms: (usize, usize),
        second_substituent: usize,
    ) -> i8 {
        let line_start = layout_output.atoms[double_bond_atoms.0].pos;
        let line_end = layout_output.atoms[double_bond_atoms.1].pos;
        side(
            line_start,
            line_end,
            layout_output.atoms[first_substituent].pos,
        ) * side(
            line_start,
            line_end,
            layout_output.atoms[second_substituent].pos,
        )
    }

    fn side(
        line_start: crate::render::Vec2,
        line_end: crate::render::Vec2,
        point: crate::render::Vec2,
    ) -> i8 {
        let cross = (line_end.x - line_start.x) * (point.y - line_start.y)
            - (line_end.y - line_start.y) * (point.x - line_start.x);
        if cross > 1e-8 {
            1
        } else {
            -1
        }
    }

    fn turn_cross(a: crate::render::Vec2, b: crate::render::Vec2, c: crate::render::Vec2) -> f64 {
        (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x)
    }

    #[test]
    fn curl_repeats_the_previous_chain_turn() {
        let normal = layout_native("CCCC").expect("normal chain failed");
        let curled = layout_native("CCC!cC").expect("curled chain failed");

        let normal_first = turn_cross(
            normal.atoms[0].pos,
            normal.atoms[1].pos,
            normal.atoms[2].pos,
        );
        let normal_second = turn_cross(
            normal.atoms[1].pos,
            normal.atoms[2].pos,
            normal.atoms[3].pos,
        );
        assert!(normal_first * normal_second < 0.0);

        let curl_first = turn_cross(
            curled.atoms[0].pos,
            curled.atoms[1].pos,
            curled.atoms[2].pos,
        );
        let curl_second = turn_cross(
            curled.atoms[1].pos,
            curled.atoms[2].pos,
            curled.atoms[3].pos,
        );
        assert!(curl_first * curl_second > 0.0);
    }

    #[test]
    fn curl_swaps_forward_branch_slots_without_overlap() {
        let layout_output = layout_native("CCC({PPh3})!cC(=O)OCC").expect("branched curl failed");
        let first = turn_cross(
            layout_output.atoms[0].pos,
            layout_output.atoms[1].pos,
            layout_output.atoms[2].pos,
        );
        let second = turn_cross(
            layout_output.atoms[1].pos,
            layout_output.atoms[2].pos,
            layout_output.atoms[4].pos,
        );
        assert!(first * second > 0.0);
        assert!(min_atom_distance(&layout_output) >= 0.5);
    }

    #[test]
    fn consecutive_curls_repeat_each_new_turn() {
        let layout_output = layout_native("CCCC!cC!cC").expect("consecutive curls failed");
        let turns = (0..3)
            .map(|i| {
                turn_cross(
                    layout_output.atoms[i + 1].pos,
                    layout_output.atoms[i + 2].pos,
                    layout_output.atoms[i + 3].pos,
                )
            })
            .collect::<Vec<_>>();
        assert!(turns.iter().all(|turn| turns[0] * turn > 0.0));
        assert!(min_atom_distance(&layout_output) >= 0.5);
    }

    #[test]
    fn crowded_substituted_curls_keep_atoms_separated() {
        let molecules = [
            "CCC(C(N)C)!cC(OC(F)C)!cCC",
            "CCC({PPh3})!cC([O-])!cC(=O)OCC",
            "CCC(CC(C)C)!cC(OC)C(NC)CC",
        ];

        for smiles in molecules {
            let layout_output =
                layout_native(smiles).unwrap_or_else(|err| panic!("{smiles}: {err}"));
            let min = min_atom_distance(&layout_output);
            assert!(min >= 0.5, "{smiles}: min atom distance {min:.3}");
        }
    }

    #[test]
    fn curl_combines_with_wedge_and_double_bonds() {
        let wedge = layout_native("CCC!c!wN").expect("curled wedge failed");
        assert_eq!(wedge.bonds[2].stereo, "wedge_up");

        let double = layout_native("CCC!c=C").expect("curled double failed");
        assert_eq!(double.bonds[2].order, 2);
        let first = turn_cross(
            double.atoms[0].pos,
            double.atoms[1].pos,
            double.atoms[2].pos,
        );
        let second = turn_cross(
            double.atoms[1].pos,
            double.atoms[2].pos,
            double.atoms[3].pos,
        );
        assert!(first * second > 0.0);
    }

    #[test]
    fn curl_requires_an_established_turn() {
        let err = layout_native("CC!cC").expect_err("early curl should fail");
        assert!(err.contains("two preceding chain bonds"));
    }

    fn side_of_line(
        start: crate::render::Vec2,
        end: crate::render::Vec2,
        point: crate::render::Vec2,
    ) -> f64 {
        (end.x - start.x) * (point.y - start.y) - (end.y - start.y) * (point.x - start.x)
    }

    #[test]
    fn curled_second_double_bond_gives_s_cis_butadiene() {
        let s_cis = layout_native("C=CC!c=C").expect("s-cis butadiene failed");
        assert_eq!(s_cis.bonds[0].order, 2);
        assert_eq!(s_cis.bonds[2].order, 2);
        let central_start = s_cis.atoms[1].pos;
        let central_end = s_cis.atoms[2].pos;
        let terminal_sides = side_of_line(central_start, central_end, s_cis.atoms[0].pos)
            * side_of_line(central_start, central_end, s_cis.atoms[3].pos);
        assert!(terminal_sides > 0.0, "termini should sit on the same side");
        assert!(min_atom_distance(&s_cis) >= 0.5);
    }

    #[test]
    fn uncurled_butadiene_stays_s_trans() {
        let s_trans = layout_native("C=CC=C").expect("butadiene failed");
        let central_start = s_trans.atoms[1].pos;
        let central_end = s_trans.atoms[2].pos;
        let terminal_sides = side_of_line(central_start, central_end, s_trans.atoms[0].pos)
            * side_of_line(central_start, central_end, s_trans.atoms[3].pos);
        assert!(terminal_sides < 0.0, "termini should sit on opposite sides");
    }

    #[test]
    fn curl_without_reference_turn_names_the_s_cis_diene() {
        let err = layout_native("C=C!cC=C").expect_err("curl without a reference turn should fail");
        assert!(err.contains("two preceding chain bonds"));
        assert!(
            err.contains("C=CC!c=C"),
            "error should show the s-cis spelling: {err}"
        );
    }

    // ── Molecular weight ──────────────────────────────────────────────────────
    //
    // Reference values are PubChem's computed molecular weights, which use the
    // IUPAC/CIAAW standard atomic weights (the same table ptable embeds via
    // PubChemElements_all.json).

    fn assert_weight(smiles: &str, expected: f64) {
        let molecular_weight = mol_weight_native(smiles).expect("mol weight failed");
        assert!(
            (molecular_weight - expected).abs() < 0.01,
            "mol_weight({smiles}) = {molecular_weight}, expected {expected}"
        );
    }

    fn assert_formula(smiles: &str, expected: &str) {
        let formula = mol_formula_native(smiles).expect("molecular formula failed");
        assert_eq!(formula, expected, "mol_formula({smiles})");
    }

    /// Extracts every SMILES string literal passed to `smiles("...")` or
    /// `mol("...")` in the visual test file, undoing Typst string escapes.
    fn test_typ_smiles_strings() -> Vec<String> {
        let src =
            std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/test.typ"))
                .expect("tests/test.typ not found");
        let src: String = src
            .lines()
            .filter(|l| !l.trim_start().starts_with("//"))
            .collect::<Vec<_>>()
            .join("\n");
        let mut found = Vec::new();
        for opener in ["smiles(\"", "mol(\""] {
            let mut rest = src.as_str();
            while let Some(pos) = rest.find(opener) {
                rest = &rest[pos + opener.len()..];
                let mut literal = String::new();
                let mut chars = rest.chars();
                while let Some(ch) = chars.next() {
                    match ch {
                        '"' => break,
                        '\\' => {
                            if let Some(esc) = chars.next() {
                                literal.push(esc);
                            }
                        }
                        _ => literal.push(ch),
                    }
                }
                found.push(literal);
            }
        }
        found.sort();
        found.dedup();
        found
    }

    #[test]
    fn every_molecule_in_test_typ_has_no_overlapping_atoms() {
        let molecules = test_typ_smiles_strings();
        assert!(
            molecules.len() > 50,
            "extraction looks broken: {molecules:?}"
        );
        let mut failures = Vec::new();
        for m in &molecules {
            match layout_native(m) {
                Ok(layout_output) => {
                    let min = min_atom_distance(&layout_output);
                    if layout_output.atoms.len() > 1 && min < 0.5 {
                        failures.push(format!("{m}: min atom distance {min:.3}"));
                    }
                }
                Err(e) => failures.push(format!("{m}: layout failed: {e}")),
            }
        }
        assert!(
            failures.is_empty(),
            "molecules with overlaps or errors:\n{}",
            failures.join("\n")
        );
    }

    // ── Aromatic ring circles ────────────────────────────────────────────────

    #[test]
    fn aromatic_input_emits_ring_circles() {
        let layout_output = layout_native("c1ccccc1").expect("benzene layout failed");
        assert_eq!(layout_output.aromatic_rings.len(), 1);
        let ring = &layout_output.aromatic_rings[0];
        // Hexagon with unit bonds: inradius ~0.866, so radius ~0.62.
        assert!((ring.radius - 0.866 * 0.72).abs() < 0.05);
        assert!(
            layout_output
                .bonds
                .iter()
                .filter(|bond| bond.aromatic)
                .count()
                == 6
        );
    }

    #[test]
    fn kekule_input_emits_no_ring_circles() {
        let layout_output = layout_native("C1=CC=CC=C1").expect("benzene layout failed");
        assert!(layout_output.aromatic_rings.is_empty());
        assert!(layout_output.bonds.iter().all(|b| !b.aromatic));
    }

    #[test]
    fn fused_aromatics_emit_one_circle_per_ring() {
        let layout_output = layout_native("c1ccc2ccccc2c1").expect("naphthalene layout failed");
        assert_eq!(layout_output.aromatic_rings.len(), 2);
    }

    #[test]
    fn aromatic_ring_with_saturated_neighbor_ring() {
        // Indane: only the aromatic ring gets a circle.
        let layout_output = layout_native("c1ccc2CCCc2c1").expect("indane layout failed");
        assert_eq!(layout_output.aromatic_rings.len(), 1);
    }

    #[test]
    fn mol_weight_water() {
        // PubChem CID 962: 18.015 g/mol
        assert_weight("O", 18.015);
    }

    #[test]
    fn mol_weight_ethanol() {
        // PubChem CID 702: 46.07 g/mol
        assert_weight("CCO", 46.069);
    }

    #[test]
    fn mol_weight_benzene_aromatic_input() {
        // PubChem CID 241: 78.11 g/mol; aromatic input exercises kekulization.
        assert_weight("c1ccccc1", 78.114);
    }

    #[test]
    fn mol_weight_glucose() {
        // PubChem CID 5793: 180.16 g/mol
        assert_weight("C(C1C(C(C(C(O1)O)O)O)O)O", 180.156);
    }

    #[test]
    fn mol_weight_pubchem_reference_set() {
        // Connectivity SMILES and reference weights queried from PubChem's
        // compound-property records; CIDs make each regression easy to audit.
        let examples = [
            ("methane", "C", 16.043, 297),
            ("ethane", "CC", 30.07, 6324),
            ("propane", "CCC", 44.10, 6334),
            ("butane", "CCCC", 58.12, 7843),
            ("isobutane", "CC(C)C", 58.12, 6360),
            ("methanol", "CO", 32.042, 887),
            ("isopropanol", "CC(C)O", 60.10, 3776),
            ("acetone", "CC(=O)C", 58.08, 180),
            ("acetic acid", "CC(=O)O", 60.05, 176),
            ("formaldehyde", "C=O", 30.026, 712),
            ("acetaldehyde", "CC=O", 44.05, 177),
            ("formic acid", "C(=O)O", 46.025, 284),
            ("phenol", "C1=CC=C(C=C1)O", 94.11, 996),
            ("toluene", "CC1=CC=CC=C1", 92.14, 1140),
            ("pyridine", "C1=CC=NC=C1", 79.10, 1049),
            ("cyclohexane", "C1CCCCC1", 84.16, 8078),
            ("naphthalene", "C1=CC=C2C=CC=CC2=C1", 128.17, 931),
            ("ethyl acetate", "CCOC(=O)C", 88.11, 8857),
            ("chloroform", "C(Cl)(Cl)Cl", 119.37, 6212),
            ("dichloromethane", "C(Cl)Cl", 84.93, 6344),
            ("hydrogen peroxide", "OO", 34.015, 784),
            ("ammonia", "N", 17.031, 222),
            ("glycine", "C(C(=O)O)N", 75.07, 750),
            ("alanine", "CC(C(=O)O)N", 89.09, 5950),
            ("urea", "C(=O)(N)N", 60.056, 1176),
            ("acetamide", "CC(=O)N", 59.07, 178),
            ("aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O", 180.16, 2244),
            ("ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", 206.28, 3672),
            ("acetaminophen", "CC(=O)NC1=CC=C(C=C1)O", 151.16, 1983),
            ("citric acid", "C(C(=O)O)C(CC(=O)O)(C(=O)O)O", 192.12, 311),
            (
                "sucrose",
                "C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O",
                342.30,
                5988,
            ),
        ];

        for (name, smiles, expected, cid) in examples {
            let molecular_weight = mol_weight_native(smiles)
                .unwrap_or_else(|error| panic!("{name} (PubChem CID {cid}, {smiles}): {error}"));
            assert!(
                (molecular_weight - expected).abs() < 0.01,
                "{name} (PubChem CID {cid}, {smiles}) = {molecular_weight}, expected {expected}"
            );
        }
    }

    #[test]
    fn mol_weight_caffeine() {
        // PubChem CID 2519: 194.19 g/mol
        assert_weight("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 194.19);
    }

    #[test]
    fn mol_weight_sodium_acetate_dot_fragments() {
        // PubChem CID 517045: 82.03 g/mol; dot-separated ion pair sums both
        // fragments (the electron mass difference of the ions is ignored, as
        // in standard formula-weight arithmetic).
        assert_weight("CC(=O)[O-].[Na+]", 82.034);
    }

    #[test]
    fn mol_weight_ammonium_bracket_h() {
        // PubChem CID 223: 18.039 g/mol; explicit bracket hydrogens counted.
        assert_weight("[NH4+]", 18.039);
    }

    #[test]
    fn mol_weight_wildcard_errors() {
        let err = mol_weight_native("*CC").expect_err("wildcard should fail");
        assert!(err.contains("wildcard"));
    }

    #[test]
    fn mol_weight_abbreviation_errors() {
        let err = mol_weight_native("{PPh3}C=O").expect_err("abbreviation should fail");
        assert!(err.contains("PPh3"));
    }

    #[test]
    fn mol_weight_isotope_errors() {
        let err = mol_weight_native("[2H]O[2H]").expect_err("isotope should fail");
        assert!(err.contains("isotope"));
    }

    #[test]
    fn molecular_formula_matches_pubchem_ethanol() {
        // PubChem CID 702: CCO -> C2H6O.
        assert_formula("CCO", "C2H6O");
    }

    #[test]
    fn molecular_formula_matches_pubchem_aromatic_benzene() {
        // PubChem CID 241: c1ccccc1 -> C6H6.
        assert_formula("c1ccccc1", "C6H6");
    }

    #[test]
    fn molecular_formula_matches_pubchem_caffeine() {
        // PubChem CID 2519: CN1C=NC2=C1C(=O)N(C(=O)N2C)C -> C8H10N4O2.
        assert_formula("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C8H10N4O2");
    }

    #[test]
    fn molecular_formula_matches_pubchem_sodium_acetate() {
        // PubChem CID 517045: CC(=O)[O-].[Na+] -> C2H3NaO2.
        assert_formula("CC(=O)[O-].[Na+]", "C2H3NaO2");
    }

    #[test]
    fn molecular_formula_matches_pubchem_ammonium_chloride() {
        // PubChem CID 25517: [NH4+].[Cl-] -> ClH4N.
        assert_formula("[NH4+].[Cl-]", "ClH4N");
    }

    #[test]
    fn molecular_formula_counts_explicit_and_implicit_hydrogens() {
        assert_formula("O", "H2O");
        assert_formula("[NH4+]", "H4N^+");
        assert_formula("C([H])([H])([H])[H]", "CH4");
    }

    #[test]
    fn molecular_formula_rejects_undefined_composition() {
        for (smiles, expected_error) in [
            ("*CC", "wildcard"),
            ("{PPh3}C=O", "PPh3"),
            ("[2H]O[2H]", "isotope"),
        ] {
            let error = mol_formula_native(smiles).expect_err("formula should fail");
            assert!(error.contains(expected_error), "{smiles}: {error}");
        }
    }
}
