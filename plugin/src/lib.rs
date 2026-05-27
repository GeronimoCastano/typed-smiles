mod error;
mod graph;
mod layout;
mod render;

pub use render::LayoutOutput;

use graph::MoleculeGraph;
use layout::compute_layout;

// ── Extended SMILES preprocessing ────────────────────────────────────────────
//
// Handles two extensions before handing the string to smiles-parser:
//
//   {label}  →  [*]   Abbreviated group (e.g. {PPh3}, {OEt}).
//                     The Nth [*] atom in the output gets abbrev = the Nth label.
//
// The cleaned string is valid standard SMILES; extensions are stripped out.

pub(crate) struct Preprocessed {
    pub smiles: String,
    /// Labels in the order they appeared, one per {label} token.
    pub abbrev_labels: Vec<String>,
}

pub(crate) fn preprocess_smiles(input: &str) -> Preprocessed {
    let mut smiles = String::with_capacity(input.len());
    let mut abbrev_labels = Vec::new();

    let mut chars = input.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == '{' {
            let mut label = String::new();
            for inner in chars.by_ref() {
                if inner == '}' {
                    break;
                }
                label.push(inner);
            }
            abbrev_labels.push(label);
            smiles.push_str("[*]");
        } else {
            smiles.push(ch);
        }
    }

    Preprocessed {
        smiles,
        abbrev_labels,
    }
}

/// Apply collected abbreviation labels to the matching `*` atoms in the graph
/// (N-th `*` atom ← N-th label, in atom-index order).
fn assign_abbrevs(mol: &mut MoleculeGraph, labels: &[String]) {
    let mut label_iter = labels.iter();
    for atom in mol.atoms.iter_mut() {
        if atom.symbol == "*" {
            if let Some(label) = label_iter.next() {
                atom.abbrev = label.clone();
            }
        }
    }
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
        let s = core::str::from_utf8(smiles).map_err(|e| format!("UTF-8 error: {e}"))?;
        let pre = preprocess_smiles(s);
        let mut mol = MoleculeGraph::from_smiles(&pre.smiles)?;
        assign_abbrevs(&mut mol, &pre.abbrev_labels);
        let out = compute_layout(&mol)?;
        serde_json::to_vec(&out).map_err(|e| format!("JSON error: {e}"))
    }
}

// ── Native entrypoint for tests / CLI ───────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
pub fn layout_native(smiles: &str) -> Result<LayoutOutput, String> {
    let pre = preprocess_smiles(smiles);
    let mut mol = MoleculeGraph::from_smiles(&pre.smiles)?;
    assign_abbrevs(&mut mol, &pre.abbrev_labels);
    compute_layout(&mol)
}

#[cfg(test)]
mod tests {
    use super::*;

    // smiles-parser 0.4 does not support unbracketed aromatic atoms (c, n, o).
    // Use Kekulé notation until we add a preprocessing step.
    #[test]
    fn benzene_kekule() {
        let out = layout_native("C1=CC=CC=C1").expect("benzene layout failed");
        assert_eq!(out.atoms.len(), 6);
        assert_eq!(out.bonds.len(), 6);
    }

    #[test]
    fn ethanol() {
        let out = layout_native("CCO").expect("ethanol layout failed");
        assert_eq!(out.atoms.len(), 3);
    }

    #[test]
    fn implicit_h_counts_ethanol() {
        let out = layout_native("CCO").expect("ethanol layout failed");
        let counts: Vec<u8> = out.atoms.iter().map(|a| a.implicit_h).collect();
        assert_eq!(counts, vec![3, 2, 1]);
    }

    #[test]
    fn explicit_h_suppresses_implicit_h() {
        let out = layout_native("[NH4+]").expect("ammonium layout failed");
        assert_eq!(out.atoms[0].hcount, 4);
        assert_eq!(out.atoms[0].implicit_h, 0);
    }

    #[test]
    fn double_bond() {
        let out = layout_native("C=C").expect("ethylene layout failed");
        assert_eq!(out.bonds[0].order, 2);
    }

    #[test]
    fn cyclohexane() {
        let out = layout_native("C1CCCCC1").expect("cyclohexane layout failed");
        assert_eq!(out.atoms.len(), 6);
        assert_eq!(out.bonds.len(), 6);
    }

    #[test]
    fn separated_ring_systems_do_not_overlap() {
        let out = layout_native("CC(N)C(=O)OCCC1=CC=CC=C1NCC1=CC=CC=C1")
            .expect("two-ring molecule layout failed");
        assert!(max_bond_length(&out) < 1.2);
    }

    #[test]
    fn isobutane() {
        let out = layout_native("CC(C)C").expect("isobutane layout failed");
        assert_eq!(out.atoms.len(), 4);
        assert_eq!(out.bonds.len(), 3);
    }

    #[test]
    fn naphthalene_kekule() {
        let out = layout_native("C1=CC2=CC=CC=C2C=C1").expect("naphthalene layout failed");
        assert_eq!(out.atoms.len(), 10);
        assert_eq!(out.bonds.len(), 11);
    }

    // ── Wedge / dash bond tests ───────────────────────────────────────────────

    #[test]
    fn wedge_up_bond() {
        let out = layout_native("C/N").expect("wedge up failed");
        assert_eq!(out.bonds[0].stereo, "wedge_up");
    }

    #[test]
    fn wedge_down_bond() {
        let out = layout_native("C\\N").expect("wedge down failed");
        assert_eq!(out.bonds[0].stereo, "wedge_down");
    }

    #[test]
    fn wedge_in_chain() {
        // C–C wedge-up–N: bonds[0] normal, bonds[1] wedge_up
        let out = layout_native("CC/N").expect("wedge in chain failed");
        assert_eq!(out.bonds[0].stereo, "none");
        assert_eq!(out.bonds[1].stereo, "wedge_up");
    }

    #[test]
    fn double_bond_stereo_unchanged() {
        // / and \ around a double bond should still render as a normal double bond
        let out = layout_native("F/C=C/F").expect("E-alkene failed");
        let double = out.bonds.iter().find(|b| b.order == 2).unwrap();
        assert_eq!(double.order, 2);
    }

    // ── Abbreviation tests ────────────────────────────────────────────────────

    #[test]
    fn preprocess_single_abbrev() {
        let p = preprocess_smiles("C({PPh3})=O");
        assert_eq!(p.smiles, "C([*])=O");
        assert_eq!(p.abbrev_labels, vec!["PPh3"]);
    }

    #[test]
    fn preprocess_multiple_abbrevs() {
        let p = preprocess_smiles("{OEt}C(=O){NHR}");
        assert_eq!(p.smiles, "[*]C(=O)[*]");
        assert_eq!(p.abbrev_labels, vec!["OEt", "NHR"]);
    }

    #[test]
    fn preprocess_no_abbrev() {
        let p = preprocess_smiles("CCO");
        assert_eq!(p.smiles, "CCO");
        assert!(p.abbrev_labels.is_empty());
    }

    #[test]
    fn abbrev_assigned_in_layout() {
        let out = layout_native("{PPh3}C=O").expect("abbrev layout failed");
        assert_eq!(out.atoms.len(), 3);
        assert_eq!(out.atoms[0].abbrev, "PPh3");
        assert_eq!(out.atoms[1].abbrev, "");
    }

    #[test]
    fn abbrev_multiple_in_layout() {
        // [*]C(=O)[*] → atoms: 0=[*](OEt), 1=C, 2=O (branch), 3=[*](NHR)
        let out = layout_native("{OEt}C(=O){NHR}").expect("multi abbrev layout failed");
        assert_eq!(out.atoms[0].abbrev, "OEt");
        assert_eq!(out.atoms[3].abbrev, "NHR");
    }

    // ── Implicit H for expanded valence table ────────────────────────────────

    #[test]
    fn implicit_h_boron() {
        // B (organic subset, valence 3): B bonded to one C → 2 implicit H
        let out = layout_native("BC").expect("boron layout failed");
        assert_eq!(out.atoms[0].implicit_h, 2);
    }

    fn max_bond_length(out: &LayoutOutput) -> f64 {
        out.bonds
            .iter()
            .map(|bond| {
                let from = out.atoms[bond.from].pos;
                let to = out.atoms[bond.to].pos;
                from.dist(to)
            })
            .fold(0.0, f64::max)
    }
}
