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
    pub abbrev_labels: Vec<AbbrevLabel>,
    /// One flag for each `/` or `\` token in `smiles`.
    /// `true` means the token came from a typed-smiles `!w`/`!h` drawing extension.
    pub forced_direction_markers: Vec<bool>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct AbbrevLabel {
    pub text: String,
    pub style: String,
}

pub(crate) fn preprocess_smiles(input: &str) -> Preprocessed {
    let mut smiles = String::with_capacity(input.len());
    let mut abbrev_labels = Vec::new();
    let mut forced_direction_markers = Vec::new();

    let mut chars = input.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == '{' {
            let mut raw_label = String::new();
            for inner in chars.by_ref() {
                if inner == '}' {
                    break;
                }
                raw_label.push(inner);
            }
            let (text, style) = raw_label
                .split_once('|')
                .map(|(text, style)| (text.to_string(), style.trim().to_string()))
                .unwrap_or_else(|| (raw_label, String::new()));
            abbrev_labels.push(AbbrevLabel { text, style });
            smiles.push_str("[*]");
        } else if ch == '!' {
            match chars.peek().copied() {
                Some('w') => {
                    chars.next();
                    smiles.push('/');
                    forced_direction_markers.push(true);
                }
                Some('h') => {
                    chars.next();
                    smiles.push('\\');
                    forced_direction_markers.push(true);
                }
                _ => {
                    smiles.push(ch);
                }
            }
        } else {
            if ch == '/' || ch == '\\' {
                forced_direction_markers.push(false);
            }
            smiles.push(ch);
        }
    }

    Preprocessed {
        smiles,
        abbrev_labels,
        forced_direction_markers,
    }
}

/// Apply collected abbreviation labels to the matching `*` atoms in the graph
/// (N-th `*` atom ← N-th label, in atom-index order).
fn assign_abbrevs(mol: &mut MoleculeGraph, labels: &[AbbrevLabel]) {
    let mut label_iter = labels.iter();
    for atom in mol.atoms.iter_mut() {
        if atom.symbol == "*" {
            if let Some(label) = label_iter.next() {
                atom.abbrev = label.text.clone();
                atom.abbrev_style = label.style.clone();
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
        let mut mol = MoleculeGraph::from_smiles(&pre.smiles, pre.forced_direction_markers)?;
        assign_abbrevs(&mut mol, &pre.abbrev_labels);
        let out = compute_layout(&mol)?;
        serde_json::to_vec(&out).map_err(|e| format!("JSON error: {e}"))
    }
}

// ── Native entrypoint for tests / CLI ───────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
pub fn layout_native(smiles: &str) -> Result<LayoutOutput, String> {
    let pre = preprocess_smiles(smiles);
    let mut mol = MoleculeGraph::from_smiles(&pre.smiles, pre.forced_direction_markers)?;
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
    fn cumulated_double_bonds_are_linear() {
        let out = layout_native("O=C=O").expect("carbon dioxide layout failed");
        assert!(atoms_are_collinear(&out, 0, 1, 2));
    }

    #[test]
    fn single_triple_chain_is_linear_at_middle_atom() {
        let out = layout_native("CC#N").expect("nitrile layout failed");
        assert!(atoms_are_collinear(&out, 0, 1, 2));
    }

    #[test]
    fn saturated_two_neighbor_atom_keeps_zigzag() {
        let out = layout_native("CCC").expect("propane layout failed");
        assert!(!atoms_are_collinear(&out, 0, 1, 2));
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
    fn steroid_ring_system_has_regular_bond_lengths() {
        let out = layout_native("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
            .expect("steroid-like molecule layout failed");
        assert!(max_bond_length(&out) < 1.25);
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

    // ── Standards-first stereochemistry and drawing extensions ───────────────

    #[test]
    fn forced_wedge_up_bond() {
        let out = layout_native("C!wN").expect("forced wedge up failed");
        assert_eq!(out.bonds[0].stereo, "wedge_up");
        assert_eq!(out.bonds[0].direction, "none");
    }

    #[test]
    fn forced_wedge_down_bond() {
        let out = layout_native("C!hN").expect("forced wedge down failed");
        assert_eq!(out.bonds[0].stereo, "wedge_down");
        assert_eq!(out.bonds[0].direction, "none");
    }

    #[test]
    fn slash_without_double_bond_is_invalid() {
        let err = layout_native("C/N").expect_err("isolated slash should fail");
        assert!(err.contains("Directional"));
    }

    #[test]
    fn forced_wedge_in_chain() {
        let out = layout_native("CC!wN").expect("forced wedge in chain failed");
        assert_eq!(out.bonds[0].stereo, "none");
        assert_eq!(out.bonds[1].stereo, "wedge_up");
    }

    #[test]
    fn alkene_directional_bonds_do_not_render_as_wedges() {
        let out = layout_native("F/C=C/F").expect("E-alkene failed");
        let double = out.bonds.iter().find(|b| b.order == 2).unwrap();
        assert_eq!(double.order, 2);
        assert!(out.bonds.iter().all(|b| b.stereo == "none"));
        assert_eq!(
            out.bonds.iter().filter(|b| b.direction != "none").count(),
            2
        );
    }

    #[test]
    fn trans_alkene_substituents_are_opposite() {
        let out = layout_native("F/C=C/F").expect("trans alkene failed");
        assert_eq!(alkene_substituent_side_product(&out), -1);
    }

    #[test]
    fn cis_alkene_substituents_are_same_side() {
        let out = layout_native("F/C=C\\F").expect("cis alkene failed");
        assert_eq!(alkene_substituent_side_product(&out), 1);
    }

    #[test]
    fn branch_direction_matches_opensmiles_examples() {
        let trans = layout_native("C(\\F)=C/F").expect("branch trans failed");
        let cis = layout_native("C(/F)=C/F").expect("branch cis failed");
        assert_eq!(alkene_substituent_side_product(&trans), -1);
        assert_eq!(alkene_substituent_side_product(&cis), 1);
    }

    #[test]
    fn conflicting_cis_trans_markers_error() {
        let err = layout_native("C/C(\\F)=C/F").expect_err("conflicting markers should fail");
        assert!(err.contains("Conflicting"));
    }

    #[test]
    fn tetrahedral_chirality_adds_rendered_stereo() {
        let out = layout_native("N[C@@H](C)C(=O)O").expect("chiral alanine failed");
        assert!(out.atoms.iter().any(|a| a.chirality == "tetra_clockwise"));
        // Exactly one wedge for the single stereocenter; direction checked below.
        let wedge_bonds = out
            .bonds
            .iter()
            .filter(|b| b.stereo != "none")
            .count();
        let wedge_h = out.atoms.iter().filter(|a| a.stereo_h != "none").count();
        assert_eq!(wedge_bonds + wedge_h, 1);
        assert!(chirality_matches_smiles("N[C@@H](C)C(=O)O"));
    }

    #[test]
    fn inverting_chirality_flips_the_wedge() {
        // Same skeleton/layout, opposite chirality token ⇒ the wedge on the chosen
        // bond must flip. Guards against regressing to a fixed @→up / @@→down map.
        let r = layout_native("N[C@@H](C)C(=O)O").expect("R failed");
        let s = layout_native("N[C@H](C)C(=O)O").expect("S failed");
        let stereo = |out: &LayoutOutput| -> String {
            out.bonds
                .iter()
                .find(|b| b.stereo != "none")
                .map(|b| b.stereo.clone())
                .or_else(|| {
                    out.atoms
                        .iter()
                        .find(|a| a.stereo_h != "none")
                        .map(|a| a.stereo_h.clone())
                })
                .unwrap_or_else(|| "none".to_string())
        };
        assert_ne!(stereo(&r), "none");
        assert_ne!(stereo(&s), "none");
        assert_ne!(stereo(&r), stereo(&s));
        assert!(chirality_matches_smiles("N[C@@H](C)C(=O)O"));
        assert!(chirality_matches_smiles("N[C@H](C)C(=O)O"));
    }

    #[test]
    fn tetrahedral_stereo_prefers_oh_substituent() {
        let out = layout_native("CC[C@@H](O)CC/C=C/CO").expect("chiral alcohol failed");
        assert!(chiral_oxygen_bond_has_stereo(&out));
    }

    #[test]
    fn chiral_alcohol_keeps_long_chain_as_continuation() {
        let out = layout_native("CC[C@@H](O)CC/C=C/CO").expect("chiral alcohol failed");
        // The long chain leaves the stereocenter (atom 2) as a straight
        // continuation; the short ethyl tail (atom 0) sits on the opposite side.
        // Direction-agnostic so it survives the conventional horizontal mirror.
        let dir = out.atoms[4].pos.x - out.atoms[2].pos.x;
        assert!(dir.abs() > 1e-6);
        assert!(
            (out.atoms[9].pos.x - out.atoms[4].pos.x) * dir > 0.0,
            "terminal alcohol should remain on the zig-zag continuation"
        );
        assert!(
            (out.atoms[0].pos.x - out.atoms[2].pos.x) * dir < 0.0,
            "short ethyl branch should sit opposite the long continuation"
        );
        assert!(out.atoms.iter().all(|atom| atom.stereo_h == "none"));
    }

    #[test]
    fn steroid_stereo_prefers_exocyclic_oh_over_ring_bond() {
        let out =
            layout_native("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
                .expect("steroid-like molecule layout failed");
        assert!(chiral_oxygen_bond_has_stereo(&out));
    }

    #[test]
    fn steroid_ring_chiral_hydrogens_are_rendered() {
        let smiles = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O";
        let out = layout_native(smiles).expect("steroid-like molecule layout failed");

        // The three ring-fusion stereocenters (no exocyclic substituent) render
        // an explicit wedge/hash hydrogen. The 17-OH carbon and the quaternary
        // C13 instead wedge their exocyclic substituent, so they get no H label.
        let stereo_h: Vec<&str> = out
            .atoms
            .iter()
            .filter_map(|atom| (atom.stereo_h != "none").then_some(atom.stereo_h.as_str()))
            .collect();
        assert_eq!(stereo_h.len(), 3);

        // Adjacent ring-fusion stereocenters must point to opposite faces.
        assert_ne!(out.atoms[5].stereo_h, "none");
        assert_ne!(out.atoms[6].stereo_h, "none");
        assert_ne!(out.atoms[5].stereo_h, out.atoms[6].stereo_h);

        // The 17-OH bond is wedged (not the hydrogen).
        assert!(chiral_oxygen_bond_has_stereo(&out));

        // Every stereocenter is depicted with the geometrically correct handedness.
        assert!(chirality_matches_smiles(smiles));
    }

    #[test]
    fn unsupported_chirality_errors() {
        let err = layout_native("C[Fe@SP1](Cl)(Br)I").expect_err("unsupported stereo should fail");
        assert!(err.contains("Only tetrahedral"));
    }

    // ── Abbreviation tests ────────────────────────────────────────────────────

    #[test]
    fn preprocess_single_abbrev() {
        let p = preprocess_smiles("C({PPh3})=O");
        assert_eq!(p.smiles, "C([*])=O");
        assert_eq!(p.abbrev_labels[0].text, "PPh3");
        assert_eq!(p.abbrev_labels[0].style, "");
    }

    #[test]
    fn preprocess_multiple_abbrevs() {
        let p = preprocess_smiles("{OEt}C(=O){NHR}");
        assert_eq!(p.smiles, "[*]C(=O)[*]");
        assert_eq!(p.abbrev_labels[0].text, "OEt");
        assert_eq!(p.abbrev_labels[1].text, "NHR");
    }

    #[test]
    fn preprocess_abbrev_style() {
        let p = preprocess_smiles("{PPh3|P}C({LG|red})=O");
        assert_eq!(p.smiles, "[*]C([*])=O");
        assert_eq!(p.abbrev_labels[0].text, "PPh3");
        assert_eq!(p.abbrev_labels[0].style, "P");
        assert_eq!(p.abbrev_labels[1].text, "LG");
        assert_eq!(p.abbrev_labels[1].style, "red");
    }

    #[test]
    fn preprocess_forced_wedge_markers() {
        let p = preprocess_smiles("C!wN!hO");
        assert_eq!(p.smiles, "C/N\\O");
        assert_eq!(p.forced_direction_markers, vec![true, true]);
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

    #[test]
    fn abbrev_style_assigned_in_layout() {
        let out = layout_native("{PPh3|P}C=O").expect("styled abbrev layout failed");
        assert_eq!(out.atoms[0].abbrev, "PPh3");
        assert_eq!(out.atoms[0].abbrev_style, "P");
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

    // ── Implicit H for expanded valence table ────────────────────────────────

    #[test]
    fn implicit_h_boron() {
        // B (organic subset, valence 3): B bonded to one C → 2 implicit H
        let out = layout_native("BC").expect("boron layout failed");
        assert_eq!(out.atoms[0].implicit_h, 2);
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

    fn atoms_are_collinear(out: &LayoutOutput, a: usize, b: usize, c: usize) -> bool {
        let a = out.atoms[a].pos;
        let b = out.atoms[b].pos;
        let c = out.atoms[c].pos;
        let cross = (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
        cross.abs() < 1e-8
    }

    /// End-to-end check: reconstructs the depicted 3D geometry from the rendered
    /// output and verifies every stereocenter's signed volume matches `@`/`@@`.
    fn chirality_matches_smiles(smiles: &str) -> bool {
        use crate::layout::{implicit_h_count, signed_volume, stereo_h_direction};

        let pre = preprocess_smiles(smiles);
        let mol = MoleculeGraph::from_smiles(&pre.smiles, pre.forced_direction_markers)
            .expect("graph build failed");
        let out = compute_layout(&mol).expect("layout failed");
        let coords: Vec<crate::render::Vec2> = out.atoms.iter().map(|a| a.pos).collect();

        for center in 0..mol.atoms.len() {
            let parity = match out.atoms[center].chirality.as_str() {
                "tetra_anti" => -1.0_f64,
                "tetra_clockwise" => 1.0_f64,
                _ => continue,
            };
            let neighbor_bonds = &mol.neighbor_bonds[center];
            let n_h = (mol.atoms[center].hcount + implicit_h_count(&mol, center)) as usize;
            if neighbor_bonds.len() + n_h != 4 || n_h > 1 {
                continue;
            }

            // Neighbor order with the implicit hydrogen inserted (same rule as the
            // renderer): after the "from" atom, or first if there is none.
            #[derive(Clone, Copy)]
            enum N {
                Bond(usize, usize),
                H,
            }
            let mut order: Vec<N> = neighbor_bonds
                .iter()
                .map(|&b| {
                    let bond = &mol.bonds[b];
                    let other = if bond.from == center { bond.to } else { bond.from };
                    N::Bond(b, other)
                })
                .collect();
            if n_h == 1 {
                let pos = if mol.has_preceding[center] { 1 } else { 0 };
                order.insert(pos.min(order.len()), N::H);
            }

            // Which bond (if any) carries the rendered wedge, and its z sign.
            let wedge_bond = neighbor_bonds
                .iter()
                .find(|&&b| out.bonds[b].stereo != "none")
                .copied();
            let h_dir = stereo_h_direction(&mol, center, &coords, None);

            let unit = |dx: f64, dy: f64| {
                let l = (dx * dx + dy * dy).sqrt();
                if l > 1e-12 {
                    (dx / l, dy / l)
                } else {
                    (dx, dy)
                }
            };
            let z_of = |s: &str| if s == "wedge_up" { 1.0 } else { -1.0 };

            let mut dirs = [[0.0_f64; 3]; 4];
            for (i, slot) in order.iter().enumerate() {
                dirs[i] = match slot {
                    N::Bond(b, other) => {
                        let (ux, uy) =
                            unit(coords[*other].x - coords[center].x, coords[*other].y - coords[center].y);
                        let z = if Some(*b) == wedge_bond {
                            z_of(&out.bonds[*b].stereo)
                        } else {
                            0.0
                        };
                        [ux, uy, z]
                    }
                    N::H => {
                        // If the H itself is wedged, use that; otherwise it sits on
                        // the face opposite the wedged heavy substituent.
                        let z = if out.atoms[center].stereo_h != "none" {
                            z_of(&out.atoms[center].stereo_h)
                        } else if let Some(b) = wedge_bond {
                            -z_of(&out.bonds[b].stereo)
                        } else {
                            0.0
                        };
                        let dir = if out.atoms[center].stereo_h != "none" {
                            out.atoms[center].stereo_h_dir
                        } else {
                            h_dir
                        };
                        [dir.x, dir.y, z]
                    }
                };
            }

            let vol = signed_volume(&dirs);
            if vol.abs() < 1e-9 || vol.signum() != parity {
                return false;
            }
        }
        true
    }

    fn chiral_oxygen_bond_has_stereo(out: &LayoutOutput) -> bool {
        out.bonds.iter().any(|bond| {
            let from = &out.atoms[bond.from];
            let to = &out.atoms[bond.to];
            bond.stereo != "none"
                && ((from.chirality != "none" && to.symbol == "O")
                    || (to.chirality != "none" && from.symbol == "O"))
        })
    }

    fn alkene_substituent_side_product(out: &LayoutOutput) -> i8 {
        let double = out.bonds.iter().find(|b| b.order == 2).unwrap();
        let a = double.from;
        let b = double.to;
        let left = out
            .bonds
            .iter()
            .find(|bond| bond.order == 1 && (bond.from == a || bond.to == a))
            .unwrap();
        let right = out
            .bonds
            .iter()
            .find(|bond| bond.order == 1 && (bond.from == b || bond.to == b))
            .unwrap();
        let left_neighbor = if left.from == a { left.to } else { left.from };
        let right_neighbor = if right.from == b {
            right.to
        } else {
            right.from
        };
        side(
            out.atoms[a].pos,
            out.atoms[b].pos,
            out.atoms[left_neighbor].pos,
        ) * side(
            out.atoms[a].pos,
            out.atoms[b].pos,
            out.atoms[right_neighbor].pos,
        )
    }

    fn side(a: crate::render::Vec2, b: crate::render::Vec2, p: crate::render::Vec2) -> i8 {
        let cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
        if cross > 1e-8 {
            1
        } else {
            -1
        }
    }
}
