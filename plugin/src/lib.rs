mod error;
mod graph;
mod kekulize;
mod layout;
mod render;

pub use render::LayoutOutput;

use graph::{ForcedBondKind, MoleculeGraph};
use layout::compute_layout;
use ptable::Element;

// ── SMILES preprocessing ─────────────────────────────────────────────────────
//
// Rewrites the input into the subset of SMILES that smiles-parser accepts and
// collects side-channel data the parser cannot carry:
//
//   {label}  →  [*]   Abbreviated group (e.g. {PPh3}, {OEt}).
//                     The Nth [*] atom in the output gets abbrev = the Nth label.
//                     A `>` marker inside the label selects the attachment glyph.
//   c, n, …  →  C, N  Unbracketed aromatic atoms are uppercased (the parser
//                     only accepts aliphatic organic-subset symbols); a marker
//                     records which atoms were aromatic so the graph stage can
//                     kekulize. Bracket aromatics ([nH], [se]) parse natively.
//
// The cleaned string is valid parser input; extensions are stripped out.

pub(crate) struct Preprocessed {
    pub smiles: String,
    /// Labels in the order they appeared, one per {label} token.
    pub abbrev_labels: Vec<AbbrevLabel>,
    /// One entry for each `/` or `\` token in `smiles`, recording whether it
    /// is plain SMILES or which `!w`/`!h`/`!s`/`!d` drawing extension it
    /// carries.
    pub forced_direction_markers: Vec<ForcedBondKind>,
    /// One flag per unbracketed organic-subset atom token in `smiles`, in
    /// writing order. `true` means the atom was written lowercase (aromatic).
    pub aromatic_atom_markers: Vec<bool>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct AbbrevLabel {
    pub text: String,
    pub style: String,
    pub anchor: usize,
    pub anchor_len: usize,
}

pub(crate) fn preprocess_smiles(input: &str) -> Result<Preprocessed, String> {
    let mut smiles = String::with_capacity(input.len());
    let mut abbrev_labels = Vec::new();
    let mut forced_direction_markers = Vec::new();
    let mut aromatic_atom_markers = Vec::new();
    let mut in_bracket = false;

    let mut chars = input.chars().peekable();
    while let Some(ch) = chars.next() {
        // Bracket atoms pass through untouched; their contents (element
        // symbols, hcount, charge) must not be mistaken for atom tokens.
        if in_bracket {
            if ch == ']' {
                in_bracket = false;
            }
            smiles.push(ch);
            continue;
        }
        if ch == '{' {
            let mut raw_label = String::new();
            for inner in chars.by_ref() {
                if inner == '}' {
                    break;
                }
                raw_label.push(inner);
            }
            abbrev_labels.push(parse_abbrev_label(&raw_label)?);
            smiles.push_str("[*]");
        } else if ch == '!' {
            match chars.peek().copied() {
                Some('w') => {
                    chars.next();
                    smiles.push('/');
                    forced_direction_markers.push(ForcedBondKind::Wedge);
                }
                Some('h') => {
                    chars.next();
                    smiles.push('\\');
                    forced_direction_markers.push(ForcedBondKind::Wedge);
                }
                Some('s') => {
                    chars.next();
                    smiles.push('/');
                    forced_direction_markers.push(ForcedBondKind::Wavy);
                }
                Some('d') => {
                    chars.next();
                    smiles.push('/');
                    forced_direction_markers.push(ForcedBondKind::Dashed);
                }
                _ => {
                    smiles.push(ch);
                }
            }
        } else if ch == '>' {
            return Err(
                "`>` is only valid inside an abbreviation label like `{>PPh3}`".to_string(),
            );
        } else {
            match ch {
                '[' => in_bracket = true,
                '/' | '\\' => forced_direction_markers.push(ForcedBondKind::Plain),
                // Aromatic organic-subset atoms: uppercase for the parser,
                // remember the aromatic designation.
                'b' | 'c' | 'n' | 'o' | 'p' | 's' => {
                    smiles.push(ch.to_ascii_uppercase());
                    aromatic_atom_markers.push(true);
                    continue;
                }
                // Aliphatic organic-subset atom starts ('l' of Cl and 'r' of
                // Br fall through to the default arm, so each two-letter
                // symbol counts once).
                'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'I' => {
                    aromatic_atom_markers.push(false);
                }
                _ => {}
            }
            smiles.push(ch);
        }
    }

    Ok(Preprocessed {
        smiles,
        abbrev_labels,
        forced_direction_markers,
        aromatic_atom_markers,
    })
}

fn parse_abbrev_label(raw_label: &str) -> Result<AbbrevLabel, String> {
    let (raw_text, style) = raw_label
        .split_once('|')
        .map(|(text, style)| (text, style.trim().to_string()))
        .unwrap_or((raw_label, String::new()));
    let marker_count = raw_text.chars().filter(|&ch| ch == '>').count();
    if marker_count > 1 {
        return Err(
            "abbreviation labels may contain at most one `>` attachment marker".to_string(),
        );
    }

    let mut text = String::with_capacity(raw_text.len());
    let mut marker = None;
    for ch in raw_text.chars() {
        if ch == '>' {
            marker = Some(text.chars().count());
        } else {
            text.push(ch);
        }
    }

    let chars: Vec<char> = text.chars().collect();
    let (anchor, anchor_len) = if let Some(pos) = marker {
        if chars.is_empty() {
            return Err("abbreviation attachment marker `>` needs a label glyph".to_string());
        }
        if pos >= chars.len() {
            return Err(
                "abbreviation attachment marker `>` must precede a label glyph".to_string(),
            );
        }
        let anchor = pos;
        let len = if chars.get(anchor).is_some_and(|ch| ch.is_ascii_uppercase())
            && chars
                .get(anchor + 1)
                .is_some_and(|ch| ch.is_ascii_lowercase())
        {
            2
        } else {
            1
        };
        (anchor, len)
    } else {
        (0, 0)
    };

    Ok(AbbrevLabel {
        text,
        style,
        anchor,
        anchor_len,
    })
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
                atom.abbrev_anchor = label.anchor;
                atom.abbrev_anchor_len = label.anchor_len;
            }
        }
    }
}

// ── Molecular weight ─────────────────────────────────────────────────────────

/// Sums standard atomic weights (PubChem / IUPAC values via ptable) over all
/// atoms, including explicit and implicit hydrogens. Errors on anything whose
/// mass is undefined: wildcard atoms, `{label}` abbreviations, and isotope
/// labels (standard atomic weights do not apply to specific nuclides).
fn compute_mol_weight(mol: &MoleculeGraph) -> Result<f64, String> {
    let h_mass = Element::Hydrogen.get_atomic_mass() as f64;
    let mut total = 0.0f64;
    for (i, atom) in mol.atoms.iter().enumerate() {
        if !atom.abbrev.is_empty() {
            return Err(format!(
                "cannot compute molecular weight: abbreviation {{{}}} has no defined composition",
                atom.abbrev
            ));
        }
        if atom.symbol == "*" {
            return Err(
                "cannot compute molecular weight: wildcard atom `*` has no mass".to_string(),
            );
        }
        if let Some(isotope) = atom.isotope {
            return Err(format!(
                "cannot compute molecular weight: isotope [{isotope}{}] needs a nuclide mass, \
                 not a standard atomic weight",
                atom.symbol
            ));
        }
        let element = Element::from_symbol(&atom.symbol).ok_or_else(|| {
            format!("cannot compute molecular weight: unknown element {}", atom.symbol)
        })?;
        total += element.get_atomic_mass() as f64;
        let hydrogens = atom.hcount as f64 + layout::implicit_h_count(mol, i) as f64;
        total += hydrogens * h_mass;
    }
    Ok(total)
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
        let pre = preprocess_smiles(s)?;
        let mut mol = MoleculeGraph::from_smiles(
            &pre.smiles,
            pre.forced_direction_markers,
            pre.aromatic_atom_markers,
        )?;
        assign_abbrevs(&mut mol, &pre.abbrev_labels);
        let out = compute_layout(&mol)?;
        serde_json::to_vec(&out).map_err(|e| format!("JSON error: {e}"))
    }

    /// Called from Typst as `smiles-plugin.mol_weight(bytes(smiles-str))`.
    /// Returns the molecular weight in g/mol as a JSON number.
    #[wasm_func]
    pub fn mol_weight(smiles: &[u8]) -> Result<Vec<u8>, String> {
        let s = core::str::from_utf8(smiles).map_err(|e| format!("UTF-8 error: {e}"))?;
        let pre = preprocess_smiles(s)?;
        let mut mol = MoleculeGraph::from_smiles(
            &pre.smiles,
            pre.forced_direction_markers,
            pre.aromatic_atom_markers,
        )?;
        assign_abbrevs(&mut mol, &pre.abbrev_labels);
        let weight = compute_mol_weight(&mol)?;
        serde_json::to_vec(&weight).map_err(|e| format!("JSON error: {e}"))
    }
}

// ── Native entrypoint for tests / CLI ───────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
pub fn layout_native(smiles: &str) -> Result<LayoutOutput, String> {
    let pre = preprocess_smiles(smiles)?;
    let mut mol = MoleculeGraph::from_smiles(
        &pre.smiles,
        pre.forced_direction_markers,
        pre.aromatic_atom_markers,
    )?;
    assign_abbrevs(&mut mol, &pre.abbrev_labels);
    compute_layout(&mol)
}

#[cfg(not(target_arch = "wasm32"))]
pub fn mol_weight_native(smiles: &str) -> Result<f64, String> {
    let pre = preprocess_smiles(smiles)?;
    let mut mol = MoleculeGraph::from_smiles(
        &pre.smiles,
        pre.forced_direction_markers,
        pre.aromatic_atom_markers,
    )?;
    assign_abbrevs(&mut mol, &pre.abbrev_labels);
    compute_mol_weight(&mol)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Smallest distance between any two distinct non-virtual atoms.
    fn min_atom_distance(out: &LayoutOutput) -> f64 {
        let mut min = f64::INFINITY;
        for i in 0..out.atoms.len() {
            for j in (i + 1)..out.atoms.len() {
                if out.atoms[i].virtual_h || out.atoms[j].virtual_h {
                    continue;
                }
                min = min.min(out.atoms[i].pos.dist(out.atoms[j].pos));
            }
        }
        min
    }

    #[test]
    fn ortho_ring_substituents_do_not_collide() {
        // Aspirin: the acetyl C=O and the carboxyl OH grow from ortho ring
        // positions toward each other and must not land on the same point.
        let out = layout_native("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap();
        assert!(
            min_atom_distance(&out) > 0.5,
            "atoms overlap: min distance {}",
            min_atom_distance(&out)
        );

        // The carboxyl carbon (atom 10: ring C9, =O 11, OH 12) must keep its
        // textbook trigonal geometry: the conflict is resolved by flipping
        // the acetyl branch aside, not by squeezing the carboxyl bonds into a
        // narrow fan or a straight line.
        let c = out.atoms[10].pos;
        let angles: Vec<f64> = [9, 11, 12]
            .iter()
            .map(|&n| {
                let p = out.atoms[n].pos;
                (p.y - c.y).atan2(p.x - c.x)
            })
            .collect();
        for i in 0..angles.len() {
            for j in (i + 1)..angles.len() {
                let mut delta = (angles[i] - angles[j]).abs();
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
        let out = layout_native("CCCC1=CC=CC=C1CCC").unwrap();
        assert!(min_atom_distance(&out) > 0.5);
    }

    #[test]
    fn benzene_kekule() {
        let out = layout_native("C1=CC=CC=C1").expect("benzene layout failed");
        assert_eq!(out.atoms.len(), 6);
        assert_eq!(out.bonds.len(), 6);
    }

    #[test]
    fn ethanol() {
        // CCO: 3 real atoms + 1 virtual H for terminal O (implicit OH)
        let out = layout_native("CCO").expect("ethanol layout failed");
        assert_eq!(out.atoms.len(), 4);
        assert!(out.atoms[3].virtual_h);
    }

    #[test]
    fn implicit_h_counts_ethanol() {
        // implicit_h counts for real atoms only (virtual H has implicit_h = 0)
        let out = layout_native("CCO").expect("ethanol layout failed");
        let counts: Vec<u8> = out.atoms.iter().map(|a| a.implicit_h).collect();
        assert_eq!(counts, vec![3, 2, 1, 0]);
    }

    #[test]
    fn explicit_h_suppresses_implicit_h() {
        let out = layout_native("[NH4+]").expect("ammonium layout failed");
        assert_eq!(out.atoms[0].hcount, 4);
        assert_eq!(out.atoms[0].implicit_h, 0);
    }

    #[test]
    fn bracket_hydrogen_atoms_fold_into_neighbor() {
        // Methane written with four explicit [H] atoms collapses into one
        // carbon carrying the folded hydrogen count; the only remaining H is the
        // virtual label placeholder, never a drawn atom.
        let out = layout_native("C([H])([H])([H])[H]").expect("methane layout failed");
        assert_eq!(out.atoms[0].symbol, "C");
        assert_eq!(out.atoms[0].hcount, 4);
        assert_eq!(out.atoms[0].implicit_h, 0);
        let drawn_h = out
            .atoms
            .iter()
            .filter(|a| a.symbol == "H" && !a.virtual_h)
            .count();
        assert_eq!(drawn_h, 0);
    }

    #[test]
    fn folded_hydrogens_drop_from_neighbor_ordering() {
        // The dichloromethyl carbon keeps only its three heavy neighbors after
        // the explicit hydrogen is folded away.
        let out = layout_native("ClC([H])Cl").expect("dichloromethane layout failed");
        let heavy = out.atoms.iter().filter(|a| a.symbol != "H").count();
        assert_eq!(heavy, 3);
        assert_eq!(out.bonds.iter().filter(|b| !b.virtual_bond).count(), 2);
    }

    #[test]
    fn isotopic_and_charged_hydrogens_are_kept() {
        // Deuterium is a real, drawn atom; only the plain [H] atoms fold away.
        let out = layout_native("[2H]C([H])([H])[H]").expect("deuteromethane layout failed");
        let drawn_h = out
            .atoms
            .iter()
            .filter(|a| a.symbol == "H" && !a.virtual_h)
            .count();
        assert_eq!(drawn_h, 1);
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
    fn symmetric_hub_is_mirror_symmetric() {
        // Tetraethylmethane's four equal arms must be drawn as a mirror-symmetric
        // figure (the "tweezers" depiction), not a rotational pinwheel. The layout
        // is built about the vertical axis, so reflecting every atom across that
        // axis must reproduce the same set of positions.
        let out = layout_native("CCC(CC)(CC)CC").expect("tetraethylmethane layout failed");
        let cx = out.atoms.iter().map(|a| a.pos.x).sum::<f64>() / out.atoms.len() as f64;
        for atom in &out.atoms {
            let mirrored_x = 2.0 * cx - atom.pos.x;
            let found = out
                .atoms
                .iter()
                .any(|b| (b.pos.x - mirrored_x).abs() < 1e-6 && (b.pos.y - atom.pos.y).abs() < 1e-6);
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
        let out =
            layout_native("BrCC(=O)OC(C)(O)OC(=O)C").expect("acetal diester layout failed");
        assert!(
            min_nonbonded_distance(&out) > 0.5,
            "branches overlap: min non-bonded distance {}",
            min_nonbonded_distance(&out)
        );
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

    // ── Dot-disconnected structures ──────────────────────────────────────────

    #[test]
    fn dot_creates_no_bond() {
        let out = layout_native("CCO.CCO").expect("two ethanols failed");
        assert_eq!(out.atoms.iter().filter(|a| !a.virtual_h).count(), 6);
        // Two C-C-O fragments: 4 real bonds, none crossing the dot.
        let real: Vec<_> = out.bonds.iter().filter(|b| !b.virtual_bond).collect();
        assert_eq!(real.len(), 4);
        assert!(real.iter().all(|b| (b.from < 3) == (b.to < 3)));
    }

    #[test]
    fn salt_fragments_are_disconnected_and_ordered() {
        // Sodium acetate: no bond may touch the sodium ion, and fragments keep
        // SMILES writing order left to right with a visible gap.
        let out = layout_native("CC(=O)[O-].[Na+]").expect("sodium acetate failed");
        let na = 4;
        assert_eq!(out.atoms[na].symbol, "Na");
        assert_eq!(out.atoms[na].charge, 1);
        assert!(out.bonds.iter().all(|b| b.from != na && b.to != na));

        let max_acetate_x = (0..4).map(|i| out.atoms[i].pos.x).fold(f64::MIN, f64::max);
        assert!(
            out.atoms[na].pos.x >= max_acetate_x + 1.0,
            "Na+ should sit clearly right of the acetate fragment"
        );
    }

    #[test]
    fn bare_ion_pair() {
        let out = layout_native("[Na+].[Cl-]").expect("sodium chloride failed");
        assert_eq!(out.atoms.len(), 2);
        assert!(out.bonds.is_empty());
        assert!(out.atoms[1].pos.x > out.atoms[0].pos.x);
    }

    #[test]
    fn three_fragments() {
        let out = layout_native("O.O.O").expect("three waters failed");
        let heavy: Vec<_> = out.atoms.iter().filter(|a| !a.virtual_h).collect();
        assert_eq!(heavy.len(), 3);
        assert!(out.bonds.iter().all(|b| b.virtual_bond));
        assert!(heavy[0].pos.x < heavy[1].pos.x && heavy[1].pos.x < heavy[2].pos.x);
    }

    #[test]
    fn ring_closure_across_dot_joins_fragments() {
        // OpenSMILES: "C1.C1" is ethane — the ring-closure bond still forms
        // even though a dot separates the digits.
        let out = layout_native("C1.C1").expect("dot ring closure failed");
        assert_eq!(out.atoms.len(), 2);
        assert_eq!(out.bonds.len(), 1);
        assert_eq!(out.atoms[0].implicit_h, 3);
    }

    #[test]
    fn aromatic_fragments_kekulize_independently() {
        let out = layout_native("c1ccccc1.c1ccccc1").expect("two benzenes failed");
        assert_eq!(out.atoms.len(), 12);
        assert_eq!(order_counts(&out), (6, 6));
        assert!(out.bonds.iter().all(|b| (b.from < 6) == (b.to < 6)));
    }

    // ── Aromatic (lowercase) SMILES and kekulization ─────────────────────────

    fn order_counts(out: &LayoutOutput) -> (usize, usize) {
        let real = out.bonds.iter().filter(|b| !b.virtual_bond);
        (
            real.clone().filter(|b| b.order == 1).count(),
            real.filter(|b| b.order == 2).count(),
        )
    }

    #[test]
    fn benzene_aromatic() {
        let out = layout_native("c1ccccc1").expect("aromatic benzene failed");
        assert_eq!(out.atoms.len(), 6);
        assert_eq!(order_counts(&out), (3, 3));
        assert!(out.atoms.iter().all(|a| a.implicit_h == 1));
    }

    #[test]
    fn aromatic_atom_indices_match_writing_order() {
        // Kekulization must not reorder atoms: index N is the Nth atom token,
        // so show-indices / highlight / arrow references keep working.
        let out = layout_native("Cc1ccncc1").expect("4-methylpyridine failed");
        let symbols: Vec<&str> = out
            .atoms
            .iter()
            .filter(|a| !a.virtual_h)
            .map(|a| a.symbol.as_str())
            .collect();
        assert_eq!(symbols, vec!["C", "C", "C", "C", "N", "C", "C"]);
    }

    #[test]
    fn pyridine_aromatic() {
        let out = layout_native("c1ccncc1").expect("pyridine failed");
        assert_eq!(out.atoms[3].symbol, "N");
        assert_eq!(out.atoms[3].implicit_h, 0);
        assert_eq!(order_counts(&out), (3, 3));
    }

    #[test]
    fn pyrrole_aromatic() {
        let out = layout_native("c1cc[nH]c1").expect("pyrrole failed");
        let n = out.atoms.iter().find(|a| a.symbol == "N").unwrap();
        assert_eq!(n.hcount, 1);
        assert_eq!(order_counts(&out), (3, 2));
    }

    #[test]
    fn furan_and_thiophene_aromatic() {
        for smiles in ["c1occc1", "c1sccc1"] {
            let out = layout_native(smiles).expect("5-ring heteroaromatic failed");
            assert_eq!(order_counts(&out), (3, 2), "wrong kekulization for {smiles}");
            let hetero = &out.atoms[1];
            assert_eq!(hetero.implicit_h, 0);
        }
    }

    #[test]
    fn imidazole_aromatic() {
        let out = layout_native("c1cnc[nH]1").expect("imidazole failed");
        assert_eq!(order_counts(&out), (3, 2));
    }

    #[test]
    fn n_methylpyrrole_aromatic() {
        // A three-connected aromatic n carries no H and no double bond.
        let out = layout_native("Cn1cccc1").expect("N-methylpyrrole failed");
        assert_eq!(out.atoms[1].symbol, "N");
        assert_eq!(out.atoms[1].implicit_h, 0);
        assert_eq!(order_counts(&out), (4, 2));
    }

    #[test]
    fn naphthalene_aromatic() {
        let out = layout_native("c1ccc2ccccc2c1").expect("aromatic naphthalene failed");
        assert_eq!(out.atoms.len(), 10);
        assert_eq!(order_counts(&out), (6, 5));
    }

    #[test]
    fn indane_mixed_aromatic_aliphatic() {
        // Spec example: aromatic ring fused to an aliphatic ring.
        let out = layout_native("c1ccc2CCCc2c1").expect("indane failed");
        assert_eq!(out.atoms.len(), 9);
        assert_eq!(order_counts(&out), (7, 3));
    }

    #[test]
    fn biphenyl_explicit_and_implicit_single_link() {
        for smiles in ["c1ccccc1-c1ccccc1", "c1ccccc1c1ccccc1"] {
            let out = layout_native(smiles).expect("biphenyl failed");
            assert_eq!(out.atoms.len(), 12);
            assert_eq!(order_counts(&out), (7, 6), "wrong kekulization for {smiles}");
            // The inter-ring bond stays single.
            let link = out
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
        let out = layout_native("O=c1cccc[nH]1").expect("2-pyridinone failed");
        assert_eq!(order_counts(&out), (4, 3));
    }

    #[test]
    fn explicit_aromatic_bond_symbol() {
        let out = layout_native("c1:c:c:c:c:c1").expect("explicit ':' benzene failed");
        assert_eq!(order_counts(&out), (3, 3));
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
        let out = layout_native("c1cc*cc1").expect("aromatic ring with wildcard failed");
        let (_, doubles) = order_counts(&out);
        assert!(doubles >= 2, "expected an alternating pattern, got {doubles} double bonds");
    }

    #[test]
    fn azulene_aromatic() {
        // Non-alternant 5-7 fused system: every carbon needs a double bond.
        let out = layout_native("c1ccc2cccc2cc1").expect("azulene failed");
        assert_eq!(out.atoms.len(), 10);
        assert_eq!(order_counts(&out), (6, 5));
    }

    #[test]
    fn caffeine_aromatic() {
        let out =
            layout_native("Cn1cnc2c1c(=O)n(C)c(=O)n2C").expect("caffeine failed");
        // Purine core: the imidazole C=N plus the C4=C5 bridge double bond,
        // and the two exocyclic carbonyls.
        assert_eq!(out.atoms.iter().filter(|a| !a.virtual_h).count(), 14);
        let (_, doubles) = order_counts(&out);
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
        assert!(err.contains("kekulize"));
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
    fn forced_wavy_bond() {
        let out = layout_native("C!sN").expect("forced wavy bond failed");
        assert_eq!(out.bonds[0].stereo, "wavy");
        assert_eq!(out.bonds[0].direction, "none");
    }

    #[test]
    fn forced_dashed_bond() {
        let out = layout_native("C!dN").expect("forced dashed bond failed");
        assert_eq!(out.bonds[0].stereo, "dashed");
        assert_eq!(out.bonds[0].direction, "none");
    }

    #[test]
    fn forced_wavy_and_dashed_in_chain() {
        let out = layout_native("CC!sO!dN").expect("wavy/dashed chain failed");
        assert_eq!(out.bonds[0].stereo, "none");
        assert_eq!(out.bonds[1].stereo, "wavy");
        assert_eq!(out.bonds[2].stereo, "dashed");
    }

    #[test]
    fn forced_wavy_does_not_disturb_real_directional_bonds() {
        // A genuine trans alkene next to a forced wavy bond: the wavy marker
        // must not consume or shift the cis/trans direction tokens.
        let out = layout_native("F/C=C/C!sN").expect("mixed directional/wavy failed");
        let wavy = out.bonds.iter().filter(|b| b.stereo == "wavy").count();
        assert_eq!(wavy, 1);
        let directional = out.bonds.iter().filter(|b| b.direction != "none").count();
        assert_eq!(directional, 2);
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

    // ── Extended chirality classes and quadruple bonds ───────────────────────

    /// Cosine of the angle neighbor-a → center → neighbor-b.
    fn bond_angle_cos(out: &LayoutOutput, center: usize, a: usize, b: usize) -> f64 {
        let c = out.atoms[center].pos;
        let (ax, ay) = (out.atoms[a].pos.x - c.x, out.atoms[a].pos.y - c.y);
        let (bx, by) = (out.atoms[b].pos.x - c.x, out.atoms[b].pos.y - c.y);
        (ax * bx + ay * by) / ((ax * ax + ay * ay).sqrt() * (bx * bx + by * by).sqrt())
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
        let z = layout_native("N[Pt@SP3](N)(Cl)Cl").expect("SP3 failed");
        assert!(bond_angle_cos(&z, 1, 3, 4).abs() < 1e-6);
        assert!((bond_angle_cos(&z, 1, 0, 4) + 1.0).abs() < 1e-6);
    }

    #[test]
    fn square_planar_all_neighbors_at_right_angles() {
        let out = layout_native("C[Fe@SP1](Cl)(Br)I").expect("SP iron failed");
        assert_eq!(out.atoms[1].chirality, "square_planar");
        for (a, b) in [(0, 2), (2, 3), (3, 4)] {
            assert!(bond_angle_cos(&out, 1, a, b).abs() < 1e-6);
        }
    }

    #[test]
    fn tb_oh_al_accepted_without_stereo_marks() {
        for smiles in [
            "S[As@TB1](F)(Cl)(Br)N",
            "C[Co@OH1](F)(Cl)(Br)(I)N",
            "NC(Br)=[C@AL1]=C(O)C",
        ] {
            let out = layout_native(smiles).expect("extended chirality should parse");
            assert!(
                out.atoms.iter().any(|a| a.chirality == "undepicted"),
                "missing undepicted center for {smiles}"
            );
            assert!(out.bonds.iter().all(|b| b.stereo == "none"));
            assert!(out.atoms.iter().all(|a| a.stereo_h == "none"));
        }
    }

    #[test]
    fn quadruple_bond_order() {
        // The classic metal-metal quadruple bond, e.g. in [Re2Cl8]2-.
        let out = layout_native("[Re]$[Re]").expect("quadruple bond failed");
        assert_eq!(out.bonds[0].order, 4);
    }

    #[test]
    fn quadruple_bond_counts_toward_valence() {
        let out = layout_native("C$C").expect("C$C failed");
        assert_eq!(out.bonds[0].order, 4);
        assert!(out.atoms.iter().all(|a| a.implicit_h == 0));
    }

    // ── Abbreviation tests ────────────────────────────────────────────────────

    #[test]
    fn preprocess_single_abbrev() {
        let p = preprocess_smiles("C({PPh3})=O").expect("preprocess failed");
        assert_eq!(p.smiles, "C([*])=O");
        assert_eq!(p.abbrev_labels[0].text, "PPh3");
        assert_eq!(p.abbrev_labels[0].style, "");
        assert_eq!(p.abbrev_labels[0].anchor_len, 0);
    }

    #[test]
    fn preprocess_multiple_abbrevs() {
        let p = preprocess_smiles("{OEt}C(=O){NHR}").expect("preprocess failed");
        assert_eq!(p.smiles, "[*]C(=O)[*]");
        assert_eq!(p.abbrev_labels[0].text, "OEt");
        assert_eq!(p.abbrev_labels[1].text, "NHR");
    }

    #[test]
    fn preprocess_abbrev_style() {
        let p = preprocess_smiles("{PPh3|P}C({LG|red})=O").expect("preprocess failed");
        assert_eq!(p.smiles, "[*]C([*])=O");
        assert_eq!(p.abbrev_labels[0].text, "PPh3");
        assert_eq!(p.abbrev_labels[0].style, "P");
        assert_eq!(p.abbrev_labels[1].text, "LG");
        assert_eq!(p.abbrev_labels[1].style, "red");
    }

    #[test]
    fn preprocess_abbrev_anchor_positions() {
        let p = preprocess_smiles("{>CAT}C({C>AT}){CA>T}").expect("preprocess failed");
        assert_eq!(p.smiles, "[*]C([*])[*]");
        assert_eq!(p.abbrev_labels[0].text, "CAT");
        assert_eq!(p.abbrev_labels[0].anchor, 0);
        assert_eq!(p.abbrev_labels[0].anchor_len, 1);
        assert_eq!(p.abbrev_labels[1].text, "CAT");
        assert_eq!(p.abbrev_labels[1].anchor, 1);
        assert_eq!(p.abbrev_labels[1].anchor_len, 1);
        assert_eq!(p.abbrev_labels[2].text, "CAT");
        assert_eq!(p.abbrev_labels[2].anchor, 2);
        assert_eq!(p.abbrev_labels[2].anchor_len, 1);
    }

    #[test]
    fn preprocess_rejects_arrow_marker_outside_abbrev() {
        assert!(preprocess_smiles("C>C").is_err());
        assert!(preprocess_smiles("{CAT>}C").is_err());
    }

    #[test]
    fn preprocess_forced_wedge_markers() {
        let p = preprocess_smiles("C!wN!hO").expect("preprocess failed");
        assert_eq!(p.smiles, "C/N\\O");
        assert_eq!(
            p.forced_direction_markers,
            vec![ForcedBondKind::Wedge, ForcedBondKind::Wedge]
        );
    }

    #[test]
    fn preprocess_no_abbrev() {
        let p = preprocess_smiles("CCO").expect("preprocess failed");
        assert_eq!(p.smiles, "CCO");
        assert!(p.abbrev_labels.is_empty());
        assert_eq!(p.aromatic_atom_markers, vec![false, false, false]);
    }

    #[test]
    fn preprocess_uppercases_aromatic_atoms() {
        let p = preprocess_smiles("Clc1ccccc1").expect("preprocess failed");
        assert_eq!(p.smiles, "ClC1CCCCC1");
        // One marker per unbracketed organic-subset atom, in writing order;
        // the two-letter Cl counts once.
        assert_eq!(
            p.aromatic_atom_markers,
            vec![false, true, true, true, true, true, true]
        );
    }

    #[test]
    fn preprocess_leaves_bracket_atoms_alone() {
        // Bracket contents are the parser's business: [nH] parses natively as
        // an aromatic atom, and the 'c' in [Sc] is not an aromatic carbon.
        let p = preprocess_smiles("c1cc[nH]c1[Sc]").expect("preprocess failed");
        assert_eq!(p.smiles, "C1CC[nH]C1[Sc]");
        assert_eq!(p.aromatic_atom_markers, vec![true, true, true, true]);
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
    fn abbrev_anchor_assigned_in_layout() {
        let out = layout_native("{>PPh3}C=O").expect("anchored abbrev layout failed");
        assert_eq!(out.atoms[0].abbrev, "PPh3");
        assert_eq!(out.atoms[0].abbrev_anchor, 0);
        assert_eq!(out.atoms[0].abbrev_anchor_len, 1);
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
        let out = layout_native("[OH-]").expect("[OH-] layout failed");
        assert_eq!(out.atoms.len(), 2);
        assert_eq!(out.atoms[0].symbol, "O");
        assert!(out.atoms[1].virtual_h);
        assert_eq!(out.atoms[1].symbol, "H");
        assert_eq!(out.bonds.len(), 1);
        assert!(out.bonds[0].virtual_bond);
    }

    #[test]
    fn ammonium_h_group_is_one_index() {
        // [NH4+]: despite hcount=4 the H-label is one glyph → exactly 1 virtual H.
        let out = layout_native("[NH4+]").expect("[NH4+] layout failed");
        assert_eq!(out.atoms.len(), 2); // 1 N + 1 virtual H group
        assert!(out.atoms[1].virtual_h);
        assert_eq!(out.bonds.len(), 1);
        assert!(out.bonds[0].virtual_bond);
    }

    #[test]
    fn virtual_h_has_valid_position() {
        let out = layout_native("[OH-]").expect("[OH-] layout failed");
        let o = out.atoms[0].pos;
        let h = out.atoms[1].pos;
        let dist = ((h.x - o.x).powi(2) + (h.y - o.y).powi(2)).sqrt();
        assert!((dist - 0.35).abs() < 1e-6, "H should be 0.35 bond lengths from O, got {dist}");
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

    /// Smallest distance between any two atoms that are not bonded to each other.
    /// A value well below the ~1.0 bond length signals two parts of the molecule
    /// being laid out on top of one another.
    fn min_nonbonded_distance(out: &LayoutOutput) -> f64 {
        let bonded: std::collections::HashSet<(usize, usize)> = out
            .bonds
            .iter()
            .map(|b| (b.from.min(b.to), b.from.max(b.to)))
            .collect();
        let mut min = f64::INFINITY;
        for i in 0..out.atoms.len() {
            for j in (i + 1)..out.atoms.len() {
                if bonded.contains(&(i, j)) {
                    continue;
                }
                min = min.min(out.atoms[i].pos.dist(out.atoms[j].pos));
            }
        }
        min
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

        let pre = preprocess_smiles(smiles).expect("preprocess failed");
        let mol = MoleculeGraph::from_smiles(
            &pre.smiles,
            pre.forced_direction_markers,
            pre.aromatic_atom_markers,
        )
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

    // ── Molecular weight ──────────────────────────────────────────────────────
    //
    // Reference values are PubChem's computed molecular weights, which use the
    // IUPAC/CIAAW standard atomic weights (the same table ptable embeds via
    // PubChemElements_all.json).

    fn assert_weight(smiles: &str, expected: f64) {
        let w = mol_weight_native(smiles).expect("mol weight failed");
        assert!(
            (w - expected).abs() < 0.01,
            "mol_weight({smiles}) = {w}, expected {expected}"
        );
    }

    /// Extracts every SMILES string literal passed to `smiles("...")` or
    /// `mol("...")` in the visual test file, undoing Typst string escapes.
    fn test_typ_smiles_strings() -> Vec<String> {
        let src = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../tests/test.typ"
        ))
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
        assert!(molecules.len() > 50, "extraction looks broken: {molecules:?}");
        let mut failures = Vec::new();
        for m in &molecules {
            match layout_native(m) {
                Ok(out) => {
                    let min = min_atom_distance(&out);
                    if out.atoms.len() > 1 && min < 0.5 {
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
        let out = layout_native("c1ccccc1").expect("benzene layout failed");
        assert_eq!(out.aromatic_rings.len(), 1);
        let ring = &out.aromatic_rings[0];
        // Hexagon with unit bonds: inradius ~0.866, so radius ~0.62.
        assert!((ring.radius - 0.866 * 0.72).abs() < 0.05);
        assert!(out.bonds.iter().filter(|b| b.aromatic).count() == 6);
    }

    #[test]
    fn kekule_input_emits_no_ring_circles() {
        let out = layout_native("C1=CC=CC=C1").expect("benzene layout failed");
        assert!(out.aromatic_rings.is_empty());
        assert!(out.bonds.iter().all(|b| !b.aromatic));
    }

    #[test]
    fn fused_aromatics_emit_one_circle_per_ring() {
        let out = layout_native("c1ccc2ccccc2c1").expect("naphthalene layout failed");
        assert_eq!(out.aromatic_rings.len(), 2);
    }

    #[test]
    fn aromatic_ring_with_saturated_neighbor_ring() {
        // Indane: only the aromatic ring gets a circle.
        let out = layout_native("c1ccc2CCCc2c1").expect("indane layout failed");
        assert_eq!(out.aromatic_rings.len(), 1);
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
}
