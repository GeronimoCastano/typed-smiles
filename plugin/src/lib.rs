mod error;
mod graph;
mod layout;
mod render;

pub use render::LayoutOutput;

use graph::MoleculeGraph;
use layout::compute_layout;

// ── WASM / Typst plugin entrypoint ──────────────────────────────────────────

#[cfg(target_arch = "wasm32")]
mod wasm_entrypoint {
    use super::*;
    use wasm_minimal_protocol::*;

    initiate_protocol!();

    /// Called from Typst as `smiles-plugin.layout(bytes(smiles-str))`.
    /// Returns JSON-encoded `LayoutOutput`.
    #[wasm_func]
    pub fn layout(smiles: &[u8]) -> Result<Vec<u8>, String> {
        let s = core::str::from_utf8(smiles).map_err(|e| format!("UTF-8 error: {e}"))?;
        let mol = MoleculeGraph::from_smiles(s)?;
        let out = compute_layout(&mol)?;
        serde_json::to_vec(&out).map_err(|e| format!("JSON error: {e}"))
    }
}

// ── Native entrypoint for tests / CLI ───────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
pub fn layout_native(smiles: &str) -> Result<LayoutOutput, String> {
    let mol = MoleculeGraph::from_smiles(smiles)?;
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
    fn isobutane() {
        let out = layout_native("CC(C)C").expect("isobutane layout failed");
        assert_eq!(out.atoms.len(), 4);
        assert_eq!(out.bonds.len(), 3);
    }

    #[test]
    fn naphthalene_kekule() {
        // Naphthalene — fused bicyclic
        let out = layout_native("C1=CC2=CC=CC=C2C=C1").expect("naphthalene layout failed");
        assert_eq!(out.atoms.len(), 10);
        assert_eq!(out.bonds.len(), 11);
    }
}
