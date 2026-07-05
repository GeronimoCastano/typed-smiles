use serde::{Deserialize, Serialize};

/// Atom position and identity sent to the Typst renderer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtomOutput {
    pub symbol: String,
    pub pos: Vec2,
    /// Explicit hydrogen count (shown as subscript when > 0)
    pub hcount: u8,
    /// Computed implicit hydrogen count (shown only when requested)
    #[serde(default)]
    pub implicit_h: u8,
    /// Formal charge (shown as superscript when != 0)
    pub charge: i8,
    /// Mass number from a bracket isotope (shown as a leading superscript when > 0).
    #[serde(default)]
    pub isotope: u16,
    /// Non-bonding electron-pair count for optional skeletal lone-pair rendering.
    #[serde(default)]
    pub lone_pairs: u8,
    /// Unit directions from this atom toward each rendered lone-pair group.
    #[serde(default)]
    pub lone_pair_dirs: Vec<Vec2>,
    /// Non-empty when this atom was created by a `{label}` abbreviation.
    /// The renderer displays this text instead of `symbol`.
    #[serde(default)]
    pub abbrev: String,
    /// Optional color/style token from `{label|style}`.
    #[serde(default)]
    pub abbrev_style: String,
    /// Character index of the abbreviation attachment glyph.
    #[serde(default)]
    pub abbrev_anchor: usize,
    /// Character length of the abbreviation attachment glyph. Zero means centered.
    #[serde(default)]
    pub abbrev_anchor_len: usize,
    /// Parsed bracket chirality, e.g. "tetra_anti" or "tetra_clockwise".
    #[serde(default)]
    pub chirality: String,
    /// Stereochemical implicit/explicit hydrogen: "none" | "wedge_up" | "wedge_down".
    #[serde(default)]
    pub stereo_h: String,
    /// Unit direction from this atom toward the rendered stereochemical hydrogen.
    #[serde(default)]
    pub stereo_h_dir: Vec2,
    /// True for H atoms synthesized from bracket hcount to make them addressable
    /// via atom() references. The Typst renderer skips drawing these atoms and their
    /// bonds; they only appear in show-indices overlays and arrow endpoints.
    #[serde(default)]
    pub virtual_h: bool,
}

/// A bond between two atoms by their index in `LayoutOutput::atoms`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BondOutput {
    pub from: usize,
    pub to: usize,
    /// 1 = single, 2 = double, 3 = triple, 4 = aromatic
    pub order: u8,
    /// Rendered stereochemistry: "none" | "wedge_up" | "wedge_down"
    pub stereo: String,
    /// True when stereo came from a typed-smiles drawing extension (`!w`/`!h`).
    #[serde(default)]
    pub forced_stereo: bool,
    /// OpenSMILES directional marker: "none" | "up" | "down"
    #[serde(default)]
    pub direction: String,
    /// For double bonds that are part of a ring: unit vector pointing from the
    /// bond midpoint toward the ring centroid (i.e., the "inside" direction).
    /// Both components are 0.0 for non-ring bonds → use symmetric offset.
    pub inner_x: f64,
    pub inner_y: f64,
    /// True for bonds connecting a heavy atom to a virtual_h atom.
    /// The Typst renderer skips drawing these bonds; they exist only so
    /// atom() references resolve to meaningful positions.
    #[serde(default)]
    pub virtual_bond: bool,
    /// True for ring bonds that were aromatic in the input. Lets the renderer
    /// offer the inscribed-circle depiction instead of alternating doubles.
    #[serde(default)]
    pub aromatic: bool,
}

/// An aromatic ring eligible for the inscribed-circle depiction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AromaticRing {
    /// Ring centroid in layout coordinates.
    pub center: Vec2,
    /// Circle radius in bond-length units.
    pub radius: f64,
}

/// 2D coordinate pair in layout-space units (1 unit = 1 bond length).
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

impl Vec2 {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn dist(self, other: Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        (dx * dx + dy * dy).sqrt()
    }
}

/// Top-level output passed back to Typst as JSON.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayoutOutput {
    pub atoms: Vec<AtomOutput>,
    pub bonds: Vec<BondOutput>,
    /// Rings whose bonds were all aromatic in the input, for the optional
    /// circle depiction.
    #[serde(default)]
    pub aromatic_rings: Vec<AromaticRing>,
    /// Bounding box dimensions in bond-length units (for auto-scaling in Typst)
    pub bbox_width: f64,
    pub bbox_height: f64,
}
