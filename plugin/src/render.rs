use serde::{Deserialize, Serialize};

/// Atom position and identity sent to the Typst renderer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtomOutput {
    pub symbol: String,
    pub pos: Vec2,
    /// Explicit hydrogen count (shown as subscript when > 0)
    pub hcount: u8,
    /// Formal charge (shown as superscript when != 0)
    pub charge: i8,
}

/// A bond between two atoms by their index in `LayoutOutput::atoms`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BondOutput {
    pub from: usize,
    pub to: usize,
    /// 1 = single, 2 = double, 3 = triple, 4 = aromatic
    pub order: u8,
    /// Stereochemistry: "none" | "wedge_up" | "wedge_down"
    pub stereo: String,
    /// For double bonds that are part of a ring: unit vector pointing from the
    /// bond midpoint toward the ring centroid (i.e., the "inside" direction).
    /// Both components are 0.0 for non-ring bonds → use symmetric offset.
    pub inner_x: f64,
    pub inner_y: f64,
}

/// 2D coordinate pair in layout-space units (1 unit = 1 bond length).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
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

    pub fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y)
    }

    pub fn scale(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s)
    }
}

/// Top-level output passed back to Typst as JSON.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayoutOutput {
    pub atoms: Vec<AtomOutput>,
    pub bonds: Vec<BondOutput>,
    /// Bounding box dimensions in bond-length units (for auto-scaling in Typst)
    pub bbox_width: f64,
    pub bbox_height: f64,
}
