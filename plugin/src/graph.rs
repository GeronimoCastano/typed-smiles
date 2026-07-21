/// Builds an internal molecular graph by walking the smiles_parser Chain AST.
///
/// We do NOT use smiles_parser::graph::MoleculeGraph because:
///   - it ignores ring-closure bonds (cyclic molecules would be wrong)
///   - it panics on N, S, halogens, and other common heteroatoms
///   - it adds explicit H atoms (which we want to keep implicit for depiction)
///
/// We walk the Chain tree recursively and handle ring-closure bonds ourselves
/// via a HashMap keyed on the SMILES ring-opening digit.
use std::collections::{HashMap, VecDeque};

use ptable::Element;
use smiles_parser::{
    chain as parse_chain, Atom as SAtom, Bond as SBond, BondOrDot, BracketAtom, Chirality, Symbol,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
}

impl BondOrder {
    pub fn as_u8(self) -> u8 {
        match self {
            Self::Single => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Quadruple => 4,
            // Aromatic orders are internal: kekulization replaces them with
            // single/double before any valence arithmetic or JSON output.
            Self::Aromatic => 1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondStereo {
    None,
    WedgeUp,
    WedgeDown,
    Wavy,
    Dashed,
}

impl BondStereo {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::WedgeUp => "wedge_up",
            Self::WedgeDown => "wedge_down",
            Self::Wavy => "wavy",
            Self::Dashed => "dashed",
        }
    }
}

/// Meaning of a `/` or `\` placeholder after preprocessing. `Directional`
/// tokens are genuine SMILES directional bonds; the others come from
/// typed-smiles drawing extensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondMarkerStyle {
    Directional,
    Plain,
    WedgeUp,
    WedgeDown,
    Wavy,
    Dashed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondMarker {
    pub style: BondMarkerStyle,
    pub order: BondOrder,
    pub curl: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondDirection {
    None,
    Up,
    Down,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondSpec {
    pub order: BondOrder,
    pub stereo: BondStereo,
    pub direction: BondDirection,
    pub forced_stereo: bool,
    pub curl: bool,
}

impl BondSpec {
    fn single() -> Self {
        Self {
            order: BondOrder::Single,
            stereo: BondStereo::None,
            direction: BondDirection::None,
            forced_stereo: false,
            curl: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomChirality {
    None,
    TetraAnti,
    TetraClockwise,
    /// `@SP1`..`@SP3`: square-planar shape class (1 = U, 2 = 4, 3 = Z).
    /// Depicted exactly, since the geometry is planar.
    SquarePlanar(u8),
    /// `@TB`, `@OH`, `@AL`: accepted; connectivity is drawn without stereo
    /// decoration.
    Undepicted,
    Unsupported,
}

impl AtomChirality {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::TetraAnti => "tetra_anti",
            Self::TetraClockwise => "tetra_clockwise",
            Self::SquarePlanar(_) => "square_planar",
            Self::Undepicted => "undepicted",
            Self::Unsupported => "unsupported",
        }
    }

    /// Whether this is a tetrahedral center whose stereo is rendered as
    /// wedge/hash bonds (possibly carried by the implicit hydrogen).
    pub fn is_tetrahedral(self) -> bool {
        matches!(self, Self::TetraAnti | Self::TetraClockwise)
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    /// Element symbol, e.g. "C", "N", "O"
    pub symbol: String,
    pub aromatic: bool,
    pub hcount: u8,
    pub has_explicit_h: bool,
    /// Mass number from a bracket isotope specification, e.g. `[2H]`.
    pub isotope: Option<u16>,
    pub charge: i8,
    pub chirality: AtomChirality,
    /// Non-empty when this atom was created by a `{label}` abbreviation substitution.
    pub abbrev: String,
    pub abbrev_style: String,
    pub abbrev_anchor: usize,
    pub abbrev_anchor_len: usize,
}

#[derive(Debug, Clone)]
pub struct Bond {
    pub from: usize,
    pub to: usize,
    pub order: BondOrder,
    pub stereo: BondStereo,
    pub direction: BondDirection,
    pub forced_stereo: bool,
    /// Repeat the preceding chain turn instead of alternating the zigzag.
    pub curl: bool,
    /// True for ring bonds that were aromatic in the input before
    /// kekulization assigned them a single or double order.
    pub aromatic: bool,
}

#[derive(Debug, Default)]
pub struct MoleculeGraph {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    /// adj[i] = [(neighbor_atom_idx, bond_idx), ...]
    pub adj: Vec<Vec<(usize, usize)>>,
    /// Bond indices per atom in SMILES writing order (ring bonds ordered by digit
    /// position). Used for OpenSMILES tetrahedral neighbor ordering.
    pub neighbor_bonds: Vec<Vec<usize>>,
    /// Whether each atom has a preceding ("from") atom to its left in the SMILES.
    pub has_preceding: Vec<bool>,
    /// Preceding atom in SMILES writing order, when one exists.
    pub preceding_atom: Vec<Option<usize>>,
}

impl MoleculeGraph {
    pub fn from_smiles(
        smiles: &str,
        forced_direction_markers: Vec<BondMarker>,
        aromatic_atom_markers: Vec<bool>,
    ) -> Result<Self, String> {
        // smiles_parser::chain takes &[u8] and returns IResult<&[u8], Chain>
        let (remaining, chain_ast) = parse_chain(smiles.as_bytes())
            .map_err(|e| format!("SMILES parse failed for {:?}: {e:?}", smiles))?;
        if !remaining.is_empty() {
            let tail = std::str::from_utf8(remaining).unwrap_or("<invalid UTF-8>");
            return Err(format!(
                "SMILES parse stopped before trailing input {:?} in {:?}",
                tail, smiles
            ));
        }

        let mut builder = GraphBuilder {
            forced_direction_markers: VecDeque::from(forced_direction_markers),
            aromatic_atom_markers: VecDeque::from(aromatic_atom_markers),
            ..GraphBuilder::default()
        };
        builder.walk_chain(&chain_ast, None, None);

        // Drop any unpatched ring-opening slot (only possible for malformed SMILES).
        let neighbor_bonds = builder
            .neighbor_slots
            .into_iter()
            .map(|slots| {
                slots
                    .into_iter()
                    .filter(|&b| b >= 0)
                    .map(|b| b as usize)
                    .collect()
            })
            .collect();

        let mut mol = MoleculeGraph {
            atoms: builder.atoms,
            bonds: builder.bonds,
            adj: builder.adj,
            neighbor_bonds,
            has_preceding: builder.has_preceding,
            preceding_atom: builder.preceding_atom,
        };
        crate::kekulize::kekulize(&mut mol, &builder.bond_implicit)?;
        Ok(mol)
    }

    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }
}

// ── Internal tree walker ──────────────────────────────────────────────────────

#[derive(Default)]
struct GraphBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    adj: Vec<Vec<(usize, usize)>>,
    /// Per-atom bond indices in SMILES writing order. A value of -1 is a
    /// placeholder reserved for a ring-opening bond, patched when the ring closes.
    neighbor_slots: Vec<Vec<i64>>,
    /// has_preceding[i] = atom `i` had a "from" atom to its left.
    has_preceding: Vec<bool>,
    preceding_atom: Vec<Option<usize>>,
    /// ring_number → (atom_idx, bond_spec_or_none, slot_index_in_opener)
    open_rings: HashMap<u8, (usize, Option<BondSpec>, usize)>,
    /// One entry for every `/` or `\` bond token seen by smiles-parser,
    /// recording whether it is genuine SMILES or a typed-smiles drawing
    /// extension (`!w`/`!h`/`!s`/`!d`) and which forced style it carries.
    forced_direction_markers: VecDeque<BondMarker>,
    /// One flag per unbracketed organic-subset atom token, in writing order.
    /// `true` means the atom was written lowercase (aromatic) and uppercased
    /// during preprocessing. Bracket atoms carry their own aromatic flag.
    aromatic_atom_markers: VecDeque<bool>,
    /// bond_implicit[i] = bond `i` had no bond symbol in the SMILES. Implicit
    /// bonds between aromatic atoms are aromatic; explicit `-` bonds are not.
    bond_implicit: Vec<bool>,
}

impl GraphBuilder {
    fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        self.adj.push(Vec::new());
        self.neighbor_slots.push(Vec::new());
        self.has_preceding.push(false);
        self.preceding_atom.push(None);
        idx
    }

    fn add_bond(&mut self, from: usize, to: usize, spec: BondSpec, implicit: bool) {
        let b_idx = self.bonds.len();
        self.bonds.push(Bond {
            from,
            to,
            order: spec.order,
            stereo: spec.stereo,
            direction: spec.direction,
            forced_stereo: spec.forced_stereo,
            curl: spec.curl,
            aromatic: false,
        });
        self.bond_implicit.push(implicit);
        self.adj[from].push((to, b_idx));
        self.adj[to].push((from, b_idx));
    }

    /// Walk a Chain node.
    ///
    /// `prev_atom`    — index of the preceding atom to bond to (or None at root).
    /// `incoming_bond` — bond order for the edge (prev → cur). Determined by the
    ///                   *parent* chain's `bond_or_dot`, not this node's own field.
    ///
    /// Key invariant: `chain.bond_or_dot` describes the bond from the CURRENT
    /// atom OUTWARD to the next atom in `chain.chain`. It must be forwarded as
    /// `incoming_bond` when recursing into `chain.chain`, not read by the child.
    fn walk_chain(
        &mut self,
        chain: &smiles_parser::Chain,
        prev_atom: Option<usize>,
        incoming: Option<BondSpec>,
    ) {
        let ba = &chain.branched_atom;

        // Add this atom. Unbracketed organic-subset atoms consume one aromatic
        // marker each; bracket atoms carry the flag from the parser instead.
        let mut atom = smiles_atom_to_atom(&ba.atom);
        if matches!(ba.atom, SAtom::AliphaticOrganic(_))
            && self.aromatic_atom_markers.pop_front().unwrap_or(false)
        {
            atom.aromatic = true;
        }

        // Fold a terminal plain hydrogen (`[H]`) into its heavy neighbor's
        // hydrogen count instead of keeping it as a drawn atom, mirroring how
        // depiction tools treat explicit hydrogens. Isotopic, charged, or
        // further-bonded hydrogens, and H2, stay as real atoms.
        if let Some(prev) = prev_atom {
            if atom.symbol == "H"
                && atom.charge == 0
                && atom.isotope.is_none()
                && ba.ring_bonds.is_empty()
                && ba.branches.is_empty()
                && chain.chain.is_none()
                && self.atoms[prev].symbol != "H"
            {
                self.atoms[prev].hcount = self.atoms[prev].hcount.saturating_add(1);
                self.atoms[prev].has_explicit_h = true;
                // Release the neighbor slot the parent reserved for this bond so
                // it is dropped from the parent's neighbor ordering.
                if let Some(slot) = self.neighbor_slots[prev].last_mut() {
                    *slot = -1;
                }
                return;
            }
        }

        let cur = self.add_atom(atom);

        // Bond to the previous atom using the bond the PARENT passed down; the
        // incoming bond is the first neighbor slot for `cur`.
        if let Some(prev) = prev_atom {
            self.has_preceding[cur] = true;
            self.preceding_atom[cur] = Some(prev);
            let b_idx = self.bonds.len();
            let implicit = incoming.is_none();
            self.add_bond(prev, cur, incoming.unwrap_or_else(BondSpec::single), implicit);
            self.neighbor_slots[cur].push(b_idx as i64);
        }

        // Process ring-closure bonds. The digit is written here, so it occupies a
        // neighbor slot at this position even though the bond is created later at
        // the closing atom (reserved as -1, patched on close).
        for rb in &ba.ring_bonds {
            let ring_num = rb.ring_number;
            let ring_spec = rb.bond.as_ref().map(|b| self.sparser_to_bond_spec(b));

            if let Some((open_atom, open_spec, open_slot)) = self.open_rings.remove(&ring_num) {
                let implicit = ring_spec.is_none() && open_spec.is_none();
                let spec = ring_spec.or(open_spec).unwrap_or_else(BondSpec::single);
                let b_idx = self.bonds.len();
                self.add_bond(open_atom, cur, spec, implicit);
                self.neighbor_slots[open_atom][open_slot] = b_idx as i64;
                self.neighbor_slots[cur].push(b_idx as i64);
            } else {
                let slot = self.neighbor_slots[cur].len();
                self.neighbor_slots[cur].push(-1);
                self.open_rings.insert(ring_num, (cur, ring_spec, slot));
            }
        }

        // Process branches — each branch starts from `cur`. A dot separator
        // means "no bond": the branch is a disconnected fragment and gets no
        // bond and no neighbor slot on `cur`.
        for branch in &ba.branches {
            if matches!(branch.bond_or_dot, Some(BondOrDot::Dot(_))) {
                self.walk_chain(&branch.chain, None, None);
                continue;
            }
            let branch_bond = branch.bond_or_dot.as_ref().and_then(|bod| match bod {
                BondOrDot::Bond(b) => Some(self.sparser_to_bond_spec(b)),
                BondOrDot::Dot(_) => None,
            });
            self.neighbor_slots[cur].push(self.bonds.len() as i64);
            self.walk_chain(&branch.chain, Some(cur), branch_bond);
        }

        // Continue the main chain, forwarding THIS node's outgoing bond. A dot
        // starts a new disconnected fragment instead of bonding to `cur`.
        if let Some(next) = &chain.chain {
            if matches!(chain.bond_or_dot, Some(BondOrDot::Dot(_))) {
                self.walk_chain(next, None, None);
                return;
            }
            let outgoing = chain.bond_or_dot.as_ref().and_then(|bod| match bod {
                BondOrDot::Bond(b) => Some(self.sparser_to_bond_spec(b)),
                BondOrDot::Dot(_) => None,
            });
            self.neighbor_slots[cur].push(self.bonds.len() as i64);
            self.walk_chain(next, Some(cur), outgoing);
        }
    }

    fn sparser_to_bond_spec(&mut self, b: &SBond) -> BondSpec {
        let marker = if matches!(b, SBond::Up | SBond::Down) {
            self.forced_direction_markers
                .pop_front()
                .unwrap_or(BondMarker {
                    style: BondMarkerStyle::Directional,
                    order: BondOrder::Single,
                    curl: false,
                })
        } else {
            return BondSpec {
                order: sparser_bond_to_order(b),
                stereo: BondStereo::None,
                direction: sparser_bond_to_direction(b),
                forced_stereo: false,
                curl: false,
            };
        };
        let (stereo, direction, forced_stereo) = match marker.style {
            BondMarkerStyle::Directional => {
                (BondStereo::None, sparser_bond_to_direction(b), false)
            }
            BondMarkerStyle::Plain => (BondStereo::None, BondDirection::None, false),
            BondMarkerStyle::WedgeUp => (BondStereo::WedgeUp, BondDirection::None, true),
            BondMarkerStyle::WedgeDown => (BondStereo::WedgeDown, BondDirection::None, true),
            BondMarkerStyle::Wavy => (BondStereo::Wavy, BondDirection::None, false),
            BondMarkerStyle::Dashed => (BondStereo::Dashed, BondDirection::None, false),
        };
        BondSpec {
            order: marker.order,
            stereo,
            direction,
            forced_stereo,
            curl: marker.curl,
        }
    }
}

// ── Conversion helpers ────────────────────────────────────────────────────────

fn smiles_atom_to_atom(atom: &SAtom) -> Atom {
    match atom {
        SAtom::AliphaticOrganic(a) => Atom {
            symbol: element_symbol(a.element),
            aromatic: false,
            hcount: 0,
            has_explicit_h: false,
            isotope: None,
            charge: 0,
            chirality: AtomChirality::None,
            abbrev: String::new(),
            abbrev_style: String::new(),
            abbrev_anchor: 0,
            abbrev_anchor_len: 0,
        },
        SAtom::Bracket(b) => bracket_to_atom(b),
        SAtom::Unknown => Atom {
            symbol: "*".to_string(),
            aromatic: false,
            hcount: 0,
            has_explicit_h: false,
            isotope: None,
            charge: 0,
            chirality: AtomChirality::None,
            abbrev: String::new(),
            abbrev_style: String::new(),
            abbrev_anchor: 0,
            abbrev_anchor_len: 0,
        },
    }
}

fn bracket_to_atom(b: &BracketAtom) -> Atom {
    let (symbol, aromatic) = match &b.symbol {
        Symbol::ElementSymbol(e) => (element_symbol(*e), false),
        Symbol::AromaticSymbol(e) => (element_symbol(*e), true),
        Symbol::Unknown => ("*".to_string(), false),
    };
    Atom {
        symbol,
        aromatic,
        hcount: b.hcount,
        has_explicit_h: true,
        isotope: b.isotope,
        charge: b.charge,
        chirality: b
            .chiral
            .map(chirality_to_atom_chirality)
            .unwrap_or(AtomChirality::None),
        abbrev: String::new(),
        abbrev_style: String::new(),
        abbrev_anchor: 0,
        abbrev_anchor_len: 0,
    }
}

fn chirality_to_atom_chirality(chirality: Chirality) -> AtomChirality {
    match chirality {
        Chirality::Anticlockwise | Chirality::Tetrahedral(1) => AtomChirality::TetraAnti,
        Chirality::Clockwise | Chirality::Tetrahedral(2) => AtomChirality::TetraClockwise,
        Chirality::SquarePlanar(class) => AtomChirality::SquarePlanar(class),
        Chirality::Allenal(_) | Chirality::TrigonalBipyramidal(_) | Chirality::Octahedral(_) => {
            AtomChirality::Undepicted
        }
        _ => AtomChirality::Unsupported,
    }
}

fn element_symbol(e: Element) -> String {
    e.get_symbol().to_string()
}

fn sparser_bond_to_order(b: &SBond) -> BondOrder {
    match b {
        SBond::Single | SBond::Up | SBond::Down => BondOrder::Single,
        SBond::Double => BondOrder::Double,
        SBond::Triple => BondOrder::Triple,
        SBond::Quadruple => BondOrder::Quadruple,
        SBond::Aromatic => BondOrder::Aromatic,
    }
}

fn sparser_bond_to_direction(b: &SBond) -> BondDirection {
    match b {
        SBond::Up => BondDirection::Up,
        SBond::Down => BondDirection::Down,
        _ => BondDirection::None,
    }
}
