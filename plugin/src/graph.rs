/// Builds an internal molecular graph by walking the smiles_parser Chain AST.
///
/// The parser's graph conversion does not preserve all information needed for
/// depiction, including ring closures, common heteroatoms, and folded explicit
/// hydrogens.
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
    /// Explicit lone-pair count from an abbreviation `lp=N` modifier. Zero when
    /// none was declared; abbreviations never infer lone pairs from their text.
    pub abbrev_lone_pairs: u8,
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
        let (remaining, chain) = parse_chain(smiles.as_bytes())
            .map_err(|error| format!("SMILES parse failed for {:?}: {error:?}", smiles))?;
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
        builder.walk_chain(&chain, None, None);

        let neighbor_bonds = builder
            .neighbor_slots
            .into_iter()
            .map(|slots| slots.into_iter().flatten().collect())
            .collect();

        let mut molecule = MoleculeGraph {
            atoms: builder.atoms,
            bonds: builder.bonds,
            adj: builder.adj,
            neighbor_bonds,
            has_preceding: builder.has_preceding,
            preceding_atom: builder.preceding_atom,
        };
        crate::kekulize::kekulize(&mut molecule, &builder.bond_implicit)?;
        Ok(molecule)
    }

    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }
}

#[derive(Default)]
struct GraphBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    adj: Vec<Vec<(usize, usize)>>,
    /// Per-atom bond indices in SMILES writing order. An empty slot reserves
    /// the writing position of a ring-opening bond until the ring closes.
    neighbor_slots: Vec<Vec<Option<usize>>>,
    /// has_preceding[i] = atom `i` had a "from" atom to its left.
    has_preceding: Vec<bool>,
    preceding_atom: Vec<Option<usize>>,
    open_rings: HashMap<u8, OpenRing>,
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

struct OpenRing {
    atom_index: usize,
    bond_specification: Option<BondSpec>,
    neighbor_slot: usize,
}

impl GraphBuilder {
    fn add_atom(&mut self, atom: Atom) -> usize {
        let atom_index = self.atoms.len();
        self.atoms.push(atom);
        self.adj.push(Vec::new());
        self.neighbor_slots.push(Vec::new());
        self.has_preceding.push(false);
        self.preceding_atom.push(None);
        atom_index
    }

    fn add_bond(&mut self, from: usize, to: usize, bond_specification: BondSpec, implicit: bool) {
        let bond_index = self.bonds.len();
        self.bonds.push(Bond {
            from,
            to,
            order: bond_specification.order,
            stereo: bond_specification.stereo,
            direction: bond_specification.direction,
            forced_stereo: bond_specification.forced_stereo,
            curl: bond_specification.curl,
            aromatic: false,
        });
        self.bond_implicit.push(implicit);
        self.adj[from].push((to, bond_index));
        self.adj[to].push((from, bond_index));
    }

    /// The parent passes its outgoing bond to the child because `bond_or_dot`
    /// belongs to the current chain node rather than its successor.
    fn walk_chain(
        &mut self,
        chain: &smiles_parser::Chain,
        preceding_atom: Option<usize>,
        incoming_bond: Option<BondSpec>,
    ) {
        let branched_atom = &chain.branched_atom;
        let mut atom = smiles_atom_to_atom(&branched_atom.atom);
        if matches!(branched_atom.atom, SAtom::AliphaticOrganic(_))
            && self.aromatic_atom_markers.pop_front().unwrap_or(false)
        {
            atom.aromatic = true;
        }

        if self.fold_terminal_hydrogen(&atom, chain, preceding_atom) {
            return;
        }

        let current_atom = self.add_atom(atom);
        if let Some(preceding_atom) = preceding_atom {
            self.connect_to_preceding_atom(current_atom, preceding_atom, incoming_bond);
        }

        self.process_ring_bonds(current_atom, &branched_atom.ring_bonds);
        self.process_branches(current_atom, &branched_atom.branches);
        self.continue_main_chain(current_atom, chain);
    }

    fn fold_terminal_hydrogen(
        &mut self,
        atom: &Atom,
        chain: &smiles_parser::Chain,
        preceding_atom: Option<usize>,
    ) -> bool {
        let Some(preceding_atom) = preceding_atom else {
            return false;
        };
        let branched_atom = &chain.branched_atom;
        let should_fold = atom.symbol == "H"
            && atom.charge == 0
            && atom.isotope.is_none()
            && branched_atom.ring_bonds.is_empty()
            && branched_atom.branches.is_empty()
            && chain.chain.is_none()
            && self.atoms[preceding_atom].symbol != "H";
        if !should_fold {
            return false;
        }

        self.atoms[preceding_atom].hcount = self.atoms[preceding_atom].hcount.saturating_add(1);
        self.atoms[preceding_atom].has_explicit_h = true;
        if let Some(neighbor_slot) = self.neighbor_slots[preceding_atom].last_mut() {
            *neighbor_slot = None;
        }
        true
    }

    fn connect_to_preceding_atom(
        &mut self,
        current_atom: usize,
        preceding_atom: usize,
        incoming_bond: Option<BondSpec>,
    ) {
        self.has_preceding[current_atom] = true;
        self.preceding_atom[current_atom] = Some(preceding_atom);
        let bond_index = self.bonds.len();
        let implicit = incoming_bond.is_none();
        self.add_bond(
            preceding_atom,
            current_atom,
            incoming_bond.unwrap_or_else(BondSpec::single),
            implicit,
        );
        self.neighbor_slots[current_atom].push(Some(bond_index));
    }

    fn process_ring_bonds(&mut self, current_atom: usize, ring_bonds: &[smiles_parser::RingBond]) {
        for ring_bond in ring_bonds {
            let ring_number = ring_bond.ring_number;
            let closing_specification = ring_bond
                .bond
                .as_ref()
                .map(|bond| self.parser_bond_to_specification(bond));

            if let Some(open_ring) = self.open_rings.remove(&ring_number) {
                self.close_ring(current_atom, open_ring, closing_specification);
            } else {
                self.open_ring(current_atom, ring_number, closing_specification);
            }
        }
    }

    fn close_ring(
        &mut self,
        current_atom: usize,
        open_ring: OpenRing,
        closing_specification: Option<BondSpec>,
    ) {
        let implicit = closing_specification.is_none() && open_ring.bond_specification.is_none();
        let bond_specification = closing_specification
            .or(open_ring.bond_specification)
            .unwrap_or_else(BondSpec::single);
        let bond_index = self.bonds.len();
        self.add_bond(
            open_ring.atom_index,
            current_atom,
            bond_specification,
            implicit,
        );
        self.neighbor_slots[open_ring.atom_index][open_ring.neighbor_slot] = Some(bond_index);
        self.neighbor_slots[current_atom].push(Some(bond_index));
    }

    fn open_ring(
        &mut self,
        current_atom: usize,
        ring_number: u8,
        bond_specification: Option<BondSpec>,
    ) {
        let neighbor_slot = self.neighbor_slots[current_atom].len();
        self.neighbor_slots[current_atom].push(None);
        self.open_rings.insert(
            ring_number,
            OpenRing {
                atom_index: current_atom,
                bond_specification,
                neighbor_slot,
            },
        );
    }

    fn process_branches(&mut self, current_atom: usize, branches: &[smiles_parser::Branch]) {
        for branch in branches {
            if matches!(branch.bond_or_dot, Some(BondOrDot::Dot(_))) {
                self.walk_chain(&branch.chain, None, None);
                continue;
            }

            let branch_bond =
                branch
                    .bond_or_dot
                    .as_ref()
                    .and_then(|bond_or_dot| match bond_or_dot {
                        BondOrDot::Bond(bond) => Some(self.parser_bond_to_specification(bond)),
                        BondOrDot::Dot(_) => None,
                    });
            self.neighbor_slots[current_atom].push(Some(self.bonds.len()));
            self.walk_chain(&branch.chain, Some(current_atom), branch_bond);
        }
    }

    fn continue_main_chain(&mut self, current_atom: usize, chain: &smiles_parser::Chain) {
        if let Some(next) = &chain.chain {
            if matches!(chain.bond_or_dot, Some(BondOrDot::Dot(_))) {
                self.walk_chain(next, None, None);
                return;
            }

            let outgoing_bond =
                chain
                    .bond_or_dot
                    .as_ref()
                    .and_then(|bond_or_dot| match bond_or_dot {
                        BondOrDot::Bond(bond) => Some(self.parser_bond_to_specification(bond)),
                        BondOrDot::Dot(_) => None,
                    });
            self.neighbor_slots[current_atom].push(Some(self.bonds.len()));
            self.walk_chain(next, Some(current_atom), outgoing_bond);
        }
    }

    fn parser_bond_to_specification(&mut self, parser_bond: &SBond) -> BondSpec {
        let marker = if matches!(parser_bond, SBond::Up | SBond::Down) {
            self.forced_direction_markers
                .pop_front()
                .unwrap_or(BondMarker {
                    style: BondMarkerStyle::Directional,
                    order: BondOrder::Single,
                    curl: false,
                })
        } else {
            return BondSpec {
                order: parser_bond_order(parser_bond),
                stereo: BondStereo::None,
                direction: parser_bond_direction(parser_bond),
                forced_stereo: false,
                curl: false,
            };
        };
        let (stereo, direction, forced_stereo) = match marker.style {
            BondMarkerStyle::Directional => {
                (BondStereo::None, parser_bond_direction(parser_bond), false)
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

fn smiles_atom_to_atom(atom: &SAtom) -> Atom {
    match atom {
        SAtom::AliphaticOrganic(organic_atom) => {
            atom_with_symbol(element_symbol(organic_atom.element))
        }
        SAtom::Bracket(bracket_atom) => bracket_to_atom(bracket_atom),
        SAtom::Unknown => atom_with_symbol("*".to_string()),
    }
}

fn atom_with_symbol(symbol: String) -> Atom {
    Atom {
        symbol,
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
        abbrev_lone_pairs: 0,
    }
}

fn bracket_to_atom(bracket_atom: &BracketAtom) -> Atom {
    let (symbol, aromatic) = match &bracket_atom.symbol {
        Symbol::ElementSymbol(element) => (element_symbol(*element), false),
        Symbol::AromaticSymbol(element) => (element_symbol(*element), true),
        Symbol::Unknown => ("*".to_string(), false),
    };
    Atom {
        symbol,
        aromatic,
        hcount: bracket_atom.hcount,
        has_explicit_h: true,
        isotope: bracket_atom.isotope,
        charge: bracket_atom.charge,
        chirality: bracket_atom
            .chiral
            .map(chirality_to_atom_chirality)
            .unwrap_or(AtomChirality::None),
        abbrev: String::new(),
        abbrev_style: String::new(),
        abbrev_anchor: 0,
        abbrev_anchor_len: 0,
        abbrev_lone_pairs: 0,
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

fn element_symbol(element: Element) -> String {
    element.get_symbol().to_string()
}

fn parser_bond_order(parser_bond: &SBond) -> BondOrder {
    match parser_bond {
        SBond::Single | SBond::Up | SBond::Down => BondOrder::Single,
        SBond::Double => BondOrder::Double,
        SBond::Triple => BondOrder::Triple,
        SBond::Quadruple => BondOrder::Quadruple,
        SBond::Aromatic => BondOrder::Aromatic,
    }
}

fn parser_bond_direction(parser_bond: &SBond) -> BondDirection {
    match parser_bond {
        SBond::Up => BondDirection::Up,
        SBond::Down => BondDirection::Down,
        _ => BondDirection::None,
    }
}
