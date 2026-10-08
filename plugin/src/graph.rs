/// The molecular graph shared by layout, kekulization, formulas, and
/// substructure matching, together with the builder that the SMILES reader
/// drives while it walks the input.
///
/// Bonds are oriented from the atom whose bond symbol describes them. A chain
/// or branch bond therefore points from the preceding atom to the new atom, and
/// a ring-closure bond points away from the end that carries its bond symbol.
/// Directional `/` and `\` symbols and the tip of a drawn wedge are read in that
/// orientation.
use std::collections::HashMap;

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

    pub fn is_wedge(self) -> bool {
        matches!(self, Self::WedgeUp | Self::WedgeDown)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondDirection {
    None,
    Up,
    Down,
}

impl BondDirection {
    /// The same directional bond read from its other end: `a/b` is `b\a`.
    pub fn reversed(self) -> Self {
        match self {
            Self::None => Self::None,
            Self::Up => Self::Down,
            Self::Down => Self::Up,
        }
    }
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
    pub fn single() -> Self {
        Self {
            order: BondOrder::Single,
            stereo: BondStereo::None,
            direction: BondDirection::None,
            forced_stereo: false,
            curl: false,
        }
    }

    pub fn with_order(order: BondOrder) -> Self {
        Self {
            order,
            ..Self::single()
        }
    }

    /// Whether this bond carries a typed-smiles drawing extension.
    fn is_drawing_extension(self) -> bool {
        self.stereo != BondStereo::None || self.curl
    }
}

/// A bond symbol as written in the input, with its 1-based character position.
#[derive(Debug, Clone, PartialEq)]
pub struct WrittenBond {
    pub specification: BondSpec,
    pub text: String,
    pub position: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomChirality {
    None,
    /// `@` or `@TH1`.
    TetraAnti,
    /// `@@` or `@TH2`.
    TetraClockwise,
    /// `@AL1`, `@AL2`: allene-like (extended tetrahedral) centers.
    Allenal(u8),
    /// `@SP1`..`@SP3`: square-planar shape class (1 = U, 2 = 4, 3 = Z).
    /// Depicted exactly, since the geometry is planar.
    SquarePlanar(u8),
    /// `@TB1`..`@TB20`.
    TrigonalBipyramidal(u8),
    /// `@OH1`..`@OH30`.
    Octahedral(u8),
}

impl AtomChirality {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::TetraAnti => "tetra_anti",
            Self::TetraClockwise => "tetra_clockwise",
            Self::Allenal(_) => "allenal",
            Self::SquarePlanar(_) => "square_planar",
            Self::TrigonalBipyramidal(_) => "trigonal_bipyramidal",
            Self::Octahedral(_) => "octahedral",
        }
    }

    /// The SMILES notation for this configuration, as used in diagnostics.
    pub fn notation(self) -> String {
        match self {
            Self::None => String::new(),
            Self::TetraAnti => "@".to_string(),
            Self::TetraClockwise => "@@".to_string(),
            Self::Allenal(class) => format!("@AL{class}"),
            Self::SquarePlanar(class) => format!("@SP{class}"),
            Self::TrigonalBipyramidal(class) => format!("@TB{class}"),
            Self::Octahedral(class) => format!("@OH{class}"),
        }
    }

    /// Whether this is a tetrahedral center whose stereo is rendered as
    /// wedge/hash bonds (possibly carried by the implicit hydrogen).
    pub fn is_tetrahedral(self) -> bool {
        matches!(self, Self::TetraAnti | Self::TetraClockwise)
    }

    /// The same tetrahedral configuration after swapping two neighbors.
    fn inverted(self) -> Self {
        match self {
            Self::TetraAnti => Self::TetraClockwise,
            Self::TetraClockwise => Self::TetraAnti,
            other => other,
        }
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
    /// OpenSMILES atom class, commonly used as a reaction atom map, e.g. the
    /// 7 in `[CH3:7]`. Zero means unmapped, as in OpenSMILES. Map numbers are
    /// user-chosen labels and never stand in for atom indices.
    pub atom_map: u32,
    /// 1-based character position of the atom in the input, for diagnostics.
    pub source_position: usize,
    /// Non-empty when this atom was created by a `{label}` abbreviation substitution.
    pub abbrev: String,
    pub abbrev_style: String,
    pub abbrev_anchor: usize,
    pub abbrev_anchor_len: usize,
    /// Explicit lone-pair count from an abbreviation `lp=N` modifier. Zero when
    /// none was declared; abbreviations never infer lone pairs from their text.
    pub abbrev_lone_pairs: u8,
    /// Page-space rendering displacement for an abbreviation, in bond lengths.
    /// The layout algorithm never reads these values.
    pub abbrev_offset_x: f64,
    pub abbrev_offset_y: f64,
}

impl Atom {
    pub fn new(symbol: &str, source_position: usize) -> Self {
        Self {
            symbol: symbol.to_string(),
            aromatic: false,
            hcount: 0,
            has_explicit_h: false,
            isotope: None,
            charge: 0,
            chirality: AtomChirality::None,
            atom_map: 0,
            source_position,
            abbrev: String::new(),
            abbrev_style: String::new(),
            abbrev_anchor: 0,
            abbrev_anchor_len: 0,
            abbrev_lone_pairs: 0,
            abbrev_offset_x: 0.0,
            abbrev_offset_y: 0.0,
        }
    }

    /// How the atom is written in SMILES, for diagnostics: `c` for an aromatic
    /// carbon, `C` otherwise.
    pub fn written_symbol(&self) -> String {
        if self.aromatic {
            self.symbol.to_lowercase()
        } else {
            self.symbol.clone()
        }
    }

    /// Whether this is a plain hydrogen atom that can be shown as part of its
    /// neighbor's hydrogen count without losing any written information.
    pub fn is_foldable_hydrogen(&self) -> bool {
        self.symbol == "H"
            && self.charge == 0
            && self.isotope.is_none()
            && self.hcount == 0
            && self.atom_map == 0
            && self.chirality == AtomChirality::None
    }
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
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Index of the bond joining two atoms, if they are bonded.
    pub(crate) fn bond_between(&self, first_atom: usize, second_atom: usize) -> Option<usize> {
        self.adj[first_atom]
            .iter()
            .find_map(|&(neighbor, bond_index)| (neighbor == second_atom).then_some(bond_index))
    }
}

/// Accumulates atoms and bonds in SMILES writing order and checks that every
/// bond joins two distinct atoms at most once.
#[derive(Default)]
pub struct GraphBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    adj: Vec<Vec<(usize, usize)>>,
    /// Per-atom bond indices in SMILES writing order. An empty slot reserves
    /// the writing position of a ring-opening bond until the ring closes.
    neighbor_slots: Vec<Vec<Option<usize>>>,
    has_preceding: Vec<bool>,
    preceding_atom: Vec<Option<usize>>,
    open_rings: HashMap<u8, OpenRing>,
    /// bond_implicit[i] = bond `i` had no bond symbol in the SMILES. Implicit
    /// bonds between aromatic atoms are aromatic; explicit `-` bonds are not.
    bond_implicit: Vec<bool>,
}

struct OpenRing {
    atom_index: usize,
    bond: Option<WrittenBond>,
    neighbor_slot: usize,
    position: usize,
}

impl GraphBuilder {
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let atom_index = self.atoms.len();
        self.atoms.push(atom);
        self.adj.push(Vec::new());
        self.neighbor_slots.push(Vec::new());
        self.has_preceding.push(false);
        self.preceding_atom.push(None);
        atom_index
    }

    /// Joins a newly written atom to the atom before it in the chain or branch.
    pub fn bond_to_preceding_atom(
        &mut self,
        preceding_atom: usize,
        current_atom: usize,
        written_bond: Option<&WrittenBond>,
    ) {
        self.has_preceding[current_atom] = true;
        self.preceding_atom[current_atom] = Some(preceding_atom);
        let bond_index = self.add_bond(
            preceding_atom,
            current_atom,
            written_bond.map_or_else(BondSpec::single, |bond| bond.specification),
            written_bond.is_none(),
        );
        self.neighbor_slots[preceding_atom].push(Some(bond_index));
        self.neighbor_slots[current_atom].push(Some(bond_index));
    }

    /// Records a plain hydrogen written as its own atom, such as the `[H]` in
    /// `C([H])O`, as part of its neighbor's hydrogen count.
    ///
    /// Tetrahedral neighbor order places a counted hydrogen immediately after
    /// the preceding atom (or first, when there is none). Moving the written
    /// hydrogen there passes over every neighbor written between those two
    /// positions; each one swaps the configuration, so an odd number of them
    /// inverts `@`/`@@` and the depicted stereocenter stays the written one.
    pub fn fold_hydrogen_into(&mut self, parent_atom: usize) {
        let written_position = self.neighbor_slots[parent_atom].len();
        let counted_position = usize::from(self.has_preceding[parent_atom]);
        let parent = &mut self.atoms[parent_atom];
        if (written_position - counted_position) % 2 == 1 {
            parent.chirality = parent.chirality.inverted();
        }
        parent.hcount = parent.hcount.saturating_add(1);
        parent.has_explicit_h = true;
    }

    /// Whether a plain hydrogen bonded to `parent_atom` may become part of the
    /// parent's hydrogen count. Non-tetrahedral configurations depend on the
    /// written position of every neighbor, so their hydrogens stay atoms.
    pub fn can_fold_hydrogen_into(&self, parent_atom: usize) -> bool {
        let parent = &self.atoms[parent_atom];
        parent.symbol != "H"
            && matches!(
                parent.chirality,
                AtomChirality::None | AtomChirality::TetraAnti | AtomChirality::TetraClockwise
            )
    }

    /// Opens or closes ring `ring_number` at `current_atom`.
    pub fn ring_bond(
        &mut self,
        current_atom: usize,
        ring_number: u8,
        written_bond: Option<WrittenBond>,
        position: usize,
    ) -> Result<(), String> {
        match self.open_rings.remove(&ring_number) {
            Some(open_ring) => {
                self.close_ring(current_atom, ring_number, open_ring, written_bond, position)
            }
            None => {
                self.open_ring(current_atom, ring_number, written_bond, position);
                Ok(())
            }
        }
    }

    fn open_ring(
        &mut self,
        current_atom: usize,
        ring_number: u8,
        bond: Option<WrittenBond>,
        position: usize,
    ) {
        let neighbor_slot = self.neighbor_slots[current_atom].len();
        self.neighbor_slots[current_atom].push(None);
        self.open_rings.insert(
            ring_number,
            OpenRing {
                atom_index: current_atom,
                bond,
                neighbor_slot,
                position,
            },
        );
    }

    fn close_ring(
        &mut self,
        current_atom: usize,
        ring_number: u8,
        open_ring: OpenRing,
        closing_bond: Option<WrittenBond>,
        position: usize,
    ) -> Result<(), String> {
        let opening_atom = open_ring.atom_index;
        if opening_atom == current_atom {
            return Err(format!(
                "ring closure {} at character {position} closes on the atom that opened it at \
                 character {}; a bond cannot join an atom to itself, so close the ring on a \
                 different atom",
                ring_label(ring_number),
                open_ring.position
            ));
        }
        if self.adj[opening_atom]
            .iter()
            .any(|&(neighbor, _)| neighbor == current_atom)
        {
            return Err(format!(
                "ring closure {} at character {position} would add a second bond between \
                 atoms that are already bonded ({} at character {} and {} at character {}); \
                 write one bond with a bond-order symbol such as `=` instead",
                ring_label(ring_number),
                self.atoms[opening_atom].written_symbol(),
                self.atoms[opening_atom].source_position,
                self.atoms[current_atom].written_symbol(),
                self.atoms[current_atom].source_position,
            ));
        }

        let (from, to, specification) = match (&open_ring.bond, &closing_bond) {
            (None, None) => (opening_atom, current_atom, BondSpec::single()),
            (Some(opening), None) => (opening_atom, current_atom, opening.specification),
            (None, Some(closing)) => (current_atom, opening_atom, closing.specification),
            (Some(opening), Some(closing)) => {
                require_matching_ring_bonds(ring_number, opening, closing)?;
                (opening_atom, current_atom, opening.specification)
            }
        };
        let implicit = open_ring.bond.is_none() && closing_bond.is_none();
        let bond_index = self.add_bond(from, to, specification, implicit);
        self.neighbor_slots[opening_atom][open_ring.neighbor_slot] = Some(bond_index);
        self.neighbor_slots[current_atom].push(Some(bond_index));
        Ok(())
    }

    fn add_bond(
        &mut self,
        from: usize,
        to: usize,
        bond_specification: BondSpec,
        implicit: bool,
    ) -> usize {
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
        bond_index
    }

    /// Completes the graph once the whole input has been read.
    pub fn finish(self) -> Result<(MoleculeGraph, Vec<bool>), String> {
        if let Some(unclosed_ring) = self
            .open_rings
            .iter()
            .min_by_key(|(_, open_ring)| open_ring.position)
        {
            let (ring_number, open_ring) = unclosed_ring;
            return Err(format!(
                "ring closure {} opened at character {} is never closed; repeat the ring number \
                 on the atom that closes the ring",
                ring_label(*ring_number),
                open_ring.position
            ));
        }

        let neighbor_bonds = self
            .neighbor_slots
            .into_iter()
            .map(|slots| slots.into_iter().flatten().collect())
            .collect();
        let molecule = MoleculeGraph {
            atoms: self.atoms,
            bonds: self.bonds,
            adj: self.adj,
            neighbor_bonds,
            has_preceding: self.has_preceding,
            preceding_atom: self.preceding_atom,
        };
        Ok((molecule, self.bond_implicit))
    }
}

/// Both ends of a ring closure may carry a bond symbol only when they describe
/// the same bond. Directional symbols are read away from the end that carries
/// them, so `/1` at one end agrees with `\1` at the other.
fn require_matching_ring_bonds(
    ring_number: u8,
    opening: &WrittenBond,
    closing: &WrittenBond,
) -> Result<(), String> {
    let opening_specification = opening.specification;
    let closing_specification = closing.specification;
    if opening_specification.is_drawing_extension() || closing_specification.is_drawing_extension()
    {
        return Err(format!(
            "ring closure {} has bond symbols on both ends (`{}` at character {} and `{}` at \
             character {}); write a drawing extension on one end of the ring closure only",
            ring_label(ring_number),
            opening.text,
            opening.position,
            closing.text,
            closing.position
        ));
    }
    let same_order = opening_specification.order == closing_specification.order;
    let same_direction =
        opening_specification.direction == closing_specification.direction.reversed();
    if same_order && same_direction {
        return Ok(());
    }

    let explanation = if same_order {
        "; a direction is read away from the end that carries it, so matching ends use \
         opposite symbols such as `/1` and `\\1`"
    } else {
        ""
    };
    Err(format!(
        "ring closure {} has conflicting bond symbols `{}` at character {} and `{}` at \
         character {}{explanation}; write the bond symbol on one end only",
        ring_label(ring_number),
        opening.text,
        opening.position,
        closing.text,
        closing.position
    ))
}

/// Ring numbers above 9 are written with a `%` prefix.
pub fn ring_label(ring_number: u8) -> String {
    if ring_number > 9 {
        format!("%{ring_number:02}")
    } else {
        ring_number.to_string()
    }
}
