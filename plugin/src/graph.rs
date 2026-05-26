/// Builds an internal molecular graph by walking the smiles_parser Chain AST.
///
/// We do NOT use smiles_parser::graph::MoleculeGraph because:
///   - it ignores ring-closure bonds (cyclic molecules would be wrong)
///   - it panics on N, S, halogens, and other common heteroatoms
///   - it adds explicit H atoms (which we want to keep implicit for depiction)
///
/// We walk the Chain tree recursively and handle ring-closure bonds ourselves
/// via a HashMap keyed on the SMILES ring-opening digit.

use std::collections::HashMap;

use smiles_parser::{
    chain as parse_chain, Atom as SAtom, Bond as SBond, BondOrDot, BracketAtom, Symbol,
};
use ptable::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondOrder {
    pub fn as_u8(self) -> u8 {
        match self {
            Self::Single => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Aromatic => 4,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    /// Element symbol, e.g. "C", "N", "O"
    pub symbol: String,
    #[allow(dead_code)]
    pub aromatic: bool,
    pub hcount: u8,
    pub charge: i8,
}

#[derive(Debug, Clone)]
pub struct Bond {
    pub from: usize,
    pub to: usize,
    pub order: BondOrder,
}

#[derive(Debug, Default)]
pub struct MoleculeGraph {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    /// adj[i] = [(neighbor_atom_idx, bond_idx), ...]
    pub adj: Vec<Vec<(usize, usize)>>,
}

impl MoleculeGraph {
    pub fn from_smiles(smiles: &str) -> Result<Self, String> {
        // smiles_parser::chain takes &[u8] and returns IResult<&[u8], Chain>
        let (_, chain_ast) = parse_chain(smiles.as_bytes())
            .map_err(|e| format!("SMILES parse failed for {:?}: {e:?}", smiles))?;

        let mut builder = GraphBuilder::default();
        builder.walk_chain(&chain_ast, None, None);

        Ok(MoleculeGraph {
            atoms: builder.atoms,
            bonds: builder.bonds,
            adj: builder.adj,
        })
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
    /// ring_number → (atom_idx, bond_order_or_none)
    open_rings: HashMap<u8, (usize, Option<BondOrder>)>,
}

impl GraphBuilder {
    fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        self.adj.push(Vec::new());
        idx
    }

    fn add_bond(&mut self, from: usize, to: usize, order: BondOrder) {
        let b_idx = self.bonds.len();
        self.bonds.push(Bond { from, to, order });
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
        incoming_bond: Option<BondOrder>,
    ) {
        let ba = &chain.branched_atom;

        // Add this atom
        let atom = smiles_atom_to_atom(&ba.atom);
        let cur = self.add_atom(atom);

        // Bond to the previous atom using the bond the PARENT passed down
        if let Some(prev) = prev_atom {
            self.add_bond(prev, cur, incoming_bond.unwrap_or(BondOrder::Single));
        }

        // Process ring-closure bonds
        for rb in &ba.ring_bonds {
            let ring_num = rb.ring_number;
            let ring_order = rb.bond.as_ref().map(sparser_bond_to_order);

            if let Some((open_atom, open_order)) = self.open_rings.remove(&ring_num) {
                // Close the ring; prefer whichever side declared an explicit bond type
                let order = ring_order.or(open_order).unwrap_or(BondOrder::Single);
                self.add_bond(open_atom, cur, order);
            } else {
                self.open_rings.insert(ring_num, (cur, ring_order));
            }
        }

        // Process branches — each branch starts from `cur`
        for branch in &ba.branches {
            let branch_bond = branch
                .bond_or_dot
                .as_ref()
                .and_then(|bod| match bod {
                    BondOrDot::Bond(b) => Some(sparser_bond_to_order(b)),
                    BondOrDot::Dot(_) => None, // disconnected salt; no bond
                });
            self.walk_chain(&branch.chain, Some(cur), branch_bond);
        }

        // Continue the main chain, forwarding THIS node's outgoing bond
        if let Some(next) = &chain.chain {
            let outgoing = chain.bond_or_dot.as_ref().and_then(|bod| match bod {
                BondOrDot::Bond(b) => Some(sparser_bond_to_order(b)),
                BondOrDot::Dot(_) => None,
            });
            self.walk_chain(next, Some(cur), outgoing);
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
            charge: 0,
        },
        SAtom::Bracket(b) => bracket_to_atom(b),
        SAtom::Unknown => Atom {
            symbol: "*".to_string(),
            aromatic: false,
            hcount: 0,
            charge: 0,
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
        charge: b.charge,
    }
}

fn element_symbol(e: Element) -> String {
    e.get_symbol().to_string()
}

fn sparser_bond_to_order(b: &SBond) -> BondOrder {
    match b {
        SBond::Single | SBond::Up | SBond::Down => BondOrder::Single,
        SBond::Double => BondOrder::Double,
        SBond::Triple | SBond::Quadruple => BondOrder::Triple,
        SBond::Aromatic => BondOrder::Aromatic,
    }
}

