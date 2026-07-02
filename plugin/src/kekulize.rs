/// Kekulization of aromatic SMILES input (OpenSMILES §3.5).
///
/// Aromatic (lowercase) notation leaves the single/double bond assignment to
/// the reader. This module converts the aromatic bonds of a parsed molecule
/// into an explicit Kekulé structure:
///
///   1. Implicit bonds (and `:` bonds) between aromatic atoms are aromatic;
///      explicit `-` bonds are not (e.g. the biphenyl linker).
///   2. Aromatic bonds outside any ring degrade to single bonds; an aromatic
///      atom left without an aromatic ring bond is an error, since aromatic
///      atoms must be part of a ring.
///   3. Every aromatic atom that still has room for one more bond — by its
///      normal valence, adjusted for charge and declared hydrogens — must
///      receive exactly one double bond. Atoms whose valence is already
///      saturated (pyrrole [nH], furan o, exocyclic C=O carbons, …)
///      contribute a lone pair or nothing instead.
///   4. A perfect matching over those atoms picks the double bonds; the
///      remaining aromatic bonds become single. If no perfect matching
///      exists, the aromatic system cannot be kekulized and the SMILES is
///      rejected.
///
/// Atom and bond indices are untouched, so writing-order references
/// (`show-indices`, highlights, mechanism arrows) stay valid for aromatic
/// and Kekulé input alike.
use crate::graph::{BondOrder, MoleculeGraph};

pub(crate) fn kekulize(mol: &mut MoleculeGraph, bond_implicit: &[bool]) -> Result<(), String> {
    let any_aromatic = mol.atoms.iter().any(|a| a.aromatic)
        || mol.bonds.iter().any(|b| b.order == BondOrder::Aromatic);
    if !any_aromatic {
        return Ok(());
    }

    mark_aromatic_bonds(mol, bond_implicit)?;
    demote_nonring_aromatic_bonds(mol);
    check_aromatic_atoms_in_rings(mol)?;
    assign_double_bonds(mol)?;

    // Whatever aromatic bonds the matching did not turn into double bonds are
    // the single bonds of the Kekulé structure.
    for bond in mol.bonds.iter_mut() {
        if bond.order == BondOrder::Aromatic {
            bond.order = BondOrder::Single;
        }
    }
    Ok(())
}

/// Wildcard atoms may sit inside an aromatic ring and bond aromatically to
/// their aromatic neighbors without being aromatic themselves.
fn aromatic_pair(mol: &MoleculeGraph, a: usize, b: usize) -> bool {
    let capable = |i: usize| mol.atoms[i].aromatic || mol.atoms[i].symbol == "*";
    capable(a) && capable(b) && (mol.atoms[a].aromatic || mol.atoms[b].aromatic)
}

fn mark_aromatic_bonds(mol: &mut MoleculeGraph, bond_implicit: &[bool]) -> Result<(), String> {
    for i in 0..mol.bonds.len() {
        let (from, to) = (mol.bonds[i].from, mol.bonds[i].to);
        match mol.bonds[i].order {
            BondOrder::Aromatic => {
                if !aromatic_pair(mol, from, to) {
                    return Err(
                        "Aromatic bond ':' must connect two aromatic atoms".to_string()
                    );
                }
            }
            BondOrder::Single => {
                if bond_implicit.get(i).copied().unwrap_or(false) && aromatic_pair(mol, from, to)
                {
                    mol.bonds[i].order = BondOrder::Aromatic;
                }
            }
            _ => {}
        }
    }
    Ok(())
}

fn demote_nonring_aromatic_bonds(mol: &mut MoleculeGraph) {
    for i in 0..mol.bonds.len() {
        if mol.bonds[i].order == BondOrder::Aromatic && !bond_in_ring(mol, i) {
            mol.bonds[i].order = BondOrder::Single;
        }
    }
}

/// A bond lies in a ring iff its endpoints stay connected without it.
fn bond_in_ring(mol: &MoleculeGraph, bond_idx: usize) -> bool {
    let (from, to) = (mol.bonds[bond_idx].from, mol.bonds[bond_idx].to);
    let mut seen = vec![false; mol.n_atoms()];
    let mut stack = vec![from];
    seen[from] = true;
    while let Some(u) = stack.pop() {
        for &(v, b) in &mol.adj[u] {
            if b == bond_idx || seen[v] {
                continue;
            }
            if v == to {
                return true;
            }
            seen[v] = true;
            stack.push(v);
        }
    }
    false
}

fn check_aromatic_atoms_in_rings(mol: &MoleculeGraph) -> Result<(), String> {
    for (idx, atom) in mol.atoms.iter().enumerate() {
        if atom.aromatic
            && !mol.adj[idx]
                .iter()
                .any(|&(_, b)| mol.bonds[b].order == BondOrder::Aromatic)
        {
            return Err(format!(
                "Aromatic atom '{}' (atom {idx}) must be part of an aromatic ring",
                atom.symbol.to_lowercase()
            ));
        }
    }
    Ok(())
}

/// Normal valences of the elements that can be written aromatic.
fn aromatic_valence(symbol: &str) -> Option<i16> {
    match symbol {
        "B" => Some(3),
        "C" => Some(4),
        "N" | "P" | "As" => Some(3),
        "O" | "S" | "Se" => Some(2),
        _ => None,
    }
}

/// Whether an aromatic atom must receive one double bond during kekulization.
///
/// Counting aromatic bonds as single, the atom needs a double bond exactly
/// when its bonds plus declared hydrogens leave room under the charge-adjusted
/// normal valence (one π electron in the OpenSMILES table); a saturated atom
/// contributes a lone pair or an empty orbital instead (zero or two).
fn needs_double_bond(mol: &MoleculeGraph, idx: usize) -> bool {
    let atom = &mol.atoms[idx];
    if !atom.aromatic {
        return false;
    }
    let Some(valence) = aromatic_valence(&atom.symbol) else {
        return false;
    };
    let sigma: i16 = mol.adj[idx]
        .iter()
        .map(|&(_, b)| match mol.bonds[b].order {
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
        })
        .sum();
    valence + atom.charge as i16 - atom.hcount as i16 - sigma >= 1
}

fn assign_double_bonds(mol: &mut MoleculeGraph) -> Result<(), String> {
    let n = mol.n_atoms();
    let needs: Vec<bool> = (0..n).map(|i| needs_double_bond(mol, i)).collect();
    // Wildcards may absorb a double bond when the alternation requires it, but
    // are never required to take one.
    let flexible: Vec<bool> = (0..n).map(|i| mol.atoms[i].symbol == "*").collect();

    let mut mate: Vec<Option<usize>> = vec![None; n];
    if !match_atoms(mol, &needs, &flexible, &mut mate, 0) {
        return Err(
            "Cannot kekulize aromatic system: no valid alternating double-bond assignment \
             exists (check hydrogen counts, e.g. pyrrole is c1cc[nH]c1)"
                .to_string(),
        );
    }

    for u in 0..n {
        if let Some(v) = mate[u] {
            if u < v {
                let bond_idx = mol.adj[u]
                    .iter()
                    .find_map(|&(w, b)| {
                        (w == v && mol.bonds[b].order == BondOrder::Aromatic).then_some(b)
                    })
                    .expect("matched atoms share an aromatic bond");
                mol.bonds[bond_idx].order = BondOrder::Double;
            }
        }
    }
    Ok(())
}

/// Backtracking search for a matching that pairs every atom in `needs` with an
/// aromatic-bonded partner. Aromatic systems are small, so exhaustive search
/// with backtracking is fast enough.
fn match_atoms(
    mol: &MoleculeGraph,
    needs: &[bool],
    flexible: &[bool],
    mate: &mut Vec<Option<usize>>,
    start: usize,
) -> bool {
    let Some(u) = (start..needs.len()).find(|&i| needs[i] && mate[i].is_none()) else {
        return true;
    };
    for &(v, b) in &mol.adj[u] {
        if mol.bonds[b].order != BondOrder::Aromatic
            || mate[v].is_some()
            || !(needs[v] || flexible[v])
        {
            continue;
        }
        mate[u] = Some(v);
        mate[v] = Some(u);
        if match_atoms(mol, needs, flexible, mate, u + 1) {
            return true;
        }
        mate[u] = None;
        mate[v] = None;
    }
    false
}
