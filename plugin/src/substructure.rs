//! A bounded SMARTS subset matcher on the depiction graph. No second SMILES
//! parser or aromaticity model: all returned indices refer to real graph atoms.
use crate::graph::{BondOrder, MoleculeGraph};
use ptable::Element;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

const MAX_PATTERN_BYTES: usize = 4096;
const MAX_QUERY_ATOMS: usize = 64;
const MAX_DEPTH: usize = 16;
const MAX_VISITS: usize = 1_000_000;
const MAX_MATCHES: usize = 4096;

#[derive(Debug, Serialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct SubstructureMatch {
    pub atoms: Vec<usize>,
    /// Query edges only, as sorted endpoint pairs (not induced-subgraph edges).
    pub bonds: Vec<(usize, usize)>,
}

#[derive(Debug)]
enum AtomQuery {
    Any,
    Element(String, bool),
    AtomicNumber(usize),
    Aromatic(bool),
    Hydrogens(usize),
    Degree(usize),
    Connectivity(usize),
    Charge(i8),
    Ring(bool),
    Recursive(Box<Query>),
    Not(Box<AtomQuery>),
    And(Box<AtomQuery>, Box<AtomQuery>),
    Or(Box<AtomQuery>, Box<AtomQuery>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BondQuery {
    Implicit,
    Single,
    Double,
    Triple,
    Aromatic,
    Any,
}

#[derive(Debug, Default)]
struct Query {
    atoms: Vec<AtomQuery>,
    bonds: Vec<(usize, usize, BondQuery)>,
}

struct Parser<'a> {
    input: &'a [u8],
    pos: usize,
    depth: usize,
}

impl Parser<'_> {
    fn peek(&self) -> Option<u8> {
        self.input.get(self.pos).copied()
    }

    fn error(&self, problem: &str) -> String {
        format!(
            "{problem} at byte {}; use the documented SMARTS subset",
            self.pos
        )
    }

    fn enter(&mut self) -> Result<(), String> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return Err(self.error("pattern nesting exceeds 16 levels"));
        }
        Ok(())
    }

    fn number(&mut self) -> Result<usize, String> {
        let start = self.pos;
        while self.peek().is_some_and(|ch| ch.is_ascii_digit()) {
            self.pos += 1;
        }
        if start == self.pos {
            return Err(self.error("expected an integer"));
        }
        std::str::from_utf8(&self.input[start..self.pos])
            .unwrap()
            .parse()
            .map_err(|_| self.error("integer is too large"))
    }

    fn bond(&mut self) -> Option<BondQuery> {
        let bond = match self.peek()? {
            b'-' => BondQuery::Single,
            b'=' => BondQuery::Double,
            b'#' => BondQuery::Triple,
            b':' => BondQuery::Aromatic,
            b'~' => BondQuery::Any,
            _ => return None,
        };
        self.pos += 1;
        Some(bond)
    }

    fn graph(&mut self, nested: bool) -> Result<Query, String> {
        let mut query = Query::default();
        let mut current = None;
        let mut pending = None;
        let mut branches = Vec::new();
        let mut rings: BTreeMap<usize, (usize, Option<BondQuery>)> = BTreeMap::new();
        // A branch must actually introduce an atom, not just a ring marker.
        let mut needs_atom = true;
        loop {
            match self.peek() {
                None => break,
                Some(b')') if branches.is_empty() && nested => break,
                Some(b'(') => {
                    let origin =
                        current.ok_or_else(|| self.error("branch needs a preceding atom"))?;
                    if pending.is_some() || needs_atom {
                        return Err(self.error("unexpected branch"));
                    }
                    self.enter()?;
                    branches.push((origin, query.atoms.len()));
                    self.pos += 1;
                    needs_atom = true;
                }
                Some(b')') => {
                    let (origin, count) =
                        branches.pop().ok_or_else(|| self.error("unmatched ')'"))?;
                    if needs_atom || pending.is_some() || query.atoms.len() == count {
                        return Err(self.error("empty or incomplete branch"));
                    }
                    self.depth -= 1;
                    current = Some(origin);
                    self.pos += 1;
                }
                Some(b'.') => {
                    if current.is_none() || pending.is_some() || needs_atom {
                        return Err(self.error("unexpected '.'"));
                    }
                    current = None;
                    self.pos += 1;
                    needs_atom = true;
                }
                Some(ch) if ch.is_ascii_digit() || ch == b'%' => {
                    let origin = current.ok_or_else(|| self.error("ring closure needs an atom"))?;
                    if needs_atom {
                        return Err(self.error("branch needs an atom before a ring closure"));
                    }
                    let ring = if ch == b'%' {
                        self.pos += 1;
                        let start = self.pos;
                        let number = self.number()?;
                        if self.pos - start != 2 {
                            return Err(self.error("'%' ring numbers need exactly two digits"));
                        }
                        number
                    } else {
                        self.pos += 1;
                        (ch - b'0') as usize
                    };
                    if let Some((other, opening)) = rings.remove(&ring) {
                        if origin == other
                            || query.bonds.iter().any(|&(a, b, _)| {
                                (a == origin && b == other) || (a == other && b == origin)
                            })
                        {
                            return Err(self.error("ring closure creates a self or duplicate bond"));
                        }
                        if opening.is_some() && pending.is_some() && opening != pending {
                            return Err(self.error("conflicting ring bond specifications"));
                        }
                        query.bonds.push((
                            other,
                            origin,
                            pending.take().or(opening).unwrap_or(BondQuery::Implicit),
                        ));
                    } else {
                        rings.insert(ring, (origin, pending.take()));
                    }
                }
                Some(b'-' | b'=' | b'#' | b':' | b'~') => {
                    if current.is_none() || pending.is_some() {
                        return Err(self.error("bond needs two atoms"));
                    }
                    pending = self.bond();
                }
                Some(_) => {
                    let atom = self.atom()?;
                    if query.atoms.len() >= MAX_QUERY_ATOMS {
                        return Err(self.error("patterns support at most 64 atoms"));
                    }
                    let index = query.atoms.len();
                    query.atoms.push(atom);
                    if let Some(from) = current {
                        query.bonds.push((
                            from,
                            index,
                            pending.take().unwrap_or(BondQuery::Implicit),
                        ));
                    }
                    current = Some(index);
                    needs_atom = false;
                }
            }
        }
        if needs_atom || pending.is_some() {
            return Err(self.error("empty pattern or missing atom after a bond or '.'"));
        }
        if !branches.is_empty() {
            return Err(self.error("unclosed '(' branch"));
        }
        if !rings.is_empty() {
            return Err(self.error("unclosed ring number; repeat it on the closing atom"));
        }
        Ok(query)
    }

    fn atom(&mut self) -> Result<AtomQuery, String> {
        match self.peek() {
            Some(b'[') => {
                self.pos += 1;
                self.enter()?;
                let query = self.low_and()?;
                if self.peek() != Some(b']') {
                    return Err(self.error("unclosed '[' or unsupported atom predicate"));
                }
                self.pos += 1;
                self.depth -= 1;
                Ok(query)
            }
            Some(b'*') => {
                self.pos += 1;
                Ok(AtomQuery::Any)
            }
            Some(b'a' | b'A') => {
                let aromatic = self.peek() == Some(b'a');
                self.pos += 1;
                Ok(AtomQuery::Aromatic(aromatic))
            }
            _ => self.element(),
        }
    }

    fn element(&mut self) -> Result<AtomQuery, String> {
        let first = self.peek().ok_or_else(|| self.error("expected an atom"))?;
        if !first.is_ascii_alphabetic() {
            return Err(self.error("expected an atom; unsupported SMARTS syntax"));
        }
        let aromatic = first.is_ascii_lowercase();
        let mut symbol = String::from(first.to_ascii_uppercase() as char);
        self.pos += 1;
        // Consume a second letter only for a valid element (Cc means C then c).
        if let Some(second) = self.peek().filter(u8::is_ascii_lowercase) {
            let candidate = format!("{symbol}{}", second as char);
            if crate::element_from_symbol(&candidate).is_some()
                && (!aromatic || matches!(candidate.as_str(), "As" | "Se"))
            {
                symbol = candidate;
                self.pos += 1;
            }
        }
        if crate::element_from_symbol(&symbol).is_none()
            || (aromatic
                && !matches!(
                    symbol.as_str(),
                    "B" | "C" | "N" | "O" | "P" | "S" | "As" | "Se"
                ))
        {
            return Err(self.error("unknown element or unsupported aromatic symbol"));
        }
        Ok(AtomQuery::Element(symbol, aromatic))
    }

    // SMARTS precedence: negation > juxtaposition/& > comma > semicolon.
    fn low_and(&mut self) -> Result<AtomQuery, String> {
        let mut expr = self.or()?;
        while self.peek() == Some(b';') {
            self.pos += 1;
            expr = AtomQuery::And(Box::new(expr), Box::new(self.or()?));
        }
        Ok(expr)
    }

    fn or(&mut self) -> Result<AtomQuery, String> {
        let mut expr = self.high_and()?;
        while self.peek() == Some(b',') {
            self.pos += 1;
            expr = AtomQuery::Or(Box::new(expr), Box::new(self.high_and()?));
        }
        Ok(expr)
    }

    fn high_and(&mut self) -> Result<AtomQuery, String> {
        let mut expr = self.primitive()?;
        while !matches!(self.peek(), None | Some(b']' | b';' | b',')) {
            if self.peek() == Some(b'&') {
                self.pos += 1;
            }
            expr = AtomQuery::And(Box::new(expr), Box::new(self.primitive()?));
        }
        Ok(expr)
    }

    fn primitive(&mut self) -> Result<AtomQuery, String> {
        // Predicate letters also begin element symbols (He, Xe, Dy, Rn).
        if let (Some(first), Some(&second)) = (self.peek(), self.input.get(self.pos + 1)) {
            if first.is_ascii_uppercase()
                && second.is_ascii_lowercase()
                && crate::element_from_symbol(&format!("{}{}", first as char, second as char))
                    .is_some()
            {
                return self.element();
            }
        }
        match self.peek() {
            Some(b'!') => {
                self.pos += 1;
                self.enter()?;
                let query = self.primitive()?;
                self.depth -= 1;
                Ok(AtomQuery::Not(Box::new(query)))
            }
            Some(b'$') => {
                self.pos += 1;
                if self.peek() != Some(b'(') {
                    return Err(self.error("recursive SMARTS must use '$(...)'"));
                }
                self.pos += 1;
                self.enter()?;
                let query = self.graph(true)?;
                if self.peek() != Some(b')') {
                    return Err(self.error("unclosed recursive SMARTS"));
                }
                self.pos += 1;
                self.depth -= 1;
                Ok(AtomQuery::Recursive(Box::new(query)))
            }
            Some(b'#') => {
                self.pos += 1;
                let number = self.number()?;
                if Element::from_atomic_number(number).is_none() {
                    return Err(self.error("atomic number must be between 1 and 118"));
                }
                Ok(AtomQuery::AtomicNumber(number))
            }
            Some(b'H') if self.input.get(self.pos.wrapping_sub(1)) == Some(&b'[')
                && !self.input.get(self.pos + 1).is_some_and(u8::is_ascii_digit) => self.element(),
            Some(b'H' | b'X' | b'D') => {
                let kind = self.peek().unwrap();
                self.pos += 1;
                let number = if self.peek().is_some_and(|ch| ch.is_ascii_digit()) {
                    self.number()?
                } else if kind == b'H' {
                    1
                } else {
                    return Err(self.error("X and D predicates need an integer"));
                };
                Ok(match kind {
                    b'H' => AtomQuery::Hydrogens(number),
                    b'X' => AtomQuery::Connectivity(number),
                    _ => AtomQuery::Degree(number),
                })
            }
            Some(b'R') => {
                self.pos += 1;
                if self.peek() == Some(b'0') {
                    self.pos += 1;
                    Ok(AtomQuery::Ring(false))
                } else if self.peek().is_some_and(|ch| ch.is_ascii_digit()) {
                    Err(self.error("ring counts are unsupported; use R or R0"))
                } else {
                    Ok(AtomQuery::Ring(true))
                }
            }
            Some(b'+' | b'-') => {
                let sign = self.peek().unwrap();
                self.pos += 1;
                let value = if self.peek().is_some_and(|ch| ch.is_ascii_digit()) {
                    self.number()?
                } else {
                    let mut count = 1;
                    while self.peek() == Some(sign) {
                        count += 1;
                        self.pos += 1;
                    }
                    count
                };
                let signed = (value as i128) * if sign == b'+' { 1 } else { -1 };
                let charge = i8::try_from(signed).map_err(|_| self.error("charge is outside -128..127"))?;
                Ok(AtomQuery::Charge(charge))
            }
            Some(b'@' | b'h' | b'r' | b'v' | b'x' | b'^' | b':') => Err(self.error("unsupported atom predicate (stereo, isotope, implicit-H, ring size/count, valence, hybridization, and atom maps are not supported)")),
            Some(ch) if ch.is_ascii_digit() => Err(self.error("isotope predicates are not supported")),
            Some(b'[') => Err(self.error("nested bracket atoms are not valid SMARTS")),
            _ => self.atom(),
        }
    }
}

fn parse_pattern(pattern: &str) -> Result<Query, String> {
    if pattern.len() > MAX_PATTERN_BYTES {
        return Err("pattern exceeds 4096 bytes; use a smaller substructure".into());
    }
    Parser {
        input: pattern.as_bytes(),
        pos: 0,
        depth: 0,
    }
    .graph(false)
}

struct Matcher<'a> {
    molecule: &'a MoleculeGraph,
    visits: usize,
}

impl Matcher<'_> {
    fn visit(&mut self) -> Result<(), String> {
        self.visits += 1;
        if self.visits > MAX_VISITS {
            return Err("SMARTS search exceeded 1000000 steps; use a more specific pattern (no partial matches returned)".into());
        }
        Ok(())
    }

    fn atom_matches(&mut self, expr: &AtomQuery, index: usize) -> Result<bool, String> {
        self.visit()?;
        let atom = &self.molecule.atoms[index];
        // Opaque labels/wildcards have no known chemistry. Only '*' may select
        // them, including through logical expressions; don't infer from label text.
        Ok(match expr {
            AtomQuery::Any => true,
            AtomQuery::Element(symbol, aromatic) => {
                atom.symbol == *symbol && atom.aromatic == *aromatic
            }
            AtomQuery::AtomicNumber(number) => crate::element_from_symbol(&atom.symbol)
                .is_some_and(|element| element.get_atomic_number() == *number),
            AtomQuery::Aromatic(aromatic) => atom.symbol != "*" && atom.aromatic == *aromatic,
            AtomQuery::Hydrogens(count) => {
                atom.symbol != "*" && self.total_hydrogens(index) == *count
            }
            AtomQuery::Degree(count) => {
                atom.symbol != "*" && self.molecule.adj[index].len() == *count
            }
            AtomQuery::Connectivity(count) => {
                atom.symbol != "*"
                    && self.molecule.adj[index].len()
                        + crate::hydrogen_count(self.molecule, index) as usize
                        == *count
            }
            AtomQuery::Charge(charge) => atom.symbol != "*" && atom.charge == *charge,
            AtomQuery::Ring(in_ring) => self.in_ring(index) == *in_ring,
            AtomQuery::Not(inner) => !self.atom_matches(inner, index)?,
            AtomQuery::And(left, right) => {
                self.atom_matches(left, index)? && self.atom_matches(right, index)?
            }
            AtomQuery::Or(left, right) => {
                self.atom_matches(left, index)? || self.atom_matches(right, index)?
            }
            AtomQuery::Recursive(query) => {
                let mut mapping = vec![None; query.atoms.len()];
                self.search(query, &mut mapping, Some(index), &mut |_| Ok(true))?
            }
        })
    }

    fn total_hydrogens(&self, index: usize) -> usize {
        crate::hydrogen_count(self.molecule, index) as usize
            + self.molecule.adj[index]
                .iter()
                .filter(|&&(neighbor, _)| self.molecule.atoms[neighbor].symbol == "H")
                .count()
    }

    fn in_ring(&self, index: usize) -> bool {
        // An atom is in a cycle when two of its neighbors remain connected
        // without passing through the atom. This avoids choosing an SSSR basis.
        for &(start, _) in &self.molecule.adj[index] {
            let mut seen = vec![false; self.molecule.n_atoms()];
            seen[index] = true;
            seen[start] = true;
            let mut stack = vec![start];
            while let Some(next) = stack.pop() {
                for &(neighbor, _) in &self.molecule.adj[next] {
                    if !seen[neighbor] {
                        if self.molecule.adj[index]
                            .iter()
                            .any(|&(end, _)| end == neighbor && end != start)
                        {
                            return true;
                        }
                        seen[neighbor] = true;
                        stack.push(neighbor);
                    }
                }
            }
        }
        false
    }

    fn bond_matches(&self, kind: BondQuery, index: usize) -> bool {
        let bond = &self.molecule.bonds[index];
        match kind {
            BondQuery::Any => true,
            BondQuery::Aromatic => bond.aromatic,
            BondQuery::Implicit => bond.aromatic || bond.order == BondOrder::Single,
            BondQuery::Single => !bond.aromatic && bond.order == BondOrder::Single,
            BondQuery::Double => !bond.aromatic && bond.order == BondOrder::Double,
            BondQuery::Triple => !bond.aromatic && bond.order == BondOrder::Triple,
        }
    }

    // The callback returns true to stop (for anchored recursive predicates).
    fn search(
        &mut self,
        query: &Query,
        mapping: &mut [Option<usize>],
        anchor: Option<usize>,
        on_match: &mut impl FnMut(&[Option<usize>]) -> Result<bool, String>,
    ) -> Result<bool, String> {
        self.visit()?;
        let Some(next) = mapping.iter().position(Option::is_none) else {
            return on_match(mapping);
        };
        // Grow along a mapped query edge whenever possible; disconnected query
        // components fall back to all target atoms.
        let connected = query.bonds.iter().find_map(|&(a, b, kind)| {
            if a == next {
                mapping[b].map(|atom| (atom, kind))
            } else if b == next {
                mapping[a].map(|atom| (atom, kind))
            } else {
                None
            }
        });
        let candidates: Vec<_> = if let Some(anchor) = anchor.filter(|_| next == 0) {
            vec![anchor]
        } else if let Some((atom, kind)) = connected {
            self.molecule.adj[atom]
                .iter()
                .filter_map(|&(neighbor, bond)| self.bond_matches(kind, bond).then_some(neighbor))
                .collect()
        } else {
            (0..self.molecule.n_atoms()).collect()
        };
        for candidate in candidates {
            self.visit()?;
            if mapping.contains(&Some(candidate))
                || !self.atom_matches(&query.atoms[next], candidate)?
            {
                continue;
            }
            let edges_match = query.bonds.iter().all(|&(a, b, kind)| {
                let other = if a == next {
                    mapping[b]
                } else if b == next {
                    mapping[a]
                } else {
                    None
                };
                other.is_none_or(|other| {
                    self.molecule.adj[candidate]
                        .iter()
                        .any(|&(neighbor, bond)| neighbor == other && self.bond_matches(kind, bond))
                })
            });
            if edges_match {
                mapping[next] = Some(candidate);
                if self.search(query, mapping, anchor, on_match)? {
                    return Ok(true);
                }
                mapping[next] = None;
            }
        }
        Ok(false)
    }
}

pub(crate) fn find_matches(
    molecule: &MoleculeGraph,
    pattern: &str,
) -> Result<Vec<SubstructureMatch>, String> {
    let query = parse_pattern(pattern)
        .map_err(|error| format!("typed-smiles: invalid SMARTS {pattern:?}: {error}"))?;
    let mut matches = BTreeSet::new();
    let mut matcher = Matcher {
        molecule,
        visits: 0,
    };
    matcher.search(&query, &mut vec![None; query.atoms.len()], None, &mut |mapping| {
        let mut atoms: Vec<_> = mapping.iter().map(|atom| atom.unwrap()).collect();
        atoms.sort_unstable();
        let mut bonds: Vec<_> = query.bonds.iter().map(|&(a, b, _)| {
            let (a, b) = (mapping[a].unwrap(), mapping[b].unwrap());
            (a.min(b), a.max(b))
        }).collect();
        bonds.sort_unstable();
        matches.insert(SubstructureMatch { atoms, bonds });
        if matches.len() > MAX_MATCHES {
            return Err("SMARTS search exceeds 4096 distinct matches; use a more specific pattern (no partial matches returned)".into());
        }
        Ok(false)
    }).map_err(|error| format!("typed-smiles: SMARTS {pattern:?}: {error}"))?;
    Ok(matches.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matches(smiles: &str, pattern: &str) -> Vec<SubstructureMatch> {
        find_matches(&crate::parse_molecule(smiles).unwrap(), pattern).unwrap()
    }

    #[test]
    fn aspirin_acid() {
        assert_eq!(
            matches("CC(=O)OC1=CC=CC=C1C(=O)O", "C(=O)[OX2H1]"),
            vec![SubstructureMatch {
                atoms: vec![10, 11, 12],
                bonds: vec![(10, 11), (10, 12)]
            }]
        );
    }

    #[test]
    fn all_occurrences_and_symmetry() {
        assert_eq!(
            matches("OCCO", "[OX2H1]")
                .iter()
                .map(|m| &m.atoms)
                .collect::<Vec<_>>(),
            vec![&vec![0], &vec![3]]
        );
        assert_eq!(matches("CC", "CC").len(), 1);
        assert_eq!(matches("c1ccccc1", "c1ccccc1").len(), 1);
        assert_eq!(matches("CCC", "CC").len(), 2);
    }

    #[test]
    fn aromatic_flags_survive_kekulization() {
        assert_eq!(matches("c1ccccc1", "c:c").len(), 6);
        assert!(matches("c1ccccc1", "c=c").is_empty());
        assert!(matches("c1ccccc1", "C").is_empty());
        assert!(matches("C1=CC=CC=C1", "c").is_empty());
        assert_eq!(matches("c1ccccc1-c2ccccc2", "c-c").len(), 1);
        assert_eq!(matches("c1cc[nH]c1", "[nH1]").len(), 1);
    }

    #[test]
    fn hydrogens_charges_and_connectivity() {
        assert_eq!(matches("CCO", "[OX2H1]").len(), 1);
        assert_eq!(matches("CC[OH]", "[OX2H1]").len(), 1);
        assert_eq!(matches("CCO[H]", "[OX2H1]").len(), 1);
        assert!(matches("CC[O-]", "[OX2H1]").is_empty());
        assert_eq!(matches("C[NH3+]", "[N+;X4;H3]").len(), 1);
        assert_eq!(matches("[NH4+]", "[N+;X4;D0]").len(), 1);
        assert_eq!(matches("[O-]", "[O-1;H0]").len(), 1);
        assert_eq!(matches("[2H]O", "[OX2H2]").len(), 1);
        assert_eq!(matches("[H]", "[H]").len(), 1);
        assert!(matches("O", "[H]").is_empty());
    }

    #[test]
    fn query_edges_only_and_distinct_edge_sets() {
        let found = matches("C1CC1", "CCC");
        assert_eq!(found.len(), 3);
        assert!(found.iter().all(|m| m.bonds.len() == 2));
        assert_eq!(matches("C1CC1", "C1CC1")[0].bonds.len(), 3);
    }

    #[test]
    fn logical_precedence_recursive_context_and_atomic_numbers() {
        assert_eq!(matches("CO.CN.C[O-]", "[O,N;+0]").len(), 2);
        assert_eq!(matches("CO.CN", "[#8,#7]").len(), 2);
        assert_eq!(
            matches("CC(=O)N.CN", "[NX3;+0;!$(N-C=O)]")[0].atoms,
            vec![5]
        );
        assert_eq!(matches("CO.CC(=O)O", "[OX2H1;$([O]-[CX4])]").len(), 1);
        assert_eq!(matches("CO", "[!#6]")[0].atoms, vec![1]);
        assert_eq!(matches("CC", "[#6&#6]").len(), 2);
    }

    #[test]
    fn rings_fragments_wildcards_and_extensions() {
        assert_eq!(matches("CC1CC1", "[R]").len(), 3);
        assert_eq!(matches("CC1CC1", "[R0]").len(), 1);
        assert_eq!(matches("C.C", "C.C").len(), 1);
        assert_eq!(matches("C%12CC%12", "C%12CC%12").len(), 1);
        assert_eq!(matches("C!wC{OH}", "C-C").len(), 1);
        assert!(matches("C{OH}", "O").is_empty());
        assert_eq!(matches("C{OH}", "*").len(), 2);
    }

    #[test]
    fn invalid_and_unsupported_patterns_are_diagnostics() {
        for pattern in [
            "",
            "C(",
            "C()",
            "C1CC",
            "C11",
            "C1C1",
            "C=",
            "C..O",
            "C.",
            ".C",
            "[",
            "[]",
            "[[O]]",
            "[O[H]]",
            "[O;]",
            "[O,]",
            "[O&]",
            "C)",
            "C==O",
            "[XeQ]",
            "[R2]",
            "[r6]",
            "[C@H]",
            "[13C]",
            "[h1]",
            "[C:1]",
            "C/C",
            "C$C",
            "[$()]",
            "[H999999999999999999999999]",
            "[C+999999999999999999999]",
        ] {
            let error =
                find_matches(&crate::parse_molecule("CCO").unwrap(), pattern).expect_err(pattern);
            assert!(error.contains("invalid SMARTS"), "{pattern}: {error}");
        }
    }

    #[test]
    fn limits_do_not_return_partial_results() {
        let molecule = crate::parse_molecule("CCO").unwrap();
        let mut matcher = Matcher {
            molecule: &molecule,
            visits: MAX_VISITS,
        };
        assert!(matcher.visit().unwrap_err().contains("no partial matches"));
        assert!(parse_pattern(&"C".repeat(MAX_QUERY_ATOMS + 1)).is_err());
        assert!(parse_pattern(&"C".repeat(MAX_PATTERN_BYTES + 1)).is_err());
        assert!(parse_pattern(&format!("[{}C]", "!".repeat(MAX_DEPTH + 1))).is_err());
        let symmetric = crate::parse_molecule("C.C.C.C.C.C.C.C.C.C").unwrap();
        assert!(find_matches(&symmetric, "C.C.C.C.C.C.C.C.C.C")
            .unwrap_err()
            .contains("no partial matches"));
    }

    fn group(name: &str) -> &'static str {
        let registry = include_str!("../../src/substructure.typ");
        let line = registry
            .lines()
            .find(|line| line.trim_start().starts_with(&format!("{name}:")))
            .unwrap();
        line.split_once(':')
            .unwrap()
            .1
            .trim()
            .trim_end_matches(',')
            .trim_matches('"')
    }

    #[test]
    fn named_group_definitions_match_examples() {
        for (name, smiles, count) in [
            ("carboxylic-acid", "CC(=O)O", 1),
            ("carboxylate", "CC(=O)[O-]", 1),
            ("alcohol", "OCCO", 2),
            ("phenol", "Oc1ccccc1", 1),
            ("amine", "CN(C)C", 1),
            ("ester", "CC(=O)OCC", 1),
            ("amide", "CC(=O)N", 1),
            ("carbonyl", "CC(=O)OC1=CC=CC=C1C(=O)O", 2),
            ("aldehyde", "CC=O", 1),
            ("aldehyde", "C=O", 1),
            ("ketone", "CC(=O)C", 1),
            ("nitrile", "CC#N", 1),
            ("ether", "COC", 1),
            ("thiol", "CCS", 1),
            ("nitro", "C[N+](=O)[O-]", 1),
            ("alkene", "C=CC=C", 2),
            ("alkyne", "CC#CC", 1),
        ] {
            assert_eq!(
                matches(smiles, group(name)).len(),
                count,
                "{name}: {smiles}"
            );
        }
    }

    #[test]
    fn named_groups_distinguish_similar_chemistry() {
        for (name, smiles) in [
            ("alcohol", "CC(=O)O"),
            ("alcohol", "Oc1ccccc1"),
            ("alcohol", "CC[O-]"),
            ("phenol", "CCO"),
            ("carboxylic-acid", "CC(=O)[O-]"),
            ("carboxylic-acid", "CC(=O)OC"),
            ("amine", "CC(=O)N"),
            ("amine", "CS(=O)(=O)NC"),
            ("amine", "C[NH3+]"),
            ("amine", "N"),
            ("amine", "n1ccccc1"),
            ("ester", "CC(=O)O"),
            ("ester", "COC"),
            ("amide", "CN"),
            ("ketone", "CC=O"),
            ("ketone", "CC(=O)O"),
            ("aldehyde", "CC(=O)C"),
            ("aldehyde", "C(=O)O"),
            ("ether", "CCO"),
            ("ether", "CC(=O)OC"),
            ("nitrile", "CC(=O)N"),
            ("alkene", "c1ccccc1"),
        ] {
            assert!(matches(smiles, group(name)).is_empty(), "{name}: {smiles}");
        }
    }

    #[test]
    fn pubchem_corpus_agrees_with_independent_match_sets() {
        let corpus: serde_json::Value = serde_json::from_str(include_str!(
            "../../tests/fixtures/substructure-pubchem.json"
        ))
        .unwrap();
        let patterns = corpus["patterns"].as_object().unwrap();
        for case in corpus["cases"].as_array().unwrap() {
            for variant in case["variants"].as_array().unwrap() {
                let smiles = variant["smiles"].as_str().unwrap();
                for (name, pattern) in patterns {
                    assert_eq!(group(name), pattern.as_str().unwrap());
                    let found = matches(smiles, group(name));
                    let actual = serde_json::to_value(found).unwrap();
                    let expected = variant["matches"]
                        .get(name)
                        .cloned()
                        .unwrap_or_else(|| serde_json::json!([]));
                    assert_eq!(
                        actual, expected,
                        "{} / {} / {name}: {smiles}",
                        case["name"], variant["notation"]
                    );
                }
            }
        }
    }
}
