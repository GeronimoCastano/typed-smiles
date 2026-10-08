/// Reads typed-smiles input into a molecular graph.
///
/// The accepted language is OpenSMILES plus the package's `{label}`
/// abbreviations and `!` drawing extensions, and two relaxations that hand
/// written SMILES commonly use: ring-closure numbers written after an atom's
/// branches (`C1=CCCC(=O)1`) and branches holding only ring closures (`C(1)`).
/// Atom order is writing order, so `show-indices`, highlights, and mechanism
/// references count atoms the way the author wrote them.
///
/// The reader is a single iterative pass: branch depth and chain length are
/// limited by memory, never by the call stack. Every diagnostic names the
/// 1-based character position in the original input.
use crate::graph::{
    Atom, AtomChirality, BondDirection, BondOrder, BondSpec, BondStereo, GraphBuilder,
    MoleculeGraph, WrittenBond,
};
use crate::label::parse_abbreviation_label;

/// Largest formal charge magnitude accepted on one atom, the range OpenSMILES
/// requires every reader to support.
const MAXIMUM_CHARGE: i8 = 15;

pub(crate) fn parse_smiles(input: &str) -> Result<MoleculeGraph, String> {
    let written = read_written_molecule(input)?;
    let mut molecule = written.molecule;
    crate::kekulize::kekulize(&mut molecule, &written.implicit_bonds)
        .map_err(|error| format!("{error}{}", written.unbracketed_element_hints))?;
    Ok(molecule)
}

/// The graph exactly as written, before aromatic bonds receive Kekulé orders.
pub(crate) struct WrittenMolecule {
    pub molecule: MoleculeGraph,
    /// Per bond: whether it was written without a bond symbol.
    pub implicit_bonds: Vec<bool>,
    unbracketed_element_hints: String,
}

pub(crate) fn read_written_molecule(input: &str) -> Result<WrittenMolecule, String> {
    let mut reader = SmilesReader::new(input);
    reader.read_all_tokens()?;
    let unbracketed_element_hints = reader.unbracketed_element_hints();
    let (molecule, implicit_bonds) = reader.builder.finish()?;
    Ok(WrittenMolecule {
        molecule,
        implicit_bonds,
        unbracketed_element_hints,
    })
}

/// The kind of the most recent token, which decides what may follow it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PreviousToken {
    Start,
    Atom,
    RingBond,
    BranchOpen,
    BranchClose,
    Bond,
    Dot,
}

struct OpenBranch {
    branch_point: usize,
    position: usize,
    holds_only_ring_bonds: bool,
}

struct PendingBond {
    bond: WrittenBond,
    follows: PreviousToken,
}

struct SmilesReader {
    characters: Vec<char>,
    cursor: usize,
    builder: GraphBuilder,
    previous_token: PreviousToken,
    /// Atom that the next atom bonds to; `None` at the start of a fragment.
    previous_atom: Option<usize>,
    pending_bond: Option<PendingBond>,
    open_branches: Vec<OpenBranch>,
    dot_position: usize,
    /// Two-letter element spellings such as `Co` that also read as an
    /// organic atom followed by an aromatic one, kept to explain later errors.
    ambiguous_element_spellings: Vec<(String, usize)>,
}

impl SmilesReader {
    fn new(input: &str) -> Self {
        Self {
            characters: input.chars().collect(),
            cursor: 0,
            builder: GraphBuilder::default(),
            previous_token: PreviousToken::Start,
            previous_atom: None,
            pending_bond: None,
            open_branches: Vec::new(),
            dot_position: 0,
            ambiguous_element_spellings: Vec::new(),
        }
    }

    fn peek(&self) -> Option<char> {
        self.characters.get(self.cursor).copied()
    }

    fn peek_at(&self, offset: usize) -> Option<char> {
        self.characters.get(self.cursor + offset).copied()
    }

    /// 1-based character position of the next unread character.
    fn position(&self) -> usize {
        self.cursor + 1
    }

    fn read_all_tokens(&mut self) -> Result<(), String> {
        if self
            .characters
            .iter()
            .all(|character| character.is_whitespace())
        {
            return Err("the expression is empty; provide at least one atom".to_string());
        }
        self.reject_unmatched_closing_delimiters()?;
        while let Some(character) = self.peek() {
            let position = self.position();
            match character {
                '[' => {
                    let atom = self.read_bracket_atom()?;
                    self.place_atom(atom)?;
                }
                '{' => {
                    let atom = self.read_abbreviation_atom()?;
                    self.place_atom(atom)?;
                }
                '0'..='9' | '%' => {
                    let ring_number = self.read_ring_number()?;
                    self.place_ring_bond(ring_number, position)?;
                }
                '-' | '=' | '#' | '$' | ':' | '/' | '\\' => {
                    self.cursor += 1;
                    self.place_bond(standard_bond(character, position))?;
                }
                '!' => {
                    let bond = self.read_drawing_extension()?;
                    self.place_bond(bond)?;
                }
                '(' => {
                    self.cursor += 1;
                    self.open_branch(position)?;
                }
                ')' => {
                    self.cursor += 1;
                    self.close_branch(position)?;
                }
                '.' => {
                    self.cursor += 1;
                    self.place_dot(position)?;
                }
                _ => {
                    let atom = self.read_unbracketed_atom()?;
                    self.place_atom(atom)?;
                }
            }
        }
        self.finish_input()
    }

    /// A stray `]` or `}` usually means a mistyped bracket atom or label, which
    /// explains the input better than whatever token happens to come first.
    fn reject_unmatched_closing_delimiters(&self) -> Result<(), String> {
        let mut open_delimiter = None;
        for (index, &character) in self.characters.iter().enumerate() {
            match (open_delimiter, character) {
                (None, '[' | '{') => open_delimiter = Some(character),
                (Some('['), ']') | (Some('{'), '}') => open_delimiter = None,
                (None, ']' | '}') => {
                    return Err(self.unexpected_character_message(character, index + 1));
                }
                _ => {}
            }
        }
        Ok(())
    }

    // ── Structure ────────────────────────────────────────────────────────────

    fn place_atom(&mut self, atom: Atom) -> Result<(), String> {
        let position = atom.source_position;
        self.reject_inside_ring_only_branch("an atom", position)?;
        let pending_bond = self.pending_bond.take().map(|pending| pending.bond);
        let preceding_atom = self.previous_atom;
        self.previous_token = PreviousToken::Atom;

        if let Some(parent_atom) = preceding_atom {
            if self.should_fold_hydrogen(&atom, parent_atom, pending_bond.as_ref()) {
                self.builder.fold_hydrogen_into(parent_atom);
                return Ok(());
            }
        }

        let atom_index = self.builder.add_atom(atom);
        if let Some(preceding_atom) = preceding_atom {
            self.builder
                .bond_to_preceding_atom(preceding_atom, atom_index, pending_bond.as_ref());
        }
        self.previous_atom = Some(atom_index);
        Ok(())
    }

    /// A plain terminal hydrogen, such as `[H]` in `C([H])O` or `C[H]`, joins
    /// its neighbor's hydrogen count so it draws like an implicit hydrogen.
    fn should_fold_hydrogen(
        &self,
        atom: &Atom,
        parent_atom: usize,
        pending_bond: Option<&WrittenBond>,
    ) -> bool {
        let bonded_by_plain_single_bond =
            pending_bond.is_none_or(|bond| bond.specification == BondSpec::single());
        let is_terminal = matches!(self.peek(), None | Some(')'));
        atom.is_foldable_hydrogen()
            && bonded_by_plain_single_bond
            && is_terminal
            && self.builder.can_fold_hydrogen_into(parent_atom)
    }

    fn place_ring_bond(&mut self, ring_number: u8, position: usize) -> Result<(), String> {
        let context = match &self.pending_bond {
            Some(pending) => pending.follows,
            None => self.previous_token,
        };
        let atom_index = match (context, self.previous_atom) {
            (
                PreviousToken::Atom | PreviousToken::RingBond | PreviousToken::BranchClose,
                Some(atom_index),
            ) => atom_index,
            (PreviousToken::BranchOpen, _) => {
                let open_branch = self
                    .open_branches
                    .last_mut()
                    .expect("a branch is open after `(`");
                open_branch.holds_only_ring_bonds = true;
                open_branch.branch_point
            }
            _ => {
                return Err(format!(
                    "ring-closure number {} at character {position} must follow an atom",
                    crate::graph::ring_label(ring_number)
                ));
            }
        };
        let bond = self.pending_bond.take().map(|pending| pending.bond);
        self.builder
            .ring_bond(atom_index, ring_number, bond, position)?;
        self.previous_token = PreviousToken::RingBond;
        Ok(())
    }

    fn place_bond(&mut self, bond: WrittenBond) -> Result<(), String> {
        match self.previous_token {
            PreviousToken::Atom
            | PreviousToken::RingBond
            | PreviousToken::BranchOpen
            | PreviousToken::BranchClose => {
                self.pending_bond = Some(PendingBond {
                    follows: self.previous_token,
                    bond,
                });
                self.previous_token = PreviousToken::Bond;
                Ok(())
            }
            PreviousToken::Bond => {
                let previous = &self
                    .pending_bond
                    .as_ref()
                    .expect("a bond is pending after a bond token")
                    .bond;
                Err(format!(
                    "bond `{}` at character {} directly follows bond `{}` at character {}; \
                     write one bond symbol between two atoms",
                    bond.text, bond.position, previous.text, previous.position
                ))
            }
            PreviousToken::Start | PreviousToken::Dot => Err(format!(
                "bond `{}` at character {} must follow an atom",
                bond.text, bond.position
            )),
        }
    }

    fn open_branch(&mut self, position: usize) -> Result<(), String> {
        self.reject_inside_ring_only_branch("a branch", position)?;
        match self.previous_token {
            PreviousToken::Atom | PreviousToken::RingBond | PreviousToken::BranchClose => {}
            PreviousToken::Bond => {
                let bond = &self
                    .pending_bond
                    .as_ref()
                    .expect("a bond is pending after a bond token")
                    .bond;
                return Err(format!(
                    "branch at character {position} follows bond `{}` at character {}; write \
                     the bond inside the branch, as in `C(=O)C`",
                    bond.text, bond.position
                ));
            }
            PreviousToken::BranchOpen => {
                return Err(format!(
                    "branch at character {position} opens directly inside another branch; \
                     start each branch with an atom, as in `C(C(C)C)`"
                ));
            }
            PreviousToken::Start | PreviousToken::Dot => {
                return Err(format!(
                    "branch at character {position} must follow the atom it branches from"
                ));
            }
        }
        let branch_point = self
            .previous_atom
            .expect("an atom precedes every branch opening");
        self.open_branches.push(OpenBranch {
            branch_point,
            position,
            holds_only_ring_bonds: false,
        });
        self.previous_token = PreviousToken::BranchOpen;
        Ok(())
    }

    fn close_branch(&mut self, position: usize) -> Result<(), String> {
        if self.open_branches.is_empty() {
            return Err(format!(
                "unmatched `)` at character {position}; every `)` must close a branch opened \
                 with `(`"
            ));
        }
        self.require_completed_token(&format!("`)` at character {position}"))?;
        if self.previous_token == PreviousToken::BranchOpen {
            return Err(format!(
                "empty branch `()` at character {}; remove it or write an atom inside it",
                position - 1
            ));
        }
        let open_branch = self
            .open_branches
            .pop()
            .expect("an open branch was checked above");
        self.previous_atom = Some(open_branch.branch_point);
        self.previous_token = PreviousToken::BranchClose;
        Ok(())
    }

    fn place_dot(&mut self, position: usize) -> Result<(), String> {
        self.reject_inside_ring_only_branch("a `.`", position)?;
        match self.previous_token {
            PreviousToken::Atom
            | PreviousToken::RingBond
            | PreviousToken::BranchOpen
            | PreviousToken::BranchClose => {}
            PreviousToken::Bond => {
                let bond = &self
                    .pending_bond
                    .as_ref()
                    .expect("a bond is pending after a bond token")
                    .bond;
                return Err(format!(
                    "bond `{}` at character {} is followed by `.` instead of an atom; remove \
                     the bond symbol or the `.`",
                    bond.text, bond.position
                ));
            }
            PreviousToken::Start | PreviousToken::Dot => {
                return Err(format!(
                    "`.` at character {position} must separate two fragments, as in `[Na+].[Cl-]`"
                ));
            }
        }
        self.previous_atom = None;
        self.previous_token = PreviousToken::Dot;
        self.dot_position = position;
        Ok(())
    }

    fn finish_input(&mut self) -> Result<(), String> {
        self.require_completed_token("the end of the SMILES")?;
        if let Some(open_branch) = self.open_branches.last() {
            return Err(format!(
                "the branch opened at character {} is never closed; add `)` after the branch",
                open_branch.position
            ));
        }
        Ok(())
    }

    /// Rejects a bond or `.` that `boundary` (a `)` or the end of the input)
    /// follows directly, since both need an atom after them.
    fn require_completed_token(&self, boundary: &str) -> Result<(), String> {
        match self.previous_token {
            PreviousToken::Bond => {
                let bond = &self
                    .pending_bond
                    .as_ref()
                    .expect("a bond is pending after a bond token")
                    .bond;
                Err(format!(
                    "bond `{}` at character {} is followed by {boundary} instead of an atom; \
                     add the bonded atom or remove the bond symbol",
                    bond.text, bond.position
                ))
            }
            PreviousToken::Dot => Err(format!(
                "`.` at character {} is followed by {boundary} instead of an atom; remove the \
                 `.` or add the next fragment",
                self.dot_position
            )),
            _ => Ok(()),
        }
    }

    fn reject_inside_ring_only_branch(&self, token: &str, position: usize) -> Result<(), String> {
        match self.open_branches.last() {
            Some(open_branch) if open_branch.holds_only_ring_bonds => Err(format!(
                "{token} at character {position} appears in the branch opened at character {}, \
                 which started with a ring-closure number; such a branch may hold only \
                 ring-closure numbers, so write atoms in a separate branch",
                open_branch.position
            )),
            _ => Ok(()),
        }
    }

    // ── Atoms ────────────────────────────────────────────────────────────────

    fn read_unbracketed_atom(&mut self) -> Result<Atom, String> {
        let position = self.position();
        let character = self.peek().expect("a character is available");
        let next_character = self.peek_at(1);

        if let Some(symbol) = two_letter_organic_symbol(character, next_character) {
            self.cursor += 2;
            return Ok(Atom::new(symbol, position));
        }
        if matches!(character, 'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I') {
            self.check_unbracketed_two_letter_element(character, next_character, position)?;
            self.cursor += 1;
            return Ok(Atom::new(&character.to_string(), position));
        }
        if matches!(character, 'b' | 'c' | 'n' | 'o' | 'p' | 's') {
            self.cursor += 1;
            let mut atom = Atom::new(&character.to_ascii_uppercase().to_string(), position);
            atom.aromatic = true;
            return Ok(atom);
        }
        if character == '*' {
            self.cursor += 1;
            return Ok(Atom::new("*", position));
        }
        Err(self.unexpected_character_message(character, position))
    }

    /// Elements outside the organic subset must be bracketed. A spelling like
    /// `Na` is rejected with that correction; `Co` is valid SMILES for carbon
    /// bonded to aromatic oxygen, so it is only remembered for later hints.
    fn check_unbracketed_two_letter_element(
        &mut self,
        character: char,
        next_character: Option<char>,
        position: usize,
    ) -> Result<(), String> {
        let Some(next_character) = next_character.filter(char::is_ascii_lowercase) else {
            return Ok(());
        };
        let symbol = format!("{character}{next_character}");
        if crate::element_from_symbol(&symbol).is_none() {
            return Ok(());
        }
        if matches!(next_character, 'b' | 'c' | 'n' | 'o' | 'p' | 's') {
            self.ambiguous_element_spellings.push((symbol, position));
            return Ok(());
        }
        Err(format!(
            "element `{symbol}` at character {position} is outside the organic subset; write \
             it in brackets, as `[{symbol}]`"
        ))
    }

    fn unexpected_character_message(&self, character: char, position: usize) -> String {
        match character {
            'H' => format!(
                "`H` at character {position} must be written in brackets; hydrogens on \
                 organic atoms are implicit (`C` is methane), or write them in the bracket \
                 atom as `[CH3]` or `[H]`"
            ),
            '@' => format!(
                "chirality `@` at character {position} must be written inside a bracket atom, \
                 as in `[C@H]` or `[C@@](F)(Cl)Br`"
            ),
            '+' => format!(
                "charge `+` at character {position} must be written inside a bracket atom, as \
                 in `[NH4+]`"
            ),
            ']' => {
                format!("unmatched `]` at character {position}; bracket atoms must start with `[`")
            }
            '}' => format!(
                "unmatched `}}` at character {position}; custom labels must start with `{{`"
            ),
            '>' | '<' if self.is_dative_bond_arrow(character, position) => format!(
                "dative bond arrow at character {position} is not OpenSMILES; write the \
                 coordination bond as a single bond `-`"
            ),
            '>' => format!(
                "`>` at character {position} is only valid inside an abbreviation label like \
                 `{{>PPh3}}`; reaction SMILES such as `A>>B` are not drawn by smiles(), so \
                 draw each species with mol() inside reaction()"
            ),
            character if character.is_whitespace() => format!(
                "SMILES cannot contain whitespace, but character {position} is {character:?}; \
                 remove it"
            ),
            character if character.is_ascii_uppercase() => {
                let symbol = match self.peek_at(1).filter(char::is_ascii_lowercase) {
                    Some(lowercase)
                        if crate::element_from_symbol(&format!("{character}{lowercase}"))
                            .is_some() =>
                    {
                        format!("{character}{lowercase}")
                    }
                    _ => character.to_string(),
                };
                if crate::element_from_symbol(&symbol).is_some() {
                    format!(
                        "element `{symbol}` at character {position} is outside the organic \
                         subset; write it in brackets, as `[{symbol}]`"
                    )
                } else {
                    format!("`{symbol}` at character {position} is not an element symbol")
                }
            }
            character => format!(
                "unexpected `{character}` at character {position}; expected an atom, bond, \
                 ring-closure number, branch, or `.`"
            ),
        }
    }

    /// Whether `character` at 1-based `position` belongs to `->` or `<-`.
    fn is_dative_bond_arrow(&self, character: char, position: usize) -> bool {
        let index = position - 1;
        match character {
            '>' => index > 0 && self.characters[index - 1] == '-',
            _ => self.characters.get(index + 1) == Some(&'-'),
        }
    }

    fn read_abbreviation_atom(&mut self) -> Result<Atom, String> {
        let position = self.position();
        self.cursor += 1;
        let mut body = String::new();
        loop {
            match self.peek() {
                Some('}') => break,
                Some(character) => {
                    body.push(character);
                    self.cursor += 1;
                }
                None => {
                    return Err(format!(
                        "unclosed custom label starting at character {position}; add `}}` \
                         after the label"
                    ));
                }
            }
        }
        self.cursor += 1;
        let label = parse_abbreviation_label(&body)
            .map_err(|error| format!("{error} (label at character {position})"))?;

        // A label stands for an unspecified group: it behaves like `[*]`.
        let mut atom = Atom::new("*", position);
        atom.has_explicit_h = true;
        atom.abbrev = label.text;
        atom.abbrev_style = label.style;
        atom.abbrev_anchor = label.anchor;
        atom.abbrev_anchor_len = label.anchor_len;
        atom.abbrev_lone_pairs = label.lone_pairs.unwrap_or(0);
        let (offset_x, offset_y) = label.offset.unwrap_or((0.0, 0.0));
        atom.abbrev_offset_x = offset_x;
        atom.abbrev_offset_y = offset_y;
        Ok(atom)
    }

    // ── Bracket atoms ────────────────────────────────────────────────────────

    /// Reads `[isotope? symbol chirality? hydrogens? charge? class?]`.
    fn read_bracket_atom(&mut self) -> Result<Atom, String> {
        let start = self.position();
        let closing = self.characters[self.cursor..]
            .iter()
            .position(|&character| character == ']')
            .map(|offset| self.cursor + offset);
        let Some(closing) = closing else {
            return Err(format!(
                "unclosed bracket atom starting at character {start}; add `]` after the atom \
                 specification"
            ));
        };
        let text: String = self.characters[self.cursor..=closing].iter().collect();
        self.cursor += 1;

        let isotope = self.read_isotope(&text)?;
        let (symbol, aromatic) = self.read_bracket_symbol(&text)?;
        let chirality = self.read_chirality(&text)?;
        let hydrogen_count = self.read_hydrogen_count();
        let charge = self.read_charge(&text)?;
        let atom_map = self.read_atom_map(&text)?;
        if self.cursor != closing {
            return Err(self.bracket_order_message(&text));
        }
        self.cursor += 1;

        let mut atom = Atom::new(&symbol, start);
        atom.aromatic = aromatic;
        atom.isotope = isotope;
        atom.chirality = chirality;
        atom.hcount = hydrogen_count;
        atom.has_explicit_h = true;
        atom.charge = charge;
        atom.atom_map = atom_map;
        Ok(atom)
    }

    fn read_digits(&mut self) -> String {
        let mut digits = String::new();
        while let Some(digit) = self.peek().filter(char::is_ascii_digit) {
            digits.push(digit);
            self.cursor += 1;
        }
        digits
    }

    fn read_isotope(&mut self, bracket_text: &str) -> Result<Option<u16>, String> {
        let position = self.position();
        let digits = self.read_digits();
        if digits.is_empty() {
            return Ok(None);
        }
        if digits.len() > 3 {
            return Err(format!(
                "isotope `{digits}` at character {position} in `{bracket_text}` has more than \
                 three digits; write a mass number such as `13` in `[13C]`"
            ));
        }
        Ok(Some(digits.parse().expect("three digits fit in u16")))
    }

    fn read_bracket_symbol(&mut self, bracket_text: &str) -> Result<(String, bool), String> {
        let position = self.position();
        let Some(first) = self.peek() else {
            unreachable!("the bracket text ends with `]`");
        };
        let second = self.peek_at(1);

        if first == '*' {
            self.cursor += 1;
            return Ok(("*".to_string(), false));
        }
        if first.is_ascii_lowercase() {
            let two_letters = second.map(|second| format!("{first}{second}"));
            if let Some(symbol) = two_letters.filter(|symbol| symbol == "se" || symbol == "as") {
                self.cursor += 2;
                return Ok((uppercase_first(&symbol), true));
            }
            if matches!(first, 'b' | 'c' | 'n' | 'o' | 'p' | 's') {
                self.cursor += 1;
                return Ok((first.to_ascii_uppercase().to_string(), true));
            }
            return Err(format!(
                "`{first}` at character {position} in `{bracket_text}` is not an element \
                 symbol; OpenSMILES aromatic symbols are b, c, n, o, p, s, se, and as"
            ));
        }
        if first.is_ascii_uppercase() {
            if let Some(second) = second.filter(char::is_ascii_lowercase) {
                let symbol = format!("{first}{second}");
                if crate::element_from_symbol(&symbol).is_some() {
                    self.cursor += 2;
                    return Ok((symbol, false));
                }
            }
            let symbol = first.to_string();
            if crate::element_from_symbol(&symbol).is_some() {
                self.cursor += 1;
                return Ok((symbol, false));
            }
            let written = match second.filter(char::is_ascii_lowercase) {
                Some(second) => format!("{first}{second}"),
                None => symbol,
            };
            return Err(format!(
                "`{written}` at character {position} in `{bracket_text}` is not an element \
                 symbol"
            ));
        }
        if first == ']' {
            return Err(format!(
                "empty bracket atom `[]` at character {}; write an element symbol or `*` \
                 inside the brackets",
                position - 1
            ));
        }
        Err(format!(
            "`{first}` at character {position} in `{bracket_text}` is not an element symbol; \
             bracket atoms are written as [isotope symbol chirality hydrogens charge class], \
             e.g. `[13CH3+:1]`"
        ))
    }

    fn read_chirality(&mut self, bracket_text: &str) -> Result<AtomChirality, String> {
        if self.peek() != Some('@') {
            return Ok(AtomChirality::None);
        }
        let position = self.position();
        self.cursor += 1;
        if self.peek() == Some('@') {
            self.cursor += 1;
            return Ok(AtomChirality::TetraClockwise);
        }

        let class_name: String = [self.peek(), self.peek_at(1)].iter().flatten().collect();
        let Some((maximum, description)) = chirality_class_range(&class_name) else {
            return Ok(AtomChirality::TetraAnti);
        };
        self.cursor += 2;
        let digits = self.read_digits();
        let class_number = digits
            .parse::<u8>()
            .ok()
            .filter(|class| (1..=maximum).contains(class));
        let Some(class_number) = class_number else {
            return Err(format!(
                "chirality `@{class_name}{digits}` at character {position} in `{bracket_text}` \
                 is not a valid {description} class; use @{class_name}1 through \
                 @{class_name}{maximum}"
            ));
        };
        Ok(match class_name.as_str() {
            "TH" if class_number == 1 => AtomChirality::TetraAnti,
            "TH" => AtomChirality::TetraClockwise,
            "AL" => AtomChirality::Allenal(class_number),
            "SP" => AtomChirality::SquarePlanar(class_number),
            "TB" => AtomChirality::TrigonalBipyramidal(class_number),
            _ => AtomChirality::Octahedral(class_number),
        })
    }

    fn read_hydrogen_count(&mut self) -> u8 {
        if self.peek() != Some('H') {
            return 0;
        }
        self.cursor += 1;
        match self.peek().and_then(|digit| digit.to_digit(10)) {
            Some(count) => {
                self.cursor += 1;
                count as u8
            }
            None => 1,
        }
    }

    /// Reads `+`, `-`, `+2`, `-15`, or the deprecated `++` and `--`.
    fn read_charge(&mut self, bracket_text: &str) -> Result<i8, String> {
        let Some(sign) = self
            .peek()
            .filter(|&character| character == '+' || character == '-')
        else {
            return Ok(0);
        };
        let position = self.position();
        self.cursor += 1;
        let sign_value: i8 = if sign == '+' { 1 } else { -1 };

        if self.peek() == Some(sign) {
            self.cursor += 1;
            if matches!(self.peek(), Some('+' | '-' | '0'..='9')) {
                return Err(invalid_charge_message(bracket_text, position));
            }
            return Ok(2 * sign_value);
        }
        let digits = self.read_digits();
        if matches!(self.peek(), Some('+' | '-')) || digits.len() > 2 {
            return Err(invalid_charge_message(bracket_text, position));
        }
        let magnitude: i8 = if digits.is_empty() {
            1
        } else {
            digits.parse().expect("two digits fit in i8")
        };
        if magnitude > MAXIMUM_CHARGE {
            return Err(format!(
                "charge `{sign}{digits}` at character {position} in `{bracket_text}` is out of \
                 range; charges between -{MAXIMUM_CHARGE} and +{MAXIMUM_CHARGE} are supported"
            ));
        }
        Ok(sign_value * magnitude)
    }

    fn read_atom_map(&mut self, bracket_text: &str) -> Result<u32, String> {
        if self.peek() != Some(':') {
            return Ok(0);
        }
        let position = self.position();
        self.cursor += 1;
        let digits = self.read_digits();
        if digits.is_empty() {
            return Err(format!(
                "atom map `:` at character {position} in `{bracket_text}` needs a number, as \
                 in `[CH3:7]`"
            ));
        }
        digits.parse().map_err(|_| {
            format!(
                "atom map `:{digits}` at character {position} in `{bracket_text}` is too \
                 large; use a number up to {}",
                u32::MAX
            )
        })
    }

    fn bracket_order_message(&self, bracket_text: &str) -> String {
        let position = self.position();
        let character = self.peek().expect("the bracket text ends with `]`");
        let hint = match character {
            'H' => "; write the hydrogen count before the charge and atom map, as in `[NH4+:1]`",
            '@' => "; write the chirality right after the element symbol, as in `[C@@H]`",
            '+' | '-' => "; write a single charge before the atom map, as in `[NH4+:1]`",
            ':' => "; write one atom map at the end, as in `[CH3:7]`",
            _ => "",
        };
        format!(
            "unexpected `{character}` at character {position} in bracket atom `{bracket_text}`; \
             bracket atoms are written as [isotope symbol chirality hydrogens charge class], \
             e.g. `[13CH3+:1]`{hint}"
        )
    }

    // ── Bonds and ring closures ─────────────────────────────────────────────

    fn read_ring_number(&mut self) -> Result<u8, String> {
        let position = self.position();
        let character = self.peek().expect("a ring-closure character is available");
        if character != '%' {
            self.cursor += 1;
            return Ok(character.to_digit(10).expect("ring digit") as u8);
        }
        let digits: String = [self.peek_at(1), self.peek_at(2)]
            .iter()
            .flatten()
            .take_while(|digit| digit.is_ascii_digit())
            .collect();
        if digits.len() != 2 {
            return Err(format!(
                "`%` at character {position} must be followed by two digits, as in `%10`; \
                 ring numbers 0 to 9 are written without `%`"
            ));
        }
        self.cursor += 3;
        Ok(digits.parse().expect("two digits fit in u8"))
    }

    /// Reads a typed-smiles drawing extension: `!w`, `!h`, `!s`, `!d`, or a
    /// curl `!c` optionally combined with one of them, then an optional bond
    /// order symbol.
    fn read_drawing_extension(&mut self) -> Result<WrittenBond, String> {
        let position = self.position();
        let start = self.cursor;
        self.cursor += 1;
        let extension = self
            .peek()
            .ok_or_else(|| format!("incomplete `!` drawing extension at character {position}"))?;
        self.cursor += 1;

        let (stereo, curl) = if extension == 'c' {
            let stereo = if self.peek() == Some('!') {
                let combined = self.peek_at(1).ok_or_else(|| {
                    format!("incomplete drawing extension after `!c` at character {position}")
                })?;
                let stereo = drawing_extension_stereo(combined).ok_or_else(|| {
                    format!(
                        "unknown drawing extension `!{combined}` after `!c` at character \
                         {position}; combine `!c` with !w, !h, !s, or !d"
                    )
                })?;
                self.cursor += 2;
                stereo
            } else {
                BondStereo::None
            };
            (stereo, true)
        } else {
            let stereo = drawing_extension_stereo(extension).ok_or_else(|| {
                format!(
                    "unknown drawing extension `!{extension}` at character {position}; use !w, \
                     !h, !s, !d, or !c"
                )
            })?;
            (stereo, false)
        };

        let order = match self.peek().and_then(bond_order_symbol) {
            Some(order) => {
                self.cursor += 1;
                order
            }
            None => BondOrder::Single,
        };
        if stereo != BondStereo::None && order != BondOrder::Single {
            return Err(format!(
                "wedge, hash, wavy, and dashed drawing extensions require a single bond \
                 (character {position})"
            ));
        }
        Ok(WrittenBond {
            specification: BondSpec {
                order,
                stereo,
                direction: BondDirection::None,
                forced_stereo: stereo.is_wedge(),
                curl,
            },
            text: self.characters[start..self.cursor].iter().collect(),
            position,
        })
    }

    /// Explains a kekulization failure caused by an element spelled without
    /// brackets, such as `Co` read as carbon and aromatic oxygen.
    fn unbracketed_element_hints(&self) -> String {
        self.ambiguous_element_spellings
            .iter()
            .map(|(symbol, position)| {
                format!(
                    "; if `{symbol}` at character {position} means the element, write it in \
                     brackets as `[{symbol}]`"
                )
            })
            .collect()
    }
}

fn two_letter_organic_symbol(first: char, second: Option<char>) -> Option<&'static str> {
    match (first, second) {
        ('C', Some('l')) => Some("Cl"),
        ('B', Some('r')) => Some("Br"),
        _ => None,
    }
}

fn standard_bond(symbol: char, position: usize) -> WrittenBond {
    let direction = match symbol {
        '/' => BondDirection::Up,
        '\\' => BondDirection::Down,
        _ => BondDirection::None,
    };
    let order = bond_order_symbol(symbol).unwrap_or(BondOrder::Single);
    WrittenBond {
        specification: BondSpec {
            direction,
            ..BondSpec::with_order(order)
        },
        text: symbol.to_string(),
        position,
    }
}

fn bond_order_symbol(symbol: char) -> Option<BondOrder> {
    match symbol {
        '-' => Some(BondOrder::Single),
        '=' => Some(BondOrder::Double),
        '#' => Some(BondOrder::Triple),
        '$' => Some(BondOrder::Quadruple),
        ':' => Some(BondOrder::Aromatic),
        _ => None,
    }
}

fn drawing_extension_stereo(extension: char) -> Option<BondStereo> {
    match extension {
        'w' => Some(BondStereo::WedgeUp),
        'h' => Some(BondStereo::WedgeDown),
        's' => Some(BondStereo::Wavy),
        'd' => Some(BondStereo::Dashed),
        _ => None,
    }
}

/// The largest class number and a description for each named chirality class.
fn chirality_class_range(class_name: &str) -> Option<(u8, &'static str)> {
    match class_name {
        "TH" => Some((2, "tetrahedral")),
        "AL" => Some((2, "allenal")),
        "SP" => Some((3, "square-planar")),
        "TB" => Some((20, "trigonal-bipyramidal")),
        "OH" => Some((30, "octahedral")),
        _ => None,
    }
}

fn invalid_charge_message(bracket_text: &str, position: usize) -> String {
    format!(
        "charge at character {position} in `{bracket_text}` is malformed; write one sign with \
         an optional count, as in `[Fe+3]` or `[O-]`"
    )
}

fn uppercase_first(symbol: &str) -> String {
    let mut characters = symbol.chars();
    characters
        .next()
        .map(|first| first.to_ascii_uppercase().to_string() + characters.as_str())
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use crate::graph::{AtomChirality, BondDirection, BondOrder};
    use crate::parse_molecule;

    fn error_for(smiles: &str) -> String {
        parse_molecule(smiles).expect_err(smiles)
    }

    #[test]
    fn ring_closures_cannot_bond_an_atom_to_itself() {
        let error = error_for("C11");
        assert!(error.contains("ring closure 1 at character 3"));
        assert!(error.contains("cannot join an atom to itself"));
    }

    #[test]
    fn ring_closures_cannot_duplicate_an_existing_bond() {
        for smiles in ["C12C12", "C1C1", "C(C1)1"] {
            let error = error_for(smiles);
            assert!(error.contains("already bonded"), "{smiles}: {error}");
        }
        assert!(parse_molecule("C1CC.C1").is_ok());
        assert!(parse_molecule("C1CC1C1CC1").is_ok());
        assert!(parse_molecule("C1.C1").is_ok());
    }

    #[test]
    fn ring_closure_bond_symbols_must_agree() {
        let error = error_for("C=1CCCCC-1");
        assert!(
            error.contains("conflicting bond symbols `=` at character 2 and `-` at character 9")
        );
        assert!(parse_molecule("C=1CCCCC=1").is_ok());
        assert!(parse_molecule("C=1CCCCC1").is_ok());
        assert!(parse_molecule("C1CCCCC=1").is_ok());

        // A direction is read away from the end that carries it.
        assert!(parse_molecule("F/C=C/1.Cl\\1").is_ok());
        let error = error_for("F/C=C/1.Cl/1");
        assert!(error.contains("opposite symbols"));
        let error = error_for("C!w1CCCCC!w1");
        assert!(error.contains("one end of the ring closure only"));
    }

    #[test]
    fn ring_closure_bonds_point_away_from_their_written_symbol() {
        let molecule = parse_molecule("F/C=C1.Cl/1").expect("closing direction failed");
        let ring_bond = &molecule.bonds[2];
        assert_eq!((ring_bond.from, ring_bond.to), (3, 2));
        assert_eq!(ring_bond.direction, BondDirection::Up);

        let molecule = parse_molecule("C1CCCCC=1").expect("closing order failed");
        let ring_bond = molecule.bonds.last().expect("ring bond");
        assert_eq!((ring_bond.from, ring_bond.to), (5, 0));
        assert_eq!(ring_bond.order, BondOrder::Double);
    }

    #[test]
    fn atom_maps_are_parsed_separately_from_indices() {
        let molecule = parse_molecule("[CH3:7]O.[NH4+:12]").expect("mapped atoms failed");
        assert_eq!(molecule.atoms[0].atom_map, 7);
        assert_eq!(molecule.atoms[0].hcount, 3);
        assert_eq!(molecule.atoms[1].atom_map, 0);
        assert_eq!(molecule.atoms[2].atom_map, 12);
        assert_eq!(molecule.atoms[2].charge, 1);
        assert_eq!(parse_molecule("[CH3:0]O").unwrap().atoms[0].atom_map, 0);
        assert_eq!(
            parse_molecule("[13C@@H:4](F)(Cl)Br").unwrap().atoms[0].atom_map,
            4
        );
    }

    #[test]
    fn mapped_hydrogens_stay_atoms() {
        let molecule = parse_molecule("C[H:3]").expect("mapped hydrogen failed");
        assert_eq!(molecule.atoms.len(), 2);
        assert_eq!(molecule.atoms[1].atom_map, 3);
        assert_eq!(parse_molecule("C[H]").unwrap().atoms.len(), 1);
    }

    #[test]
    fn malformed_atom_maps_are_rejected() {
        assert!(error_for("[CH3:]O").contains("needs a number"));
        assert!(error_for("[CH3:7:8]O").contains("write one atom map at the end"));
        assert!(error_for("[C:7H3]O").contains("hydrogen count before the charge and atom map"));
        assert!(error_for("[C:99999999999]").contains("too large"));
    }

    #[test]
    fn dangling_bonds_and_dots_are_rejected() {
        assert!(error_for("CC=").contains("followed by the end of the SMILES"));
        assert!(error_for("C(C=)C").contains("followed by `)` at character 5"));
        assert!(error_for("C.").contains("`.` at character 2"));
        assert!(error_for(".C").contains("must separate two fragments"));
        assert!(error_for("C..C").contains("must separate two fragments"));
        assert!(error_for("C==C").contains("directly follows bond"));
        assert!(error_for("C=(C)C").contains("write the bond inside the branch"));
        assert!(parse_molecule("C(.C)C").is_ok());
    }

    #[test]
    fn branch_errors_name_their_positions() {
        assert!(error_for("C(C").contains("opened at character 2 is never closed"));
        assert!(error_for("C)C").contains("unmatched `)` at character 2"));
        assert!(error_for("C()C").contains("empty branch"));
        assert!(error_for("C((C))C").contains("opens directly inside another branch"));
        assert!(error_for("(C)C").contains("must follow the atom it branches from"));
        assert!(error_for("C(1C)CCC1").contains("may hold only ring-closure numbers"));
        assert!(parse_molecule("C(1)CCCC1").is_ok());
    }

    #[test]
    fn charges_follow_opensmiles() {
        assert_eq!(parse_molecule("[O--]").unwrap().atoms[0].charge, -2);
        assert_eq!(parse_molecule("[Fe+3]").unwrap().atoms[0].charge, 3);
        assert_eq!(parse_molecule("[C+15]").unwrap().atoms[0].charge, 15);
        for smiles in ["[C+-]", "[Fe+++]", "[C++2]", "[C+2+]", "[C+100]"] {
            assert!(error_for(smiles).contains("malformed"), "{smiles}");
        }
        assert!(error_for("[C+16]").contains("between -15 and +15 are supported"));
    }

    #[test]
    fn bracket_atom_errors_explain_the_property_order() {
        assert!(error_for("[Xx]").contains("`Xx` at character 2"));
        assert!(error_for("[]").contains("empty bracket atom"));
        assert!(error_for("[1000C]").contains("more than three digits"));
        assert!(error_for("[C@TH3](F)(Cl)(Br)I").contains("use @TH1 through @TH2"));
        assert!(error_for("[C@OH31]").contains("use @OH1 through @OH30"));
        assert!(error_for("[C+H]").contains("hydrogen count before the charge"));
        assert!(error_for("[NH2").contains("unclosed bracket atom"));
    }

    #[test]
    fn every_chirality_class_number_parses() {
        for class in 1..=30 {
            let smiles = format!("C[Co@OH{class}](F)(Cl)(Br)(I)N");
            let molecule = parse_molecule(&smiles).expect(&smiles);
            assert_eq!(
                molecule.atoms[1].chirality,
                AtomChirality::Octahedral(class)
            );
        }
        for class in 1..=20 {
            let smiles = format!("S[As@TB{class}](F)(Cl)(Br)N");
            assert!(parse_molecule(&smiles).is_ok(), "{smiles}");
        }
        assert_eq!(
            parse_molecule("[C@TH2](F)(Cl)(Br)I").unwrap().atoms[0].chirality,
            AtomChirality::TetraClockwise
        );
    }

    #[test]
    fn unbracketed_atom_errors_suggest_corrections() {
        assert!(error_for("CH3").contains("`[CH3]`"));
        assert!(error_for("NaCl").contains("`[Na]`"));
        assert!(error_for("CK").contains("`[K]`"));
        assert!(error_for("C@H").contains("inside a bracket atom"));
        assert!(error_for("C+").contains("`[NH4+]`"));
        assert!(error_for("CC>>CC").contains("reaction SMILES"));
        assert!(error_for("[NH3]->[Pt]").contains("dative bond arrow at character 7"));
        assert!(error_for("[Pt]<-[NH3]").contains("dative bond arrow at character 5"));
        assert!(error_for("C C").contains("whitespace"));
        assert!(error_for("Co").contains("write it in brackets as `[Co]`"));
        assert!(error_for("OMe}").contains("unmatched `}` at character 4"));
        assert!(error_for("C%1CC").contains("two digits"));
    }

    #[test]
    fn unknown_symbols_are_errors_rather_than_crashes() {
        // Text that sorts before every element symbol still yields a diagnostic.
        for smiles in ["A", "BAD", "[A]", "[Aa]", "C{X|#12}", "C{X|A}", "C{X|@}"] {
            assert!(parse_molecule(smiles).is_err(), "{smiles}");
        }
        assert!(error_for("BAD").contains("`A` at character 2 is not an element symbol"));
        assert!(error_for("C{X|#12}").contains("unknown abbreviation style `#12`"));
        assert!(crate::element_from_symbol("Ac").is_some());
        assert!(crate::element_from_symbol("Og").is_some());
        assert!(crate::element_from_symbol("H").is_some());
        assert!(crate::element_from_symbol("").is_none());
    }

    #[test]
    fn folded_hydrogens_keep_the_written_configuration() {
        // Each pair writes the same enantiomer with the hydrogen in different places.
        for (folded, implicit) in [
            ("F[C@](Cl)([H])Br", "F[C@@H](Cl)Br"),
            ("F[C@](Cl)(Br)[H]", "F[C@H](Cl)Br"),
            ("[C@](F)(Cl)(Br)[H]", "[C@@H](F)(Cl)Br"),
            ("[C@]([H])(F)(Cl)Br", "[C@H](F)(Cl)Br"),
        ] {
            let folded_molecule = parse_molecule(folded).expect(folded);
            let implicit_molecule = parse_molecule(implicit).expect(implicit);
            assert_eq!(
                folded_molecule.atoms[folded_molecule.atoms.len() - 4].chirality,
                implicit_molecule.atoms[implicit_molecule.atoms.len() - 4].chirality,
                "{folded} vs {implicit}"
            );
        }
        // A hydrogen on a square-planar center keeps its written position.
        let molecule = parse_molecule("Cl[Pt@SP1](Cl)(N)[H]").expect("square planar hydride");
        assert_eq!(molecule.atoms.len(), 5);
    }

    #[test]
    fn long_chains_and_deep_branches_parse_without_recursion() {
        let chain = "C".repeat(20_000);
        assert_eq!(parse_molecule(&chain).unwrap().atoms.len(), 20_000);
        let depth = 5_000;
        let nested = "C(".repeat(depth) + "C" + &")".repeat(depth);
        assert_eq!(parse_molecule(&nested).unwrap().atoms.len(), depth + 1);
    }
}
