mod error;
mod graph;
mod kekulize;
mod layout;
mod render;

pub use render::LayoutOutput;

use graph::{BondMarker, BondMarkerStyle, BondOrder, MoleculeGraph};
use layout::compute_layout;
use ptable::Element;
use std::iter::Peekable;
use std::str::Chars;

// ── SMILES preprocessing ─────────────────────────────────────────────────────
//
// Rewrites the input into the subset of SMILES that smiles-parser accepts and
// collects side-channel data the parser cannot carry:
//
//   {label}  →  [*]   Abbreviated group (e.g. {PPh3}, {OEt}).
//                     The Nth [*] atom in the output gets abbrev = the Nth label.
//                     A `>` marker inside the label selects the attachment glyph.
//   c, n, …  →  C, N  Unbracketed aromatic atoms are uppercased (the parser
//                     only accepts aliphatic organic-subset symbols); a marker
//                     records which atoms were aromatic so the graph stage can
//                     kekulize. Bracket aromatics ([nH], [se]) parse natively.
//
// The cleaned string is valid parser input; extensions are stripped out.

pub(crate) struct PreprocessedSmiles {
    pub smiles: String,
    /// Labels in the order they appeared, one per {label} token.
    pub abbrev_labels: Vec<AbbreviationLabel>,
    /// One entry for each `/` or `\` token in `smiles`, recording whether it
    /// is plain SMILES or which typed-smiles bond/layout extension it carries.
    pub forced_direction_markers: Vec<BondMarker>,
    /// One flag per unbracketed organic-subset atom token in `smiles`, in
    /// writing order. `true` means the atom was written lowercase (aromatic).
    pub aromatic_atom_markers: Vec<bool>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AbbreviationLabel {
    pub text: String,
    pub style: String,
    pub anchor: usize,
    pub anchor_len: usize,
    /// Explicit non-bonding electron-pair count from an inline `lp=N` modifier.
    /// `None` means the label declared no lone pairs.
    pub lone_pairs: Option<u8>,
    /// Optional page-space displacement in bond-length units. Layout coordinates
    /// remain unchanged; the Typst renderer applies this after molecular rotation.
    pub offset: Option<(f64, f64)>,
}

pub(crate) fn preprocess_smiles(input: &str) -> Result<PreprocessedSmiles, String> {
    let mut parser_compatible_smiles = String::with_capacity(input.len());
    let mut abbreviation_labels = Vec::new();
    let mut bond_markers = Vec::new();
    let mut aromatic_markers = Vec::new();
    let mut in_bracket = false;

    let mut characters = input.chars().peekable();
    while let Some(character) = characters.next() {
        // Bracket contents pass through as a unit because their letters do not
        // represent independent organic-subset atom tokens.
        if in_bracket {
            if character == ']' {
                in_bracket = false;
            }
            parser_compatible_smiles.push(character);
            continue;
        }

        match character {
            '{' => {
                let raw_label = collect_abbreviation_body(&mut characters)?;
                abbreviation_labels.push(parse_abbreviation_label(&raw_label)?);
                parser_compatible_smiles.push_str("[*]");
            }
            '!' => {
                bond_markers.push(parse_drawing_extension(&mut characters)?);
                parser_compatible_smiles.push('/');
            }
            '>' => {
                return Err(
                    "`>` is only valid inside an abbreviation label like `{>PPh3}`".to_string(),
                );
            }
            '[' => {
                in_bracket = true;
                parser_compatible_smiles.push(character);
            }
            ']' => {
                return Err("unmatched `]`; bracket atoms must start with `[`".to_string());
            }
            '}' => {
                return Err("unmatched `}`; custom labels must start with `{`".to_string());
            }
            '/' | '\\' => {
                bond_markers.push(BondMarker {
                    style: BondMarkerStyle::Directional,
                    order: BondOrder::Single,
                    curl: false,
                });
                parser_compatible_smiles.push(character);
            }
            'b' | 'c' | 'n' | 'o' | 'p' | 's' => {
                parser_compatible_smiles.push(character.to_ascii_uppercase());
                aromatic_markers.push(true);
            }
            'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'I' => {
                aromatic_markers.push(false);
                parser_compatible_smiles.push(character);
            }
            _ => parser_compatible_smiles.push(character),
        }
    }
    if in_bracket {
        return Err("unclosed bracket atom; add `]` after the atom specification".to_string());
    }

    let smiles = normalize_post_branch_ring_bonds(&parser_compatible_smiles)?;

    Ok(PreprocessedSmiles {
        smiles,
        abbrev_labels: abbreviation_labels,
        forced_direction_markers: bond_markers,
        aromatic_atom_markers: aromatic_markers,
    })
}

fn collect_abbreviation_body(characters: &mut Peekable<Chars<'_>>) -> Result<String, String> {
    let mut label = String::new();
    for character in characters.by_ref() {
        if character == '}' {
            return Ok(label);
        }
        label.push(character);
    }
    Err("unclosed custom label; add `}` after the label".to_string())
}

fn parse_drawing_extension(characters: &mut Peekable<Chars<'_>>) -> Result<BondMarker, String> {
    let extension = characters
        .next()
        .ok_or_else(|| "incomplete `!` drawing extension".to_string())?;
    let (style, curl) = parse_drawing_style(extension, characters)?;
    let order = consume_optional_bond_order(characters);

    if style != BondMarkerStyle::Plain && order != BondOrder::Single {
        return Err(
            "wedge, hash, wavy, and dashed drawing extensions require a single bond".to_string(),
        );
    }

    Ok(BondMarker { style, order, curl })
}

fn parse_drawing_style(
    extension: char,
    characters: &mut Peekable<Chars<'_>>,
) -> Result<(BondMarkerStyle, bool), String> {
    if extension == 'c' {
        let style = consume_optional_curl_style(characters)?;
        return Ok((style, true));
    }

    let style = bond_marker_style(extension)
        .ok_or_else(|| format!("unknown drawing extension `!{extension}`"))?;
    Ok((style, false))
}

fn consume_optional_curl_style(
    characters: &mut Peekable<Chars<'_>>,
) -> Result<BondMarkerStyle, String> {
    if characters.peek() != Some(&'!') {
        return Ok(BondMarkerStyle::Plain);
    }

    characters.next();
    let extension = characters
        .next()
        .ok_or_else(|| "incomplete drawing extension after `!c`".to_string())?;
    bond_marker_style(extension)
        .ok_or_else(|| format!("unknown drawing extension `!{extension}` after `!c`"))
}

fn bond_marker_style(extension: char) -> Option<BondMarkerStyle> {
    match extension {
        'w' => Some(BondMarkerStyle::WedgeUp),
        'h' => Some(BondMarkerStyle::WedgeDown),
        's' => Some(BondMarkerStyle::Wavy),
        'd' => Some(BondMarkerStyle::Dashed),
        _ => None,
    }
}

fn consume_optional_bond_order(characters: &mut Peekable<Chars<'_>>) -> BondOrder {
    let order = match characters.peek() {
        Some('-') => BondOrder::Single,
        Some('=') => BondOrder::Double,
        Some('#') => BondOrder::Triple,
        Some('$') => BondOrder::Quadruple,
        Some(':') => BondOrder::Aromatic,
        _ => return BondOrder::Single,
    };
    characters.next();
    order
}

fn normalize_post_branch_ring_bonds(input: &str) -> Result<String, String> {
    // Ring closures are branch-like atom modifiers in SMILES, but the parser
    // crate only accepts them before parenthesized branches on the same atom.
    let bytes = input.as_bytes();
    let mut normalized = String::with_capacity(input.len());
    let mut cursor = 0;

    while cursor < bytes.len() {
        let Some(atom_end) = atom_token_end(input, cursor)? else {
            normalized.push(bytes[cursor] as char);
            cursor += 1;
            continue;
        };

        let atom_token = &input[cursor..atom_end];
        let mut suffix_cursor = atom_end;
        let mut ring_bonds = String::new();
        let mut branches = Vec::new();
        let mut has_non_ring_branch = false;

        loop {
            if let Some(end) = ring_bond_token_end(input, suffix_cursor, !has_non_ring_branch) {
                ring_bonds.push_str(&input[suffix_cursor..end]);
                suffix_cursor = end;
                continue;
            }

            if bytes.get(suffix_cursor) == Some(&b'(') {
                let branch_end = matching_paren_end(input, suffix_cursor)?;
                let branch_contents =
                    normalize_post_branch_ring_bonds(&input[suffix_cursor + 1..branch_end - 1])?;
                if ring_bond_sequence_end(&branch_contents, false) == Some(branch_contents.len()) {
                    ring_bonds.push_str(&branch_contents);
                } else {
                    branches.push(format!("({branch_contents})"));
                    has_non_ring_branch = true;
                }
                suffix_cursor = branch_end;
                continue;
            }

            break;
        }

        normalized.push_str(atom_token);
        normalized.push_str(&ring_bonds);
        for branch in branches {
            normalized.push_str(&branch);
        }
        cursor = suffix_cursor;
    }

    Ok(normalized)
}

fn atom_token_end(input: &str, start: usize) -> Result<Option<usize>, String> {
    let bytes = input.as_bytes();
    let Some(&first_byte) = bytes.get(start) else {
        return Ok(None);
    };

    if first_byte == b'[' {
        let mut cursor = start + 1;
        while cursor < bytes.len() {
            if bytes[cursor] == b']' {
                return Ok(Some(cursor + 1));
            }
            cursor += 1;
        }
        return Err("unclosed bracket atom".to_string());
    }

    if start + 1 < bytes.len()
        && ((bytes[start] == b'C' && bytes[start + 1] == b'l')
            || (bytes[start] == b'B' && bytes[start + 1] == b'r'))
    {
        return Ok(Some(start + 2));
    }

    if matches!(
        first_byte,
        b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I' | b'*'
    ) {
        return Ok(Some(start + 1));
    }

    Ok(None)
}

fn matching_paren_end(input: &str, start: usize) -> Result<usize, String> {
    let bytes = input.as_bytes();
    let mut depth = 0usize;
    let mut cursor = start;

    while cursor < bytes.len() {
        match bytes[cursor] {
            b'[' => {
                cursor += 1;
                while cursor < bytes.len() && bytes[cursor] != b']' {
                    cursor += 1;
                }
                if cursor == bytes.len() {
                    return Err("unclosed bracket atom".to_string());
                }
            }
            b'(' => depth += 1,
            b')' => {
                depth = depth.saturating_sub(1);
                if depth == 0 {
                    return Ok(cursor + 1);
                }
            }
            _ => {}
        }
        cursor += 1;
    }

    Err("unclosed branch".to_string())
}

fn ring_bond_sequence_end(input: &str, allow_directional: bool) -> Option<usize> {
    let mut cursor = 0;
    let mut found_ring_bond = false;
    while let Some(end) = ring_bond_token_end(input, cursor, allow_directional) {
        found_ring_bond = true;
        cursor = end;
    }
    found_ring_bond.then_some(cursor)
}

fn ring_bond_token_end(input: &str, start: usize, allow_directional: bool) -> Option<usize> {
    let bytes = input.as_bytes();
    if start >= bytes.len() {
        return None;
    }

    let mut cursor = start;
    if matches!(bytes[cursor], b'-' | b'=' | b'#' | b'$' | b':')
        || (allow_directional && matches!(bytes[cursor], b'/' | b'\\'))
    {
        cursor += 1;
        if cursor >= bytes.len() {
            return None;
        }
    }

    if bytes[cursor].is_ascii_digit() {
        return Some(cursor + 1);
    }

    if bytes[cursor] == b'%'
        && cursor + 2 < bytes.len()
        && bytes[cursor + 1].is_ascii_digit()
        && bytes[cursor + 2].is_ascii_digit()
    {
        return Some(cursor + 3);
    }

    None
}

/// Parses one `{...}` abbreviation body into its displayed text, optional
/// attachment marker, optional style, and named rendering modifiers.
///
/// The body is split on `|` into fields:
///   - the first field is the displayed label, optionally carrying one `>`
///     attachment marker before the anchor glyph;
///   - an optional plain style field (a color or element token);
///   - named `lp=N` and `offset=(x,y)` modifiers.
///
/// The plain style field, when present, must precede any named modifier. A
/// second field that begins with a named modifier therefore means the label has
/// no style.
fn parse_abbreviation_label(raw_label: &str) -> Result<AbbreviationLabel, String> {
    let mut fields = raw_label.split('|');
    let raw_text = fields.next().unwrap_or("").trim();

    let mut style = String::new();
    let mut style_set = false;
    let mut lone_pairs: Option<u8> = None;
    let mut offset: Option<(f64, f64)> = None;
    let mut seen_named = false;
    for field in fields {
        let trimmed = field.trim();
        if let Some((raw_modifier_name, raw_modifier_value)) = trimmed.split_once('=') {
            let modifier_name = raw_modifier_name.trim();
            let modifier_value = raw_modifier_value.trim();
            match modifier_name {
                "lp" => {
                    if lone_pairs.is_some() {
                        return Err(
                            "abbreviation label has more than one `lp=` modifier".to_string()
                        );
                    }
                    let count: u8 = modifier_value.parse().map_err(|_| {
                        format!(
                            "abbreviation `lp=` needs an integer from 1 to 4, got `{modifier_value}`"
                        )
                    })?;
                    if !(1..=4).contains(&count) {
                        return Err(format!(
                            "abbreviation `lp=` must be from 1 to 4, got {count}"
                        ));
                    }
                    lone_pairs = Some(count);
                }
                "offset" => {
                    if offset.is_some() {
                        return Err(
                            "abbreviation label has more than one `offset=` modifier".to_string()
                        );
                    }
                    offset = Some(parse_abbreviation_offset(modifier_value)?);
                }
                other => {
                    return Err(format!("unknown abbreviation modifier `{other}=`"));
                }
            }
            seen_named = true;
        } else {
            if seen_named {
                return Err(
                    "abbreviation style must come before named modifiers like `lp=`".to_string(),
                );
            }
            if style_set {
                return Err("abbreviation label has more than one style field".to_string());
            }
            style = trimmed.to_string();
            style_set = true;
        }
    }
    validate_abbreviation_style(&style)?;

    let marker_count = raw_text.chars().filter(|&ch| ch == '>').count();
    if marker_count > 1 {
        return Err(
            "abbreviation labels may contain at most one `>` attachment marker".to_string(),
        );
    }

    let mut text = String::with_capacity(raw_text.len());
    let mut attachment_marker = None;
    for character in raw_text.chars() {
        if character == '>' {
            attachment_marker = Some(text.chars().count());
        } else {
            text.push(character);
        }
    }

    let label_characters: Vec<char> = text.chars().collect();
    let (anchor, anchor_len) = if let Some(marker_position) = attachment_marker {
        if label_characters.is_empty() {
            return Err("abbreviation attachment marker `>` needs a label glyph".to_string());
        }
        if marker_position >= label_characters.len() {
            return Err(
                "abbreviation attachment marker `>` must precede a label glyph".to_string(),
            );
        }
        let anchor = marker_position;
        let anchor_len = if label_characters
            .get(anchor)
            .is_some_and(|character| character.is_ascii_uppercase())
            && label_characters
                .get(anchor + 1)
                .is_some_and(|character| character.is_ascii_lowercase())
        {
            2
        } else {
            1
        };
        (anchor, anchor_len)
    } else {
        (0, 0)
    };

    Ok(AbbreviationLabel {
        text,
        style,
        anchor,
        anchor_len,
        lone_pairs,
        offset,
    })
}

fn validate_abbreviation_style(style: &str) -> Result<(), String> {
    if style.is_empty() {
        return Ok(());
    }
    const NAMED_COLORS: &[&str] = &[
        "red", "blue", "green", "black", "gray", "grey", "silver", "white", "orange", "yellow",
        "brown", "pink", "purple", "cyan", "lime", "teal", "maroon", "navy",
    ];
    let valid_hex = style.len() == 7
        && style.starts_with('#')
        && style[1..].bytes().all(|byte| byte.is_ascii_hexdigit());
    if valid_hex || NAMED_COLORS.contains(&style) || Element::from_symbol(style).is_some() {
        return Ok(());
    }
    Err(format!(
        "unknown abbreviation style `{style}`; use an element symbol, a supported \
         color name, or a #RRGGBB color"
    ))
}

fn parse_abbreviation_offset(raw_offset: &str) -> Result<(f64, f64), String> {
    let value = raw_offset.trim();
    let Some(components) = value
        .strip_prefix('(')
        .and_then(|inner| inner.strip_suffix(')'))
    else {
        return Err(format!(
            "abbreviation `offset=` needs `(x, y)`, got `{raw_offset}`"
        ));
    };

    let mut coordinates = components.split(',').map(str::trim);
    let Some(raw_x) = coordinates.next() else {
        return Err(format!(
            "abbreviation `offset=` needs two numbers, got `{raw_offset}`"
        ));
    };
    let Some(raw_y) = coordinates.next() else {
        return Err(format!(
            "abbreviation `offset=` needs two numbers, got `{raw_offset}`"
        ));
    };
    if coordinates.next().is_some() || raw_x.is_empty() || raw_y.is_empty() {
        return Err(format!(
            "abbreviation `offset=` needs exactly two numbers, got `{raw_offset}`"
        ));
    }

    let parse_coordinate = |coordinate: &str| -> Result<f64, String> {
        let parsed = coordinate.parse::<f64>().map_err(|_| {
            format!("abbreviation `offset=` coordinates must be numbers, got `{coordinate}`")
        })?;
        if !parsed.is_finite() {
            return Err(format!(
                "abbreviation `offset=` coordinates must be finite, got `{coordinate}`"
            ));
        }
        Ok(parsed)
    };

    Ok((parse_coordinate(raw_x)?, parse_coordinate(raw_y)?))
}

/// Apply collected abbreviation labels to the matching `*` atoms in the graph
/// (N-th `*` atom ← N-th label, in atom-index order).
fn assign_abbreviation_labels(molecule: &mut MoleculeGraph, labels: &[AbbreviationLabel]) {
    let mut label_iter = labels.iter();
    for atom in molecule.atoms.iter_mut().filter(|atom| atom.symbol == "*") {
        let Some(label) = label_iter.next() else {
            break;
        };
        atom.abbrev = label.text.clone();
        atom.abbrev_style = label.style.clone();
        atom.abbrev_anchor = label.anchor;
        atom.abbrev_anchor_len = label.anchor_len;
        atom.abbrev_lone_pairs = label.lone_pairs.unwrap_or(0);
        let (offset_x, offset_y) = label.offset.unwrap_or((0.0, 0.0));
        atom.abbrev_offset_x = offset_x;
        atom.abbrev_offset_y = offset_y;
    }
}

fn parse_molecule(smiles: &str) -> Result<MoleculeGraph, String> {
    if smiles.trim().is_empty() {
        return Err(
            "invalid SMILES: the expression is empty; provide at least one atom".to_string(),
        );
    }
    let preprocessed =
        preprocess_smiles(smiles).map_err(|error| format!("invalid SMILES {smiles:?}: {error}"))?;
    let mut molecule = MoleculeGraph::from_smiles(
        &preprocessed.smiles,
        preprocessed.forced_direction_markers,
        preprocessed.aromatic_atom_markers,
    )
    .map_err(|error| format!("invalid SMILES {smiles:?}: {error}"))?;
    assign_abbreviation_labels(&mut molecule, &preprocessed.abbrev_labels);
    Ok(molecule)
}

fn atomic_mass(atom: &graph::Atom) -> Result<f64, String> {
    if !atom.abbrev.is_empty() {
        return Err(format!(
            "cannot compute molecular weight: abbreviation {{{}}} has no defined composition",
            atom.abbrev
        ));
    }
    if atom.symbol == "*" {
        return Err("cannot compute molecular weight: wildcard atom `*` has no mass".to_string());
    }
    if let Some(isotope) = atom.isotope {
        return Err(format!(
            "cannot compute molecular weight: isotope [{isotope}{}] needs a nuclide mass, \
             not a standard atomic weight",
            atom.symbol
        ));
    }

    Element::from_symbol(&atom.symbol)
        .map(|element| element.get_atomic_mass() as f64)
        .ok_or_else(|| {
            format!(
                "cannot compute molecular weight: unknown element {}",
                atom.symbol
            )
        })
}

fn hydrogen_count(molecule: &MoleculeGraph, atom_index: usize) -> u8 {
    molecule.atoms[atom_index].hcount + layout::implicit_h_count(molecule, atom_index)
}

fn compute_molecular_weight(molecule: &MoleculeGraph) -> Result<f64, String> {
    let hydrogen_mass = Element::Hydrogen.get_atomic_mass() as f64;
    let mut molecular_weight = 0.0;

    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        molecular_weight += atomic_mass(atom)?;
        molecular_weight += f64::from(hydrogen_count(molecule, atom_index)) * hydrogen_mass;
    }

    Ok(molecular_weight)
}

// ── WASM / Typst plugin entrypoint ──────────────────────────────────────────

#[cfg(target_arch = "wasm32")]
mod wasm_entrypoint {
    use super::*;
    use wasm_minimal_protocol::*;

    initiate_protocol!();

    /// Called from Typst as `smiles-plugin.layout(bytes(smiles-str))`.
    /// Returns JSON-encoded `LayoutOutput`.
    /// Accepts extended SMILES with `{label}` abbreviation syntax.
    #[wasm_func]
    pub fn layout(smiles: &[u8]) -> Result<Vec<u8>, String> {
        let smiles =
            core::str::from_utf8(smiles).map_err(|error| format!("UTF-8 error: {error}"))?;
        let molecule = parse_molecule(smiles)?;
        let layout = compute_layout(&molecule)?;
        serde_json::to_vec(&layout).map_err(|error| format!("JSON error: {error}"))
    }

    /// Called from Typst as `smiles-plugin.mol_weight(bytes(smiles-str))`.
    /// Returns the molecular weight in g/mol as a JSON number.
    #[wasm_func]
    pub fn mol_weight(smiles: &[u8]) -> Result<Vec<u8>, String> {
        let smiles =
            core::str::from_utf8(smiles).map_err(|error| format!("UTF-8 error: {error}"))?;
        let molecule = parse_molecule(smiles)?;
        let molecular_weight = compute_molecular_weight(&molecule)?;
        serde_json::to_vec(&molecular_weight).map_err(|error| format!("JSON error: {error}"))
    }
}

// ── Native entrypoint for tests / CLI ───────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
pub fn layout_native(smiles: &str) -> Result<LayoutOutput, String> {
    let molecule = parse_molecule(smiles)?;
    compute_layout(&molecule)
}

#[cfg(not(target_arch = "wasm32"))]
pub fn mol_weight_native(smiles: &str) -> Result<f64, String> {
    let molecule = parse_molecule(smiles)?;
    compute_molecular_weight(&molecule)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Smallest distance between any two distinct non-virtual atoms.
    fn min_atom_distance(layout_output: &LayoutOutput) -> f64 {
        let mut min = f64::INFINITY;
        for i in 0..layout_output.atoms.len() {
            for j in (i + 1)..layout_output.atoms.len() {
                if layout_output.atoms[i].virtual_h || layout_output.atoms[j].virtual_h {
                    continue;
                }
                min = min.min(layout_output.atoms[i].pos.dist(layout_output.atoms[j].pos));
            }
        }
        min
    }

    #[test]
    fn ortho_ring_substituents_do_not_collide() {
        // Aspirin: the acetyl C=O and the carboxyl OH grow from ortho ring
        // positions toward each other and must not land on the same point.
        let layout_output = layout_native("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap();
        assert!(
            min_atom_distance(&layout_output) > 0.5,
            "atoms overlap: min distance {}",
            min_atom_distance(&layout_output)
        );

        // The carboxyl carbon (atom 10: ring C9, =O 11, OH 12) must keep its
        // textbook trigonal geometry: the conflict is resolved by flipping
        // the acetyl branch aside, not by squeezing the carboxyl bonds into a
        // narrow fan or a straight line.
        let carboxyl_center = layout_output.atoms[10].pos;
        let angles: Vec<f64> = [9, 11, 12]
            .iter()
            .map(|&neighbor_index| {
                let neighbor_position = layout_output.atoms[neighbor_index].pos;
                (neighbor_position.y - carboxyl_center.y)
                    .atan2(neighbor_position.x - carboxyl_center.x)
            })
            .collect();
        for first_angle_index in 0..angles.len() {
            for second_angle_index in (first_angle_index + 1)..angles.len() {
                let mut delta = (angles[first_angle_index] - angles[second_angle_index]).abs();
                if delta > std::f64::consts::PI {
                    delta = 2.0 * std::f64::consts::PI - delta;
                }
                let deg = delta.to_degrees();
                assert!(
                    (deg - 120.0).abs() < 1.0,
                    "carboxyl bond angle {deg:.1} degrees is not the ideal 120"
                );
            }
        }
    }

    #[test]
    fn crowded_branch_chains_do_not_collide() {
        // Two ortho substituent chains long enough to sweep past each other.
        let layout_output = layout_native("CCCC1=CC=CC=C1CCC").unwrap();
        assert!(min_atom_distance(&layout_output) > 0.5);
    }

    #[test]
    fn spiro_rings_share_one_atom_without_overlap() {
        // 1,6-dioxaspiro[4.4]nonane: two five-membered rings joined at a single
        // spiro atom. Each ring must close as a regular pentagon rather than
        // unravel into a chain, so no atoms collide.
        let layout_output = layout_native("CC[C@H](O1)CC[C@@]12CCCO2").unwrap();
        assert!(
            min_atom_distance(&layout_output) > 0.5,
            "spiro atoms overlap: min distance {}",
            min_atom_distance(&layout_output)
        );
        // The two ring oxygens (atoms 3 and 10) both neighbor the spiro carbon
        // (atom 6) and must sit one bond length away from it.
        let spiro = layout_output.atoms[6].pos;
        for oxygen_index in [3usize, 10] {
            let distance = layout_output.atoms[oxygen_index].pos.dist(spiro);
            assert!(
                (distance - 1.0).abs() < 0.05,
                "ring O-spiro bond length {} for atom {}",
                distance,
                oxygen_index
            );
        }
    }

    #[test]
    fn benzene_kekule() {
        let layout_output = layout_native("C1=CC=CC=C1").expect("benzene layout failed");
        assert_eq!(layout_output.atoms.len(), 6);
        assert_eq!(layout_output.bonds.len(), 6);
    }

    #[test]
    fn ethanol() {
        // CCO: 3 real atoms + 1 virtual H for terminal O (implicit OH)
        let layout_output = layout_native("CCO").expect("ethanol layout failed");
        assert_eq!(layout_output.atoms.len(), 4);
        assert!(layout_output.atoms[3].virtual_h);
    }

    #[test]
    fn implicit_h_counts_ethanol() {
        // implicit_h counts for real atoms only (virtual H has implicit_h = 0)
        let layout_output = layout_native("CCO").expect("ethanol layout failed");
        let counts: Vec<u8> = layout_output
            .atoms
            .iter()
            .map(|atom| atom.implicit_h)
            .collect();
        assert_eq!(counts, vec![3, 2, 1, 0]);
    }

    #[test]
    fn explicit_h_suppresses_implicit_h() {
        let layout_output = layout_native("[NH4+]").expect("ammonium layout failed");
        assert_eq!(layout_output.atoms[0].hcount, 4);
        assert_eq!(layout_output.atoms[0].implicit_h, 0);
    }

    #[test]
    fn bracket_hydrogen_atoms_fold_into_neighbor() {
        // Methane written with four explicit [H] atoms collapses into one
        // carbon carrying the folded hydrogen count; the only remaining H is the
        // virtual label placeholder, never a drawn atom.
        let layout_output = layout_native("C([H])([H])([H])[H]").expect("methane layout failed");
        assert_eq!(layout_output.atoms[0].symbol, "C");
        assert_eq!(layout_output.atoms[0].hcount, 4);
        assert_eq!(layout_output.atoms[0].implicit_h, 0);
        let drawn_hydrogen_count = layout_output
            .atoms
            .iter()
            .filter(|atom| atom.symbol == "H" && !atom.virtual_h)
            .count();
        assert_eq!(drawn_hydrogen_count, 0);
    }

    #[test]
    fn folded_hydrogens_drop_from_neighbor_ordering() {
        // The dichloromethyl carbon keeps only its three heavy neighbors after
        // the explicit hydrogen is folded away.
        let layout_output = layout_native("ClC([H])Cl").expect("dichloromethane layout failed");
        let heavy_atom_count = layout_output
            .atoms
            .iter()
            .filter(|atom| atom.symbol != "H")
            .count();
        assert_eq!(heavy_atom_count, 3);
        assert_eq!(
            layout_output
                .bonds
                .iter()
                .filter(|bond| !bond.virtual_bond)
                .count(),
            2
        );
    }

    #[test]
    fn isotopic_and_charged_hydrogens_are_kept() {
        // Deuterium is a real, drawn atom; only the plain [H] atoms fold away.
        let layout_output =
            layout_native("[2H]C([H])([H])[H]").expect("deuteromethane layout failed");
        let drawn_hydrogen_count = layout_output
            .atoms
            .iter()
            .filter(|atom| atom.symbol == "H" && !atom.virtual_h)
            .count();
        assert_eq!(drawn_hydrogen_count, 1);
    }

    #[test]
    fn double_bond() {
        let layout_output = layout_native("C=C").expect("ethylene layout failed");
        assert_eq!(layout_output.bonds[0].order, 2);
    }

    #[test]
    fn cumulated_double_bonds_are_linear() {
        let layout_output = layout_native("O=C=O").expect("carbon dioxide layout failed");
        assert!(atoms_are_collinear(&layout_output, 0, 1, 2));
    }

    #[test]
    fn single_triple_chain_is_linear_at_middle_atom() {
        let layout_output = layout_native("CC#N").expect("nitrile layout failed");
        assert!(atoms_are_collinear(&layout_output, 0, 1, 2));
    }

    #[test]
    fn saturated_two_neighbor_atom_keeps_zigzag() {
        let layout_output = layout_native("CCC").expect("propane layout failed");
        assert!(!atoms_are_collinear(&layout_output, 0, 1, 2));
    }

    #[test]
    fn cyclohexane() {
        let layout_output = layout_native("C1CCCCC1").expect("cyclohexane layout failed");
        assert_eq!(layout_output.atoms.len(), 6);
        assert_eq!(layout_output.bonds.len(), 6);
    }

    #[test]
    fn separated_ring_systems_do_not_overlap() {
        let layout_output = layout_native("CC(N)C(=O)OCCC1=CC=CC=C1NCC1=CC=CC=C1")
            .expect("two-ring molecule layout failed");
        assert!(max_bond_length(&layout_output) < 1.2);
    }

    #[test]
    fn symmetric_hub_is_mirror_symmetric() {
        // Tetraethylmethane's four equal arms must be drawn as a mirror-symmetric
        // figure (the "tweezers" depiction), not a rotational pinwheel. The layout
        // is built about the vertical axis, so reflecting every atom across that
        // axis must reproduce the same set of positions.
        let layout_output =
            layout_native("CCC(CC)(CC)CC").expect("tetraethylmethane layout failed");
        let cx = layout_output
            .atoms
            .iter()
            .map(|atom| atom.pos.x)
            .sum::<f64>()
            / layout_output.atoms.len() as f64;
        for atom in &layout_output.atoms {
            let mirrored_x = 2.0 * cx - atom.pos.x;
            let found = layout_output.atoms.iter().any(|b| {
                (b.pos.x - mirrored_x).abs() < 1e-6 && (b.pos.y - atom.pos.y).abs() < 1e-6
            });
            assert!(
                found,
                "no mirror partner for atom at ({:.3}, {:.3})",
                atom.pos.x, atom.pos.y
            );
        }
    }

    #[test]
    fn branch_point_routes_largest_subtree_straight_ahead() {
        // The central acetal carbon has four substituents; the long ester chain
        // must continue away from the rest of the molecule instead of folding
        // back over the other ester group.
        let layout_output =
            layout_native("BrCC(=O)OC(C)(O)OC(=O)C").expect("acetal diester layout failed");
        assert!(
            min_nonbonded_distance(&layout_output) > 0.5,
            "branches overlap: min non-bonded distance {}",
            min_nonbonded_distance(&layout_output)
        );
    }

    #[test]
    fn steroid_ring_system_has_regular_bond_lengths() {
        let layout_output =
            layout_native("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
                .expect("steroid-like molecule layout failed");
        assert!(max_bond_length(&layout_output) < 1.25);
    }

    #[test]
    fn isobutane() {
        let layout_output = layout_native("CC(C)C").expect("isobutane layout failed");
        assert_eq!(layout_output.atoms.len(), 4);
        assert_eq!(layout_output.bonds.len(), 3);
    }

    #[test]
    fn naphthalene_kekule() {
        let layout_output =
            layout_native("C1=CC2=CC=CC=C2C=C1").expect("naphthalene layout failed");
        assert_eq!(layout_output.atoms.len(), 10);
        assert_eq!(layout_output.bonds.len(), 11);
    }

    // ── Dot-disconnected structures ──────────────────────────────────────────

    #[test]
    fn dot_creates_no_bond() {
        let layout_output = layout_native("CCO.CCO").expect("two ethanols failed");
        assert_eq!(
            layout_output.atoms.iter().filter(|a| !a.virtual_h).count(),
            6
        );
        // Two C-C-O fragments: 4 real bonds, none crossing the dot.
        let real: Vec<_> = layout_output
            .bonds
            .iter()
            .filter(|b| !b.virtual_bond)
            .collect();
        assert_eq!(real.len(), 4);
        assert!(real.iter().all(|b| (b.from < 3) == (b.to < 3)));
    }

    #[test]
    fn salt_fragments_are_disconnected_and_ordered() {
        // Sodium acetate: no bond may touch the sodium ion, and fragments keep
        // SMILES writing order left to right with a visible gap.
        let layout_output = layout_native("CC(=O)[O-].[Na+]").expect("sodium acetate failed");
        let sodium_index = 4;
        assert_eq!(layout_output.atoms[sodium_index].symbol, "Na");
        assert_eq!(layout_output.atoms[sodium_index].charge, 1);
        assert!(layout_output
            .bonds
            .iter()
            .all(|bond| bond.from != sodium_index && bond.to != sodium_index));

        let max_acetate_x = (0..4)
            .map(|atom_index| layout_output.atoms[atom_index].pos.x)
            .fold(f64::MIN, f64::max);
        assert!(
            layout_output.atoms[sodium_index].pos.x >= max_acetate_x + 1.0,
            "Na+ should sit clearly right of the acetate fragment"
        );
    }

    #[test]
    fn bare_ion_pair() {
        let layout_output = layout_native("[Na+].[Cl-]").expect("sodium chloride failed");
        assert_eq!(layout_output.atoms.len(), 2);
        assert!(layout_output.bonds.is_empty());
        assert!(layout_output.atoms[1].pos.x > layout_output.atoms[0].pos.x);
    }

    #[test]
    fn three_fragments() {
        let layout_output = layout_native("O.O.O").expect("three waters failed");
        let heavy: Vec<_> = layout_output
            .atoms
            .iter()
            .filter(|a| !a.virtual_h)
            .collect();
        assert_eq!(heavy.len(), 3);
        assert!(layout_output.bonds.iter().all(|bond| bond.virtual_bond));
        assert!(heavy[0].pos.x < heavy[1].pos.x && heavy[1].pos.x < heavy[2].pos.x);
    }

    #[test]
    fn ring_closure_across_dot_joins_fragments() {
        // OpenSMILES: "C1.C1" is ethane — the ring-closure bond still forms
        // even though a dot separates the digits.
        let layout_output = layout_native("C1.C1").expect("dot ring closure failed");
        assert_eq!(layout_output.atoms.len(), 2);
        assert_eq!(layout_output.bonds.len(), 1);
        assert_eq!(layout_output.atoms[0].implicit_h, 3);
    }

    #[test]
    fn aromatic_fragments_kekulize_independently() {
        let layout_output = layout_native("c1ccccc1.c1ccccc1").expect("two benzenes failed");
        assert_eq!(layout_output.atoms.len(), 12);
        assert_eq!(order_counts(&layout_output), (6, 6));
        assert!(layout_output
            .bonds
            .iter()
            .all(|b| (b.from < 6) == (b.to < 6)));
    }

    // ── Aromatic (lowercase) SMILES and kekulization ─────────────────────────

    fn order_counts(layout_output: &LayoutOutput) -> (usize, usize) {
        let real = layout_output.bonds.iter().filter(|b| !b.virtual_bond);
        (
            real.clone().filter(|bond| bond.order == 1).count(),
            real.filter(|bond| bond.order == 2).count(),
        )
    }

    #[test]
    fn benzene_aromatic() {
        let layout_output = layout_native("c1ccccc1").expect("aromatic benzene failed");
        assert_eq!(layout_output.atoms.len(), 6);
        assert_eq!(order_counts(&layout_output), (3, 3));
        assert!(layout_output.atoms.iter().all(|atom| atom.implicit_h == 1));
    }

    #[test]
    fn aromatic_atom_indices_match_writing_order() {
        // Kekulization must not reorder atoms: index N is the Nth atom token,
        // so show-indices / highlight / arrow references keep working.
        let layout_output = layout_native("Cc1ccncc1").expect("4-methylpyridine failed");
        let symbols: Vec<&str> = layout_output
            .atoms
            .iter()
            .filter(|a| !a.virtual_h)
            .map(|atom| atom.symbol.as_str())
            .collect();
        assert_eq!(symbols, vec!["C", "C", "C", "C", "N", "C", "C"]);
    }

    #[test]
    fn pyridine_aromatic() {
        let layout_output = layout_native("c1ccncc1").expect("pyridine failed");
        assert_eq!(layout_output.atoms[3].symbol, "N");
        assert_eq!(layout_output.atoms[3].implicit_h, 0);
        assert_eq!(order_counts(&layout_output), (3, 3));
    }

    #[test]
    fn pyrrole_aromatic() {
        let layout_output = layout_native("c1cc[nH]c1").expect("pyrrole failed");
        let nitrogen = layout_output
            .atoms
            .iter()
            .find(|atom| atom.symbol == "N")
            .unwrap();
        assert_eq!(nitrogen.hcount, 1);
        assert_eq!(order_counts(&layout_output), (3, 2));
    }

    #[test]
    fn furan_and_thiophene_aromatic() {
        for smiles in ["c1occc1", "c1sccc1"] {
            let layout_output = layout_native(smiles).expect("5-ring heteroaromatic failed");
            assert_eq!(
                order_counts(&layout_output),
                (3, 2),
                "wrong kekulization for {smiles}"
            );
            let hetero = &layout_output.atoms[1];
            assert_eq!(hetero.implicit_h, 0);
        }
    }

    #[test]
    fn imidazole_aromatic() {
        let layout_output = layout_native("c1cnc[nH]1").expect("imidazole failed");
        assert_eq!(order_counts(&layout_output), (3, 2));
    }

    #[test]
    fn n_methylpyrrole_aromatic() {
        // A three-connected aromatic n carries no H and no double bond.
        let layout_output = layout_native("Cn1cccc1").expect("N-methylpyrrole failed");
        assert_eq!(layout_output.atoms[1].symbol, "N");
        assert_eq!(layout_output.atoms[1].implicit_h, 0);
        assert_eq!(order_counts(&layout_output), (4, 2));
    }

    #[test]
    fn naphthalene_aromatic() {
        let layout_output = layout_native("c1ccc2ccccc2c1").expect("aromatic naphthalene failed");
        assert_eq!(layout_output.atoms.len(), 10);
        assert_eq!(order_counts(&layout_output), (6, 5));
    }

    #[test]
    fn indane_mixed_aromatic_aliphatic() {
        // Spec example: aromatic ring fused to an aliphatic ring.
        let layout_output = layout_native("c1ccc2CCCc2c1").expect("indane failed");
        assert_eq!(layout_output.atoms.len(), 9);
        assert_eq!(order_counts(&layout_output), (7, 3));
    }

    #[test]
    fn biphenyl_explicit_and_implicit_single_link() {
        for smiles in ["c1ccccc1-c1ccccc1", "c1ccccc1c1ccccc1"] {
            let layout_output = layout_native(smiles).expect("biphenyl failed");
            assert_eq!(layout_output.atoms.len(), 12);
            assert_eq!(
                order_counts(&layout_output),
                (7, 6),
                "wrong kekulization for {smiles}"
            );
            // The inter-ring bond stays single.
            let link = layout_output
                .bonds
                .iter()
                .find(|b| (b.from < 6) != (b.to < 6))
                .unwrap();
            assert_eq!(link.order, 1);
        }
    }

    #[test]
    fn pyridinone_exocyclic_double_bond() {
        // The carbonyl carbon already has its double bond outside the ring and
        // must not receive another one during kekulization.
        let layout_output = layout_native("O=c1cccc[nH]1").expect("2-pyridinone failed");
        assert_eq!(order_counts(&layout_output), (4, 3));
    }

    #[test]
    fn explicit_aromatic_bond_symbol() {
        let layout_output = layout_native("c1:c:c:c:c:c1").expect("explicit ':' benzene failed");
        assert_eq!(order_counts(&layout_output), (3, 3));
    }

    #[test]
    fn charged_aromatic_rings() {
        // Pyridinium: the protonated nitrogen still takes part in a double bond.
        let pyridinium = layout_native("c1cc[nH+]cc1").expect("pyridinium failed");
        assert_eq!(order_counts(&pyridinium), (3, 3));

        // Pyrylium: positively charged oxygen participates in a double bond.
        let pyrylium = layout_native("c1cc[o+]cc1").expect("pyrylium failed");
        assert_eq!(order_counts(&pyrylium), (3, 3));
    }

    #[test]
    fn wildcard_in_aromatic_ring() {
        let layout_output = layout_native("c1cc*cc1").expect("aromatic ring with wildcard failed");
        let (_, doubles) = order_counts(&layout_output);
        assert!(
            doubles >= 2,
            "expected an alternating pattern, got {doubles} double bonds"
        );
    }

    #[test]
    fn azulene_aromatic() {
        // Non-alternant 5-7 fused system: every carbon needs a double bond.
        let layout_output = layout_native("c1ccc2cccc2cc1").expect("azulene failed");
        assert_eq!(layout_output.atoms.len(), 10);
        assert_eq!(order_counts(&layout_output), (6, 5));
    }

    #[test]
    fn caffeine_aromatic() {
        let layout_output = layout_native("Cn1cnc2c1c(=O)n(C)c(=O)n2C").expect("caffeine failed");
        // Purine core: the imidazole C=N plus the C4=C5 bridge double bond,
        // and the two exocyclic carbonyls.
        assert_eq!(
            layout_output.atoms.iter().filter(|a| !a.virtual_h).count(),
            14
        );
        let (_, doubles) = order_counts(&layout_output);
        assert_eq!(doubles, 4);
    }

    #[test]
    fn aromatic_bond_symbol_requires_aromatic_atoms() {
        let err = layout_native("C:C").expect_err("':' between aliphatic atoms should fail");
        assert!(err.contains("aromatic"));
    }

    #[test]
    fn aromatic_atom_outside_ring_errors() {
        let err = layout_native("Cc").expect_err("acyclic aromatic atom should fail");
        assert!(err.contains("ring"));
    }

    #[test]
    fn unkekulizable_ring_errors() {
        // Pyrrole written without its hydrogen has five atoms all demanding a
        // double bond; no perfect matching exists.
        let err = layout_native("c1ccnc1").expect_err("H-less pyrrole should fail");
        assert!(err.contains("kekulize"));
    }

    // ── Standards-first stereochemistry and drawing extensions ───────────────

    #[test]
    fn forced_wedge_up_bond() {
        let layout_output = layout_native("C!wN").expect("forced wedge up failed");
        assert_eq!(layout_output.bonds[0].stereo, "wedge_up");
        assert_eq!(layout_output.bonds[0].direction, "none");
        assert!(layout_output.bonds[0].forced_stereo);
    }

    #[test]
    fn forced_wedge_down_bond() {
        let layout_output = layout_native("C!hN").expect("forced wedge down failed");
        assert_eq!(layout_output.bonds[0].stereo, "wedge_down");
        assert_eq!(layout_output.bonds[0].direction, "none");
        assert!(layout_output.bonds[0].forced_stereo);
    }

    #[test]
    fn forced_wavy_bond() {
        let layout_output = layout_native("C!sN").expect("forced wavy bond failed");
        assert_eq!(layout_output.bonds[0].stereo, "wavy");
        assert_eq!(layout_output.bonds[0].direction, "none");
        assert!(!layout_output.bonds[0].forced_stereo);
    }

    #[test]
    fn forced_dashed_bond() {
        let layout_output = layout_native("C!dN").expect("forced dashed bond failed");
        assert_eq!(layout_output.bonds[0].stereo, "dashed");
        assert_eq!(layout_output.bonds[0].direction, "none");
        assert!(!layout_output.bonds[0].forced_stereo);
    }

    #[test]
    fn forced_wavy_and_dashed_in_chain() {
        let layout_output = layout_native("CC!sO!dN").expect("wavy/dashed chain failed");
        assert_eq!(layout_output.bonds[0].stereo, "none");
        assert_eq!(layout_output.bonds[1].stereo, "wavy");
        assert_eq!(layout_output.bonds[2].stereo, "dashed");
    }

    #[test]
    fn forced_wavy_does_not_disturb_real_directional_bonds() {
        // A genuine trans alkene next to a forced wavy bond: the wavy marker
        // must not consume or shift the cis/trans direction tokens.
        let layout_output = layout_native("F/C=C/C!sN").expect("mixed directional/wavy failed");
        let wavy = layout_output
            .bonds
            .iter()
            .filter(|bond| bond.stereo == "wavy")
            .count();
        assert_eq!(wavy, 1);
        let directional = layout_output
            .bonds
            .iter()
            .filter(|bond| bond.direction != "none")
            .count();
        assert_eq!(directional, 2);
    }

    #[test]
    fn slash_without_double_bond_is_invalid() {
        let err = layout_native("C/N").expect_err("isolated slash should fail");
        assert!(err.contains("Directional"));
    }

    #[test]
    fn forced_wedge_in_chain() {
        let layout_output = layout_native("CC!wN").expect("forced wedge in chain failed");
        assert_eq!(layout_output.bonds[0].stereo, "none");
        assert_eq!(layout_output.bonds[1].stereo, "wedge_up");
        assert!(layout_output.bonds[1].forced_stereo);
    }

    #[test]
    fn alkene_directional_bonds_do_not_render_as_wedges() {
        let layout_output = layout_native("F/C=C/F").expect("E-alkene failed");
        let double = layout_output
            .bonds
            .iter()
            .find(|bond| bond.order == 2)
            .unwrap();
        assert_eq!(double.order, 2);
        assert!(layout_output.bonds.iter().all(|bond| bond.stereo == "none"));
        assert_eq!(
            layout_output
                .bonds
                .iter()
                .filter(|bond| bond.direction != "none")
                .count(),
            2
        );
    }

    #[test]
    fn trans_alkene_substituents_are_opposite() {
        let layout_output = layout_native("F/C=C/F").expect("trans alkene failed");
        assert_eq!(alkene_substituent_side_product(&layout_output), -1);
    }

    #[test]
    fn cis_alkene_substituents_are_same_side() {
        let layout_output = layout_native("F/C=C\\F").expect("cis alkene failed");
        assert_eq!(alkene_substituent_side_product(&layout_output), 1);
    }

    #[test]
    fn branch_direction_matches_opensmiles_examples() {
        let trans = layout_native("C(\\F)=C/F").expect("branch trans failed");
        let cis = layout_native("C(/F)=C/F").expect("branch cis failed");
        assert_eq!(alkene_substituent_side_product(&trans), -1);
        assert_eq!(alkene_substituent_side_product(&cis), 1);
    }

    #[test]
    fn pyrethroid_like_smiles_with_multiple_alkene_markers_parses() {
        let smiles = "CC1=C(C(=O)C[C@@H]1OC(=O)[C@@H]2[C@H](C2(C)C)/C=C(\\C)/C(=O)OC)C/C=C\\C=C";
        let layout_output = layout_native(smiles).expect("pyrethroid-like molecule failed");
        assert_eq!(
            layout_output
                .atoms
                .iter()
                .filter(|atom| !atom.virtual_h)
                .count(),
            27
        );
        assert_eq!(
            layout_output
                .bonds
                .iter()
                .filter(|bond| bond.direction != "none")
                .count(),
            5
        );
    }

    #[test]
    fn conflicting_cis_trans_markers_error() {
        let err = layout_native("C/C(\\F)=C/F").expect_err("conflicting markers should fail");
        assert!(err.contains("Conflicting"));
    }

    #[test]
    fn tetrahedral_chirality_adds_rendered_stereo() {
        let layout_output = layout_native("N[C@@H](C)C(=O)O").expect("chiral alanine failed");
        assert!(layout_output
            .atoms
            .iter()
            .any(|atom| atom.chirality == "tetra_clockwise"));
        // Exactly one wedge for the single stereocenter; direction checked below.
        let wedge_bonds = layout_output
            .bonds
            .iter()
            .filter(|bond| bond.stereo != "none")
            .count();
        let wedge_h = layout_output
            .atoms
            .iter()
            .filter(|atom| atom.stereo_h != "none")
            .count();
        assert_eq!(wedge_bonds + wedge_h, 1);
        assert!(chirality_matches_smiles("N[C@@H](C)C(=O)O"));
    }

    #[test]
    fn inverting_chirality_flips_the_wedge() {
        // Same skeleton/layout, opposite chirality token ⇒ the wedge on the chosen
        // bond must flip. Guards against regressing to a fixed @→up / @@→down map.
        let clockwise_layout = layout_native("N[C@@H](C)C(=O)O").expect("R failed");
        let anticlockwise_layout = layout_native("N[C@H](C)C(=O)O").expect("S failed");
        let stereo = |layout_output: &LayoutOutput| -> String {
            layout_output
                .bonds
                .iter()
                .find(|bond| bond.stereo != "none")
                .map(|bond| bond.stereo.clone())
                .or_else(|| {
                    layout_output
                        .atoms
                        .iter()
                        .find(|atom| atom.stereo_h != "none")
                        .map(|atom| atom.stereo_h.clone())
                })
                .unwrap_or_else(|| "none".to_string())
        };
        assert_ne!(stereo(&clockwise_layout), "none");
        assert_ne!(stereo(&anticlockwise_layout), "none");
        assert_ne!(stereo(&clockwise_layout), stereo(&anticlockwise_layout));
        assert!(chirality_matches_smiles("N[C@@H](C)C(=O)O"));
        assert!(chirality_matches_smiles("N[C@H](C)C(=O)O"));
    }

    #[test]
    fn tetrahedral_stereo_prefers_oh_substituent() {
        let layout_output = layout_native("CC[C@@H](O)CC/C=C/CO").expect("chiral alcohol failed");
        assert!(chiral_oxygen_bond_has_stereo(&layout_output));
    }

    #[test]
    fn chiral_alcohol_keeps_long_chain_as_continuation() {
        let layout_output = layout_native("CC[C@@H](O)CC/C=C/CO").expect("chiral alcohol failed");
        // The long chain leaves the stereocenter (atom 2) as a straight
        // continuation; the short ethyl tail (atom 0) sits on the opposite side.
        // Direction-agnostic so it survives the conventional horizontal mirror.
        let continuation_direction = layout_output.atoms[4].pos.x - layout_output.atoms[2].pos.x;
        assert!(continuation_direction.abs() > 1e-6);
        assert!(
            (layout_output.atoms[9].pos.x - layout_output.atoms[4].pos.x) * continuation_direction
                > 0.0,
            "terminal alcohol should remain on the zig-zag continuation"
        );
        assert!(
            (layout_output.atoms[0].pos.x - layout_output.atoms[2].pos.x) * continuation_direction
                < 0.0,
            "short ethyl branch should sit opposite the long continuation"
        );
        assert!(layout_output
            .atoms
            .iter()
            .all(|atom| atom.stereo_h == "none"));
    }

    #[test]
    fn steroid_stereo_prefers_exocyclic_oh_over_ring_bond() {
        let layout_output =
            layout_native("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O")
                .expect("steroid-like molecule layout failed");
        assert!(chiral_oxygen_bond_has_stereo(&layout_output));
    }

    #[test]
    fn steroid_ring_chiral_hydrogens_are_rendered() {
        let smiles = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O";
        let layout_output = layout_native(smiles).expect("steroid-like molecule layout failed");

        // The three ring-fusion stereocenters (no exocyclic substituent) render
        // an explicit wedge/hash hydrogen. The 17-OH carbon and the quaternary
        // C13 instead wedge their exocyclic substituent, so they get no H label.
        let stereo_h: Vec<&str> = layout_output
            .atoms
            .iter()
            .filter_map(|atom| (atom.stereo_h != "none").then_some(atom.stereo_h.as_str()))
            .collect();
        assert_eq!(stereo_h.len(), 3);

        // Adjacent ring-fusion stereocenters must point to opposite faces.
        assert_ne!(layout_output.atoms[5].stereo_h, "none");
        assert_ne!(layout_output.atoms[6].stereo_h, "none");
        assert_ne!(
            layout_output.atoms[5].stereo_h,
            layout_output.atoms[6].stereo_h
        );

        // The 17-OH bond is wedged (not the hydrogen).
        assert!(chiral_oxygen_bond_has_stereo(&layout_output));

        // Every stereocenter is depicted with the geometrically correct handedness.
        assert!(chirality_matches_smiles(smiles));
    }

    // ── Extended chirality classes and quadruple bonds ───────────────────────

    /// Cosine of the angle neighbor-a → center → neighbor-b.
    fn bond_angle_cos(
        layout_output: &LayoutOutput,
        center: usize,
        first_neighbor: usize,
        second_neighbor: usize,
    ) -> f64 {
        let center_position = layout_output.atoms[center].pos;
        let (first_offset_x, first_offset_y) = (
            layout_output.atoms[first_neighbor].pos.x - center_position.x,
            layout_output.atoms[first_neighbor].pos.y - center_position.y,
        );
        let (second_offset_x, second_offset_y) = (
            layout_output.atoms[second_neighbor].pos.x - center_position.x,
            layout_output.atoms[second_neighbor].pos.y - center_position.y,
        );
        (first_offset_x * second_offset_x + first_offset_y * second_offset_y)
            / ((first_offset_x * first_offset_x + first_offset_y * first_offset_y).sqrt()
                * (second_offset_x * second_offset_x + second_offset_y * second_offset_y).sqrt())
    }

    #[test]
    fn square_planar_cis_and_trans() {
        // Neighbors in writing order: N(0), N(2), Cl(3), Cl(4) around Pt(1).
        // @SP1 (U shape) puts the two Cl on adjacent corners: cis, 90° apart.
        let cis = layout_native("N[Pt@SP1](N)(Cl)Cl").expect("cisplatin failed");
        assert!(bond_angle_cos(&cis, 1, 3, 4).abs() < 1e-6);
        // @SP2 (4 shape) pairs Cl trans to Cl: 180° apart.
        let trans = layout_native("N[Pt@SP2](N)(Cl)Cl").expect("transplatin failed");
        assert!((bond_angle_cos(&trans, 1, 3, 4) + 1.0).abs() < 1e-6);
        // @SP3 (Z shape) pairs N1 trans to Cl4, so the Cl pair is cis again.
        let third_configuration = layout_native("N[Pt@SP3](N)(Cl)Cl").expect("SP3 failed");
        assert!(bond_angle_cos(&third_configuration, 1, 3, 4).abs() < 1e-6);
        assert!((bond_angle_cos(&third_configuration, 1, 0, 4) + 1.0).abs() < 1e-6);
    }

    #[test]
    fn square_planar_all_neighbors_at_right_angles() {
        let layout_output = layout_native("C[Fe@SP1](Cl)(Br)I").expect("SP iron failed");
        assert_eq!(layout_output.atoms[1].chirality, "square_planar");
        for (a, b) in [(0, 2), (2, 3), (3, 4)] {
            assert!(bond_angle_cos(&layout_output, 1, a, b).abs() < 1e-6);
        }
    }

    #[test]
    fn tb_oh_al_accepted_without_stereo_marks() {
        for smiles in [
            "S[As@TB1](F)(Cl)(Br)N",
            "C[Co@OH1](F)(Cl)(Br)(I)N",
            "NC(Br)=[C@AL1]=C(O)C",
        ] {
            let layout_output = layout_native(smiles).expect("extended chirality should parse");
            assert!(
                layout_output
                    .atoms
                    .iter()
                    .any(|atom| atom.chirality == "undepicted"),
                "missing undepicted center for {smiles}"
            );
            assert!(layout_output.bonds.iter().all(|bond| bond.stereo == "none"));
            assert!(layout_output
                .atoms
                .iter()
                .all(|atom| atom.stereo_h == "none"));
        }
    }

    #[test]
    fn quadruple_bond_order() {
        // The classic metal-metal quadruple bond, e.g. in [Re2Cl8]2-.
        let layout_output = layout_native("[Re]$[Re]").expect("quadruple bond failed");
        assert_eq!(layout_output.bonds[0].order, 4);
    }

    #[test]
    fn quadruple_bond_counts_toward_valence() {
        let layout_output = layout_native("C$C").expect("C$C failed");
        assert_eq!(layout_output.bonds[0].order, 4);
        assert!(layout_output.atoms.iter().all(|atom| atom.implicit_h == 0));
    }

    // ── Abbreviation tests ────────────────────────────────────────────────────

    #[test]
    fn preprocess_single_abbrev() {
        let preprocessed = preprocess_smiles("C({PPh3})=O").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "C([*])=O");
        assert_eq!(preprocessed.abbrev_labels[0].text, "PPh3");
        assert_eq!(preprocessed.abbrev_labels[0].style, "");
        assert_eq!(preprocessed.abbrev_labels[0].anchor_len, 0);
    }

    #[test]
    fn preprocess_multiple_abbrevs() {
        let preprocessed = preprocess_smiles("{OEt}C(=O){NHR}").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "[*]C(=O)[*]");
        assert_eq!(preprocessed.abbrev_labels[0].text, "OEt");
        assert_eq!(preprocessed.abbrev_labels[1].text, "NHR");
    }

    #[test]
    fn preprocess_abbrev_style() {
        let preprocessed = preprocess_smiles("{PPh3|P}C({LG|red})=O").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "[*]C([*])=O");
        assert_eq!(preprocessed.abbrev_labels[0].text, "PPh3");
        assert_eq!(preprocessed.abbrev_labels[0].style, "P");
        assert_eq!(preprocessed.abbrev_labels[1].text, "LG");
        assert_eq!(preprocessed.abbrev_labels[1].style, "red");
    }

    #[test]
    fn preprocess_abbrev_anchor_positions() {
        let preprocessed = preprocess_smiles("{>CAT}C({C>AT}){CA>T}").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "[*]C([*])[*]");
        assert_eq!(preprocessed.abbrev_labels[0].text, "CAT");
        assert_eq!(preprocessed.abbrev_labels[0].anchor, 0);
        assert_eq!(preprocessed.abbrev_labels[0].anchor_len, 1);
        assert_eq!(preprocessed.abbrev_labels[1].text, "CAT");
        assert_eq!(preprocessed.abbrev_labels[1].anchor, 1);
        assert_eq!(preprocessed.abbrev_labels[1].anchor_len, 1);
        assert_eq!(preprocessed.abbrev_labels[2].text, "CAT");
        assert_eq!(preprocessed.abbrev_labels[2].anchor, 2);
        assert_eq!(preprocessed.abbrev_labels[2].anchor_len, 1);
    }

    #[test]
    fn preprocess_rejects_arrow_marker_outside_abbrev() {
        assert!(preprocess_smiles("C>C").is_err());
        assert!(preprocess_smiles("{CAT>}C").is_err());
    }

    #[test]
    fn malformed_delimiters_report_actionable_errors() {
        let unclosed_label = layout_native("C{OH").expect_err("unclosed label should fail");
        assert!(unclosed_label.contains("unclosed custom label"));
        assert!(unclosed_label.contains("add `}`"));

        let unclosed_bracket = layout_native("C[NH2").expect_err("unclosed bracket should fail");
        assert!(unclosed_bracket.contains("unclosed bracket atom"));
        assert!(unclosed_bracket.contains("add `]`"));

        let unmatched_label_end =
            layout_native("COH}").expect_err("unmatched label end should fail");
        assert!(unmatched_label_end.contains("unmatched `}`"));
    }

    #[test]
    fn empty_smiles_reports_actionable_error() {
        let error = layout_native("   ").expect_err("empty SMILES should fail");
        assert!(error.contains("expression is empty"));
        assert!(error.contains("at least one atom"));
    }

    #[test]
    fn unclosed_ring_reports_ring_number_and_correction() {
        let error = layout_native("C1CC").expect_err("unclosed ring should fail");
        assert!(error.contains("ring closure 1"));
        assert!(error.contains("never closed"));
        assert!(error.contains("repeat each ring number"));
    }

    #[test]
    fn unknown_abbreviation_style_reports_supported_forms() {
        let error = layout_native("C{LG|chartreuse}").expect_err("unknown style should fail");
        assert!(error.contains("unknown abbreviation style `chartreuse`"));
        assert!(error.contains("element symbol"));
        assert!(error.contains("#RRGGBB"));
    }

    #[test]
    fn preprocess_forced_wedge_markers() {
        let preprocessed = preprocess_smiles("C!wN!hO").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "C/N/O");
        assert_eq!(
            preprocessed.forced_direction_markers,
            vec![
                BondMarker {
                    style: BondMarkerStyle::WedgeUp,
                    order: BondOrder::Single,
                    curl: false,
                },
                BondMarker {
                    style: BondMarkerStyle::WedgeDown,
                    order: BondOrder::Single,
                    curl: false,
                },
            ]
        );
    }

    #[test]
    fn preprocess_curl_markers_and_combinations() {
        let preprocessed = preprocess_smiles("CCC!cC").expect("plain curl failed");
        assert_eq!(preprocessed.smiles, "CCC/C");
        assert_eq!(
            preprocessed.forced_direction_markers,
            vec![BondMarker {
                style: BondMarkerStyle::Plain,
                order: BondOrder::Single,
                curl: true,
            }]
        );

        let preprocessed = preprocess_smiles("CCC!c!wC").expect("curl wedge failed");
        assert_eq!(preprocessed.smiles, "CCC/C");
        assert_eq!(
            preprocessed.forced_direction_markers[0].style,
            BondMarkerStyle::WedgeUp
        );
        assert!(preprocessed.forced_direction_markers[0].curl);

        let preprocessed = preprocess_smiles("CCC!c=C").expect("curl double failed");
        assert_eq!(preprocessed.smiles, "CCC/C");
        assert_eq!(
            preprocessed.forced_direction_markers[0].order,
            BondOrder::Double
        );
        assert!(preprocessed.forced_direction_markers[0].curl);
    }

    #[test]
    fn preprocess_moves_ring_closures_written_after_branches() {
        let preprocessed = preprocess_smiles("C1=CCCC(=O)1").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "C1=CCCC1(=O)");

        let preprocessed = preprocess_smiles("C(=O)(O)1N").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "C1(=O)(O)N");
    }

    #[test]
    fn preprocess_no_abbrev() {
        let preprocessed = preprocess_smiles("CCO").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "CCO");
        assert!(preprocessed.abbrev_labels.is_empty());
        assert_eq!(
            preprocessed.aromatic_atom_markers,
            vec![false, false, false]
        );
    }

    #[test]
    fn preprocess_uppercases_aromatic_atoms() {
        let preprocessed = preprocess_smiles("Clc1ccccc1").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "ClC1CCCCC1");
        // One marker per unbracketed organic-subset atom, in writing order;
        // the two-letter Cl counts once.
        assert_eq!(
            preprocessed.aromatic_atom_markers,
            vec![false, true, true, true, true, true, true]
        );
    }

    #[test]
    fn preprocess_leaves_bracket_atoms_alone() {
        // Bracket contents are the parser's business: [nH] parses natively as
        // an aromatic atom, and the 'c' in [Sc] is not an aromatic carbon.
        let preprocessed = preprocess_smiles("c1cc[nH]c1[Sc]").expect("preprocess failed");
        assert_eq!(preprocessed.smiles, "C1CC[nH]C1[Sc]");
        assert_eq!(
            preprocessed.aromatic_atom_markers,
            vec![true, true, true, true]
        );
    }

    #[test]
    fn abbrev_assigned_in_layout() {
        let layout_output = layout_native("{PPh3}C=O").expect("abbrev layout failed");
        assert_eq!(layout_output.atoms.len(), 3);
        assert_eq!(layout_output.atoms[0].abbrev, "PPh3");
        assert_eq!(layout_output.atoms[1].abbrev, "");
    }

    #[test]
    fn abbrev_multiple_in_layout() {
        // [*]C(=O)[*] → atoms: 0=[*](OEt), 1=C, 2=O (branch), 3=[*](NHR)
        let layout_output = layout_native("{OEt}C(=O){NHR}").expect("multi abbrev layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "OEt");
        assert_eq!(layout_output.atoms[3].abbrev, "NHR");
    }

    #[test]
    fn abbrev_style_assigned_in_layout() {
        let layout_output = layout_native("{PPh3|P}C=O").expect("styled abbrev layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "PPh3");
        assert_eq!(layout_output.atoms[0].abbrev_style, "P");
    }

    #[test]
    fn abbrev_anchor_assigned_in_layout() {
        let layout_output = layout_native("{>PPh3}C=O").expect("anchored abbrev layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "PPh3");
        assert_eq!(layout_output.atoms[0].abbrev_anchor, 0);
        assert_eq!(layout_output.atoms[0].abbrev_anchor_len, 1);
    }

    // ── Abbreviation lone-pair modifier (`lp=N`) ─────────────────────────────

    #[test]
    fn abbrev_lp_parses_all_counts() {
        for n in 1u8..=4 {
            let raw = format!("Cl|Cl|lp={n}");
            let label = parse_abbreviation_label(&raw).expect("lp parse failed");
            assert_eq!(label.text, "Cl");
            assert_eq!(label.style, "Cl");
            assert_eq!(label.lone_pairs, Some(n));
        }
    }

    #[test]
    fn abbrev_lp_without_style() {
        let label = parse_abbreviation_label("OR|lp=2").expect("lp without style failed");
        assert_eq!(label.text, "OR");
        assert_eq!(label.style, "");
        assert_eq!(label.lone_pairs, Some(2));
    }

    #[test]
    fn abbrev_lp_with_anchor_marker() {
        let label = parse_abbreviation_label(">PPh_3|P|lp=1").expect("anchored lp failed");
        assert_eq!(label.text, "PPh_3");
        assert_eq!(label.anchor, 0);
        assert_eq!(label.anchor_len, 1);
        assert_eq!(label.style, "P");
        assert_eq!(label.lone_pairs, Some(1));
    }

    #[test]
    fn abbrev_lp_with_anchor_no_style() {
        let label = parse_abbreviation_label(">OR|lp=2").expect("anchored lp no style failed");
        assert_eq!(label.text, "OR");
        assert_eq!(label.anchor_len, 1);
        assert_eq!(label.style, "");
        assert_eq!(label.lone_pairs, Some(2));
    }

    #[test]
    fn abbrev_existing_syntax_preserved() {
        assert_eq!(parse_abbreviation_label("PPh3").unwrap().lone_pairs, None);
        assert_eq!(parse_abbreviation_label("PPh3|P").unwrap().style, "P");
        assert_eq!(parse_abbreviation_label("PPh3|P").unwrap().lone_pairs, None);
        let anchored = parse_abbreviation_label(">PPh3").unwrap();
        assert_eq!(anchored.anchor_len, 1);
        assert_eq!(anchored.lone_pairs, None);
        let anchored_styled = parse_abbreviation_label(">PPh3|red").unwrap();
        assert_eq!(anchored_styled.style, "red");
        assert_eq!(anchored_styled.lone_pairs, None);
    }

    #[test]
    fn abbrev_lp_rejects_out_of_range_and_malformed() {
        assert!(parse_abbreviation_label("Cl|Cl|lp=0").is_err());
        assert!(parse_abbreviation_label("Cl|Cl|lp=5").is_err());
        assert!(parse_abbreviation_label("Cl|Cl|lp=-1").is_err());
        assert!(parse_abbreviation_label("Cl|Cl|lp=2.5").is_err());
        assert!(parse_abbreviation_label("Cl|Cl|lp=two").is_err());
        assert!(parse_abbreviation_label("Cl|Cl|lp=").is_err());
        // Duplicate lp modifiers.
        assert!(parse_abbreviation_label("Cl|Cl|lp=1|lp=2").is_err());
        // Unknown named modifier.
        assert!(parse_abbreviation_label("Cl|Cl|foo=1").is_err());
        // Style after a named modifier.
        assert!(parse_abbreviation_label("Cl|lp=1|Cl").is_err());
    }

    #[test]
    fn abbrev_lp_count_reaches_layout_output() {
        let layout_output = layout_native("{>Cl|Cl|lp=3}C").expect("lp layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "Cl");
        assert_eq!(layout_output.atoms[0].lone_pairs, 3);
        // A fallback direction record is emitted per declared pair so lp()
        // references stay resolvable.
        assert_eq!(layout_output.atoms[0].lone_pair_dirs.len(), 3);
    }

    // ── Abbreviation rendering offset (`offset=(x,y)`) ─────────────────────

    #[test]
    fn abbrev_offset_parses_with_style_and_lone_pairs() {
        let label = parse_abbreviation_label("H|grey|lp=1|offset=(0.1, -0.2)")
            .expect("offset parse failed");
        assert_eq!(label.text, "H");
        assert_eq!(label.style, "grey");
        assert_eq!(label.lone_pairs, Some(1));
        assert_eq!(label.offset, Some((0.1, -0.2)));

        let reversed = parse_abbreviation_label("H|offset=(-.25, .4)|lp=2")
            .expect("reordered modifiers failed");
        assert_eq!(reversed.style, "");
        assert_eq!(reversed.lone_pairs, Some(2));
        assert_eq!(reversed.offset, Some((-0.25, 0.4)));
    }

    #[test]
    fn abbrev_offset_rejects_malformed_values() {
        assert!(parse_abbreviation_label("H|offset=0.1,0.2").is_err());
        assert!(parse_abbreviation_label("H|offset=(0.1)").is_err());
        assert!(parse_abbreviation_label("H|offset=(0.1,0.2,0.3)").is_err());
        assert!(parse_abbreviation_label("H|offset=(left,0.2)").is_err());
        assert!(parse_abbreviation_label("H|offset=(NaN,0.2)").is_err());
        assert!(parse_abbreviation_label("H|offset=(0.1,0.2)|offset=(0.3,0.4)").is_err());
    }

    #[test]
    fn abbrev_offset_reaches_output_without_changing_layout_coordinates() {
        let baseline = layout_native("C{H}").expect("baseline layout failed");
        let displaced = layout_native("C{H|offset=(0.3,-0.2)}").expect("offset layout failed");

        assert_eq!(baseline.atoms[0].pos.x, displaced.atoms[0].pos.x);
        assert_eq!(baseline.atoms[0].pos.y, displaced.atoms[0].pos.y);
        assert_eq!(baseline.atoms[1].pos.x, displaced.atoms[1].pos.x);
        assert_eq!(baseline.atoms[1].pos.y, displaced.atoms[1].pos.y);
        assert_eq!(baseline.bbox_width, displaced.bbox_width);
        assert_eq!(baseline.bbox_height, displaced.bbox_height);
        assert_eq!(displaced.atoms[1].abbrev_offset_x, 0.3);
        assert_eq!(displaced.atoms[1].abbrev_offset_y, -0.2);
    }

    #[test]
    fn abbrev_without_lp_has_no_pairs() {
        let layout_output = layout_native("{PPh3}C=O").expect("no-lp layout failed");
        assert_eq!(layout_output.atoms[0].abbrev, "PPh3");
        assert_eq!(layout_output.atoms[0].lone_pairs, 0);
        assert_eq!(layout_output.atoms[0].lone_pair_dirs.len(), 0);
    }

    #[test]
    fn bracket_atom_is_not_literal_label() {
        let bracket = layout_native("[N]").expect("bracket nitrogen failed");
        let label = layout_native("{N}").expect("literal label failed");
        assert_eq!(bracket.atoms[0].symbol, "N");
        assert_eq!(bracket.atoms[0].abbrev, "");
        assert_eq!(label.atoms[0].symbol, "*");
        assert_eq!(label.atoms[0].abbrev, "N");
    }

    // ── Virtual H atoms for bracket-notation hydrogens ───────────────────────

    #[test]
    fn bracket_h_group_is_one_addressable_index() {
        // [OH-]: one O heavy atom + one virtual H group → 2 atoms total.
        let layout_output = layout_native("[OH-]").expect("[OH-] layout failed");
        assert_eq!(layout_output.atoms.len(), 2);
        assert_eq!(layout_output.atoms[0].symbol, "O");
        assert!(layout_output.atoms[1].virtual_h);
        assert_eq!(layout_output.atoms[1].symbol, "H");
        assert_eq!(layout_output.bonds.len(), 1);
        assert!(layout_output.bonds[0].virtual_bond);
    }

    #[test]
    fn ammonium_h_group_is_one_index() {
        // [NH4+]: despite hcount=4 the H-label is one glyph → exactly 1 virtual H.
        let layout_output = layout_native("[NH4+]").expect("[NH4+] layout failed");
        assert_eq!(layout_output.atoms.len(), 2); // 1 N + 1 virtual H group
        assert!(layout_output.atoms[1].virtual_h);
        assert_eq!(layout_output.bonds.len(), 1);
        assert!(layout_output.bonds[0].virtual_bond);
    }

    #[test]
    fn virtual_h_has_valid_position() {
        let layout_output = layout_native("[OH-]").expect("[OH-] layout failed");
        let oxygen_position = layout_output.atoms[0].pos;
        let hydrogen_position = layout_output.atoms[1].pos;
        let distance = ((hydrogen_position.x - oxygen_position.x).powi(2)
            + (hydrogen_position.y - oxygen_position.y).powi(2))
        .sqrt();
        assert!(
            (distance - 0.35).abs() < 1e-6,
            "H should be 0.35 bond lengths from O, got {distance}"
        );
    }

    // ── Implicit H for expanded valence table ────────────────────────────────

    #[test]
    fn implicit_h_boron() {
        // B (organic subset, valence 3): B bonded to one C → 2 implicit H
        let layout_output = layout_native("BC").expect("boron layout failed");
        assert_eq!(layout_output.atoms[0].implicit_h, 2);
    }

    #[test]
    fn lone_pair_counts_common_organic_atoms() {
        let alcohol = layout_native("CCO").expect("alcohol layout failed");
        assert_eq!(alcohol.atoms[2].lone_pairs, 2);
        assert_eq!(alcohol.atoms[2].lone_pair_dirs.len(), 2);

        let amine = layout_native("CCN").expect("amine layout failed");
        assert_eq!(amine.atoms[2].lone_pairs, 1);

        let chloride = layout_native("CCl").expect("chloride layout failed");
        assert_eq!(chloride.atoms[1].lone_pairs, 3);
    }

    #[test]
    fn lone_pair_counts_respect_charge_and_explicit_hydrogens() {
        let ammonium = layout_native("[NH4+]").expect("ammonium layout failed");
        assert_eq!(ammonium.atoms[0].lone_pairs, 0);

        let formate = layout_native("[O-]C=O").expect("formate layout failed");
        assert_eq!(formate.atoms[0].lone_pairs, 3);
        assert_eq!(formate.atoms[2].lone_pairs, 2);
    }

    fn max_bond_length(layout_output: &LayoutOutput) -> f64 {
        layout_output
            .bonds
            .iter()
            .map(|bond| {
                let from = layout_output.atoms[bond.from].pos;
                let to = layout_output.atoms[bond.to].pos;
                from.dist(to)
            })
            .fold(0.0, f64::max)
    }

    /// Smallest distance between any two atoms that are not bonded to each other.
    /// A value well below the ~1.0 bond length signals two parts of the molecule
    /// being laid out on top of one another.
    fn min_nonbonded_distance(layout_output: &LayoutOutput) -> f64 {
        let bonded: std::collections::HashSet<(usize, usize)> = layout_output
            .bonds
            .iter()
            .map(|b| (b.from.min(b.to), b.from.max(b.to)))
            .collect();
        let mut min = f64::INFINITY;
        for i in 0..layout_output.atoms.len() {
            for j in (i + 1)..layout_output.atoms.len() {
                if bonded.contains(&(i, j)) {
                    continue;
                }
                min = min.min(layout_output.atoms[i].pos.dist(layout_output.atoms[j].pos));
            }
        }
        min
    }

    fn atoms_are_collinear(
        layout_output: &LayoutOutput,
        first_atom: usize,
        middle_atom: usize,
        last_atom: usize,
    ) -> bool {
        let first_position = layout_output.atoms[first_atom].pos;
        let middle_position = layout_output.atoms[middle_atom].pos;
        let last_position = layout_output.atoms[last_atom].pos;
        let cross = (middle_position.x - first_position.x) * (last_position.y - middle_position.y)
            - (middle_position.y - first_position.y) * (last_position.x - middle_position.x);
        cross.abs() < 1e-8
    }

    /// End-to-end check: reconstructs the depicted 3D geometry from the rendered
    /// output and verifies every stereocenter's signed volume matches `@`/`@@`.
    fn chirality_matches_smiles(smiles: &str) -> bool {
        use crate::layout::{implicit_h_count, signed_volume, stereochemical_hydrogen_direction};

        let preprocessed = preprocess_smiles(smiles).expect("preprocess failed");
        let molecule = MoleculeGraph::from_smiles(
            &preprocessed.smiles,
            preprocessed.forced_direction_markers,
            preprocessed.aromatic_atom_markers,
        )
        .expect("graph build failed");
        let layout_output = compute_layout(&molecule).expect("layout failed");
        let coordinates: Vec<crate::render::Vec2> =
            layout_output.atoms.iter().map(|atom| atom.pos).collect();

        for center in 0..molecule.atoms.len() {
            let parity = match layout_output.atoms[center].chirality.as_str() {
                "tetra_anti" => -1.0_f64,
                "tetra_clockwise" => 1.0_f64,
                _ => continue,
            };
            let neighbor_bonds = &molecule.neighbor_bonds[center];
            let hydrogen_count =
                (molecule.atoms[center].hcount + implicit_h_count(&molecule, center)) as usize;
            if neighbor_bonds.len() + hydrogen_count != 4 || hydrogen_count > 1 {
                continue;
            }

            // Neighbor order with the implicit hydrogen inserted (same rule as the
            // renderer): after the "from" atom, or first if there is none.
            #[derive(Clone, Copy)]
            enum TetrahedralNeighbor {
                Bond(usize, usize),
                Hydrogen,
            }
            let mut neighbor_order: Vec<TetrahedralNeighbor> = neighbor_bonds
                .iter()
                .map(|&bond_index| {
                    let bond = &molecule.bonds[bond_index];
                    let other = if bond.from == center {
                        bond.to
                    } else {
                        bond.from
                    };
                    TetrahedralNeighbor::Bond(bond_index, other)
                })
                .collect();
            if hydrogen_count == 1 {
                let insertion_position = if molecule.has_preceding[center] { 1 } else { 0 };
                neighbor_order.insert(
                    insertion_position.min(neighbor_order.len()),
                    TetrahedralNeighbor::Hydrogen,
                );
            }

            // Which bond (if any) carries the rendered wedge, and its z sign.
            let wedge_bond = neighbor_bonds
                .iter()
                .find(|&&bond_index| layout_output.bonds[bond_index].stereo != "none")
                .copied();
            let hydrogen_direction =
                stereochemical_hydrogen_direction(&molecule, center, &coordinates, None);

            let normalize = |offset_x: f64, offset_y: f64| {
                let length = (offset_x * offset_x + offset_y * offset_y).sqrt();
                if length > 1e-12 {
                    (offset_x / length, offset_y / length)
                } else {
                    (offset_x, offset_y)
                }
            };
            let stereo_depth = |stereo: &str| if stereo == "wedge_up" { 1.0 } else { -1.0 };

            let mut directions = [[0.0_f64; 3]; 4];
            for (neighbor_index, neighbor) in neighbor_order.iter().enumerate() {
                directions[neighbor_index] = match neighbor {
                    TetrahedralNeighbor::Bond(bond_index, other) => {
                        let (direction_x, direction_y) = normalize(
                            coordinates[*other].x - coordinates[center].x,
                            coordinates[*other].y - coordinates[center].y,
                        );
                        let depth = if Some(*bond_index) == wedge_bond {
                            stereo_depth(&layout_output.bonds[*bond_index].stereo)
                        } else {
                            0.0
                        };
                        [direction_x, direction_y, depth]
                    }
                    TetrahedralNeighbor::Hydrogen => {
                        // If the H itself is wedged, use that; otherwise it sits on
                        // the face opposite the wedged heavy substituent.
                        let depth = if layout_output.atoms[center].stereo_h != "none" {
                            stereo_depth(&layout_output.atoms[center].stereo_h)
                        } else if let Some(bond_index) = wedge_bond {
                            -stereo_depth(&layout_output.bonds[bond_index].stereo)
                        } else {
                            0.0
                        };
                        let direction = if layout_output.atoms[center].stereo_h != "none" {
                            layout_output.atoms[center].stereo_h_dir
                        } else {
                            hydrogen_direction
                        };
                        [direction.x, direction.y, depth]
                    }
                };
            }

            let volume = signed_volume(&directions);
            if volume.abs() < 1e-9 || volume.signum() != parity {
                return false;
            }
        }
        true
    }

    fn chiral_oxygen_bond_has_stereo(layout_output: &LayoutOutput) -> bool {
        layout_output.bonds.iter().any(|bond| {
            let from = &layout_output.atoms[bond.from];
            let to = &layout_output.atoms[bond.to];
            bond.stereo != "none"
                && ((from.chirality != "none" && to.symbol == "O")
                    || (to.chirality != "none" && from.symbol == "O"))
        })
    }

    fn alkene_substituent_side_product(layout_output: &LayoutOutput) -> i8 {
        let double = layout_output
            .bonds
            .iter()
            .find(|bond| bond.order == 2)
            .unwrap();
        let first_alkene_atom = double.from;
        let second_alkene_atom = double.to;
        let left = layout_output
            .bonds
            .iter()
            .find(|bond| {
                bond.order == 1 && (bond.from == first_alkene_atom || bond.to == first_alkene_atom)
            })
            .unwrap();
        let right = layout_output
            .bonds
            .iter()
            .find(|bond| {
                bond.order == 1
                    && (bond.from == second_alkene_atom || bond.to == second_alkene_atom)
            })
            .unwrap();
        let left_neighbor = if left.from == first_alkene_atom {
            left.to
        } else {
            left.from
        };
        let right_neighbor = if right.from == second_alkene_atom {
            right.to
        } else {
            right.from
        };
        side(
            layout_output.atoms[first_alkene_atom].pos,
            layout_output.atoms[second_alkene_atom].pos,
            layout_output.atoms[left_neighbor].pos,
        ) * side(
            layout_output.atoms[first_alkene_atom].pos,
            layout_output.atoms[second_alkene_atom].pos,
            layout_output.atoms[right_neighbor].pos,
        )
    }

    fn side(
        line_start: crate::render::Vec2,
        line_end: crate::render::Vec2,
        point: crate::render::Vec2,
    ) -> i8 {
        let cross = (line_end.x - line_start.x) * (point.y - line_start.y)
            - (line_end.y - line_start.y) * (point.x - line_start.x);
        if cross > 1e-8 {
            1
        } else {
            -1
        }
    }

    fn turn_cross(a: crate::render::Vec2, b: crate::render::Vec2, c: crate::render::Vec2) -> f64 {
        (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x)
    }

    #[test]
    fn curl_repeats_the_previous_chain_turn() {
        let normal = layout_native("CCCC").expect("normal chain failed");
        let curled = layout_native("CCC!cC").expect("curled chain failed");

        let normal_first = turn_cross(
            normal.atoms[0].pos,
            normal.atoms[1].pos,
            normal.atoms[2].pos,
        );
        let normal_second = turn_cross(
            normal.atoms[1].pos,
            normal.atoms[2].pos,
            normal.atoms[3].pos,
        );
        assert!(normal_first * normal_second < 0.0);

        let curl_first = turn_cross(
            curled.atoms[0].pos,
            curled.atoms[1].pos,
            curled.atoms[2].pos,
        );
        let curl_second = turn_cross(
            curled.atoms[1].pos,
            curled.atoms[2].pos,
            curled.atoms[3].pos,
        );
        assert!(curl_first * curl_second > 0.0);
    }

    #[test]
    fn curl_swaps_forward_branch_slots_without_overlap() {
        let layout_output = layout_native("CCC({PPh3})!cC(=O)OCC").expect("branched curl failed");
        let first = turn_cross(
            layout_output.atoms[0].pos,
            layout_output.atoms[1].pos,
            layout_output.atoms[2].pos,
        );
        let second = turn_cross(
            layout_output.atoms[1].pos,
            layout_output.atoms[2].pos,
            layout_output.atoms[4].pos,
        );
        assert!(first * second > 0.0);
        assert!(min_atom_distance(&layout_output) >= 0.5);
    }

    #[test]
    fn consecutive_curls_repeat_each_new_turn() {
        let layout_output = layout_native("CCCC!cC!cC").expect("consecutive curls failed");
        let turns = (0..3)
            .map(|i| {
                turn_cross(
                    layout_output.atoms[i + 1].pos,
                    layout_output.atoms[i + 2].pos,
                    layout_output.atoms[i + 3].pos,
                )
            })
            .collect::<Vec<_>>();
        assert!(turns.iter().all(|turn| turns[0] * turn > 0.0));
        assert!(min_atom_distance(&layout_output) >= 0.5);
    }

    #[test]
    fn crowded_substituted_curls_keep_atoms_separated() {
        let molecules = [
            "CCC(C(N)C)!cC(OC(F)C)!cCC",
            "CCC({PPh3})!cC([O-])!cC(=O)OCC",
            "CCC(CC(C)C)!cC(OC)C(NC)CC",
        ];

        for smiles in molecules {
            let layout_output =
                layout_native(smiles).unwrap_or_else(|err| panic!("{smiles}: {err}"));
            let min = min_atom_distance(&layout_output);
            assert!(min >= 0.5, "{smiles}: min atom distance {min:.3}");
        }
    }

    #[test]
    fn curl_combines_with_wedge_and_double_bonds() {
        let wedge = layout_native("CCC!c!wN").expect("curled wedge failed");
        assert_eq!(wedge.bonds[2].stereo, "wedge_up");

        let double = layout_native("CCC!c=C").expect("curled double failed");
        assert_eq!(double.bonds[2].order, 2);
        let first = turn_cross(
            double.atoms[0].pos,
            double.atoms[1].pos,
            double.atoms[2].pos,
        );
        let second = turn_cross(
            double.atoms[1].pos,
            double.atoms[2].pos,
            double.atoms[3].pos,
        );
        assert!(first * second > 0.0);
    }

    #[test]
    fn curl_requires_an_established_turn() {
        let err = layout_native("CC!cC").expect_err("early curl should fail");
        assert!(err.contains("two preceding chain bonds"));
    }

    // ── Molecular weight ──────────────────────────────────────────────────────
    //
    // Reference values are PubChem's computed molecular weights, which use the
    // IUPAC/CIAAW standard atomic weights (the same table ptable embeds via
    // PubChemElements_all.json).

    fn assert_weight(smiles: &str, expected: f64) {
        let molecular_weight = mol_weight_native(smiles).expect("mol weight failed");
        assert!(
            (molecular_weight - expected).abs() < 0.01,
            "mol_weight({smiles}) = {molecular_weight}, expected {expected}"
        );
    }

    /// Extracts every SMILES string literal passed to `smiles("...")` or
    /// `mol("...")` in the visual test file, undoing Typst string escapes.
    fn test_typ_smiles_strings() -> Vec<String> {
        let src =
            std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/test.typ"))
                .expect("tests/test.typ not found");
        let src: String = src
            .lines()
            .filter(|l| !l.trim_start().starts_with("//"))
            .collect::<Vec<_>>()
            .join("\n");
        let mut found = Vec::new();
        for opener in ["smiles(\"", "mol(\""] {
            let mut rest = src.as_str();
            while let Some(pos) = rest.find(opener) {
                rest = &rest[pos + opener.len()..];
                let mut literal = String::new();
                let mut chars = rest.chars();
                while let Some(ch) = chars.next() {
                    match ch {
                        '"' => break,
                        '\\' => {
                            if let Some(esc) = chars.next() {
                                literal.push(esc);
                            }
                        }
                        _ => literal.push(ch),
                    }
                }
                found.push(literal);
            }
        }
        found.sort();
        found.dedup();
        found
    }

    #[test]
    fn every_molecule_in_test_typ_has_no_overlapping_atoms() {
        let molecules = test_typ_smiles_strings();
        assert!(
            molecules.len() > 50,
            "extraction looks broken: {molecules:?}"
        );
        let mut failures = Vec::new();
        for m in &molecules {
            match layout_native(m) {
                Ok(layout_output) => {
                    let min = min_atom_distance(&layout_output);
                    if layout_output.atoms.len() > 1 && min < 0.5 {
                        failures.push(format!("{m}: min atom distance {min:.3}"));
                    }
                }
                Err(e) => failures.push(format!("{m}: layout failed: {e}")),
            }
        }
        assert!(
            failures.is_empty(),
            "molecules with overlaps or errors:\n{}",
            failures.join("\n")
        );
    }

    // ── Aromatic ring circles ────────────────────────────────────────────────

    #[test]
    fn aromatic_input_emits_ring_circles() {
        let layout_output = layout_native("c1ccccc1").expect("benzene layout failed");
        assert_eq!(layout_output.aromatic_rings.len(), 1);
        let ring = &layout_output.aromatic_rings[0];
        // Hexagon with unit bonds: inradius ~0.866, so radius ~0.62.
        assert!((ring.radius - 0.866 * 0.72).abs() < 0.05);
        assert!(
            layout_output
                .bonds
                .iter()
                .filter(|bond| bond.aromatic)
                .count()
                == 6
        );
    }

    #[test]
    fn kekule_input_emits_no_ring_circles() {
        let layout_output = layout_native("C1=CC=CC=C1").expect("benzene layout failed");
        assert!(layout_output.aromatic_rings.is_empty());
        assert!(layout_output.bonds.iter().all(|b| !b.aromatic));
    }

    #[test]
    fn fused_aromatics_emit_one_circle_per_ring() {
        let layout_output = layout_native("c1ccc2ccccc2c1").expect("naphthalene layout failed");
        assert_eq!(layout_output.aromatic_rings.len(), 2);
    }

    #[test]
    fn aromatic_ring_with_saturated_neighbor_ring() {
        // Indane: only the aromatic ring gets a circle.
        let layout_output = layout_native("c1ccc2CCCc2c1").expect("indane layout failed");
        assert_eq!(layout_output.aromatic_rings.len(), 1);
    }

    #[test]
    fn mol_weight_water() {
        // PubChem CID 962: 18.015 g/mol
        assert_weight("O", 18.015);
    }

    #[test]
    fn mol_weight_ethanol() {
        // PubChem CID 702: 46.07 g/mol
        assert_weight("CCO", 46.069);
    }

    #[test]
    fn mol_weight_benzene_aromatic_input() {
        // PubChem CID 241: 78.11 g/mol; aromatic input exercises kekulization.
        assert_weight("c1ccccc1", 78.114);
    }

    #[test]
    fn mol_weight_glucose() {
        // PubChem CID 5793: 180.16 g/mol
        assert_weight("C(C1C(C(C(C(O1)O)O)O)O)O", 180.156);
    }

    #[test]
    fn mol_weight_caffeine() {
        // PubChem CID 2519: 194.19 g/mol
        assert_weight("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 194.19);
    }

    #[test]
    fn mol_weight_sodium_acetate_dot_fragments() {
        // PubChem CID 517045: 82.03 g/mol; dot-separated ion pair sums both
        // fragments (the electron mass difference of the ions is ignored, as
        // in standard formula-weight arithmetic).
        assert_weight("CC(=O)[O-].[Na+]", 82.034);
    }

    #[test]
    fn mol_weight_ammonium_bracket_h() {
        // PubChem CID 223: 18.039 g/mol; explicit bracket hydrogens counted.
        assert_weight("[NH4+]", 18.039);
    }

    #[test]
    fn mol_weight_wildcard_errors() {
        let err = mol_weight_native("*CC").expect_err("wildcard should fail");
        assert!(err.contains("wildcard"));
    }

    #[test]
    fn mol_weight_abbreviation_errors() {
        let err = mol_weight_native("{PPh3}C=O").expect_err("abbreviation should fail");
        assert!(err.contains("PPh3"));
    }

    #[test]
    fn mol_weight_isotope_errors() {
        let err = mol_weight_native("[2H]O[2H]").expect_err("isotope should fail");
        assert!(err.contains("isotope"));
    }
}
