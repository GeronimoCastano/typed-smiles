/// Parses the body of a `{label}` abbreviation: the displayed text with an
/// optional `>` attachment marker, an optional color or element style, and
/// the named `lp=N` and `offset=(x,y)` rendering modifiers.

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
pub(crate) fn parse_abbreviation_label(raw_label: &str) -> Result<AbbreviationLabel, String> {
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
    if valid_hex || NAMED_COLORS.contains(&style) || crate::element_from_symbol(style).is_some() {
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
