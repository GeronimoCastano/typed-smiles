use std::fmt;

#[derive(Debug)]
pub enum SmilesError {
    Parse(String),
    Layout(String),
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(msg) => write!(f, "parse error: {msg}"),
            Self::Layout(msg) => write!(f, "layout error: {msg}"),
        }
    }
}

impl From<SmilesError> for String {
    fn from(e: SmilesError) -> Self {
        e.to_string()
    }
}
