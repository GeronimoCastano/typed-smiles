use std::fmt;

#[derive(Debug)]
pub enum SmilesError {
    Parse(String),
    Layout(String),
}

impl fmt::Display for SmilesError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(message) => write!(formatter, "parse error: {message}"),
            Self::Layout(message) => write!(formatter, "layout error: {message}"),
        }
    }
}

impl From<SmilesError> for String {
    fn from(error: SmilesError) -> Self {
        error.to_string()
    }
}
