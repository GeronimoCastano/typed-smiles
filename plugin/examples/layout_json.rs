//! Prints one JSON object per input line for independent verification tools.
//!
//! Each stdin line holds a SMILES string, optionally preceded by an identifier
//! and a tab. The output line contains either the layout returned to Typst or
//! the diagnostic the package would show.

use std::io::{self, BufRead, Write};

fn main() {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut output = stdout.lock();
    for line in stdin.lock().lines() {
        let line = line.expect("stdin is readable");
        let (identifier, smiles) = line.split_once('\t').unwrap_or(("", line.as_str()));
        let record = match typst_smiles_plugin::layout_native(smiles) {
            Ok(layout) => {
                serde_json::json!({ "id": identifier, "smiles": smiles, "layout": layout })
            }
            Err(error) => serde_json::json!({ "id": identifier, "smiles": smiles, "error": error }),
        };
        writeln!(output, "{record}").expect("stdout is writable");
    }
}
