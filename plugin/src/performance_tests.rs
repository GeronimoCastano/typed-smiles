//! Measures each stage of the SMILES pipeline on representative and stress
//! molecules. Run with:
//!
//! `cargo test --release --lib measure_pipeline_stages -- --ignored --nocapture`
//!
//! The first run of each molecule is reported separately from the median of
//! later runs, which see warm caches and a warm allocator.

use std::time::{Duration, Instant};

use crate::kekulize::kekulize;
use crate::layout::compute_layout;
use crate::parser::read_written_molecule;

type StageSelector = fn(&StageTimes) -> Duration;

struct StageTimes {
    parse: Duration,
    kekulize: Duration,
    layout: Duration,
    serialize: Duration,
}

fn time_stages(smiles: &str) -> StageTimes {
    let start = Instant::now();
    let written = read_written_molecule(smiles).expect("benchmark molecules parse");
    let parse = start.elapsed();

    let mut molecule = written.molecule;
    let start = Instant::now();
    kekulize(&mut molecule, &written.implicit_bonds).expect("benchmark molecules kekulize");
    let kekulize_time = start.elapsed();

    let start = Instant::now();
    let layout_output = compute_layout(&molecule).expect("benchmark molecules lay out");
    let layout = start.elapsed();

    let start = Instant::now();
    let json = serde_json::to_vec(&layout_output).expect("layouts serialize");
    let serialize = start.elapsed();
    assert!(!json.is_empty());

    StageTimes {
        parse,
        kekulize: kekulize_time,
        layout,
        serialize,
    }
}

fn median(mut durations: Vec<Duration>) -> Duration {
    durations.sort_unstable();
    durations[durations.len() / 2]
}

fn microseconds(duration: Duration) -> String {
    format!("{:.1}", duration.as_secs_f64() * 1e6)
}

pub(crate) fn benchmark_molecules() -> Vec<(&'static str, String)> {
    vec![
        ("ethanol", "CCO".to_string()),
        ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C".to_string()),
        (
            "cholesterol",
            "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
                .to_string(),
        ),
        (
            "strychnine",
            "O=C1C[C@@H]2OCC=C3CN4CC[C@@]56[C@H]4C[C@H]3[C@H]2[C@H]6N1c1ccccc15".to_string(),
        ),
        ("coronene", "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61".to_string()),
        (
            "fullerene",
            "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41"
                .to_string(),
        ),
        ("glycine-30-mer", format!("N{}O", "CC(=O)N".repeat(29) + "CC(=O)")),
        ("macrocycle-100", format!("C1{}C1", "C".repeat(98))),
        ("chain-1000", "C".repeat(1000)),
        ("nested-branches-300", "C(".repeat(300) + "C" + &")".repeat(300)),
    ]
}

#[test]
#[ignore = "measurement, not a correctness check"]
fn measure_pipeline_stages() {
    const WARM_RUNS: usize = 25;
    println!("| molecule | atoms | stage | cold µs | warm median µs |\n|---|---:|---|---:|---:|");
    for (name, smiles) in benchmark_molecules() {
        let atom_count = crate::parse_molecule(&smiles).unwrap().n_atoms();
        let cold = time_stages(&smiles);
        let warm: Vec<StageTimes> = (0..WARM_RUNS).map(|_| time_stages(&smiles)).collect();
        let stages: [(&str, Duration, StageSelector); 4] = [
            ("parse", cold.parse, |times| times.parse),
            ("kekulize", cold.kekulize, |times| times.kekulize),
            ("layout", cold.layout, |times| times.layout),
            ("serialize", cold.serialize, |times| times.serialize),
        ];
        for (stage, cold_time, select) in stages {
            let warm_median = median(warm.iter().map(select).collect());
            println!(
                "| {name} | {atom_count} | {stage} | {} | {} |",
                microseconds(cold_time),
                microseconds(warm_median)
            );
        }
    }
}
