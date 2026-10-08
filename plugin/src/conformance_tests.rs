//! Checks the categorized corpus in tests/fixtures/smiles-conformance.json.
//!
//! Parsing expectations and layout quality are separate: a molecule can parse
//! correctly and still draw with overlapping atoms. tests/verify-conformance.py
//! records RDKit's verdicts and the layout classes this test compares against.

use crate::layout_native;
use crate::render::LayoutOutput;

/// Geometric problems in a drawing, named as tests/verify-conformance.py names them.
fn layout_problems(layout_output: &LayoutOutput) -> Vec<&'static str> {
    let atoms: Vec<usize> = (0..layout_output.atoms.len())
        .filter(|&index| !layout_output.atoms[index].virtual_h)
        .collect();
    let bonds: Vec<(usize, usize)> = layout_output
        .bonds
        .iter()
        .filter(|bond| !bond.virtual_bond)
        .map(|bond| (bond.from, bond.to))
        .collect();
    let position = |index: usize| layout_output.atoms[index].pos;
    let mut problems = Vec::new();

    if bonds.iter().any(|&(from, to)| {
        let length = position(from).distance_to(position(to));
        !(0.75..=1.35).contains(&length)
    }) {
        problems.push("distorted-bonds");
    }
    let is_bonded = |first: usize, second: usize| {
        bonds
            .iter()
            .any(|&(from, to)| (from, to) == (first, second) || (from, to) == (second, first))
    };
    let atoms_overlap = atoms.iter().enumerate().any(|(offset, &first)| {
        atoms[offset + 1..].iter().any(|&second| {
            !is_bonded(first, second) && position(first).distance_to(position(second)) < 0.5
        })
    });
    if atoms_overlap {
        problems.push("overlapping-atoms");
    }
    let bonds_cross = bonds.iter().enumerate().any(|(offset, &(a, b))| {
        bonds[offset + 1..].iter().any(|&(c, d)| {
            ![c, d].contains(&a)
                && ![c, d].contains(&b)
                && segments_cross(position(a), position(b), position(c), position(d))
        })
    });
    if bonds_cross {
        problems.push("crossing-bonds");
    }
    problems.sort_unstable();
    problems
}

fn segments_cross(
    a: crate::render::Vec2,
    b: crate::render::Vec2,
    c: crate::render::Vec2,
    d: crate::render::Vec2,
) -> bool {
    let orientation = |p: crate::render::Vec2, q: crate::render::Vec2, r: crate::render::Vec2| {
        (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x)
    };
    orientation(a, b, c) * orientation(a, b, d) < -1e-9
        && orientation(c, d, a) * orientation(c, d, b) < -1e-9
}

#[test]
fn conformance_corpus_matches_reviewed_expectations() {
    let corpus: serde_json::Value =
        serde_json::from_str(include_str!("../../tests/fixtures/smiles-conformance.json"))
            .expect("the conformance corpus is valid JSON");
    let cases = corpus["cases"].as_array().expect("cases array");
    assert!(cases.len() >= 100);

    for case in cases {
        let label = format!("{}/{}", case["category"], case["name"]);
        let smiles = case["smiles"].as_str().expect("smiles string");
        let result = layout_native(smiles);
        match case["expect"].as_str() {
            Some("reject") => {
                let error = result.expect_err(&label);
                let fragment = case["error"]
                    .as_str()
                    .expect("rejections name an error fragment");
                assert!(error.contains(fragment), "{label}: {error}");
            }
            Some("accept") => {
                let layout_output = result.unwrap_or_else(|error| panic!("{label}: {error}"));
                let expects_undepicted = case["undepicted"].as_bool().unwrap_or(false);
                assert_eq!(
                    !layout_output.undepicted_stereo.is_empty(),
                    expects_undepicted,
                    "{label}: {:?}",
                    layout_output.undepicted_stereo
                );
                let problems = layout_problems(&layout_output).join(",");
                let recorded = case["layout"]
                    .as_str()
                    .expect("accepted cases record a layout");
                if recorded == "clean" {
                    assert!(problems.is_empty(), "{label} now draws with {problems}");
                } else {
                    assert_eq!(
                        problems, recorded,
                        "{label}: update the corpus with --write"
                    );
                }
            }
            other => panic!("{label}: unknown expectation {other:?}"),
        }
    }
}
