//! Geometric quality checks for generated 2D layouts.
//!
//! A molecule that parses correctly can still be drawn badly. These checks
//! measure drawing quality on their own terms, separate from parsing, so a
//! layout regression is reported as a layout failure rather than hidden
//! behind a successful parse.

use crate::geometry::{point_to_segment_distance, segments_cross};
use crate::render::LayoutOutput;

/// Drawing-quality measurements for one layout, in bond-length units.
#[derive(Debug, Clone, Copy)]
pub(crate) struct LayoutQuality {
    /// Largest relative deviation of a drawn bond from the unit bond length.
    pub worst_bond_length_error: f64,
    /// Smallest distance between two atoms that are not bonded to each other.
    pub closest_nonbonded_atoms: f64,
    /// Smallest distance from an atom to a bond it does not belong to.
    pub closest_atom_to_foreign_bond: f64,
    /// Number of pairs of bonds without a shared atom whose segments cross.
    pub bond_crossings: usize,
}

/// Thresholds that separate an acceptable drawing from a layout failure.
#[derive(Debug, Clone, Copy)]
pub(crate) struct QualityLimits {
    pub max_bond_length_error: f64,
    pub min_nonbonded_distance: f64,
    pub min_atom_to_bond_distance: f64,
    pub max_bond_crossings: usize,
}

/// Limits for a drawing that reads as a clean textbook structure.
pub(crate) const CLEAN_DRAWING: QualityLimits = QualityLimits {
    max_bond_length_error: 0.25,
    min_nonbonded_distance: 0.5,
    min_atom_to_bond_distance: 0.3,
    max_bond_crossings: 0,
};

impl LayoutQuality {
    pub(crate) fn violations(&self, limits: QualityLimits) -> Vec<String> {
        let mut violations = Vec::new();
        if self.worst_bond_length_error > limits.max_bond_length_error {
            violations.push(format!(
                "bond length deviates by {:.2} (limit {:.2})",
                self.worst_bond_length_error, limits.max_bond_length_error
            ));
        }
        if self.closest_nonbonded_atoms < limits.min_nonbonded_distance {
            violations.push(format!(
                "nonbonded atoms {:.2} apart (limit {:.2})",
                self.closest_nonbonded_atoms, limits.min_nonbonded_distance
            ));
        }
        if self.closest_atom_to_foreign_bond < limits.min_atom_to_bond_distance {
            violations.push(format!(
                "atom {:.2} from a foreign bond (limit {:.2})",
                self.closest_atom_to_foreign_bond, limits.min_atom_to_bond_distance
            ));
        }
        if self.bond_crossings > limits.max_bond_crossings {
            violations.push(format!(
                "{} crossing bond pairs (limit {})",
                self.bond_crossings, limits.max_bond_crossings
            ));
        }
        violations
    }
}

/// Measures a rendered layout, ignoring the virtual hydrogens that exist only
/// as reference targets.
pub(crate) fn measure_layout_quality(layout: &LayoutOutput) -> LayoutQuality {
    let drawn_bonds: Vec<(usize, usize)> = layout
        .bonds
        .iter()
        .filter(|bond| !bond.virtual_bond)
        .map(|bond| (bond.from, bond.to))
        .collect();
    let drawn_atoms: Vec<usize> = (0..layout.atoms.len())
        .filter(|&atom_index| !layout.atoms[atom_index].virtual_h)
        .collect();
    let position = |atom_index: usize| layout.atoms[atom_index].pos;
    let bonded = |first_atom: usize, second_atom: usize| {
        drawn_bonds.iter().any(|&(from, to)| {
            (from == first_atom && to == second_atom) || (from == second_atom && to == first_atom)
        })
    };

    let worst_bond_length_error = drawn_bonds
        .iter()
        .map(|&(from, to)| (position(from).distance_to(position(to)) - 1.0).abs())
        .fold(0.0, f64::max);

    let mut closest_nonbonded_atoms = f64::INFINITY;
    for (list_index, &first_atom) in drawn_atoms.iter().enumerate() {
        for &second_atom in &drawn_atoms[list_index + 1..] {
            if !bonded(first_atom, second_atom) {
                closest_nonbonded_atoms = closest_nonbonded_atoms
                    .min(position(first_atom).distance_to(position(second_atom)));
            }
        }
    }

    let mut closest_atom_to_foreign_bond = f64::INFINITY;
    for &atom_index in &drawn_atoms {
        for &(from, to) in &drawn_bonds {
            if atom_index == from || atom_index == to {
                continue;
            }
            closest_atom_to_foreign_bond = closest_atom_to_foreign_bond.min(
                point_to_segment_distance(position(atom_index), position(from), position(to)),
            );
        }
    }

    let mut bond_crossings = 0;
    for (list_index, &(first_from, first_to)) in drawn_bonds.iter().enumerate() {
        for &(second_from, second_to) in &drawn_bonds[list_index + 1..] {
            let share_atom = first_from == second_from
                || first_from == second_to
                || first_to == second_from
                || first_to == second_to;
            if !share_atom
                && segments_cross(
                    position(first_from),
                    position(first_to),
                    position(second_from),
                    position(second_to),
                )
            {
                bond_crossings += 1;
            }
        }
    }

    LayoutQuality {
        worst_bond_length_error,
        closest_nonbonded_atoms,
        closest_atom_to_foreign_bond,
        bond_crossings,
    }
}

/// Limits for the perspective views used for cages and for bicyclics whose
/// one-atom bridge carries substituents. Foreshortened bonds and bonds
/// passing behind others are part of those drawings.
pub(crate) const PERSPECTIVE_DRAWING: QualityLimits = QualityLimits {
    max_bond_length_error: 0.32,
    min_nonbonded_distance: 0.5,
    min_atom_to_bond_distance: 0.3,
    max_bond_crossings: 2,
};

/// Twistane has no template, so arcs and relaxation leave stretched bonds
/// and one crossing.
const TWISTANE_ALLOWANCE: QualityLimits = QualityLimits {
    max_bond_length_error: 0.6,
    max_bond_crossings: 1,
    ..CLEAN_DRAWING
};

/// The ethanamine bridge of morphine is drawn across one ring bond, as in
/// most textbook drawings.
const MORPHINE_ALLOWANCE: QualityLimits = QualityLimits {
    max_bond_crossings: 1,
    ..CLEAN_DRAWING
};

/// Triptycene's three benzene blades cannot all lie flat; one overlaps the
/// others in every 2D drawing.
const TRIPTYCENE_ALLOWANCE: QualityLimits = QualityLimits {
    max_bond_crossings: 3,
    ..CLEAN_DRAWING
};

/// One molecule of the complex-ring corpus with the drawing quality its
/// layout must reach.
pub(crate) struct CorpusMolecule {
    pub category: &'static str,
    pub name: &'static str,
    pub smiles: &'static str,
    pub limits: QualityLimits,
}

/// Molecules whose ring systems are hard to draw, grouped by the kind of
/// difficulty they present: bridged bicyclics, cages, macrocycles, and
/// crowded fused systems.
pub(crate) const COMPLEX_RING_CORPUS: &[CorpusMolecule] = &[
    CorpusMolecule {
        category: "bridged",
        name: "norbornane",
        smiles: "C1CC2CCC1C2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "norbornene",
        smiles: "C1C2CCC1C=C2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "camphor",
        smiles: "CC1(C)C2CCC1(C)C(=O)C2",
        limits: PERSPECTIVE_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "bicyclo[2.2.2]octane",
        smiles: "C1CC2CCC1CC2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "DABCO",
        smiles: "C1CN2CCN1CC2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "quinuclidine",
        smiles: "C1CN2CCC1CC2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "bicyclo[3.3.1]nonane",
        smiles: "C1CC2CCCC(C1)C2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "bicyclo[1.1.1]pentane",
        smiles: "C1C2CC1C2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "alpha-pinene",
        smiles: "CC1=CCC2CC1C2(C)C",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "1,4-cineole",
        smiles: "CC12CCC(CC1)C(C)(C)O2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "tropanol",
        smiles: "CN1C2CCC1CC(O)C2",
        limits: PERSPECTIVE_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "cocaine",
        smiles: "CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)c1ccccc1)C(=O)OC",
        limits: PERSPECTIVE_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "quinine",
        smiles: "COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "cage",
        name: "adamantane",
        smiles: "C1C2CC3CC1CC(C2)C3",
        limits: PERSPECTIVE_DRAWING,
    },
    CorpusMolecule {
        category: "cage",
        name: "cubane",
        smiles: "C12C3C4C1C5C2C3C45",
        limits: PERSPECTIVE_DRAWING,
    },
    CorpusMolecule {
        category: "cage",
        name: "twistane",
        smiles: "C1CC2CC3CCC2CC13",
        limits: TWISTANE_ALLOWANCE,
    },
    CorpusMolecule {
        category: "bridged",
        name: "morphine",
        smiles: "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
        limits: MORPHINE_ALLOWANCE,
    },
    CorpusMolecule {
        category: "bridged",
        name: "strychnine",
        smiles: "O=C1C[C@H]2OCC=C3CN4CC[C@]56[C@@H]4C[C@H]3[C@H]2[C@H]6N1c1ccccc15",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "artemisinin",
        smiles: "C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "[2.2]paracyclophane",
        smiles: "c1cc2ccc1CCc1ccc(cc1)CC2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "bridged",
        name: "triptycene",
        smiles: "c1ccc2c(c1)C1c3ccccc3C2c2ccccc12",
        limits: TRIPTYCENE_ALLOWANCE,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "cyclododecane",
        smiles: "C1CCCCCCCCCCC1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "cyclohexadecane",
        smiles: "C1CCCCCCCCCCCCCCC1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "18-crown-6",
        smiles: "C1COCCOCCOCCOCCOCCO1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "exaltolide",
        smiles: "O=C1CCCCCCCCCCCCCCO1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "muscone-like",
        smiles: "CC1CCCCCCCCCCCCC1=O",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "erythronolide",
        smiles: "CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)(C)O)C)C)O)C",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "porphine",
        smiles: "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "benzo-15-crown-5",
        smiles: "c1ccc2c(c1)OCCOCCOCCOCCO2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "cyclic hexaglycine",
        smiles: "O=C1CNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CN1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "epothilone-like",
        smiles: "CC1=CC=CC(=O)OC(C)CC(O)C(C)C(=O)C(C)(C)C(O)CC1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "coronene",
        smiles: "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "triphenylene",
        smiles: "c1ccc2c(c1)c1ccccc1c1ccccc21",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "pyrene",
        smiles: "c1cc2ccc3cccc4ccc(c1)c2c34",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "perylene",
        smiles: "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "hexahelicene",
        smiles: "c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c21",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "corannulene",
        smiles: "c1cc2ccc3ccc4ccc5ccc1c1c2c3c4c51",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "1,8-dimethylnaphthalene",
        smiles: "Cc1cccc2cccc(C)c12",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "hexaphenylbenzene",
        smiles: "c1ccc(cc1)-c1c(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c(-c2ccccc2)c1-c1ccccc1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "colchicine",
        smiles: "CC(=O)N[C@H]1CCc2cc(OC)c(OC)c(OC)c2-c2ccc(OC)c(=O)cc21",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "estradiol",
        smiles: "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "carbamazepine",
        smiles: "NC(=O)N1c2ccccc2C=Cc2ccccc21",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "bicyclo[4.1.0]heptane",
        smiles: "C1CCC2CC2C1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "benzocyclobutene",
        smiles: "C1Cc2ccccc21",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "benzocyclooctatetraene",
        smiles: "c1ccc2c(c1)C=CC=CC=C2",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "dibenzosuberone",
        smiles: "O=C1c2ccccc2CCc2ccccc12",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "ergoline acid",
        smiles: "CN1C[C@@H](C=C2[C@H]1Cc1c[nH]c3cccc2c13)C(=O)O",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "fused",
        name: "reserpine core",
        smiles: "COC(=O)[C@H]1[C@@H](OC)[C@@H](O)C[C@@H]2CN3CCc4c([nH]c5ccccc45)[C@H]3C[C@@H]21",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "(E)-cyclododecene",
        smiles: "C1CCCCC/C=C/CCCC1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "(Z)-cyclododecene",
        smiles: "C1CCCCC/C=C\\CCCC1",
        limits: CLEAN_DRAWING,
    },
    CorpusMolecule {
        category: "macrocycle",
        name: "(E)-cyclononene",
        smiles: "C1CCC/C=C/CCC1",
        limits: CLEAN_DRAWING,
    },
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layout_native;

    #[test]
    fn complex_ring_corpus_meets_its_quality_limits() {
        let failures: Vec<String> = COMPLEX_RING_CORPUS
            .iter()
            .filter_map(|molecule| {
                let layout = match layout_native(molecule.smiles) {
                    Ok(layout) => layout,
                    Err(error) => {
                        return Some(format!("{}: layout failed: {error}", molecule.name))
                    }
                };
                let violations = measure_layout_quality(&layout).violations(molecule.limits);
                (!violations.is_empty()).then(|| {
                    format!(
                        "{} ({}): {}",
                        molecule.name,
                        molecule.smiles,
                        violations.join("; ")
                    )
                })
            })
            .collect();
        assert!(
            failures.is_empty(),
            "layout failures:\n{}",
            failures.join("\n")
        );
    }

    #[test]
    #[ignore = "prints a quality report for manual inspection"]
    fn report_complex_ring_layout_quality() {
        for molecule in COMPLEX_RING_CORPUS {
            let layout = match layout_native(molecule.smiles) {
                Ok(layout) => layout,
                Err(error) => {
                    println!(
                        "{:<11} {:<26} PARSE ERROR {error}",
                        molecule.category, molecule.name
                    );
                    continue;
                }
            };
            let quality = measure_layout_quality(&layout);
            let strict_violations = quality.violations(CLEAN_DRAWING);
            println!(
                "{:<11} {:<26} bond {:>5.2} atoms {:>5.2} atom-bond {:>5.2} crossings {:>3}  {}",
                molecule.category,
                molecule.name,
                quality.worst_bond_length_error,
                quality.closest_nonbonded_atoms,
                quality.closest_atom_to_foreign_bond,
                quality.bond_crossings,
                if strict_violations.is_empty() {
                    "clean"
                } else {
                    "not clean"
                },
            );
        }
    }
}
