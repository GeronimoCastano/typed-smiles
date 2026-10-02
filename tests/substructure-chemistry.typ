#import "../src/lib.typ": smiles, substructure-matches, functional-groups
#import "../src/substructure.typ": _substructure-highlights

#let corpus = json("fixtures/substructure-pubchem.json")

// RDKit supplies the exact match sets; reviewed counts check the intended
// chemical classes separately. Run these checks through the shipped WASM API.
#for case in corpus.cases {
  for variant in case.variants {
    for (group, pattern) in functional-groups {
      let expected = variant.matches.at(group, default: ())
      let actual = substructure-matches(variant.smiles, pattern)
      let description = case.name + " / " + variant.notation + " / " + group
      assert.eq(actual, expected, message: description)
      if variant.notation == "aromatic" {
        assert.eq(actual.len(), case.at("reviewed-counts").at(group, default: 0), message: description)
      }
      for include-atoms in (true, false) {
        let annotations = _substructure-highlights(variant.smiles,
          highlight-groups: (group: group, include-atoms: include-atoms),
          highlight-unmatched: "ignore")
        assert.eq(annotations.len(), expected.len(), message: description)
        for (annotation, match) in annotations.zip(expected) {
          assert.eq(annotation.include-atoms, include-atoms)
          let bonds = annotation.ref.filter(ref => ref.__ref__ == "bond")
            .map(ref => (ref.i, ref.j))
          let atoms = annotation.ref.filter(ref => ref.__ref__ == "atom")
            .map(ref => ref.index)
          if include-atoms { atoms += bonds.flatten() }
          let expected-atoms = if include-atoms {
            match.atoms
          } else {
            match.atoms.filter(index => index not in match.bonds.flatten())
          }
          assert.eq(atoms.dedup().sorted(), expected-atoms, message: description)
          assert.eq(bonds, match.bonds, message: description)
        }
      }
    }
  }
}

// Every named group has a positive control in the independent corpus.
#for group in functional-groups.keys() {
  assert(corpus.cases.any(case => case.at("reviewed-counts").at(group, default: 0) > 0))
}

#set text(font: "New Computer Modern", size: 11pt)
#set page(margin: 2cm)

= PubChem functional-group verification

The pinned corpus checks 63 molecules against all 16 groups, with exact atom
and bond sets from RDKit. These examples use lowercase aromatic notation.

#let sample(name, groups) = {
  let case = corpus.cases.find(case => case.name == name)
  [*#name* \ #text(size: 8pt, groups.join(", ")) \
    #smiles(case.variants.first().smiles, highlight-groups: groups,
      highlight-unmatched: "ignore", scale: 0.68,
      highlight-colors: if "thiol" in groups { (rgb("#BBE1FA"),) } else { auto })]
}

#grid(columns: (1fr, 1fr, 1fr), gutter: 14pt,
  sample("aspirin", ("carboxylic-acid", "ester")),
  sample("acetaminophen", ("phenol", "amide")),
  sample("caffeine", ("carbonyl",)),
  sample("glycerol", ("alcohol",)),
  sample("arginine", ("amine",)),
  sample("nicotine", ("amine",)),
  sample("sodium acetate", ("carboxylate",)),
  sample("benzaldehyde", ("aldehyde",)),
  sample("pyruvic acid", ("ketone",)),
  sample("acrylonitrile", ("alkene", "nitrile")),
  sample("acetylene", ("alkyne",)),
  sample("cysteine", ("thiol",)),
  sample("vanillic acid", ("ether",)),
  sample("nitrobenzene", ("nitro",)),
  sample("urea", ("amide",)),
)

== Negative controls

No shading is expected for these requests; carbonyl groups remain selectable
in the anhydride and carbamate.

#grid(columns: (1fr, 1fr, 1fr), gutter: 14pt,
  sample("nitroglycerin", ("nitro",)),
  sample("acetic anhydride", ("ester",)),
  sample("ethyl carbamate", ("ester", "amine")),
  sample("guanidine", ("amine",)),
  sample("cyanamide", ("amine", "nitrile")),
  sample("carbonic acid", ("carboxylic-acid",)),
)
