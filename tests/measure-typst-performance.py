#!/usr/bin/env python3
"""Measures the shipped WASM plugin and Typst rendering per molecule.

Each measurement is the median wall time of fresh `typst compile` processes:

  * cold plugin call: one layout call, minus a document that only imports the
    package (includes loading and instantiating the WASM module);
  * warm plugin call: the extra cost of each further layout call in the same
    document;
  * rendering: the extra cost of drawing a molecule with smiles() over only
    computing its layout.

Typst caches identical plugin calls, so repeated calls use distinct variants
that add a one-atom isotope-labelled `[kLi+]` fragment. Native stage timings come from
`cargo test --release --lib measure_pipeline_stages -- --ignored --nocapture`.
"""

import statistics
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPEATS = 5
CALLS = 20
MOLECULES = {
    'ethanol': 'CCO',
    'caffeine': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
    'cholesterol': 'CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C',
    'strychnine': 'O=C1C[C@@H]2OCC=C3CN4CC[C@@]56[C@H]4C[C@H]3[C@H]2[C@H]6N1c1ccccc15',
    'fullerene': ('c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9'
                  'c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41'),
    'glycine-30-mer': 'N' + 'CC(=O)N' * 29 + 'CC(=O)O',
    'chain-300': 'C' * 300,
}


def compile_seconds(body, directory):
    source = directory / 'measure.typ'
    source.write_text('#import "/src/lib.typ": *\n#import "/src/chemistry.typ": _smiles-plugin\n' + body)
    durations = []
    for _ in range(REPEATS):
        start = time.perf_counter()
        subprocess.run(['typst', 'compile', '--root', str(ROOT), str(source), str(directory / 'out.pdf')],
                       check=True, capture_output=True)
        durations.append(time.perf_counter() - start)
    return statistics.median(durations)


def variant(smiles, index):
    return f'[{index + 1}Li+].{smiles}'


def typst_string(text):
    return '"' + text.replace('\\', '\\\\').replace('"', '\\"') + '"'


def main():
    with tempfile.TemporaryDirectory(dir=ROOT / 'tests') as temporary:
        directory = Path(temporary)
        baseline = compile_seconds('', directory)
        print(f'Typst {subprocess.run(["typst", "--version"], capture_output=True, text=True).stdout.strip()}; '
              f'empty document {baseline * 1000:.0f} ms (median of {REPEATS})\n')
        print('| molecule | cold plugin call ms | warm plugin call ms | rendering ms |')
        print('|---|---:|---:|---:|')
        for name, smiles in MOLECULES.items():
            variants = [typst_string(variant(smiles, index)) for index in range(CALLS)]
            one_call = compile_seconds(f'#json(_smiles-plugin.layout(bytes({variants[0]}))).atoms.len()', directory)
            calls = compile_seconds(''.join(f'#json(_smiles-plugin.layout(bytes({v}))).atoms.len()\n' for v in variants), directory)
            drawings = compile_seconds(''.join(f'#smiles({v})\n' for v in variants), directory)
            cold = one_call - baseline
            warm = (calls - one_call) / (CALLS - 1)
            rendering = (drawings - calls) / CALLS
            print(f'| {name} | {cold * 1000:.1f} | {warm * 1000:.2f} | {rendering * 1000:.1f} |')


if __name__ == '__main__':
    main()
