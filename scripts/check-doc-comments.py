#!/usr/bin/env python3
"""Check that `///` doc comments above public Typst functions compile as markup.

tinymist renders each `///` block as Typst markup in hovers and completions,
where `#name` starts a code expression. This script reproduces that by
compiling every block on its own as plain markup, so a broken reference fails
here instead of only in the editor.

Run from anywhere with:

    python3 -I scripts/check-doc-comments.py

The script exits with status 1 when any block fails to compile.
"""

import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCE_ROOT = REPO_ROOT / "src"

FUNCTION_DEFINITION = re.compile(r"^#let\s+([A-Za-z_][\w-]*)\s*\(")
DOC_LINE = re.compile(r"^\s*///")
DOC_PREFIX = re.compile(r"^\s*/// ?")
TYPST_ERROR = re.compile(r"^.+?:(\d+):\d+: (error: .*)$")


@dataclass(frozen=True)
class DocBlock:
    source: Path
    first_line: int
    function_name: str
    text: str


def find_doc_blocks(source: Path) -> list[DocBlock]:
    """Return the doc block directly above each top-level function definition."""
    lines = source.read_text(encoding="utf-8").splitlines()
    blocks = []
    for definition_index, line in enumerate(lines):
        match = FUNCTION_DEFINITION.match(line)
        if match is None:
            continue
        block_start = definition_index
        while block_start > 0 and DOC_LINE.match(lines[block_start - 1]):
            block_start -= 1
        if block_start == definition_index:
            continue
        doc_lines = [
            DOC_PREFIX.sub("", doc_line)
            for doc_line in lines[block_start:definition_index]
        ]
        blocks.append(
            DocBlock(
                source=source,
                first_line=block_start + 1,
                function_name=match.group(1),
                text="\n".join(doc_lines) + "\n",
            )
        )
    return blocks


def compile_block(typst: str, block: DocBlock, work_dir: Path, number: int) -> list[tuple[int, str]]:
    """Compile one block as plain markup and return (source line, message) errors."""
    block_path = work_dir / f"block-{number}.typ"
    block_path.write_text(block.text, encoding="utf-8")
    result = subprocess.run(
        [
            typst,
            "compile",
            "--diagnostic-format",
            "short",
            str(block_path),
            str(work_dir / f"block-{number}.pdf"),
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        return []

    errors = []
    for output_line in result.stderr.splitlines():
        match = TYPST_ERROR.match(output_line)
        if match is None:
            continue
        block_line = int(match.group(1))
        source_line = block.first_line + block_line - 1
        errors.append((source_line, match.group(2)))
    if errors:
        return errors
    return [(block.first_line, result.stderr.strip() or "error: typst failed")]


def main() -> int:
    typst = shutil.which("typst")
    if typst is None:
        print("error: typst was not found on PATH", file=sys.stderr)
        return 2

    blocks = [
        block
        for source in sorted(SOURCE_ROOT.rglob("*.typ"))
        for block in find_doc_blocks(source)
    ]

    failed_count = 0
    with tempfile.TemporaryDirectory(prefix="typed-smiles-doc-check-") as work_name:
        work_dir = Path(work_name)
        for number, block in enumerate(blocks):
            errors = compile_block(typst, block, work_dir, number)
            if not errors:
                continue
            failed_count += 1
            relative_source = block.source.relative_to(REPO_ROOT)
            for source_line, message in errors:
                print(
                    f"{relative_source}:{source_line}: doc block of "
                    f"`{block.function_name}`: {message}"
                )

    passed_count = len(blocks) - failed_count
    print(f"{passed_count} of {len(blocks)} doc blocks compile as markup")
    return 1 if failed_count else 0


if __name__ == "__main__":
    sys.exit(main())
