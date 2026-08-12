#!/usr/bin/env python3

"""Create a temporary source tree with apostrophes stripped from Fortran comments.

The build uses gfortran with -cpp, which can mis-handle apostrophes that appear in
Fortran comment lines. This helper copies the build-relevant source tree to a
temporary location and removes apostrophes from comment lines in copied .f90 files.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


FORTRAN_SUFFIXES = {".f90", ".f95", ".f03", ".f08"}


def strip_comment_apostrophes(text: str) -> str:
    stripped_lines = []
    for line in text.splitlines(keepends=True):
        if line.lstrip().startswith("!"):
            stripped_lines.append(line.replace("'", ""))
        else:
            stripped_lines.append(line)
    return "".join(stripped_lines)


def copy_tree(source_root: Path, target_root: Path) -> None:
    target_root.mkdir(parents=True, exist_ok=True)

    for entry_name in ("CMakeLists.txt", "cmake", "src", "include"):
        source_path = source_root / entry_name
        if not source_path.exists():
            continue

        destination_path = target_root / entry_name
        if source_path.is_dir():
            shutil.copytree(source_path, destination_path, copy_function=shutil.copy2)
        else:
            destination_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_path, destination_path)

    for fortran_path in target_root.rglob("*"):
        if not fortran_path.is_file() or fortran_path.suffix.lower() not in FORTRAN_SUFFIXES:
            continue

        original_text = fortran_path.read_text(encoding="utf-8")
        stripped_text = strip_comment_apostrophes(original_text)
        if stripped_text != original_text:
            fortran_path.write_text(stripped_text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path, help="Path to the repository root")
    parser.add_argument(
        "--target-root",
        type=Path,
        help="Destination directory for the stripped source tree. Defaults to a temporary directory.",
    )
    args = parser.parse_args()

    source_root = args.source_root.resolve()
    if not source_root.is_dir():
        raise SystemExit(f"source root does not exist or is not a directory: {source_root}")

    if args.target_root:
        target_root = args.target_root.resolve()
    else:
        build_root = source_root / "build"
        build_root.mkdir(parents=True, exist_ok=True)
        target_root = build_root / "stripped-source"

    if target_root.exists():
        shutil.rmtree(target_root)

    copy_tree(source_root, target_root)
    print(target_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())