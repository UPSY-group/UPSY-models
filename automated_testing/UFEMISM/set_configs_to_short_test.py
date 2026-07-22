#!/usr/bin/env python3
"""Shorten UFEMISM test configs by reducing the run time by 10 years."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


KEY_PATTERN = re.compile(
    r"^(?P<prefix>\s*end_time_of_run_config\s*=\s*)(?P<value>.*?)(?P<suffix>\s*(?:!.*)?)$"
)


def update_config_file(config_path: Path) -> bool:
    """Set end_time_of_run_config to start_time_of_run_config + 10.0 in one config file."""
    lines = config_path.read_text(encoding="utf-8").splitlines()
    start_time: float | None = None
    updated = False

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            continue

        if start_time is None and "start_time_of_run_config" in stripped and "=" in stripped:
            _, value_part = stripped.split("=", 1)
            value_text = value_part.split("!", 1)[0].strip().strip("'").strip('"')
            value_text = value_text.replace("D", "E").replace("d", "e")
            start_time = float(value_text)

    if start_time is None:
        raise ValueError(f"start_time_of_run_config not found in {config_path}")

    end_time = start_time + 10.0

    for index, line in enumerate(lines):
        match = KEY_PATTERN.match(line)
        if not match:
            continue

        lines[index] = f"{match.group('prefix')}{end_time:.1f}{match.group('suffix')}"
        updated = True
        break

    if not updated:
        raise ValueError(f"end_time_of_run_config not found in {config_path}")

    config_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return True


def iter_config_files(test_dir: Path) -> list[Path]:
    """Find all config_*.cfg files below one UFEMISM test directory."""
    return sorted(path for path in test_dir.rglob("config_*.cfg") if path.is_file())


def main() -> None:
    parser = argparse.ArgumentParser(description="Shorten UFEMISM benchmark config files")
    parser.add_argument("test_dir", help="Path to one UFEMISM test directory")
    args = parser.parse_args()

    test_dir = Path(args.test_dir)
    if not test_dir.is_dir():
        raise NotADirectoryError(f"Test directory not found: {test_dir}")

    config_files = iter_config_files(test_dir)
    if not config_files:
        raise FileNotFoundError(f"No config_*.cfg files found under {test_dir}")

    for config_file in config_files:
        update_config_file(config_file)


if __name__ == "__main__":
    main()