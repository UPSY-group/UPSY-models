#!/usr/bin/env python3
"""Build a static HTML site from benchmark figure artifacts."""

from __future__ import annotations

import argparse
import os
import shutil
from html import escape
from pathlib import Path


def copy_artifact_images(artifacts_dir: Path, site_images_dir: Path) -> None:
    """Copy benchmark PNG files grouped by artifact folder into site/images."""
    site_images_dir.mkdir(parents=True, exist_ok=True)

    artifact_dirs = sorted(path for path in artifacts_dir.iterdir() if path.is_dir()) if artifacts_dir.exists() else []
    copied_any = False

    for artifact_dir in artifact_dirs:
        png_paths = sorted(path for path in artifact_dir.iterdir() if path.suffix.lower() == ".png")
        if not png_paths:
            continue

        destination_dir = site_images_dir / artifact_dir.name
        destination_dir.mkdir(parents=True, exist_ok=True)

        for png_path in png_paths:
            shutil.copy2(png_path, destination_dir / png_path.name)
            copied_any = True

    if not copied_any:
        (site_images_dir.parent / "no-artifacts.txt").write_text(
            "No benchmark artifacts were downloaded.\n", encoding="utf-8"
        )


def build_index_html(site_dir: Path) -> None:
    """Write the benchmark gallery index.html based on copied images."""
    images_dir = site_dir / "images"
    artifact_dirs = sorted(path for path in images_dir.iterdir() if path.is_dir()) if images_dir.exists() else []

    run_number = os.environ.get("GITHUB_RUN_NUMBER", "unknown")
    run_id = os.environ.get("GITHUB_RUN_ID", "unknown")
    sha = os.environ.get("GITHUB_SHA", "unknown")
    repository = os.environ.get("GITHUB_REPOSITORY", "unknown")

    parts = [
        "<!doctype html>",
        '<html lang="en">',
        "<head>",
        '  <meta charset="utf-8">',
        '  <meta name="viewport" content="width=device-width, initial-scale=1">',
        "  <title>UFEMISM Weekly Benchmarks</title>",
        "  <style>",
        "    :root { color-scheme: light; }",
        "    body { font-family: Georgia, serif; margin: 2rem auto; max-width: 1100px; padding: 0 1rem 3rem; background: #f4f1e8; color: #1f2a2e; }",
        '    h1, h2 { font-family: "Avenir Next", "Segoe UI", sans-serif; }',
        "    h1 { margin-bottom: 0.25rem; }",
        "    .meta { color: #46555a; margin-bottom: 2rem; }",
        "    .artifact { margin: 2rem 0; padding: 1.25rem; background: #fffdf8; border: 1px solid #d7d0bf; border-radius: 10px; box-shadow: 0 8px 24px rgba(31, 42, 46, 0.08); }",
        "    .artifact img { display: block; width: 100%; height: auto; margin: 1rem 0; border: 1px solid #d7d0bf; border-radius: 6px; background: white; }",
        "    .empty { padding: 1rem; background: #fff4df; border-left: 4px solid #c97912; }",
        "    code { font-size: 0.95em; }",
        "  </style>",
        "</head>",
        "<body>",
        "  <h1>UFEMISM Weekly Benchmarks</h1>",
        f'  <p class="meta">Repository: <code>{escape(repository)}</code><br>Run number: <code>{escape(run_number)}</code><br>Run ID: <code>{escape(run_id)}</code><br>Commit: <code>{escape(sha[:12])}</code></p>',
    ]

    if not artifact_dirs:
        parts.append('  <div class="empty">No benchmark figures were available for this run.</div>')

    for artifact_dir in artifact_dirs:
        images = sorted(path for path in artifact_dir.iterdir() if path.suffix.lower() == ".png")
        parts.append('  <section class="artifact">')
        parts.append(f"    <h2>{escape(artifact_dir.name)}</h2>")
        if not images:
            parts.append("    <p>No PNG files were published for this benchmark.</p>")
        for image in images:
            rel_path = image.relative_to(site_dir).as_posix()
            parts.append(
                f'    <figure><img src="{escape(rel_path)}" alt="{escape(image.stem)}"><figcaption>{escape(image.name)}</figcaption></figure>'
            )
        parts.append("  </section>")

    parts.extend(["</body>", "</html>"])
    (site_dir / "index.html").write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build UFEMISM weekly benchmark static site")
    parser.add_argument("--artifacts-dir", required=True, help="Directory containing downloaded benchmark artifact folders")
    parser.add_argument("--site-dir", required=True, help="Output site directory")
    args = parser.parse_args()

    artifacts_dir = Path(args.artifacts_dir)
    site_dir = Path(args.site_dir)
    site_images_dir = site_dir / "images"

    site_dir.mkdir(parents=True, exist_ok=True)
    copy_artifact_images(artifacts_dir, site_images_dir)
    build_index_html(site_dir)


if __name__ == "__main__":
    main()