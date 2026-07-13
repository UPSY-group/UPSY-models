#!/usr/bin/env python3
"""Build a static HTML site from benchmark figure artifacts."""

from __future__ import annotations

import argparse
import os
import shutil
from html import escape
from pathlib import Path


UFEMISM_TESTS_DIR = Path(__file__).resolve().parent / "UFEMISM"


def known_test_dir_names() -> list[str]:
    """Return known UFEMISM integrated test directory names, longest first."""
    if not UFEMISM_TESTS_DIR.exists():
        return []
    return sorted(
        (path.name for path in UFEMISM_TESTS_DIR.iterdir() if path.is_dir()),
        key=len,
        reverse=True,
    )


def infer_test_dir_name(text: str, test_dir_names: list[str]) -> str | None:
    """Infer UFEMISM test directory name present in text."""
    for test_dir_name in test_dir_names:
        if test_dir_name in text:
            return test_dir_name
    return None


def parse_benchmark_info(info_path: Path) -> tuple[str | None, list[str], list[str], str]:
    """Parse benchmark metadata text into title, keywords, URLs, and description."""
    title: str | None = None
    keywords: list[str] = []
    urls: list[str] = []
    description_lines: list[str] = []

    if not info_path.exists():
        return title, keywords, urls, ""

    section = ""
    for raw_line in info_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            if section == "description":
                description_lines.append("")
            continue
        if line.startswith("#"):
            continue

        lower = line.lower().rstrip(":")
        if lower.startswith("title"):
            if ":" in line:
                title_value = line.split(":", 1)[1].strip()
                if title_value:
                    title = title_value
            continue
        if lower == "keywords":
            section = "keywords"
            continue
        if lower == "urls":
            section = "urls"
            continue
        if lower == "description":
            section = "description"
            continue

        value = line[1:].strip() if line.startswith("-") else line
        if section == "keywords":
            keywords.append(value)
        elif section == "urls":
            urls.append(value)
        elif section == "description":
            description_lines.append(raw_line.strip())

    description = "\n".join(description_lines).strip()
    return title, keywords, urls, description


def get_benchmark_info(artifact_dir_name: str) -> tuple[str | None, list[str], list[str], str]:
    """Resolve and load benchmark metadata for a given artifact directory name."""
    test_dir_name = artifact_dir_name
    prefix = "benchmark_figures_"
    if artifact_dir_name.startswith(prefix):
        test_dir_name = artifact_dir_name[len(prefix) :]

    info_path = UFEMISM_TESTS_DIR / test_dir_name / "test_metadata.txt"
    return parse_benchmark_info(info_path)


def render_artifact_section(
    parts: list[str],
    section_title: str,
    images: list[Path],
    site_dir: Path,
    benchmark_key: str,
) -> None:
    """Render one benchmark section with images and optional metadata."""
    title, keywords, urls, description = get_benchmark_info(benchmark_key)
    display_title = title if title else section_title
    parts.append('  <section class="artifact">')
    parts.append(f"    <h2>{escape(display_title)}</h2>")
    if not images:
        parts.append("    <p>No PNG files were published for this benchmark.</p>")
    for image in images:
        rel_path = image.relative_to(site_dir).as_posix()
        parts.append(
            f'    <figure><img src="{escape(rel_path)}" alt="{escape(image.stem)}"><figcaption>{escape(image.name)}</figcaption></figure>'
        )

    if keywords or urls or description:
        parts.append('    <div class="benchmark-meta">')

        if keywords:
            parts.append('      <h3 class="meta-title">Keywords</h3>')
            parts.append('      <ul class="keyword-list">')
            for keyword in keywords:
                parts.append(f'        <li class="keyword-chip">{escape(keyword)}</li>')
            parts.append("      </ul>")

        if urls:
            parts.append('      <h3 class="meta-title">Publications</h3>')
            parts.append('      <ul class="url-list">')
            for url in urls:
                safe_url = escape(url)
                parts.append(f'        <li><a href="{safe_url}" target="_blank" rel="noopener noreferrer">{safe_url}</a></li>')
            parts.append("      </ul>")

        if description:
            parts.append('      <h3 class="meta-title">Description</h3>')
            for paragraph in description.split("\n\n"):
                text = paragraph.strip()
                if text:
                    parts.append(f'      <p class="description">{escape(text)}</p>')

        parts.append("    </div>")

    parts.append("  </section>")


def copy_artifact_images(artifacts_dir: Path, site_images_dir: Path) -> None:
    """Copy benchmark PNG files into site/images from folders or a flat PNG directory."""
    site_images_dir.mkdir(parents=True, exist_ok=True)

    artifact_dirs = sorted(path for path in artifacts_dir.iterdir() if path.is_dir()) if artifacts_dir.exists() else []
    flat_png_paths = sorted(path for path in artifacts_dir.iterdir() if path.suffix.lower() == ".png") if artifacts_dir.exists() else []
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

    if flat_png_paths:
        destination_dir = site_images_dir / artifacts_dir.name
        destination_dir.mkdir(parents=True, exist_ok=True)
        for png_path in flat_png_paths:
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
        "    .artifact figure { margin: 0 0 1rem; }",
        "    .artifact figcaption { color: #4d5659; font-size: 0.95rem; margin-top: 0.4rem; }",
        "    .benchmark-meta { margin-top: 0.75rem; padding-top: 0.5rem; border-top: 1px dashed #d7d0bf; }",
        "    .meta-title { margin: 0.9rem 0 0.4rem; font-family: \"Avenir Next\", \"Segoe UI\", sans-serif; font-size: 0.95rem; color: #2e3a3f; letter-spacing: 0.01em; text-transform: uppercase; }",
        "    .keyword-list { display: flex; flex-wrap: wrap; gap: 0.5rem; margin: 0; padding: 0; list-style: none; }",
        "    .keyword-chip { display: inline-block; padding: 0.25rem 0.65rem; border-radius: 999px; background: #d8efe3; color: #174936; border: 1px solid #9ecdb6; font-size: 0.9rem; font-weight: 600; }",
        "    .url-list { margin: 0; padding-left: 1.2rem; }",
        "    .url-list li { margin: 0.2rem 0; }",
        "    .url-list a { color: #0f4b78; text-decoration-thickness: 1px; }",
        "    .description { margin: 0.45rem 0 0; line-height: 1.5; }",
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

    test_dir_names = known_test_dir_names()

    for artifact_dir in artifact_dirs:
        images = sorted(path for path in artifact_dir.iterdir() if path.suffix.lower() == ".png")

        # Local runs may put all test figures in a single flat folder (e.g. images/figures).
        if artifact_dir.name == "figures":
            grouped_images: dict[str, list[Path]] = {}
            ungrouped: list[Path] = []

            for image in images:
                inferred_test = infer_test_dir_name(image.stem, test_dir_names)
                if inferred_test is None:
                    ungrouped.append(image)
                    continue
                grouped_images.setdefault(inferred_test, []).append(image)

            for test_name in sorted(grouped_images):
                render_artifact_section(parts, test_name, grouped_images[test_name], site_dir, test_name)

            if ungrouped:
                render_artifact_section(parts, artifact_dir.name, ungrouped, site_dir, artifact_dir.name)
            continue

        inferred_test = infer_test_dir_name(artifact_dir.name, test_dir_names)
        section_title = inferred_test if inferred_test else artifact_dir.name
        benchmark_key = inferred_test if inferred_test else artifact_dir.name
        render_artifact_section(parts, section_title, images, site_dir, benchmark_key)

    parts.extend(["</body>", "</html>"])
    (site_dir / "index.html").write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build UFEMISM weekly benchmark static site")
    parser.add_argument(
        "--artifacts-dir",
        required=True,
        help="Directory containing benchmark artifact folders or local PNG files",
    )
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