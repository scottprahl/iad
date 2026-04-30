#!/usr/bin/env python3
"""Remove generated files without relying on broad git clean patterns.

By default this removes only build outputs and known generated intermediates.
Use --include-tracked-generated to also remove committed generated deliverables
that can be rebuilt from the CWEB and documentation sources.
"""

from __future__ import annotations

import argparse
import fnmatch
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

SAFE_PATTERNS = [
    ".pytest_cache",
    ".jupyter",
    ".ipynb_checkpoints",
    "ad",
    "iad",
    "mc_lost",
    "ad.exe",
    "iad.exe",
    "libiad.a",
    "libiad.so",
    "libiad.dylib",
    "libiad.dll",
    "src/*.o",
    "src/ad",
    "src/iad",
    "src/ad.exe",
    "src/iad.exe",
    "src/cone_test",
    "src/layer_test",
    "src/mc_test",
    "src/oblique_test",
    "src/mc_lost",
    "src/adchapter",
    "src/libiad.a",
    "src/libiad.so",
    "src/libiad.dylib",
    "src/libiad.dll",
    "src/libiad.h",
    "src/libiad.tmp",
    "src/libiad.tmp2",
    "src/lib_ad.h",
    "src/lib_iad.h",
    "src/*.aux",
    "src/*.dvi",
    "src/*.idx",
    "src/*.log",
    "src/*.pdf",
    "src/*.ref",
    "src/*.scn",
    "src/*.sref",
    "src/*.tex",
    "src/*.toc",
    "docs/manual.aux",
    "docs/manual.bbl",
    "docs/manual.blg",
    "docs/manual.log",
    "docs/manual.out",
    "docs/manual.toc",
    "docs/jupyter",
    "docs/tutorial/*.aux",
    "docs/tutorial/*.fdb_latexmk",
    "docs/tutorial/*.fls",
    "docs/tutorial/*.log",
    "docs/tutorial/*.out",
    "docs/tutorial/*.synctex.gz",
    "docs/tutorial/*.toc",
    "tests/.jupyter",
    "tests/.ipynb_checkpoints",
    "tests/*.txt",
    "tests/*.svg",
    "tests/*.ipynb",
    "tests/*.grid",
    "tests/*notebook",
]

TRACKED_GENERATED_EXTRAS = [
    "docs/ad_src.pdf",
    "docs/iad_src.pdf",
    "docs/manual.pdf",
]

ARCHIVE_PATTERNS = [
    "archives/iad-*.zip",
    "archives/iad-win-*.zip",
    "archives/iad-latest.zip",
    "archives/iad-win-latest.zip",
    "iad-*",
    "iad-win-*",
]

NEVER_REMOVE = {
    ".",
    ".git",
    "CLAUDE.md",
    "PLAN.md",
    "Makefile",
    "src/Makefile",
}


def git_files(args: list[str]) -> set[Path]:
    output = subprocess.check_output(["git", *args], cwd=ROOT, text=True)
    return {Path(line) for line in output.splitlines() if line}


def tracked_generated_from_cweb() -> set[Path]:
    generated: set[Path] = set()
    for source in sorted((ROOT / "src").glob("*.w")):
        stem = source.with_suffix("")
        for suffix in (".c", ".h"):
            candidate = stem.with_suffix(suffix).relative_to(ROOT)
            if candidate.exists():
                generated.add(candidate)
    # These CWEB roots do not produce committed public source files.
    generated.discard(Path("src/ad.c"))
    generated.discard(Path("src/ad.h"))
    generated.discard(Path("src/iad.c"))
    generated.discard(Path("src/iad.h"))
    generated.discard(Path("src/iad_main_mus.c"))
    generated.discard(Path("src/iad_main_mus.h"))
    return generated


def matches(path: Path, patterns: list[str]) -> bool:
    text = path.as_posix()
    for pattern in patterns:
        if "/" not in pattern:
            if "/" not in text and fnmatch.fnmatchcase(text, pattern):
                return True
            continue

        path_parts = text.split("/")
        pattern_parts = pattern.split("/")
        if len(path_parts) != len(pattern_parts):
            continue
        if all(fnmatch.fnmatchcase(part, pat) for part, pat in zip(path_parts, pattern_parts)):
            return True
    return False


def collect_candidates(include_tracked_generated: bool, include_archives: bool) -> list[Path]:
    tracked = git_files(["ls-files"])
    all_paths = [
        path.relative_to(ROOT)
        for path in ROOT.rglob("*")
        if ".git" not in path.relative_to(ROOT).parts
    ]

    candidates: set[Path] = set()
    for path in all_paths:
        if path.as_posix() in NEVER_REMOVE:
            continue
        if path in tracked and not include_tracked_generated:
            continue
        if matches(path, SAFE_PATTERNS):
            candidates.add(path)
        if include_archives and matches(path, ARCHIVE_PATTERNS):
            candidates.add(path)

    if include_tracked_generated:
        candidates.update(tracked_generated_from_cweb())
        candidates.update(Path(path) for path in TRACKED_GENERATED_EXTRAS)

    existing = [path for path in candidates if (ROOT / path).exists()]
    existing.sort(key=lambda item: (len(item.parts), item.as_posix()))
    return existing


def remove_path(path: Path) -> None:
    full_path = ROOT / path
    if full_path.is_dir():
        shutil.rmtree(full_path)
    else:
        full_path.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="remove files")
    parser.add_argument(
        "--include-tracked-generated",
        action="store_true",
        help="also remove committed generated CWEB/doc deliverables",
    )
    parser.add_argument(
        "--include-archives",
        action="store_true",
        help="also remove generated release archives and staging directories",
    )
    args = parser.parse_args()

    candidates = collect_candidates(args.include_tracked_generated, args.include_archives)
    action = "remove" if args.execute else "would remove"
    for path in candidates:
        print(f"{action} {path.as_posix()}")
    if args.execute:
        for path in sorted(candidates, key=lambda item: len(item.parts), reverse=True):
            remove_path(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
