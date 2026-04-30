#!/usr/bin/env python3
"""Small numeric comparator for iad/ad CLI regression tests."""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values = [float(value) for value in NUMBER_RE.findall(stripped)]
        if values:
            rows.append(values)
    return rows


def extract_last(path: Path) -> int:
    rows = numeric_rows(path)
    if not rows:
        return 1
    print(" ".join(f"{value:.10g}" for value in rows[-1]))
    return 0


def comparable_actual(expected: list[float], actual: list[float]) -> list[float]:
    if len(actual) == len(expected):
        return actual
    if len(actual) == len(expected) + 1 and abs(actual[0] - round(actual[0])) < 1e-12:
        return actual[1:]
    if len(actual) > len(expected):
        return actual[: len(expected)]
    raise ValueError(
        f"actual row has {len(actual)} values, expected at least {len(expected)}"
    )


def compare(expected_path: Path, actual_path: Path, tolerance: float) -> int:
    expected_rows = numeric_rows(expected_path)
    actual_rows = numeric_rows(actual_path)

    if len(expected_rows) != len(actual_rows):
        print(
            f"row count mismatch: expected {len(expected_rows)}, got {len(actual_rows)}",
            file=sys.stderr,
        )
        return 1

    for row_number, (expected, actual) in enumerate(zip(expected_rows, actual_rows), 1):
        try:
            actual_values = comparable_actual(expected, actual)
        except ValueError as error:
            print(f"row {row_number}: {error}", file=sys.stderr)
            return 1

        for column, (left, right) in enumerate(zip(expected, actual_values), 1):
            if not math.isfinite(right) or abs(left - right) > tolerance:
                print(
                    "row %d column %d mismatch: expected %.10g, got %.10g"
                    % (row_number, column, left, right),
                    file=sys.stderr,
                )
                return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract-last", metavar="FILE")
    parser.add_argument("--tolerance", type=float, default=1e-4)
    parser.add_argument("expected", nargs="?")
    parser.add_argument("actual", nargs="?")
    args = parser.parse_args()

    if args.extract_last:
        return extract_last(Path(args.extract_last))

    if not args.expected or not args.actual:
        parser.error("expected and actual files are required unless --extract-last is used")

    return compare(Path(args.expected), Path(args.actual), args.tolerance)


if __name__ == "__main__":
    raise SystemExit(main())
