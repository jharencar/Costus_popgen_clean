#!/usr/bin/env python3
"""Remove redundant rows from combined_sample_manifest.tsv (one row per sample_id)."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO / "sample_manifest" / "combined_sample_manifest.tsv"


def deduplicate(path: Path, backup: bool = True) -> tuple[int, int]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames or "sample_id" not in fieldnames:
            raise SystemExit(f"{path} must have a sample_id column")
        rows = list(reader)

    if backup:
        shutil.copy2(path, path.with_name(path.stem + "_before_dedup.tsv"))

    seen: set[str] = set()
    deduped = []
    for row in rows:
        sid = row["sample_id"]
        if sid in seen:
            continue
        seen.add(sid)
        deduped.append(row)

    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(deduped)

    return len(rows), len(deduped)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--no-backup", action="store_true")
    args = parser.parse_args()
    before, after = deduplicate(args.manifest, backup=not args.no_backup)
    print(f"{args.manifest}: {before} -> {after} rows (removed {before - after})")


if __name__ == "__main__":
    main()
