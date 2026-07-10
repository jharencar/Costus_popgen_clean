#!/usr/bin/env python3
"""Build combined wild + QTL mapping sample manifest for publication."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_WILD = REPO / "sample_manifest" / "wild_sample_metadata.tsv"
DEFAULT_QTL = REPO / "sample_manifest" / "qtl_mapping_sample_metadata.tsv"
DEFAULT_OUT = REPO / "sample_manifest" / "combined_sample_manifest.tsv"

OUTPUT_COLUMNS = [
    "cohort", "sample_name", "sample_id", "organism", "replicate",
    "development_stage", "ecotype", "tissue", "collection_date",
    "geographic_location", "lat_lon", "field_collection_id",
    "species_field_sheet", "site", "latitude", "longitude", "elevation_m",
    "greenhouse_cross_id", "mother_id", "father_id", "hyb_generation",
    "library_ID", "library_strategy", "library_source", "library_selection",
    "library_layout", "platform", "instrument_model", "fastq_r1", "fastq_r2",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def normalize_row(row: dict[str, str], cohort: str) -> dict[str, str]:
    out = {col: "" for col in OUTPUT_COLUMNS}
    out["cohort"] = cohort
    for key, val in row.items():
        if key in out:
            out[key] = val or ""
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wild", type=Path, default=DEFAULT_WILD)
    parser.add_argument("--qtl", type=Path, default=DEFAULT_QTL)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    wild_rows = [normalize_row(r, "wild_popgen") for r in read_tsv(args.wild)]
    qtl_rows = [normalize_row(r, "qtl_mapping") for r in read_tsv(args.qtl)]
    combined = wild_rows + qtl_rows
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(combined)
    print(f"Wrote {len(combined)} rows ({len(wild_rows)} wild + {len(qtl_rows)} QTL) -> {args.out}")


if __name__ == "__main__":
    main()
