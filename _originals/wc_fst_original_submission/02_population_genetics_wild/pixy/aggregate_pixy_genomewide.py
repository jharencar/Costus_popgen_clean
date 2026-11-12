#!/usr/bin/env python3
"""Aggregate pixy windowed output into genome-wide pi, dxy, and Weir & Cockerham Fst.

Usage:
    aggregate_pixy_genomewide.py <pixy_out_dir> [out_tsv]

The genome-wide value is NOT the average of per-window statistics. Instead, per the
pixy documentation, sum the raw counts/components across ALL windows (and all
chromosomes) and then recompute the statistic:

    pi, dxy  = sum(count_diffs) / sum(count_comparisons)
    wc_fst   = sum(wc_fst_a) / (sum(wc_fst_a) + sum(wc_fst_b) + sum(wc_fst_c))

The Fst calculation requires that pixy was run with --fst_components (so that the
wc_fst_a / wc_fst_b / wc_fst_c columns are present in the *_fst.txt files).
"""
import sys
import os
import glob
import csv


def to_float(x):
    """Parse a pixy numeric cell; treat 'NA'/blank/non-numeric as 0 (no contribution)."""
    try:
        return float(x)
    except (ValueError, TypeError):
        return 0.0


def read_rows(pattern):
    rows = []
    for path in sorted(glob.glob(pattern)):
        with open(path) as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    return rows


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: aggregate_pixy_genomewide.py <pixy_out_dir> [out_tsv]")
    outdir = sys.argv[1]
    out_tsv = sys.argv[2] if len(sys.argv) > 2 else os.path.join(outdir, "genomewide_summary.tsv")

    results = []  # (statistic, pop1, pop2, value)

    # ---- pi: per population ----
    pi = {}
    for row in read_rows(os.path.join(outdir, "*_pi.txt")):
        d, c = pi.get(row["pop"], (0.0, 0.0))
        pi[row["pop"]] = (d + to_float(row["count_diffs"]),
                          c + to_float(row["count_comparisons"]))
    for pop, (d, c) in sorted(pi.items()):
        results.append(("pi", pop, "", d / c if c else float("nan")))

    # ---- dxy: per population pair ----
    dxy = {}
    for row in read_rows(os.path.join(outdir, "*_dxy.txt")):
        key = (row["pop1"], row["pop2"])
        d, c = dxy.get(key, (0.0, 0.0))
        dxy[key] = (d + to_float(row["count_diffs"]),
                    c + to_float(row["count_comparisons"]))
    for (p1, p2), (d, c) in sorted(dxy.items()):
        results.append(("dxy", p1, p2, d / c if c else float("nan")))

    # ---- Weir & Cockerham Fst: per population pair (needs --fst_components) ----
    fst = {}
    missing_components = False
    for row in read_rows(os.path.join(outdir, "*_fst.txt")):
        if "wc_fst_a" not in row:
            missing_components = True
            continue
        key = (row["pop1"], row["pop2"])
        a, b, c = fst.get(key, (0.0, 0.0, 0.0))
        fst[key] = (a + to_float(row["wc_fst_a"]),
                    b + to_float(row["wc_fst_b"]),
                    c + to_float(row["wc_fst_c"]))
    for (p1, p2), (a, b, c) in sorted(fst.items()):
        denom = a + b + c
        results.append(("wc_fst", p1, p2, a / denom if denom else float("nan")))

    if missing_components and not fst:
        sys.stderr.write(
            "WARNING: no wc_fst_a/b/c columns found in *_fst.txt. "
            "Re-run pixy with --fst_components to get a genome-wide Fst.\n"
        )

    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["statistic", "pop1", "pop2", "genome_wide_value"])
        for stat, p1, p2, val in results:
            w.writerow([stat, p1, p2, f"{val:.8g}"])

    # also echo to stdout for the slurm log
    print("statistic\tpop1\tpop2\tgenome_wide_value")
    for stat, p1, p2, val in results:
        print(f"{stat}\t{p1}\t{p2}\t{val:.8g}")


if __name__ == "__main__":
    main()
