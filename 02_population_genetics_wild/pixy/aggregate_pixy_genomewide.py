#!/usr/bin/env python3
"""Aggregate pixy windowed output into genome-wide pi, dxy, and Fst.

Usage:
    aggregate_pixy_genomewide.py <pixy_out_dir> [out_tsv]

The genome-wide value is NOT the average of per-window statistics. Sum the raw
counts or estimator components across all windows and chromosomes, then divide
once (ratio of averages; Bhatia et al. 2013):

    pi, dxy    = sum(count_diffs) / sum(count_comparisons)
    hudson_fst = sum(hudson_fst_num) / sum(hudson_fst_den)
    wc_fst     = sum(wc_fst_a) / (sum(wc_fst_a) + sum(wc_fst_b) + sum(wc_fst_c))

Hudson columns are used when present (pixy --fst_type hudson --fst_components).
Weir–Cockerham components are used only when the Hudson columns are absent.
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

    # ---- Fst: Hudson ratio of sums, else Weir & Cockerham ----
    fst_rows = read_rows(os.path.join(outdir, "*_fst.txt"))
    has_hudson = any("hudson_fst_num" in row and "hudson_fst_den" in row for row in fst_rows)
    has_wc = any("wc_fst_a" in row for row in fst_rows)

    if has_hudson:
        hudson = {}
        for row in fst_rows:
            key = (row["pop1"], row["pop2"])
            num, den = hudson.get(key, (0.0, 0.0))
            hudson[key] = (num + to_float(row.get("hudson_fst_num")),
                           den + to_float(row.get("hudson_fst_den")))
        for (p1, p2), (num, den) in sorted(hudson.items()):
            results.append(("hudson_fst", p1, p2, num / den if den else float("nan")))
    elif has_wc:
        fst = {}
        for row in fst_rows:
            key = (row["pop1"], row["pop2"])
            a, b, c = fst.get(key, (0.0, 0.0, 0.0))
            fst[key] = (a + to_float(row["wc_fst_a"]),
                        b + to_float(row["wc_fst_b"]),
                        c + to_float(row["wc_fst_c"]))
        for (p1, p2), (a, b, c) in sorted(fst.items()):
            denom = a + b + c
            results.append(("wc_fst", p1, p2, a / denom if denom else float("nan")))
    elif fst_rows:
        sys.stderr.write(
            "WARNING: *_fst.txt has neither hudson_fst_num/hudson_fst_den nor "
            "wc_fst_a/b/c. Re-run pixy with --fst_type hudson --fst_components.\n"
        )

    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["statistic", "pop1", "pop2", "genome_wide_value"])
        for stat, p1, p2, val in results:
            w.writerow([stat, p1, p2, f"{val:.8g}"])

    print("statistic\tpop1\tpop2\tgenome_wide_value")
    for stat, p1, p2, val in results:
        print(f"{stat}\t{p1}\t{p2}\t{val:.8g}")


if __name__ == "__main__":
    main()
