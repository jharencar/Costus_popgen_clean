#!/bin/bash
# Build a genmap mappability mask for VCF filtering (manuscript parameters).
#
# Defaults match popgen_manu_methods.txt / NOTES.md:
#   k-mer length 150, merge distance 100 bp, exclude mappability < 0.7
#
# Usage:
#   export REF_FASTA=/path/to/C.lasius_genome.fasta
#   export GENMAP_OUTDIR=/path/to/output_dir   # optional; default: cwd
#   bash run_genmap_mask.sh
#
# Requires: genmap, bedtools, sort, awk
# Optional SLURM: submit this script with appropriate --cpus-per-task / --mem.

set -euo pipefail

REF_FASTA="${REF_FASTA:?Set REF_FASTA to the reference FASTA (C. lasius / GCA_027563935.1)}"
OUTDIR="${GENMAP_OUTDIR:-.}"
THREADS="${THREADS:-10}"

mappability_k="${MAPPABILITY_K:-150}"
mappability_merge="${MAPPABILITY_MERGE:-100}"
mappability_min="${MAPPABILITY_MIN:-0.7}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

INDEX_DIR="${OUTDIR}/assembly.genmap_index"
MAP_PREFIX="${OUTDIR}/assembly.genmap${mappability_k}_out"
BG_SORTED="${MAP_PREFIX}.s.bg"
TEMP_MASK="${OUTDIR}/temporary_mask.bed"
NON_CALLABLE="${OUTDIR}/non_callable_sites.bed"

echo "[genmap] index → ${INDEX_DIR}"
genmap index -F "${REF_FASTA}" -I "${INDEX_DIR}"

echo "[genmap] map -K ${mappability_k} -E 0"
genmap map -K "${mappability_k}" -E 0 -I "${INDEX_DIR}" \
  -O "${MAP_PREFIX}" -T "${THREADS}" \
  -bg -v &> "${OUTDIR}/assembly.genmap.${mappability_k}.log"

echo "[genmap] sort bedGraph"
sort -k1,1 -k2,2n "${MAP_PREFIX}.bedgraph" > "${BG_SORTED}"

# Low-mappability windows (score < threshold) → mask
echo "[genmap] windows with mappability < ${mappability_min} → ${TEMP_MASK}"
awk -v min_val="${mappability_min}" 'BEGIN{OFS="\t"; FS="\t"} { if ($4 < min_val) print $1,$2,$3 }' \
  "${BG_SORTED}" > "${TEMP_MASK}"

echo "[genmap] merge within ${mappability_merge} bp → ${NON_CALLABLE}"
bedtools sort -i "${TEMP_MASK}" | \
  bedtools merge -d "${mappability_merge}" -i - > "${NON_CALLABLE}"

echo "[genmap] done. Exclude mask in bcftools with: -T ^${NON_CALLABLE}"
