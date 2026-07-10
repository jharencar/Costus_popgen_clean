#!/usr/bin/env bash
# Run while dissert is mounted (default /Volumes/dissert; override with COSTUS_DISS).
# Usage: bash inventory_dissert.sh > dissert_inventory_$(date +%Y%m%d).txt

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dissert_common.sh
source "${SCRIPT_DIR}/dissert_common.sh"

echo "=== Dissert inventory $(date -Iseconds) ==="
echo "Mount: ${DISS}"
df -h "${DISS}" 2>/dev/null || { echo "ERROR: ${DISS} not mounted"; exit 1; }
echo

echo "=== Top-level sizes ==="
du -sh "${DISS}"/*/ 2>/dev/null | sort -h
echo

if [[ -d "${ARCH}" ]]; then
  echo "=== Costus_archive sizes ==="
  du -sh "${ARCH}"/*/ 2>/dev/null | sort -h
  echo
fi

echo "=== Key publication artifacts (presence + size) ==="
check_glob() {
  local label="$1"
  shift
  local paths=("$@")
  local n=0 total
  while IFS= read -r f; do
    [[ -n "$f" ]] && ((n++)) || true
  done < <(find "${paths[@]}" -type f 2>/dev/null)
  if [[ "$n" -gt 0 ]]; then
    total=$(find "${paths[@]}" -type f -print0 2>/dev/null | xargs -0 du -ch 2>/dev/null | tail -1 | awk '{print $1}')
    echo "OK  ${label}: ${n} files, ${total:-?}"
  else
    echo "MISSING  ${label}"
  fi
}

# VCFs (dissert or archive)
for pat in \
  "BAVLasi_LASIref_filtered_scaf_*.vcf.gz" \
  "renamed_AVonly_all_sites_mappability_maxDepth_filtered_chr*.vcf.gz"
do
  echo -n "${pat}: "
  find "${DISS}" -name "${pat}" 2>/dev/null | wc -l | tr -d ' '
done

check_glob "Dsuite Fstats" "${DISS}" -path "*/07_dsuite/*" -name "Fstats*" 2>/dev/null
check_glob "Dsuite Fstats (legacy humm_backup)" "${DISS}/humm_backup" -name "Fstats*" 2>/dev/null
check_glob "pixy output dirs" "${DISS}" -path "*/06_pixy/output*" 2>/dev/null
check_glob "QTL posteriors" "${DISS}" -path "*/03_aHMM_posteriors/*" -name "*.posterior" 2>/dev/null
check_glob "lepmap3" "${DISS}" -path "*/04_lepmap3/*" 2>/dev/null
check_glob "genoprob.rds" "${DISS}" -name "*genoprob*.rds" 2>/dev/null
check_glob "scantwo perm rds" "${DISS}" -path "*scantwo*" -name "*.rds" 2>/dev/null
check_glob "ANGSD abba_baba" "${DISS}" -path "*/abba_baba/*" -name "*.abbababa2" 2>/dev/null
check_glob "wild aHMM LASIref2" "${DISS}" -path "*AHV_LASIref2*" 2>/dev/null
check_glob "mapping pop FASTQs" "${DISS}/Costus_mapping_pop_raw_for_NCBI" -name "*.fastq.gz" 2>/dev/null

echo
echo "=== Mapping pop FASTQ count ==="
find "${DISS}/Costus_mapping_pop_raw_for_NCBI" -name "*.fastq.gz" 2>/dev/null | wc -l

echo
echo "=== Uncompressed legacy VCFs (space recovery candidates) ==="
find "${DISS}" -name "renamed_AHV_mappability_maxDepth_filtered_scaf*.vcf" -type f 2>/dev/null | while read -r f; do
  du -sh "$f"
done

echo
echo "=== humm_backup top-level (if still present) ==="
if [[ -d "${DISS}/humm_backup" ]]; then
  du -sh "${DISS}/humm_backup"/* 2>/dev/null | sort -h | tail -25
fi
