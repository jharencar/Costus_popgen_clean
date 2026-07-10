#!/usr/bin/env bash
# Remove known redundant / broken paths on dissert to free space.
# Default: dry-run. Pass --execute to delete.
#
# Usage:
#   bash cleanup_dissert_redundant.sh
#   bash cleanup_dissert_redundant.sh --execute

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dissert_common.sh
source "${SCRIPT_DIR}/dissert_common.sh"

EXEC=false
[[ "${1:-}" == "--execute" ]] && EXEC=true

MAPPING="${DISS}/Costus_mapping_pop_raw_for_NCBI"

rm_path() {
  local p="$1"
  local note="$2"
  [[ -e "$p" ]] || { echo "skip (missing): $p"; return 0; }
  if $EXEC; then
    echo "DELETE: $p  ($note)"
    rm -rf "$p"
  else
    echo "[dry-run] DELETE: $p  ($note)"
    du -sh "$p" 2>/dev/null || true
  fi
}

echo "=== Dissert cleanup (execute=$EXEC) ==="
df -h "$DISS" | tail -1
echo

# 818 broken symlinks -> old /Volumes/Radishes/ path; zero bytes of real data
rm_path "${MAPPING}/fastqs_for_upload" \
  "818 broken symlinks to /Volumes/Radishes/; canonical FASTQs are top-level Illumina names"

# Optional: duplicate lane kept separately during NCBI prep (22_033 alternate lane)
# NOT at top level — safe to delete if you do not need the alternate lane
# rm_path "${MAPPING}/_extras/duplicate_lanes" \
#   "212 MB alternate lane for 22_033; not in main mapping set"

# Optional: ~19 GB uncompressed VCFs (gzipped/filtered sets are in Costus_archive/02_VCFs/)
rm_path "${DISS}/Costus_archive/02_VCFs/legacy_uncompressed_AHV" \
  "~19 GB uncompressed renamed_AHV_mappability_*.vcf — verify gz sets first"

# Legacy clumped FASTQs unrelated to main wild/mapping sets
rm_path "${DISS}/_extras/_extras" \
  "1.9 GB SA49794 clumped .fq.gz — not part of main sample manifest"

echo
if ! $EXEC; then
  echo "Re-run with --execute to apply. Review commented duplicate_lanes block first."
fi
