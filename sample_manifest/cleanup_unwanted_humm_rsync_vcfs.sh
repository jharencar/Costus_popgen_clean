#!/usr/bin/env bash
# Remove VCFs accidentally pulled into Costus_archive/from_hummingbird/ trees.
# Default: dry-run. Pass --execute to delete.
#
# Usage:
#   bash cleanup_unwanted_humm_rsync_vcfs.sh
#   bash cleanup_unwanted_humm_rsync_vcfs.sh --execute

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dissert_common.sh
source "${SCRIPT_DIR}/dissert_common.sh"

EXEC=false
[[ "${1:-}" == "--execute" ]] && EXEC=true

echo "=== Remove unwanted VCFs from Hummingbird rsync dirs (execute=$EXEC) ==="

search_roots=(
  "${ARCH}/07_dsuite/from_hummingbird"
  "${ARCH}/06_pixy/from_hummingbird"
  "${ARCH}/08_ANGSD/abba_baba/from_hummingbird"
  "${ARCH}/05_rqtl/scantwo.perms/from_hummingbird"
)

patterns=(
  -iname '*OG*'
  -iname '*unfiltered*'
  -name '*.vcf'
  -name '*.vcf.gz'
  -name '*.vcf.gz.tbi'
)

count=0
for root in "${search_roots[@]}"; do
  [[ -d "$root" ]] || continue
  while IFS= read -r -d '' f; do
    ((count++)) || true
    if $EXEC; then
      echo "DELETE: $f"
      rm -f "$f"
    else
      echo "[dry-run] DELETE: $f ($(du -sh "$f" 2>/dev/null | awk '{print $1}'))"
    fi
  done < <(find "$root" \( "${patterns[0]}" -o "${patterns[1]}" -o "${patterns[2]}" -o "${patterns[3]}" -o "${patterns[4]}" \) -print0 2>/dev/null)
done

echo
echo "Matched ${count} file(s)."
if ! $EXEC; then
  echo "Re-run with --execute to delete."
fi
