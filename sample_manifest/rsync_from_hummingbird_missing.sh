#!/usr/bin/env bash
# Pull only what is still missing on $COSTUS_ARCHIVE (default /Volumes/dissert/Costus_archive) after reorganize.
# Default: dry-run (-n). Pass --execute to transfer.
#
# VCF policy: never pull VCFs from Hummingbird in this script — canonical filtered
# sets are already in Costus_archive/02_VCFs/. Excludes OG-labeled and all .vcf* files.
#
# Usage:
#   bash rsync_from_hummingbird_missing.sh
#   bash rsync_from_hummingbird_missing.sh --execute

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dissert_common.sh
source "${SCRIPT_DIR}/dissert_common.sh"

RSYNC_EXTRA="-n"
[[ "${1:-}" == "--execute" ]] && RSYNC_EXTRA=""

DEST="$ARCH"
HUMM="${HUMM_USER:-jharenca}@hb.ucsc.edu"
REMOTE_BASE="/hb/groups/kay_lab/Julia_popgen"
LOGDIR="${DEST}/manifest/rsync_logs"
mkdir -p "$LOGDIR"

# Keep SSH alive on flaky VPN (code 255 / unexpected end of file)
export RSYNC_RSH="${RSYNC_RSH:-ssh -o ServerAliveInterval=30 -o ServerAliveCountMax=10 -o TCPKeepAlive=yes}"

# exFAT dissert cannot store Unix owner/perms; omit to avoid exit code 23
RSYNC_OPTS=(
  -avP
  --partial
  --no-perms
  --no-owner
  --no-group
  --omit-dir-times
  --modify-window=2
)

# Applied to every rsync — no VCFs, no OG-labeled paths, no trash/old trees
RSYNC_EXCLUDES=(
  --exclude='trash/'
  --exclude='**/trash/'
  --exclude='old_*'
  --exclude='*_old_*'
  --exclude='*OG*'
  --exclude='*og_*'
  --exclude='*unfiltered*'
  --exclude='*.vcf'
  --exclude='*.vcf.gz'
  --exclude='*.vcf.gz.tbi'
  --exclude='*.vcf.biallelic*'
)

rsync_pull() {
  local label="$1"
  local remote="$2"
  local local_path="$3"
  shift 3
  local extra_args=("$@")
  local log="${LOGDIR}/${label}_$(date +%Y%m%d_%H%M%S).log"
  local rc=0
  mkdir -p "$local_path"
  echo "=== ${label} ==="
  echo "  remote: ${HUMM}:${remote}"
  echo "  local:  ${local_path}"
  # shellcheck disable=SC2086
  rsync "${RSYNC_OPTS[@]}" ${RSYNC_EXTRA} \
    "${RSYNC_EXCLUDES[@]}" \
    "${extra_args[@]}" \
    "${HUMM}:${remote}" "${local_path}/" 2>&1 | tee "$log" || rc=$?
  if [[ "$rc" -eq 23 ]]; then
    if grep -qiE 'unexpected end of file|code 255|io_read_|vanished|No space left' "$log"; then
      echo "ERROR: connection dropped or I/O failure — re-run when VPN is stable. See $log"
      exit 23
    fi
    if grep -qiE 'failed to set permissions|failed to set times' "$log"; then
      echo "WARNING: exFAT cannot set Unix permissions — data likely transferred; re-run with updated script."
    else
      echo "WARNING: rsync exit 23 — see $log"
    fi
    rc=0
  elif [[ "$rc" -ne 0 ]]; then
    echo "ERROR: rsync failed with exit $rc — see $log"
    exit "$rc"
  fi
  echo
}

[[ -d "$DEST" ]] || { echo "ERROR: $DEST not found. Run reorganize_dissert_archive.sh --execute first."; exit 1; }

avail="$(df_free_gb "$DEST")"
echo "=== Hummingbird gap-fill rsync (execute=$([[ -z "$RSYNC_EXTRA" ]] && echo true || echo false)) ==="
echo "Free space on dissert: ${avail} GB"
echo "VCF excludes: all .vcf*, *OG*, *unfiltered*, trash/, old_*"
echo

# 1. Dsuite — results + trio files (excludes VCFs via RSYNC_EXCLUDES)
rsync_pull "07_dsuite_full" "${REMOTE_BASE}/Dsuite" "${DEST}/07_dsuite/from_hummingbird"

# 2. Pixy — output folders only
rsync_pull "06_pixy_full" "${REMOTE_BASE}/pixy" "${DEST}/06_pixy/from_hummingbird" \
  --exclude='ALLEref/' \
  --exclude='output_0.2missingness/'

# 3. Scantwo perms — RDS only
rsync_pull "05_rqtl_scantwo_perms" "${REMOTE_BASE}/rqtl/scantwo.perms.done" "${DEST}/05_rqtl/scantwo.perms/from_hummingbird" \
  --include='*.rds' \
  --include='*/' \
  --exclude='*'

# 4. Wild aHMM — already on dissert from Mac copy
repo_wild="${DEST}/03_aHMM_posteriors/wild_AHV_LASIref2/AHV_LASIref2_260210"
if [[ -d "$repo_wild" ]]; then
  echo "=== wild aHMM: local copy present ($(du -sh "$repo_wild" | awk '{print $1}')) — skipping cluster pull ==="
  echo
fi

# 5. ANGSD abba_baba
rsync_pull "08_ANGSD_abba_baba" "${REMOTE_BASE}/ANGSD/abba_baba" "${DEST}/08_ANGSD/abba_baba/from_hummingbird" \
  --exclude='slurm_output/' \
  --exclude='older/'

echo "=== Skipped (already on dissert) ==="
echo "  - VCFs: 02_VCFs/BAVLasi_LASIref_filtered_scaf + renamed_AVonly_pixy"
echo "  - Mapping FASTQs: 09_raw_sequencing/mapping_pop_FASTQs"
echo "  - QTL aHMM posteriors + lepmap3: from alle_vill_QTL Mac copy"
echo
if [[ -n "$RSYNC_EXTRA" ]]; then
  echo "Dry run complete. Re-run with --execute to transfer."
fi
