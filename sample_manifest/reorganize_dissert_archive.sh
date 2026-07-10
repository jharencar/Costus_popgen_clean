#!/usr/bin/env bash
# Reorganize scattered Costus data on dissert into Costus_archive/ (default /Volumes/dissert).
# Default is dry-run; pass --execute to move files (same volume = no extra space).
#
# Usage:
#   bash reorganize_dissert_archive.sh           # preview
#   bash reorganize_dissert_archive.sh --execute # apply moves

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dissert_common.sh
source "${SCRIPT_DIR}/dissert_common.sh"
HUMM="${DISS}/humm_backup"
EXEC=false
[[ "${1:-}" == "--execute" ]] && EXEC=true

run() {
  if $EXEC; then
    echo "+ $*"
    eval "$@"
  else
    echo "[dry-run] $*"
  fi
}

mkdir_p() {
  if $EXEC; then mkdir -p "$1"; else echo "[dry-run] mkdir -p $1"; fi
}

mv_if() {
  local src="$1" dest="$2"
  [[ -e "$src" ]] || return 0
  mkdir_p "$(dirname "$dest")"
  if $EXEC; then mv "$src" "$dest"; else echo "[dry-run] mv $src -> $dest"; fi
}

ln_if() {
  local target="$1" link="$2"
  [[ -e "$target" ]] || return 0
  if [[ -e "$link" || -L "$link" ]]; then
    echo "skip link (exists): $link"
    return 0
  fi
  run "ln -s \"$target\" \"$link\""
}

[[ -d "$DISS" ]] || { echo "ERROR: $DISS not mounted"; exit 1; }

echo "=== Reorganize dissert -> ${ARCH} (execute=$EXEC) ==="

# --- directory layout ---
for d in manifest 01_reference 02_VCFs/BAVLasi_LASIref_filtered_scaf \
  02_VCFs/renamed_AVonly_pixy 02_VCFs/other_filtered \
  02_VCFs/legacy_uncompressed_AHV 03_aHMM_posteriors/QTL_grouped_run \
  03_aHMM_posteriors/wild_AHV_LASIref2 04_lepmap3 05_rqtl/scantwo.perms \
  06_pixy 07_dsuite 08_ANGSD/abba_baba 08_ANGSD/other \
  09_raw_sequencing 10_scripts_legacy; do
  mkdir_p "${ARCH}/${d}"
done

# --- VCFs from flat humm_backup root ---
if [[ -d "$HUMM" ]]; then
  for i in $(seq 1 9); do
    for f in \
      "BAVLasi_LASIref_filtered_scaf_${i}.vcf.gz" \
      "BAVLasi_LASIref_filtered_scaf_${i}.vcf.gz.tbi" \
      "BAVLasi_LASIref_filtered_scaf_${i}.ldepth.mean" \
      "renamed_AVonly_all_sites_mappability_maxDepth_filtered_chr${i}.vcf.gz" \
      "renamed_AVonly_all_sites_mappability_maxDepth_filtered_chr${i}.vcf.gz.tbi"
    do
      mv_if "${HUMM}/${f}" "${ARCH}/02_VCFs/$( [[ "$f" == BAVLasi* ]] && echo BAVLasi_LASIref_filtered_scaf || echo renamed_AVonly_pixy )/${f}"
    done
  done

  # Other filtered VCF sets at humm_backup root
  for f in "${HUMM}"/renamed_filtered_*.vcf.gz "${HUMM}"/renamed_filtered_no.80.per.spp.filt_*.vcf.gz; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")
    mv_if "$f" "${ARCH}/02_VCFs/other_filtered/${base}"
    [[ -e "${f}.tbi" ]] && mv_if "${f}.tbi" "${ARCH}/02_VCFs/other_filtered/${base}.tbi"
  done

  # Large uncompressed legacy VCFs (free ~19 GB after gzip counterparts verified)
  for f in "${HUMM}"/renamed_AHV_mappability_maxDepth_filtered_scaf*.vcf; do
    [[ -e "$f" ]] || continue
    mv_if "$f" "${ARCH}/02_VCFs/legacy_uncompressed_AHV/$(basename "$f")"
  done

  # Named subdirs
  mv_if "${HUMM}/ANGSD" "${ARCH}/08_ANGSD/_from_humm_backup_ANGSD"
  mv_if "${HUMM}/rqtl" "${ARCH}/05_rqtl/_from_humm_backup_rqtl"
  mv_if "${HUMM}/Costus_abba.baba_redo_stuff" "${ARCH}/08_ANGSD/abba_baba_redo_stuff"
  mv_if "${HUMM}/VCF_scripts" "${ARCH}/10_scripts_legacy/VCF_scripts"
  mv_if "${HUMM}/ch4_scripts_11_2024" "${ARCH}/10_scripts_legacy/ch4_scripts_11_2024"
  mv_if "${HUMM}/pixy_whole_genome" "${ARCH}/10_scripts_legacy/pixy_whole_genome"

  # Loose shell/R scripts at humm_backup root
  shopt -s nullglob
  for f in "${HUMM}"/*.sh "${HUMM}"/*.R "${HUMM}"/*.out "${HUMM}"/*.err "${HUMM}"/*.Rout "${HUMM}"/*.csv; do
    [[ -e "$f" ]] || continue
    mv_if "$f" "${ARCH}/10_scripts_legacy/$(basename "$f")"
  done
  shopt -u nullglob
fi

# --- symlinks to large raw data (leave in place) ---
ln_if "${DISS}/Costus_mapping_pop_raw_for_NCBI" "${ARCH}/09_raw_sequencing/mapping_pop_FASTQs"
ln_if "${DISS}/All_wild_samps_raw_fastqs" "${ARCH}/09_raw_sequencing/wild_raw_FASTQs"
ln_if "${DISS}/OG_Combined_cleaned_bams" "${ARCH}/09_raw_sequencing/wild_cleaned_bams"
ln_if "${DISS}/_extras" "${ARCH}/09_raw_sequencing/sample_metadata_extras"

# --- copy from Mac (not Hummingbird) ---
ALLE="/Users/admin-jgharenc/Documents/Github/alle_vill_QTL"
REPO="/Users/admin-jgharenc/Documents/Github/Costus_popgen_clean"

copy_local() {
  local src="$1" dest="$2"
  [[ -e "$src" ]] || { echo "skip missing: $src"; return 0; }
  mkdir_p "$(dirname "$dest")"
  if $EXEC; then
    echo "+ rsync -a \"$src\" \"$(dirname "$dest")/\""
    rsync -a "$src" "$(dirname "$dest")/"
  else
    echo "[dry-run] rsync -a $src -> $(dirname "$dest")/"
  fi
}

echo
echo "=== Copy from local Mac (alle_vill_QTL + repo) ==="
copy_local "${ALLE}/lepmap3" "${ARCH}/04_lepmap3/lepmap3"
copy_local "${ALLE}/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO" \
  "${ARCH}/03_aHMM_posteriors/QTL_grouped_run/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO"
copy_local "${ALLE}/scantwo.perms" "${ARCH}/05_rqtl/scantwo.perms/from_alle_vill_QTL"
for rds in "${ALLE}"/aHMM_251111_*_genoprob.rds; do
  [[ -e "$rds" ]] || continue
  copy_local "$rds" "${ARCH}/05_rqtl/$(basename "$rds")"
done
copy_local "${REPO}/aHMM/AHV_LASIref2_260210" "${ARCH}/03_aHMM_posteriors/wild_AHV_LASIref2/AHV_LASIref2_260210"
copy_local "${REPO}/pixy/output_251124_LASIref" "${ARCH}/06_pixy/output_251124_LASIref"
copy_local "${REPO}/Dsuite/Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120" \
  "${ARCH}/07_dsuite/Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"

# manifest copies
if $EXEC; then
  cp -f "$(dirname "$0")/humm_inventory_20260710.txt" "${ARCH}/manifest/" 2>/dev/null || true
  bash "$(dirname "$0")/inventory_dissert.sh" > "${ARCH}/manifest/dissert_inventory_$(date +%Y%m%d).txt" || true
fi

echo
echo "=== Done. Next: bash rsync_from_hummingbird_missing.sh (see HUMMINGBIRD_DATA_ARCHIVE.md) ==="
if ! $EXEC; then
  echo "Re-run with --execute to apply."
fi
