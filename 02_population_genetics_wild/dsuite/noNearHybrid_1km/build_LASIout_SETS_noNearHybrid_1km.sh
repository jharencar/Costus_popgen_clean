#!/bin/bash
# Build a Dsuite SETS file for the supplementary no-near-hybrid analysis.
#
# Keeps:
#   - all BRAC samples as brac
#   - all LASI samples as Outgroup
#   - allenii / villosissimus IDs listed in alle_vill_include_IDs.txt
# Marks every other VCF sample as xxx (including hybrids and samples within
# 1 km of a hybrid).
#
# Usage:
#   bash build_LASIout_SETS_noNearHybrid_1km.sh /path/to/BAVLasi_..._scaf_1.vcf.gz
#
# Writes LASIout_SETS_noNearHybrid_1km.txt in this directory.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INCLUDE_FILE="${INCLUDE_FILE:-${SCRIPT_DIR}/alle_vill_include_IDs.txt}"
OUT_SETS="${OUT_SETS:-${SCRIPT_DIR}/LASIout_SETS_noNearHybrid_1km.txt}"
VCF="${1:?Usage: $0 /path/to/example.vcf.gz}"

if [[ ! -f "${VCF}" ]]; then
  echo "VCF not found: ${VCF}" >&2
  exit 1
fi
if [[ ! -f "${INCLUDE_FILE}" ]]; then
  echo "Include list not found: ${INCLUDE_FILE}" >&2
  exit 1
fi

module load bcftools 2>/dev/null || true
if ! command -v bcftools >/dev/null 2>&1; then
  echo "bcftools not found" >&2
  exit 1
fi

TMP="$(mktemp)"
trap 'rm -f "${TMP}"' EXIT

# Base ID -> intended pop for kept allenii/vill (ALLE/VILL only).
declare -A KEEP
while read -r id spp; do
  [[ -z "${id:-}" || "${id}" =~ ^# ]] && continue
  KEEP["${id}"]="${spp}"
done < "${INCLUDE_FILE}"

bcftools query -l "${VCF}" > "${TMP}"

n_brac=0
n_alle=0
n_vill=0
n_out=0
n_xxx=0

: > "${OUT_SETS}"
while read -r samp; do
  [[ -z "${samp}" ]] && continue
  base="${samp%%_*}"
  # Strip trailing ab merges if still present in some VCFs (19.606ab -> 19.606).
  base_noab="${base%ab}"
  spp_tag="${samp##*_}"

  set_name="xxx"
  if [[ "${spp_tag}" == "BRAC" ]]; then
    set_name="brac"
    n_brac=$((n_brac + 1))
  elif [[ "${spp_tag}" == "LASI" || "${spp_tag}" == "Outgroup" ]]; then
    set_name="Outgroup"
    n_out=$((n_out + 1))
  elif [[ -n "${KEEP[${base_noab}]+x}" ]]; then
    case "${KEEP[${base_noab}]}" in
      ALLE) set_name="alle"; n_alle=$((n_alle + 1)) ;;
      VILL) set_name="vill"; n_vill=$((n_vill + 1)) ;;
      *) set_name="xxx"; n_xxx=$((n_xxx + 1)) ;;
    esac
  else
    n_xxx=$((n_xxx + 1))
  fi
  printf '%s\t%s\n' "${samp}" "${set_name}" >> "${OUT_SETS}"
done < "${TMP}"

echo "Wrote ${OUT_SETS}"
echo "  brac=${n_brac}  alle=${n_alle}  vill=${n_vill}  Outgroup=${n_out}  xxx=${n_xxx}"
if [[ "${n_alle}" -eq 0 || "${n_vill}" -eq 0 || "${n_brac}" -eq 0 || "${n_out}" -eq 0 ]]; then
  echo "ERROR: expected non-zero counts for brac/alle/vill/Outgroup" >&2
  exit 1
fi
if [[ "${n_alle}" -ne 12 || "${n_vill}" -ne 12 ]]; then
  echo "WARNING: expected 12 alle and 12 vill from the 1 km filter (got alle=${n_alle} vill=${n_vill})" >&2
fi
