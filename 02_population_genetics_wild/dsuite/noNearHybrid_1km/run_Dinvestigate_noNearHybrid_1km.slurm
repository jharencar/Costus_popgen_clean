#!/bin/bash
# Supplementary Dsuite Dinvestigate — exclude allenii/vill within 1 km of hybrids.
# Original w50s25 LASIout run is unchanged; this writes a separate output directory.
#
# Prerequisites (once per machine):
#   1. Stage BAVLasi VCFs (see README.md)
#   2. Install Dsuite: bash install_dsuite_farm.sh
#   3. Build SETS:
#        bash build_LASIout_SETS_noNearHybrid_1km.sh \
#          ${VCF_DIR}/${VCF_PREFIX}_1.vcf.gz
#
# Submit from this directory:
#   sbatch --export=ALL,PATH="$HOME/bin:$PATH" run_Dinvestigate_noNearHybrid_1km.slurm
#
# Environment overrides:
#   VCF_DIR, VCF_PREFIX, SETS_FILE, TRIOS_FILE, OUTDIR, WINDOW, DSUITE_BIN

#SBATCH --job-name=Dinv_noNearHyb
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=24:00:00
#SBATCH -c 1
#SBATCH --mem=120G
#SBATCH --array=1-9
#SBATCH --output=/home/jgharenc/Costus_scripts/slurm_log/Dinv_noNearHyb_%A_%a.out
#SBATCH --error=/home/jgharenc/Costus_scripts/slurm_log/Dinv_noNearHyb_%A_%a.err

set -euo pipefail

export PATH="${HOME}/bin:${PATH}"

# Prefer submit dir: Slurm copies the batch script under /var/spool/slurmd/, so
# dirname(BASH_SOURCE[0]) is wrong for locating SETS/trios next to the script.
WORKDIR="${WORKDIR:-${SLURM_SUBMIT_DIR:-${HOME}/Costus_scripts/dsuite/noNearHybrid_1km}}"
VCF_DIR="${VCF_DIR:-/quobyte/jrigrp/juliah/Costus_data_tmp/VCFs/BAVLasi_LASIref_filtered_scaf}"
VCF_PREFIX="${VCF_PREFIX:-BAVLasi_LASIref_filtered_scaf}"
SETS_FILE="${SETS_FILE:-${WORKDIR}/LASIout_SETS_noNearHybrid_1km.txt}"
TRIOS_FILE="${TRIOS_FILE:-${WORKDIR}/test_trios.txt}"
OUTDIR="${OUTDIR:-/quobyte/jrigrp/juliah/Costus_data_tmp/dsuite_noNearHybrid_1km_w50s25}"
WINDOW="${WINDOW:-50,25}"
DSUITE_BIN="${DSUITE_BIN:-Dsuite}"
CHR="${SLURM_ARRAY_TASK_ID:?Set SLURM_ARRAY_TASK_ID or use --array}"

if ! command -v "${DSUITE_BIN}" >/dev/null 2>&1; then
  echo "Dsuite not found on PATH (${PATH}). Run install_dsuite_farm.sh first." >&2
  exit 1
fi

VCF="${VCF_DIR}/${VCF_PREFIX}_${CHR}.vcf.gz"
if [[ ! -f "${VCF}" ]]; then
  echo "VCF not found: ${VCF}" >&2
  echo "Stage BAVLasi VCFs first (README.md)." >&2
  exit 1
fi
if [[ ! -f "${SETS_FILE}" ]]; then
  echo "SETS file not found: ${SETS_FILE}" >&2
  echo "Run build_LASIout_SETS_noNearHybrid_1km.sh first." >&2
  exit 1
fi
if [[ ! -f "${TRIOS_FILE}" ]]; then
  echo "Trios file not found: ${TRIOS_FILE}" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"
mkdir -p /home/jgharenc/Costus_scripts/slurm_log
cd "${OUTDIR}"

RUN_NAME="LASIout_noNearHybrid_1km_chr${CHR}"
echo "Dsuite Dinvestigate -w ${WINDOW} -n ${RUN_NAME} ${VCF} ${SETS_FILE} ${TRIOS_FILE}"
"${DSUITE_BIN}" Dinvestigate -w "${WINDOW}" -n "${RUN_NAME}" "${VCF}" "${SETS_FILE}" "${TRIOS_FILE}"
echo "scaffold ${CHR} done -> ${OUTDIR}"
