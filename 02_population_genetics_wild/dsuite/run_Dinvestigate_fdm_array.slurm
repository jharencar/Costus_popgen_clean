#!/bin/bash
# Farm / local — Dsuite Dinvestigate f_dM (50-variant windows, 25-variant step).
# Original (Hummingbird): LASIout_Dinvestigate_array.sh in _originals/
#
# Submit on Farm (example):
#   export DSUITE_WORKDIR=... VCF_DIR=... SETS_FILE=... TRIOS_FILE=...
#   sbatch --array=1-9 run_Dinvestigate_fdm_array.slurm
#
# Environment:
#   DSUITE_WORKDIR   directory with SETS and trios files (default: cwd)
#   VCF_DIR          filtered BAVLasi LASIref VCFs per scaffold
#   VCF_PREFIX       prefix before _scaf_N.vcf.gz (default: BAVLasi_LASIref_filtered_scaf)
#   SETS_FILE        population assignments (default: LASIout_SETS.txt)
#   TRIOS_FILE       trio definitions (default: test_trios.txt)
#   WINDOW           Dsuite -w argument (default: 50,25)
#   DSUITE_BIN       path to Dsuite executable (default: Dsuite)

#SBATCH --job-name=Dinvestigate_fdm
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=24:00:00
#SBATCH -c 1
#SBATCH --mem=120G
#SBATCH --array=1-9

set -euo pipefail

DSUITE_WORKDIR="${DSUITE_WORKDIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
VCF_DIR="${VCF_DIR:-${DSUITE_WORKDIR}/VCFs}"
VCF_PREFIX="${VCF_PREFIX:-BAVLasi_LASIref_filtered_scaf}"
SETS_FILE="${SETS_FILE:-LASIout_SETS.txt}"
TRIOS_FILE="${TRIOS_FILE:-test_trios.txt}"
WINDOW="${WINDOW:-50,25}"
DSUITE_BIN="${DSUITE_BIN:-Dsuite}"
CHR="${SLURM_ARRAY_TASK_ID:?Set SLURM_ARRAY_TASK_ID or use --array}"

cd "${DSUITE_WORKDIR}"
VCF="${VCF_DIR}/${VCF_PREFIX}_${CHR}.vcf.gz"
RUN_NAME="LASIout_chr${CHR}"

echo "Dsuite Dinvestigate -w ${WINDOW} -n ${RUN_NAME} ${VCF} ${SETS_FILE} ${TRIOS_FILE}"
"${DSUITE_BIN}" Dinvestigate -w "${WINDOW}" -n "${RUN_NAME}" "${VCF}" "${SETS_FILE}" "${TRIOS_FILE}"
echo "scaffold ${CHR} done"
