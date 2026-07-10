#!/bin/bash
# PIXY — pi, Fst, Dxy in 10 kb windows (one chromosome per array task).
# Original (Hummingbird): run_pixy.sh in _originals/
#
# Environment:
#   PIXY_WORKDIR     output parent directory (default: cwd)
#   VCF_DIR          directory with per-chr VCFs
#   VCF_PREFIX       prefix before _chrN.vcf.gz
#   POPFILE          populations file for pixy
#   WINDOW_SIZE      bp (default: 10000)
#   N_CORES          pixy --n_cores (default: SLURM_CPUS_PER_TASK or 8)
#   CONDA_ENV        conda env with pixy installed

#SBATCH --job-name=pixy
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=12:00:00
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --array=1-9

set -euo pipefail

PIXY_WORKDIR="${PIXY_WORKDIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
VCF_DIR="${VCF_DIR:-${PIXY_WORKDIR}/VCFs}"
VCF_PREFIX="${VCF_PREFIX:-renamed_AVonly_all_sites_mappability_maxDepth_filtered_chr}"
POPFILE="${POPFILE:-AV_sampleIDs_popfile.txt}"
WINDOW_SIZE="${WINDOW_SIZE:-10000}"
N_CORES="${N_CORES:-${SLURM_CPUS_PER_TASK:-8}}"
CONDA_ENV="${CONDA_ENV:-}"
CHR="${SLURM_ARRAY_TASK_ID:?Set SLURM_ARRAY_TASK_ID or use --array}"

if [[ -n "${CONDA_ENV}" ]]; then
  module load conda 2>/dev/null || true
  eval "$(conda shell.bash hook)"
  conda activate "${CONDA_ENV}"
fi

cd "${PIXY_WORKDIR}"
VCF="${VCF_DIR}/${VCF_PREFIX}${CHR}.vcf.gz"
OUT="${PIXY_WORKDIR}/output_chr${CHR}"

pixy --stats pi fst dxy \
  --vcf "${VCF}" \
  --populations "${POPFILE}" \
  --window_size "${WINDOW_SIZE}" \
  --n_cores "${N_CORES}" \
  --output_folder "${OUT}" \
  --output_prefix "chr${CHR}"

echo "chr${CHR} done -> ${OUT}"
