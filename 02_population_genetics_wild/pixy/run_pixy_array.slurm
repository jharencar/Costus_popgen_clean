#!/bin/bash
# PIXY — Hudson Fst, pi, and Dxy (one chromosome per array task).
#
# Submit once per window size. Output directories stay separate.
#   sbatch --export=ALL,WINDOW_SIZE=10000  run_pixy_array.slurm
#   sbatch --export=ALL,WINDOW_SIZE=500000 run_pixy_array.slurm
#
# Each window's Fst is pixy's ratio of summed Hudson components
# (Bhatia et al. 2013). --fst_components keeps hudson_fst_num and
# hudson_fst_den so later cM / QTL summaries can sum those columns
# instead of averaging window Fst values.
#
# The 50 kb Weir–Cockerham job (run_pixy_avonly.sh) is retired.
# Snapshot of that workflow: _originals/wc_fst_original_submission/.
#
# Environment:
#   WINDOW_SIZE   10000 or 500000 (default: 10000)
#   OUTDIR        overrides the default directory for that window size
#   VCF_DIR       directory with per-chr VCFs
#   VCF_PREFIX    prefix before the chromosome number and .vcf.gz
#   POPFILE       pixy populations file
#   N_CORES       pixy --n_cores (default: SLURM_CPUS_PER_TASK or 8)
#   CONDA_ENV     conda env with pixy >= 2.0 (default: pixy)

#SBATCH --job-name=pixy_hudson
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=12:00:00
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --array=1-9

set -euo pipefail

# sbatch inherits the submit shell. Farm's conda module aborts if conda is already active.
unset CONDA_EXE CONDA_PREFIX CONDA_PYTHON_EXE CONDA_SHLVL CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER || true

WINDOW_SIZE="${WINDOW_SIZE:-10000}"
case "${WINDOW_SIZE}" in
  10000)  DEFAULT_OUT=/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_10kb ;;
  500000) DEFAULT_OUT=/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_500kb ;;
  *)
    echo "WINDOW_SIZE must be 10000 or 500000 (got ${WINDOW_SIZE})" >&2
    exit 1
    ;;
esac

VCF_DIR="${VCF_DIR:-/quobyte/jrigrp/juliah/Costus_data_tmp/VCFs}"
VCF_PREFIX="${VCF_PREFIX:-renamed_AVonly_all_sites_mappability_maxDepth_filtered_chr}"
POPFILE="${POPFILE:-/home/jgharenc/Costus_scripts/pixy_whole_genome/avonly_populations.txt}"
OUTDIR="${OUTDIR:-${DEFAULT_OUT}}"
N_CORES="${N_CORES:-${SLURM_CPUS_PER_TASK:-8}}"
CONDA_ENV="${CONDA_ENV:-pixy}"
CHR="${SLURM_ARRAY_TASK_ID:?Set SLURM_ARRAY_TASK_ID or use --array}"

module load conda/latest
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

pixy_help="$(pixy --help 2>&1 || true)"
if ! grep -q -- '--fst_type' <<<"${pixy_help}"; then
  echo "This pixy build has no --fst_type. Install pixy >= 2.0 for the Hudson estimator." >&2
  exit 1
fi
if ! grep -q -- '--fst_components' <<<"${pixy_help}"; then
  echo "This pixy build has no --fst_components. Install pixy >= 2.0 so Fst can be aggregated as a ratio of sums." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"
VCF="${VCF_DIR}/${VCF_PREFIX}${CHR}.vcf.gz"
if [[ ! -f "${VCF}" ]]; then
  echo "VCF not found: ${VCF}" >&2
  exit 1
fi
if [[ ! -f "${POPFILE}" ]]; then
  echo "Populations file not found: ${POPFILE}" >&2
  exit 1
fi

echo "=== pixy Hudson chr${CHR} window ${WINDOW_SIZE} start $(date) ==="
pixy --stats pi fst dxy \
  --fst_type hudson \
  --fst_components \
  --vcf "${VCF}" \
  --populations "${POPFILE}" \
  --window_size "${WINDOW_SIZE}" \
  --n_cores "${N_CORES}" \
  --output_folder "${OUTDIR}" \
  --output_prefix "chr${CHR}"
echo "=== chr${CHR} done -> ${OUTDIR} $(date) ==="
