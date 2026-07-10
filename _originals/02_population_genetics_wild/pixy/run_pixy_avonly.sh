#!/bin/bash
#SBATCH --job-name=pixy_AVonly
#SBATCH -D /home/jgharenc/Costus_scripts/pixy_whole_genome
#SBATCH -o /home/jgharenc/Costus_scripts/slurm_log/pixy_AVonly.%A_%a.txt
#SBATCH -e /home/jgharenc/Costus_scripts/slurm_log/pixy_AVonly.%A_%a.txt
#SBATCH -A jrigrp
#SBATCH -p high
#SBATCH --time=12:00:00
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --array=1-9

# Genome-wide pi / dxy / Fst for the AVonly all-sites data (ALLE vs VILL).
#
# ARRAY JOB: one chromosome per array task (chr1-9), all running in parallel.
# Aggregation is a separate step (all chromosomes must finish first). Submit both with
# a dependency so aggregation runs automatically once every task succeeds:
#
#   aid=$(sbatch --parsable run_pixy_avonly.sh)
#   sbatch --dependency=afterok:${aid} aggregate_pixy.sh
#
# (Or, since aggregation is instant and pure-Python, just run it by hand afterwards:
#   python3 aggregate_pixy_genomewide.py /quobyte/jrigrp/juliah/Costus_data_tmp/pixy )

set -euo pipefail

# ------------------------- paths (edit if needed) -------------------------
VCF_DIR=/quobyte/jrigrp/juliah/Costus_data_tmp/VCFs
POPFILE=/home/jgharenc/Costus_scripts/pixy_whole_genome/avonly_populations.txt
OUTDIR=/quobyte/jrigrp/juliah/Costus_data_tmp/pixy
VCF_PREFIX=renamed_AVonly_all_sites_mappability_maxDepth_filtered_chr
WINDOW=50000          # only affects internal parallelism, not the aggregated result
NCORES=8              # keep in sync with #SBATCH -c
# --------------------------------------------------------------------------

# sbatch inherits the submit shell's environment; if conda was active when you
# submitted, CONDA_EXE is set and the farm conda/latest module aborts. Clear it first.
unset CONDA_EXE CONDA_PREFIX CONDA_PYTHON_EXE CONDA_SHLVL CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER

module load conda/latest
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pixy

mkdir -p "${OUTDIR}"

i=${SLURM_ARRAY_TASK_ID}
done_flag="${OUTDIR}/avonly_chr${i}.done"

# Resume support: the .done sentinel is written only after pixy exits successfully,
# so a task interrupted before completion is re-run on resubmission rather than skipped.
if [[ -f "${done_flag}" ]]; then
  echo "=== chr${i} already complete ($(cat "${done_flag}")), nothing to do ==="
  exit 0
fi

echo "=== pixy chr${i} start $(date) ==="
pixy --stats pi fst dxy \
  --vcf "${VCF_DIR}/${VCF_PREFIX}${i}.vcf.gz" \
  --populations "${POPFILE}" \
  --window_size "${WINDOW}" \
  --fst_components \
  --n_cores "${NCORES}" \
  --output_folder "${OUTDIR}" \
  --output_prefix "avonly_chr${i}"
date > "${done_flag}"
echo "=== pixy chr${i} done  $(date) ==="
