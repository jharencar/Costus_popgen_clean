#!/bin/bash
#SBATCH --job-name=rqtl_scantwo
#SBATCH -D /home/jgharenc/Costus_scripts
#SBATCH -o /home/jgharenc/Costus_scripts/slurm_log/rqtl_scantwo.%A_%a.txt
#SBATCH -e /home/jgharenc/Costus_scripts/slurm_log/rqtl_scantwo.%A_%a.txt
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=72:00:00
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --array=0-4

set -euo pipefail
TRAITS=(SD GR tough Nreabs Chlreabs)
export TRAIT="${TRAITS[$SLURM_ARRAY_TASK_ID]}"
WORKDIR="${SLURM_SUBMIT_DIR:-/home/jgharenc/Costus_scripts}"
cd "${WORKDIR}"
export WORKDIR N_PERM="${N_PERM:-1000}" N_CORES="${SLURM_CPUS_PER_TASK:-8}"
export CONDA_ENV="${CONDA_ENV:-rqtl}"
export CROSS_RDS="${CROSS_RDS:-${WORKDIR}/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_least_dups_est.rf.rds}"
export OUTDIR="${OUTDIR:-${WORKDIR}/scantwo.perms}"
module load conda && eval "$(conda shell.bash hook)" && conda activate "${CONDA_ENV}"
for p in run_scantwo_per_trait.R 03_QTL_mapping/rqtl/run_scantwo_per_trait.R; do
  [[ -f "$p" ]] && Rscript "$p" && exit 0
done
echo "ERROR: run_scantwo_per_trait.R not found" >&2; exit 1
