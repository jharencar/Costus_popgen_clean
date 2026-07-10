#!/bin/bash
#SBATCH --job-name=WatersMass_QTL
#SBATCH -D /home/jgharenc/Costus_scripts
#SBATCH -o /home/jgharenc/Costus_scripts/slurm_log/WatersMass_QTL.%j.txt
#SBATCH -e /home/jgharenc/Costus_scripts/slurm_log/WatersMass_QTL.%j.txt
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=48:00:00
#SBATCH -c 8
#SBATCH --mem=100G

# Farm (UC Davis) — Waters mass-alignment chemistry QTL.
# Hummingbird/kay_lab variant used: partition 512x64, /hb/groups/kay_lab/R_QTL (see git history).
#
# First pass (collaborator): N_PERM=100  (~hours on ~146 traits)
# Full manuscript run:        N_PERM=1000 (~24-48 h)
# Optional multi-QTL:         RUN_MULTI_QTL=1 (10% FDR traits only)
#
# Submit from ~/Costus_scripts:
#   mkdir -p slurm_log
#   sbatch run_WatersMass_chemistry_QTL.slurm
#   N_PERM=1000 sbatch run_WatersMass_chemistry_QTL.slurm
#
# Override paths if needed:
#   GENOPROB_RDS=/path/to/most_dups_genoprob.rds CONDA_ENV=my_r_env sbatch ...
#   WARNING_LOG=/path/to/WatersMass_chemistry_warnings.txt sbatch ...

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-/home/jgharenc/Costus_scripts}"
cd "${WORKDIR}"

export N_PERM="${N_PERM:-100}"
export RUN_MULTI_QTL="${RUN_MULTI_QTL:-0}"
export GENOPROB_RDS="${GENOPROB_RDS:-${WORKDIR}/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds}"
export CONDA_ENV="${CONDA_ENV:-rqtl}"

module load conda
eval "$(conda shell.bash hook)"
if [[ -n "${CONDA_ENV}" ]]; then
  conda activate "${CONDA_ENV}"
else
  # Set CONDA_ENV to an environment with qtl, dplyr, tidyr, readxl (e.g. clone kay_lab R_QTL locally).
  echo "WARNING: CONDA_ENV not set; using base/default conda env. Activate qtl env if Rscript fails." >&2
fi

run_rscript() {
  local script=$1
  if [[ -f "${script}" ]]; then
    Rscript "${script}"
  elif [[ -f "QTL/${script}" ]]; then
    Rscript "QTL/${script}"
  else
    echo "ERROR: ${script} not found in ${WORKDIR} or ${WORKDIR}/QTL/" >&2
    exit 1
  fi
}

echo "WORKDIR=${WORKDIR}"
echo "GENOPROB_RDS=${GENOPROB_RDS}"
echo "N_PERM=${N_PERM} RUN_MULTI_QTL=${RUN_MULTI_QTL}"
echo "Rscript: $(command -v Rscript)"

if [[ ! -f "${GENOPROB_RDS}" ]]; then
  echo "ERROR: GENOPROB_RDS not found: ${GENOPROB_RDS}" >&2
  echo "Copy most_dups_genoprob.rds from alle_vill_QTL (~858 MB). est.rf is NOT sufficient." >&2
  exit 1
fi

run_rscript prepare_WatersMass_chemistry.R
run_rscript WatersMass_chemistry_QTL.R
