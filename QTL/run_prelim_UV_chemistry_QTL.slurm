#!/bin/bash
#SBATCH --partition 512x64
#SBATCH --qos 512x64
#SBATCH --account 512x64
#SBATCH --job-name=QTL_UV_chem
#SBATCH -e slurm_files/QTL_UV_chem.err
#SBATCH --output=slurm_files/QTL_UV_chem.out
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8

# Run prelim UV chemistry QTL (scanone + 100 permutations per peak + plots).
# load R QTL
module load miniconda3
conda activate /hb/groups/kay_lab/R_QTL

cd /hb/groups/kay_lab/Julia_popgen/rqtl

PROJECT_DIR="${SLURM_SUBMIT_DIR:-.}"
cd "$PROJECT_DIR"

Rscript QTL/prelim_UV_chemistry_QTL_setup.R
