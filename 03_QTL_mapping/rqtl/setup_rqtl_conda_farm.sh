#!/bin/bash
# Create conda env for R/qtl on UC Davis Farm.
# Run interactively on a login node (not via sbatch):
#   bash setup_rqtl_conda_farm.sh
#
# Optional: ENV_NAME=my_rqtl bash setup_rqtl_conda_farm.sh

set -euo pipefail

ENV_NAME="${ENV_NAME:-rqtl}"

module load conda
eval "$(conda shell.bash hook)"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Environment '${ENV_NAME}' already exists."
  echo "To recreate: conda env remove -n ${ENV_NAME}"
  echo "To use:      conda activate ${ENV_NAME}"
  exit 0
fi

echo "Creating conda env '${ENV_NAME}' (conda-forge R + qtl + tidyverse helpers)..."

# conda-forge has r-qtl; bioconda not required for basic mapping.
conda create -n "${ENV_NAME}" -y -c conda-forge \
  r-base \
  r-qtl \
  r-dplyr \
  r-tidyr \
  r-readxl

conda activate "${ENV_NAME}"

echo "Testing packages..."
Rscript -e '
pkgs <- c("qtl", "dplyr", "tidyr", "readxl")
ok <- sapply(pkgs, requireNamespace, quietly = TRUE)
print(ok)
if (!all(ok)) stop("Missing package(s)")
cat("qtl version:", as.character(packageVersion("qtl")), "\n")
'

echo ""
echo "Success. Use on Farm:"
echo "  module load conda && eval \"\$(conda shell.bash hook)\" && conda activate ${ENV_NAME}"
echo "  CONDA_ENV=${ENV_NAME} sbatch slurm/run_scantwo_all_traits.slurm"
