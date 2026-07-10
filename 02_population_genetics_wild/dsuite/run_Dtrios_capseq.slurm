#!/bin/bash
# Dsuite Dtrios — Ecuador capture ABBA-BABA test (methods).
# Original (Hummingbird): Dtrios_capseq_LASIout.sh in _originals/
#
# Environment:
#   DSUITE_WORKDIR   working directory (default: cwd)
#   VCF_FILE         input VCF
#   SETS_FILE        population sets (default: capseq_LASIout_SETS.txt)
#   TREE_FILE        Newick tree (default: tree.nwk)
#   JK_BLOCKS        jackknife blocks (default: 20)
#   OUT_PREFIX       --out-prefix (default: capseq_LASIout)
#   DSUITE_BIN       path to Dsuite (default: Dsuite)

#SBATCH --job-name=capseq_Dtrios
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=4:00:00
#SBATCH -c 4
#SBATCH --mem=16G

set -euo pipefail

DSUITE_WORKDIR="${DSUITE_WORKDIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
VCF_FILE="${VCF_FILE:?Set VCF_FILE to capture-seq VCF}"
SETS_FILE="${SETS_FILE:-capseq_LASIout_SETS.txt}"
TREE_FILE="${TREE_FILE:-tree.nwk}"
JK_BLOCKS="${JK_BLOCKS:-20}"
OUT_PREFIX="${OUT_PREFIX:-capseq_LASIout}"
DSUITE_BIN="${DSUITE_BIN:-Dsuite}"

cd "${DSUITE_WORKDIR}"
echo "${DSUITE_BIN} Dtrios -k ${JK_BLOCKS} -t ${TREE_FILE} --out-prefix=${OUT_PREFIX} ${VCF_FILE} ${SETS_FILE}"
"${DSUITE_BIN}" Dtrios -k "${JK_BLOCKS}" -t "${TREE_FILE}" \
  --out-prefix="${OUT_PREFIX}" --run-name TEST \
  "${VCF_FILE}" "${SETS_FILE}"
