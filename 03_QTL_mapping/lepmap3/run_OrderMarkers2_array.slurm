#!/bin/bash
# Lep-MAP3 OrderMarkers2 per chromosome (manual map).
# Original (Hummingbird): lepmap_OrderMarkers2_manMap.sh in _originals/
#
# Environment:
#   LEPMAP_WORKDIR     directory with .call and manual map (default: cwd)
#   LEPMAP_JAR_DIR     Lep-MAP3 bin directory
#   DATA_CALL          genotype call file (default: refV2_thinned_50k_data.call)
#   MANUAL_MAP         manual map file (default: manual_map3.txt)
#   JAVA_MEM           -Xmx for java (default: 80g on Hummingbird; tune for Farm)

#SBATCH --job-name=lepmap_OrderMarkers2
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=24:00:00
#SBATCH -c 2
#SBATCH --mem=80G
#SBATCH --array=1-9

set -euo pipefail

LEPMAP_WORKDIR="${LEPMAP_WORKDIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
LEPMAP_JAR_DIR="${LEPMAP_JAR_DIR:?Set LEPMAP_JAR_DIR to Lep-MAP3 bin/}"
DATA_CALL="${DATA_CALL:-refV2_thinned_50k_data.call}"
MANUAL_MAP="${MANUAL_MAP:-manual_map3.txt}"
CHR="${SLURM_ARRAY_TASK_ID:?Set SLURM_ARRAY_TASK_ID or use --array}"

module load java 2>/dev/null || true
cd "${LEPMAP_WORKDIR}"

java -cp "${LEPMAP_JAR_DIR}" OrderMarkers2 \
  map="${MANUAL_MAP}" \
  data="${DATA_CALL}" \
  scale=2M/N \
  sexAveraged=1 \
  outputPhasedData=0 \
  chromosome="${CHR}" \
  > "order${CHR}_LodLim4.txt"

echo "OrderMarkers2 chr${CHR} done"
