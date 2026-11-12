#!/bin/bash
# Wait for Farm subset-pixy jobs, sync outputs, remake combo + circular-shift.
# Safe to re-run. From repo root:
#   bash 04_introgression_vs_QTL/noNearHybrid_1km/finish_subset_pixy_pipeline.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

JOB10="${JOB10:-40563512}"
JOB500="${JOB500:-40563513}"
OUT10_FARM=/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_10kb_noNearHybrid_1km
OUT500_FARM=/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_500kb_noNearHybrid_1km
OUT10_LOCAL="$ROOT/pixy/pixy_hudson_10kb_noNearHybrid_1km"
OUT500_LOCAL="$ROOT/pixy/pixy_hudson_500kb_noNearHybrid_1km"

echo "Waiting for Farm jobs $JOB10 (10kb) and $JOB500 (500kb)..."
while true; do
  states=$(ssh -o BatchMode=yes -o ConnectTimeout=30 farm \
    "sacct -j ${JOB10},${JOB500} -n -X -o State -P" | sort -u | tr '\n' ' ')
  echo "  states: $states"
  if echo "$states" | grep -Eq 'PENDING|RUNNING|CONFIGURING|COMPLETING'; then
    sleep 60
    continue
  fi
  break
done

ssh -o BatchMode=yes farm "sacct -j ${JOB10},${JOB500} -X --format=JobID,State,ExitCode,Elapsed -P"
fail=$(ssh -o BatchMode=yes farm \
  "sacct -j ${JOB10},${JOB500} -n -X -o State -P" | grep -cv COMPLETED || true)
if [[ "$fail" -ne 0 ]]; then
  echo "ERROR: not all tasks COMPLETED. Check ~/Costus_scripts/pixy_noNearHybrid_1km/slurm-*.out on Farm." >&2
  exit 1
fi

mkdir -p "$OUT10_LOCAL" "$OUT500_LOCAL"
rsync -avP "farm:${OUT10_FARM}/" "$OUT10_LOCAL/"
rsync -avP "farm:${OUT500_FARM}/" "$OUT500_LOCAL/"

echo "Remaking combo 4-panel + circular-shift with subset pixy + subset Dsuite..."
Rscript "$ROOT/04_introgression_vs_QTL/noNearHybrid_1km/remake_combo_and_circular_shift.R"
echo "Done. Outputs in 04_introgression_vs_QTL/noNearHybrid_1km/"
