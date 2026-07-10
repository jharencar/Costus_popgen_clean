#!/usr/bin/env bash
# Run ON HUMMINGBIRD (login node) to list dirs and sizes before archival rsync.
# Usage: bash inventory_hummingbird.sh > humm_inventory_$(date +%Y%m%d).txt

set -uo pipefail

echo "=== Hummingbird data inventory $(date -Iseconds) ==="
echo

PATHS=(
  /hb/groups/kay_lab/Julia_popgen
  /hb/scratch/jharenca/QTL_mapping
  /hb/scratch/jharenca/vill_alle_plate1
  /hb/groups/kay_lab/Davis_data/plates_8.9.10/Data
  /hb/home/jharenca/pop_gen
)

for p in "${PATHS[@]}"; do
  echo "----------------------------------------"
  echo "PATH: $p"
  if [[ -d "$p" ]]; then
    du -sh "$p" 2>/dev/null || true
    echo "Subdirs (maxdepth 2):"
    find "$p" -maxdepth 2 -type d 2>/dev/null | sort | head -80
    echo "(... truncated if >80 subdirs)"
  else
    echo "  MISSING or not accessible"
  fi
  echo
done

echo "=== Julia_popgen top-level sizes ==="
du -sh /hb/groups/kay_lab/Julia_popgen/* 2>/dev/null | sort -h || true

echo
echo "=== Key file globs (count) ==="
for g in \
  "/hb/groups/kay_lab/Julia_popgen/VCFs/*.vcf.gz" \
  "/hb/groups/kay_lab/Julia_popgen/rqtl/*.rds" \
  "/hb/scratch/jharenca/QTL_mapping/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO/no_thin_cated_posteriors_minDif0.7_grp_allF2s_CO0.05/*posterior"
do
  echo -n "$g: "
  ls -1 $g 2>/dev/null | wc -l || echo 0
done
