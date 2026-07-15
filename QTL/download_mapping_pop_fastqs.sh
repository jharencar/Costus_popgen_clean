#!/usr/bin/env bash
# Copy mapping-population FASTQs from Davis plates 8/9/10 by lab ID prefix.
# Expects filenames like 22_033_S127_L001_R1_001.fastq.gz

set -euo pipefail

FASTQ_ROOT="/hb/groups/kay_lab/Davis_data/plates_8.9.10/Data"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ID_FILE="${SCRIPT_DIR}/mapping_pop_lab_ids.txt"

DRY_RUN=0
OUT_DIR=""

usage() {
  echo "Usage: $0 [--dry-run] OUT_DIR" >&2
  echo "  OUT_DIR   Destination directory (created if missing)" >&2
  echo "  --dry-run List files that would be copied; no writes" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage ;;
    *)
      if [[ -n "$OUT_DIR" ]]; then
        usage
      fi
      OUT_DIR="$1"
      shift
      ;;
  esac
done

[[ -n "$OUT_DIR" ]] || usage
[[ -f "$ID_FILE" ]] || { echo "Missing ID file: $ID_FILE" >&2; exit 1; }
[[ -d "$FASTQ_ROOT" ]] || { echo "FASTQ root not found: $FASTQ_ROOT" >&2; exit 1; }

if [[ "$DRY_RUN" -eq 0 ]]; then
  mkdir -p "$OUT_DIR"
fi

total_files=0
missing_ids=()
wrong_count_ids=()

while IFS= read -r lab_id || [[ -n "$lab_id" ]]; do
  [[ -z "$lab_id" ]] && continue
  [[ "$lab_id" =~ ^# ]] && continue

  mapfile -t matches < <(find "$FASTQ_ROOT" -name "${lab_id}_S*_R*.fastq.gz" 2>/dev/null | sort)

  n=${#matches[@]}
  if [[ "$n" -eq 0 ]]; then
    missing_ids+=("$lab_id")
    continue
  fi
  if [[ "$n" -ne 2 ]]; then
    wrong_count_ids+=("$lab_id ($n files)")
  fi

  for f in "${matches[@]}"; do
    if [[ "$DRY_RUN" -eq 1 ]]; then
      echo "$f"
    else
      cp -n "$f" "$OUT_DIR/"
    fi
    total_files=$((total_files + 1))
  done
done < "$ID_FILE"

echo "IDs in list: $(grep -cve '^[[:space:]]*$' -e '^#' "$ID_FILE" || true)" >&2
echo "FASTQ files ${DRY_RUN:+would be }found: $total_files" >&2

if [[ ${#missing_ids[@]} -gt 0 ]]; then
  echo "Missing IDs (0 files): ${#missing_ids[@]}" >&2
  printf '  %s\n' "${missing_ids[@]}" >&2
fi

if [[ ${#wrong_count_ids[@]} -gt 0 ]]; then
  echo "Unexpected file count (not 2): ${#wrong_count_ids[@]}" >&2
  printf '  %s\n' "${wrong_count_ids[@]}" >&2
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "(dry-run: no files copied)" >&2
fi
