# Dissert mount paths (macOS default).
# Source from other scripts: source "$(dirname "$0")/dissert_common.sh"
#
# Override on Linux: export COSTUS_DISS=/run/media/julia/dissert

DISS="${COSTUS_DISS:-/Volumes/dissert}"
ARCH="${COSTUS_ARCHIVE:-${DISS}/Costus_archive}"

# Free space in GB (macOS df -g; GNU df fallback on Linux)
df_free_gb() {
  local path="$1"
  if df -g "$path" >/dev/null 2>&1; then
    df -g "$path" | awk 'NR==2 {print $4}'
  elif df -BG "$path" >/dev/null 2>&1; then
    df -BG "$path" | awk 'NR==2 {gsub(/G/,""); print $1}'
  else
    df -Pk "$path" | awk 'NR==2 {printf "%.0f\n", $4/1024/1024}'
  fi
}
