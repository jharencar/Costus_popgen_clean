#!/usr/bin/env bash
# Compress vector-friendly PDFs for a smaller manuscript build.
# Figure 1 (photo-heavy) is rasterized to 300 DPI JPEG instead of PDF compression.
#
# Usage:
#   ./scripts/compress_final_figures_for_manuscript.sh
#   ./scripts/compress_final_figures_for_manuscript.sh "/path/to/Final_figures"
#
# Requires: ghostscript (brew install ghostscript)

set -euo pipefail

FINAL_DIR="${1:-$HOME/Library/CloudStorage/GoogleDrive-jgharencar@ucdavis.edu/My Drive/Costus/PopGen/Figures&Tables/Final_figures}"
FIGURE1_BASENAME="figure_1_ecogeo_iso_summary"
FIGURE1_PDF="${FIGURE1_BASENAME}.pdf"
JPEG_DPI="${JPEG_DPI:-300}"
JPEG_QUALITY="${JPEG_QUALITY:-92}"
# /printer ≈ 300 dpi embedded images; /ebook ≈ 150 dpi (smaller). Use /ebook for max shrink.
PDF_SETTINGS="${PDF_SETTINGS:-/printer}"

GS="$(command -v gs || true)"
if [[ -z "$GS" && -x /opt/homebrew/bin/gs ]]; then
  GS="/opt/homebrew/bin/gs"
fi
if [[ -z "$GS" ]]; then
  echo "ghostscript (gs) not found. Install with: brew install ghostscript" >&2
  exit 1
fi

if [[ ! -d "$FINAL_DIR" ]]; then
  echo "Directory not found: $FINAL_DIR" >&2
  exit 1
fi

BACKUP_DIR="${FINAL_DIR}/originals_uncompressed"
mkdir -p "$BACKUP_DIR"

compress_pdf() {
  local in_pdf="$1"
  local base out_pdf tmp
  base="$(basename "$in_pdf")"
  out_pdf="${FINAL_DIR}/${base}"
  tmp="$(mktemp "${TMPDIR:-/tmp}/gs_compress.XXXXXX.pdf")"

  if [[ ! -f "${BACKUP_DIR}/${base}" ]]; then
    cp "$in_pdf" "${BACKUP_DIR}/${base}"
    echo "  backed up -> originals_uncompressed/${base}"
  fi

  "$GS" -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 \
    -dPDFSETTINGS="${PDF_SETTINGS}" \
    -dNOPAUSE -dQUIET -dBATCH \
    -sOutputFile="$tmp" "$in_pdf"

  mv "$tmp" "$out_pdf"
  echo "  compressed $(basename "$in_pdf")"
}

echo "Final figures directory: $FINAL_DIR"
echo "Ghostscript: $GS (PDF_SETTINGS=${PDF_SETTINGS}, JPEG ${JPEG_DPI} dpi)"
echo ""

# --- Figure 1: high-DPI JPEG (photos), not PDF compression ---
fig1="${FINAL_DIR}/${FIGURE1_PDF}"
if [[ -f "$fig1" ]]; then
  if [[ ! -f "${BACKUP_DIR}/${FIGURE1_PDF}" ]]; then
    cp "$fig1" "${BACKUP_DIR}/${FIGURE1_PDF}"
    echo "Backed up ${FIGURE1_PDF}"
  fi
  out_jpg="${FINAL_DIR}/${FIGURE1_BASENAME}_${JPEG_DPI}dpi.jpg"
  "$GS" -dNOPAUSE -dBATCH -sDEVICE=jpeg \
    -r"${JPEG_DPI}" -dJPEGQ="${JPEG_QUALITY}" \
    -sOutputFile="$out_jpg" "$fig1"
  echo "Wrote ${out_jpg} (${JPEG_DPI} dpi, quality ${JPEG_QUALITY})"
  echo "  (keep or remove ${FIGURE1_PDF} in Word/LaTeX as you prefer)"
else
  echo "Warning: missing ${FIGURE1_PDF}" >&2
fi
echo ""

# --- All other PDFs: lossy re-encode via pdfwrite ---
echo "Compressing PDFs (excluding ${FIGURE1_PDF}):"
shopt -s nullglob
for pdf in "$FINAL_DIR"/*.pdf; do
  base="$(basename "$pdf")"
  [[ "$base" == "$FIGURE1_PDF" ]] && continue
  before="$(stat -f%z "$pdf" 2>/dev/null || stat -c%s "$pdf")"
  compress_pdf "$pdf"
  after="$(stat -f%z "$pdf" 2>/dev/null || stat -c%s "$pdf")"
  pct=$((100 - (after * 100 / before)))
  echo "    ${base}: $(numfmt --to=iec "$before" 2>/dev/null || echo "${before}B") -> $(numfmt --to=iec "$after" 2>/dev/null || echo "${after}B") (~${pct}% smaller)"
done
echo ""
echo "Done. Uncompressed originals: ${BACKUP_DIR}"
