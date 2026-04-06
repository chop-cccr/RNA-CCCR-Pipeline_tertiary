#!/usr/bin/env bash
# Copy the entire /mnt/isilon/reference tree into a destination (default: ./reference).
# MODE=copy|hardlink (default: copy); DRY=true for preview.
set -euo pipefail

SRC="/mnt/isilon/cccr_bfx/CCCR_Pipelines/RNA-Seq_Pipeline_HPC/reference/"        # <- hard-coded source
DST="${1:-reference}"              # default destination dir in the current folder
MODE="${MODE:-copy}"               # copy | hardlink
DRY="${DRY:-false}"                # true | false

[[ -d "$SRC" ]] || { echo "ERROR: Source not found: $SRC"; exit 2; }
mkdir -p "$DST"

echo "Source : $SRC"
echo "Dest   : $DST"
echo "Mode   : $MODE"
echo "Dry    : $DRY"
echo

if [[ "$DRY" == "true" ]]; then
  echo "[DRY] Would copy entire directory tree: $SRC/ -> $DST/"
  exit 0
fi

case "$MODE" in
  hardlink)
    # Try hard links (same filesystem). If that fails, fall back to real copy.
    if cp -al "$SRC/." "$DST/" 2>/dev/null; then
      echo "Hardlinked tree successfully."
    else
      echo "Hardlink failed; falling back to real copy..."
      cp -a "$SRC/." "$DST/"
      echo "Copied tree successfully."
    fi
    ;;
  copy)
    cp -a "$SRC/." "$DST/"
    echo "Copied tree successfully."
    ;;
  *)
    echo "ERROR: Unknown MODE '$MODE' (use 'copy' or 'hardlink')"
    exit 3
    ;;
esac

