#!/usr/bin/env bash
# Stage 09c - extract unique transcript IDs from each exon GTF.
#
# Reads every *.gtf in --exons-dir and writes <basename>_IDs.txt beside it.
#
# Usage:  bash scripts/bash/09c_extract_transcript_ids.sh [--exons-dir DIR] [--log-dir DIR]
#
# Adapted from archive/transcriptome_metadata_original/extract_transcript_ids.sh
#
# Next: 09d_subset_merged_gtfs.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

EXONS_DIR="${EXONS_DIR:-results/09b_exons}"
LOG_DIR="${LOG_DIR:-results/logs}"

usage() {
  cat <<'EOF'
Usage: 09c_extract_transcript_ids.sh [--exons-dir DIR] [--log-dir DIR]

Reads every *.gtf in --exons-dir and writes <basename>_IDs.txt beside it.

  --exons-dir DIR   default results/09b_exons   (env EXONS_DIR)
  --log-dir   DIR   default results/logs        (env LOG_DIR)
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --exons-dir) EXONS_DIR="$2"; shift 2 ;;
    --log-dir)   LOG_DIR="$2";   shift 2 ;;
    --help|-h)   usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [ ! -d "$EXONS_DIR" ]; then
  echo "[ERROR] exons directory not found: $EXONS_DIR" >&2
  echo "[ERROR] run scripts/bash/09b_feelnc_split_exon_gtfs.sh first." >&2
  exit 1
fi

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/transcript_extraction_$(date +%Y%m%d_%H%M%S).log"

log() { echo "$1" | tee -a "$LOG_FILE"; }

log "=== Transcript ID extraction (stage 09c) ==="
log "Started at: $(date)"
log "Exons dir : $EXONS_DIR"
log ""

shopt -s nullglob
gtfs=("$EXONS_DIR"/*.gtf)
if [ ${#gtfs[@]} -eq 0 ]; then
  echo "[ERROR] no .gtf files in $EXONS_DIR" >&2
  exit 1
fi

for gtf_file in "${gtfs[@]}"; do
  base="$(basename "${gtf_file%.gtf}")"
  out="$EXONS_DIR/${base}_IDs.txt"
  log "Processing: $(basename "$gtf_file")"

  grep -o 'transcript_id "[^"]*"' "$gtf_file" \
    | sed 's/transcript_id "//g' \
    | sed 's/"//g' \
    | sort -u > "$out"

  count=$(wc -l < "$out")
  if [ "$count" -eq 0 ]; then
    echo "[ERROR] no transcript_id found in $gtf_file — refusing to continue with an empty ID list." >&2
    exit 1
  fi
  log "  -> $count unique transcript IDs in $(basename "$out")"
done

log ""
log "Completed at: $(date)"
echo "[INFO] stage 09c complete. Log: $LOG_FILE"
echo "[INFO] Next: scripts/bash/09d_subset_merged_gtfs.sh"
