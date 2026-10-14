#!/usr/bin/env bash
# Stage 09d - pull full records out of the merged annotated GTF.
#
# For each <base>_IDs.txt from stage 09c, writes <base>_tx.gtf containing the
# transcript and exon lines for those IDs, with the full gffcompare attribute
# set (xloc, class_code, cmp_ref, tss_id) that R01 needs.
#
# Usage:  bash scripts/bash/09d_subset_merged_gtfs.sh [--exons-dir DIR]
#                [--merged-dir DIR] [--merged-gtf NAME] [--transcripts-dir DIR]
#
# Adapted from archive/transcriptome_metadata_original/subset_merged_gtfs.sh
#
# Next: 09e_filter_out_lncRNAs.sh

set -euo pipefail

EXONS_DIR="${EXONS_DIR:-results/09b_exons}"
MERGED_DIR="${MERGED_DIR:-results/07_gffcompare}"
TRANSCRIPTS_DIR="${TRANSCRIPTS_DIR:-results/09d_transcripts}"
LOG_DIR="${LOG_DIR:-results/logs}"
MERGED_GTF_NAME="${MERGED_GTF_NAME:-merged.annotated.gtf}"

usage() {
  cat <<'EOF'
Usage: 09d_subset_merged_gtfs.sh [options]

  --exons-dir DIR        stage 09c output     (default results/09b_exons)
  --merged-dir DIR       gffcompare output    (default results/07_gffcompare)
  --merged-gtf NAME      merged GTF basename  (default merged.annotated.gtf)
  --transcripts-dir DIR  output               (default results/09d_transcripts)
  --log-dir DIR          default results/logs

For each <base>_IDs.txt in --exons-dir, writes <base>_tx.gtf into
--transcripts-dir, containing the full records for those transcript IDs.
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --exons-dir)       EXONS_DIR="$2";       shift 2 ;;
    --merged-dir)      MERGED_DIR="$2";      shift 2 ;;
    --merged-gtf)      MERGED_GTF_NAME="$2"; shift 2 ;;
    --transcripts-dir) TRANSCRIPTS_DIR="$2"; shift 2 ;;
    --log-dir)         LOG_DIR="$2";         shift 2 ;;
    --help|-h)         usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

MERGED_GTF="$MERGED_DIR/$MERGED_GTF_NAME"
if [ ! -f "$MERGED_GTF" ]; then
  echo "[ERROR] merged annotated GTF not found: $MERGED_GTF" >&2
  exit 1
fi

mkdir -p "$TRANSCRIPTS_DIR" "$LOG_DIR"
LOG_FILE="$LOG_DIR/gtf_subset_$(date +%Y%m%d_%H%M%S).log"
log() { echo "$1" | tee -a "$LOG_FILE"; }

log "=== GTF subsetting (stage 09d) ==="
log "Started at : $(date)"
log "Merged GTF : $MERGED_GTF"
log "Output dir : $TRANSCRIPTS_DIR"
log ""

shopt -s nullglob
id_files=("$EXONS_DIR"/*_IDs.txt)
if [ ${#id_files[@]} -eq 0 ]; then
  echo "[ERROR] no *_IDs.txt in $EXONS_DIR — run stage 09c first." >&2
  exit 1
fi

for id_file in "${id_files[@]}"; do
  base="$(basename "$id_file" _IDs.txt)"
  out_gtf="$TRANSCRIPTS_DIR/${base}_tx.gtf"
  n_ids=$(wc -l < "$id_file")

  log "Processing : $base"
  log "  IDs      : $n_ids"

  pattern_file="$(mktemp)"
  # shellcheck disable=SC2001
  sed 's/.*/transcript_id "&"/' "$id_file" > "$pattern_file"

  grep -F -f "$pattern_file" "$MERGED_GTF" > "$out_gtf" || true
  rm -f "$pattern_file"

  n_lines=$(wc -l < "$out_gtf")
  n_tx=$(grep -o 'transcript_id "[^"]*"' "$out_gtf" | sort -u | wc -l)
  log "  -> $(basename "$out_gtf"): $n_lines lines, $n_tx transcripts"

  if [ "$n_lines" -eq 0 ]; then
    echo "[ERROR] $out_gtf is empty. The transcript IDs in $id_file do not appear in $MERGED_GTF." >&2
    exit 1
  fi
  if [ "$n_tx" -ne "$n_ids" ]; then
    log "  WARNING: recovered $n_tx of $n_ids transcripts."
  fi
  log ""
done

log "Completed at: $(date)"
echo "[INFO] stage 09d complete. Log: $LOG_FILE"
echo "[INFO] Next: scripts/bash/09e_filter_out_lncRNAs.sh"
