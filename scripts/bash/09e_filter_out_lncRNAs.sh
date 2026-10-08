#!/usr/bin/env bash
# Stage 09e: produce the lncRNA-depleted merged GTF.
#
# Ported from provenance/babi_preprocessing/filter_out_lncRNAs.sh (B. Borges).
# This is the step that produces `filtered_merged_<id>.annotated.gtf`, which
# R01 reads as its `other_rna_gtf` input. It takes the merged annotated GTF and
# removes every record belonging to a transcript in the lncRNA set.
#
# WHAT CHANGED
# ------------
#  * DIR=/home/vmaracaja/shared-projects/sebauz_sr is now --work-dir. Note the
#    original used a DIFFERENT server root from stages 09c/09d
#    (SUZ_sr_2025), so the chain as received never ran end-to-end in one place.
#  * The original wrote <ident>_transcript_ids.txt beside its inputs and those
#    intermediates were never preserved; they now go to --transcripts-dir.
#  * The four hard-coded idents (two HISAT2 + two STAR) are now derived from
#    whatever lncRNA_tx GTFs stage 09d produced. STAR is not carried forward.
#  * set -euo pipefail added.
#
# GAWK IS REQUIRED
# ----------------
# The filter uses GNU awk's THREE-ARGUMENT match(), i.e.
#     match($0, /transcript_id "([^"]+)"/, arr)
# mawk and BSD awk accept only two arguments: against them this either errors
# or silently never populates arr[1], which would pass every line through and
# produce an unfiltered "filtered" GTF. The check below is not optional.

set -euo pipefail

TRANSCRIPTS_DIR="${TRANSCRIPTS_DIR:-results/09d_transcripts}"
MERGED_DIR="${MERGED_DIR:-results/07_gffcompare}"
MERGED_GTF_NAME="${MERGED_GTF_NAME:-merged.annotated.gtf}"
OUT_DIR="${OUT_DIR:-$TRANSCRIPTS_DIR}"

usage() {
  cat <<'EOF'
Usage: 09e_filter_out_lncRNAs.sh [options]

  --transcripts-dir DIR  stage 09d output   (default results/09d_transcripts)
  --merged-dir DIR       gffcompare output  (default results/07_gffcompare)
  --merged-gtf NAME      merged GTF name    (default merged.annotated.gtf)
  --out-dir DIR          default: same as --transcripts-dir

For each <id>.annotated.lncRNA_tx.gtf, writes
filtered_merged_<id>.annotated.gtf containing every record of the merged GTF
whose transcript is NOT in that lncRNA set.
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --transcripts-dir) TRANSCRIPTS_DIR="$2"; OUT_DIR="${OUT_DIR:-$2}"; shift 2 ;;
    --merged-dir)      MERGED_DIR="$2";      shift 2 ;;
    --merged-gtf)      MERGED_GTF_NAME="$2"; shift 2 ;;
    --out-dir)         OUT_DIR="$2";         shift 2 ;;
    --help|-h)         usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# --- hard requirement: GNU awk with 3-argument match() ----------------------
AWK_BIN="${AWK_BIN:-}"
if [ -z "$AWK_BIN" ]; then
  for cand in gawk awk; do
    if command -v "$cand" >/dev/null 2>&1 &&
       echo 'x' | "$cand" '{ if (match($0, /(x)/, a)) exit 0; exit 1 }' >/dev/null 2>&1; then
      AWK_BIN="$cand"; break
    fi
  done
fi
if [ -z "$AWK_BIN" ]; then
  cat >&2 <<'EOF'
[ERROR] GNU awk is required by this stage.

This filter uses the three-argument form of match():
    match($0, /transcript_id "([^"]+)"/, arr)
which only GNU awk implements. Under mawk or BSD awk the capture group is
never populated, every line passes the filter, and the output is an
UNFILTERED GTF that looks correct. That would silently corrupt R01.

Install gawk (it is pinned in envs/environment_qc_align.yml) or set AWK_BIN
to a GNU awk binary.
EOF
  exit 1
fi
echo "[INFO] using GNU awk: $AWK_BIN ($("$AWK_BIN" --version 2>/dev/null | head -1))"

MERGED_GTF="$MERGED_DIR/$MERGED_GTF_NAME"
[ -f "$MERGED_GTF" ] || { echo "[ERROR] merged GTF not found: $MERGED_GTF" >&2; exit 1; }
mkdir -p "$OUT_DIR"

shopt -s nullglob
lnc_files=("$TRANSCRIPTS_DIR"/*.annotated.lncRNA_tx.gtf)
if [ ${#lnc_files[@]} -eq 0 ]; then
  echo "[ERROR] no *.annotated.lncRNA_tx.gtf in $TRANSCRIPTS_DIR — run stage 09d first." >&2
  exit 1
fi

for tx_file in "${lnc_files[@]}"; do
  ident="$(basename "$tx_file" .annotated.lncRNA_tx.gtf)"
  id_list="$OUT_DIR/${ident}_transcript_ids.txt"
  output_file="$OUT_DIR/filtered_merged_${ident}.annotated.gtf"

  "$AWK_BIN" -F 'transcript_id "' '{ if (NF>1) { split($2,a,"\""); print a[1] } }' "$tx_file" \
    | sort -u > "$id_list"
  n_ids=$(wc -l < "$id_list")
  [ "$n_ids" -gt 0 ] || { echo "[ERROR] no transcript IDs extracted from $tx_file" >&2; exit 1; }

  "$AWK_BIN" -v idfile="$id_list" '
    BEGIN { while ((getline id < idfile) > 0) ids[id] = 1; close(idfile) }
    {
      if (match($0, /transcript_id "([^"]+)"/, arr)) {
        if (!(arr[1] in ids)) print $0
      } else {
        print $0
      }
    }
  ' "$MERGED_GTF" > "$output_file"

  n_in=$(wc -l < "$MERGED_GTF")
  n_out=$(wc -l < "$output_file")
  echo "[INFO] $ident: removed $n_ids lncRNA transcripts; $n_in -> $n_out lines"

  if [ "$n_out" -ge "$n_in" ]; then
    echo "[ERROR] $output_file is not smaller than the input. The filter did nothing —" >&2
    echo "[ERROR] this is the classic symptom of a non-GNU awk. Aborting." >&2
    exit 1
  fi
done

echo "[INFO] stage 09e complete -> $OUT_DIR"
echo "[INFO] Next: Rscript scripts/R/01_transcriptome_metadata_annotation.R"
