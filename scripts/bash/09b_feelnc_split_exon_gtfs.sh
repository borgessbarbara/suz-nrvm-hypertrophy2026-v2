#!/usr/bin/env bash
# Stage 09b - materialise the exons/ layout from the FEELnc output.
#
# Writes, per dataset:
#   <prefix>.annotated.gtf             the codpot candidate set
#   <prefix>.annotated.lncRNA.gtf      coding-potential lncRNA calls
#   <prefix>.annotated.otherRNAs.gtf   the complement
#
# Usage:  bash scripts/bash/09b_feelnc_split_exon_gtfs.sh [--help]
# Env:    FEELNC_OUTPUT, EXONS_DIR, PREFIX, INCLUDE_NOORF_IN_OTHER
#
# RECONSTRUCTED: the original command for this step was never recorded.
# Verify with scripts/checks/check_feelnc_split_counts.sh, which compares the
# transcript counts against the historical outputs.
# See docs/CHANGES_FROM_ORIGINALS.md and docs/OPEN_ISSUES.md.
#
# Next: 09c_extract_transcript_ids.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

FEELNC_OUTPUT="${FEELNC_OUTPUT:-results/09_feelnc_classify}"
EXONS_DIR="${EXONS_DIR:-results/09b_exons}"
PREFIX="${PREFIX:-CNE_6_24_HISAT2}"
# Keep .noORF transcripts with the otherRNAs (they are not lncRNA calls).
# Set to 0 to reproduce the historical split exactly, where 81 transcripts
# appear in the candidate set but in neither branch.
INCLUDE_NOORF_IN_OTHER="${INCLUDE_NOORF_IN_OTHER:-0}"

usage() {
  cat <<'EOF'
Usage: 09b_feelnc_split_exon_gtfs.sh

Environment variables:
  FEELNC_OUTPUT            stage 09 output dir   (default results/09_feelnc_classify)
  EXONS_DIR                output dir            (default results/09b_exons)
  PREFIX                   dataset prefix        (default CNE_6_24_HISAT2)
  INCLUDE_NOORF_IN_OTHER   1 to fold .noORF transcripts into otherRNAs
                           (default 0, reproducing the historical split)

Outputs:
  $EXONS_DIR/<PREFIX>.annotated.gtf
  $EXONS_DIR/<PREFIX>.annotated.lncRNA.gtf
  $EXONS_DIR/<PREFIX>.annotated.otherRNAs.gtf
EOF
}
[ "${1:-}" = "--help" ] && { usage; exit 0; }

CODPOT_DIR="$FEELNC_OUTPUT/codpot"
CANDIDATES="$FEELNC_OUTPUT/${PREFIX}.annotated.gtf"
CODPOT_LNC="$CODPOT_DIR/${PREFIX}.lncRNA.gtf"
CODPOT_MRNA="$CODPOT_DIR/${PREFIX}.mRNA.gtf"
CODPOT_NOORF="$CODPOT_DIR/${PREFIX}.noORF.gtf"

for f in "$CANDIDATES" "$CODPOT_LNC" "$CODPOT_MRNA"; do
  if [ ! -f "$f" ]; then
    echo "[ERROR] stage 09 output missing: $f" >&2
    echo "[ERROR] run scripts/bash/09_feelnc_classify.sh first." >&2
    exit 1
  fi
done

mkdir -p "$EXONS_DIR"

OUT_ALL="$EXONS_DIR/${PREFIX}.annotated.gtf"
OUT_LNC="$EXONS_DIR/${PREFIX}.annotated.lncRNA.gtf"
OUT_OTH="$EXONS_DIR/${PREFIX}.annotated.otherRNAs.gtf"

echo "[INFO] stage 09b  (RECONSTRUCTED — see header)"
cp "$CANDIDATES" "$OUT_ALL"
cp "$CODPOT_LNC" "$OUT_LNC"

if [ "$INCLUDE_NOORF_IN_OTHER" = "1" ] && [ -f "$CODPOT_NOORF" ]; then
  echo "[INFO] folding .noORF transcripts into otherRNAs"
  cat "$CODPOT_MRNA" "$CODPOT_NOORF" > "$OUT_OTH"
else
  cp "$CODPOT_MRNA" "$OUT_OTH"
fi

count_tx() {
  grep -o 'transcript_id "[^"]*"' "$1" 2>/dev/null \
    | sed 's/transcript_id "//;s/"//' | sort -u | wc -l
}

N_ALL=$(count_tx "$OUT_ALL")
N_LNC=$(count_tx "$OUT_LNC")
N_OTH=$(count_tx "$OUT_OTH")

echo "[INFO] unique transcripts"
printf '[INFO]   %-14s %8s   (historical: 46030)\n' "annotated"  "$N_ALL"
printf '[INFO]   %-14s %8s   (historical: 42127)\n' "lncRNA"     "$N_LNC"
printf '[INFO]   %-14s %8s   (historical:  3822)\n' "otherRNAs"  "$N_OTH"

# Disjointness is a hard invariant: the historical files share no transcript.
OVERLAP=$(comm -12 \
  <(grep -o 'transcript_id "[^"]*"' "$OUT_LNC" | sed 's/transcript_id "//;s/"//' | sort -u) \
  <(grep -o 'transcript_id "[^"]*"' "$OUT_OTH" | sed 's/transcript_id "//;s/"//' | sort -u) \
  | wc -l)
if [ "$OVERLAP" -ne 0 ]; then
  echo "[ERROR] lncRNA and otherRNAs share $OVERLAP transcripts; they must be disjoint." >&2
  exit 1
fi
echo "[INFO]   lncRNA/otherRNAs disjoint: yes"

cat > "$EXONS_DIR/RECONSTRUCTION_NOTE.txt" <<EOF
These files were produced by scripts/bash/09b_feelnc_split_exon_gtfs.sh, a
RECONSTRUCTION of a step whose original command was not recovered.

Historical reference counts (from the original transcriptome-metadata analysis (exons/), CNE_6_24_HISAT2):
  annotated  46030 unique transcripts
  lncRNA     42127
  otherRNAs   3822

This run:
  annotated  $N_ALL
  lncRNA     $N_LNC
  otherRNAs  $N_OTH

Verify with scripts/checks/check_feelnc_split_counts.sh before using these
outputs for anything that reaches the manuscript.
EOF

echo "[INFO] stage 09b complete -> $EXONS_DIR"
echo "[INFO] Next: scripts/bash/09c_extract_transcript_ids.sh"
