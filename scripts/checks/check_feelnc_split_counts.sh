#!/usr/bin/env bash
# Regression check for the RECONSTRUCTED step 09b.
#
# 09b_feelnc_split_exon_gtfs.sh rebuilds a step whose original command was
# never recovered. The only evidence we have that the reconstruction is right
# is that it reproduces the transcript counts of the historical outputs.
#
# Historical reference (the original transcriptome-metadata analysis (exons/), CNE_6_24_HISAT2):
#     <id>.annotated.gtf            46030 unique transcripts
#     <id>.annotated.lncRNA.gtf     42127
#     <id>.annotated.otherRNAs.gtf   3822
# with the invariants:
#     lncRNA and otherRNAs are DISJOINT
#     lncRNA fully CONTAINS the FEELnc classifier's lncRNA set (42040)
#     46030 - 42127 - 3822 = 81 transcripts in neither split (the noORF bucket)
#
# Run it against the historical directory to establish the baseline:
#   check_feelnc_split_counts.sh --dir "/path/to/transcriptome-metadata/exons" \
#                                --prefix CNE_6_24_HISAT2
# or against a pipeline run:
#   check_feelnc_split_counts.sh --dir results/09b_exons

set -uo pipefail

DIR="results/09b_exons"
PREFIX="CNE_6_24_HISAT2"
EXP_ALL=46030; EXP_LNC=42127; EXP_OTH=3822; TOL=100
FEELNC=""

usage() {
  cat <<'EOF'
Usage: check_feelnc_split_counts.sh [options]

  --dir DIR             directory holding the three exon GTFs [results/09b_exons]
  --prefix NAME         dataset prefix [CNE_6_24_HISAT2]
  --feelnc-classes F    optional lncRNA_classes_*.txt, to verify containment
  --expect-all N        [46030]
  --expect-lncrna N     [42127]
  --expect-other N      [3822]
  --tolerance N         absolute tolerance per count [100]
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --dir) DIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --feelnc-classes) FEELNC="$2"; shift 2 ;;
    --expect-all) EXP_ALL="$2"; shift 2 ;;
    --expect-lncrna) EXP_LNC="$2"; shift 2 ;;
    --expect-other) EXP_OTH="$2"; shift 2 ;;
    --tolerance) TOL="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

F_ALL="$DIR/${PREFIX}.annotated.gtf"
F_LNC="$DIR/${PREFIX}.annotated.lncRNA.gtf"
F_OTH="$DIR/${PREFIX}.annotated.otherRNAs.gtf"

echo "== FEELnc exon-split regression check =="
echo "   dir    : $DIR"
echo "   prefix : $PREFIX"

missing=0
for f in "$F_ALL" "$F_LNC" "$F_OTH"; do
  [ -f "$f" ] || { echo "  MISSING: $f"; missing=1; }
done
if [ "$missing" -ne 0 ]; then
  echo "SKIP: stage 09b has not been run (or --dir is wrong)."
  exit 0
fi

ids() { grep -o 'transcript_id "[^"]*"' "$1" | sed 's/transcript_id "//;s/"//' | sort -u; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
ids "$F_ALL" > "$T/all"; ids "$F_LNC" > "$T/lnc"; ids "$F_OTH" > "$T/oth"
N_ALL=$(wc -l < "$T/all"); N_LNC=$(wc -l < "$T/lnc"); N_OTH=$(wc -l < "$T/oth")

status=0
cmp_count() {
  local label="$1" got="$2" exp="$3"
  local d=$(( got > exp ? got - exp : exp - got ))
  if [ "$d" -le "$TOL" ]; then
    printf "  ok   %-10s %7d  (expected %7d, delta %d)\n" "$label" "$got" "$exp" "$d"
  else
    printf "  FAIL %-10s %7d  (expected %7d, delta %d > tolerance %d)\n" \
           "$label" "$got" "$exp" "$d" "$TOL"
    status=1
  fi
}
cmp_count annotated "$N_ALL" "$EXP_ALL"
cmp_count lncRNA    "$N_LNC" "$EXP_LNC"
cmp_count otherRNAs "$N_OTH" "$EXP_OTH"

# invariant 1: the two splits are disjoint
OVERLAP=$(comm -12 "$T/lnc" "$T/oth" | wc -l)
if [ "$OVERLAP" -eq 0 ]; then
  echo "  ok   disjoint   lncRNA and otherRNAs share no transcript"
else
  echo "  FAIL disjoint   lncRNA and otherRNAs share $OVERLAP transcripts"
  status=1
fi

# invariant 2: both splits are subsets of the candidate set
OUTSIDE=$(cat "$T/lnc" "$T/oth" | sort -u | comm -23 - "$T/all" | wc -l)
if [ "$OUTSIDE" -eq 0 ]; then
  echo "  ok   subset     both splits are contained in the candidate set"
else
  echo "  FAIL subset     $OUTSIDE transcripts are in a split but not in the candidate set"
  status=1
fi

NEITHER=$(( N_ALL - N_LNC - N_OTH ))
echo "  info noORF      $NEITHER transcripts in the candidate set but in neither split (historical: 81)"

# invariant 3: the classifier's lncRNA set must be fully contained
if [ -n "$FEELNC" ] && [ -f "$FEELNC" ]; then
  awk -F'\t' 'NR>1{print $3}' "$FEELNC" | sort -u > "$T/feelnc"
  N_FE=$(wc -l < "$T/feelnc")
  NOTIN=$(comm -13 "$T/lnc" "$T/feelnc" | wc -l)
  if [ "$NOTIN" -eq 0 ]; then
    echo "  ok   contains   all $N_FE classifier lncRNA transcripts are in the lncRNA split"
  else
    echo "  FAIL contains   $NOTIN of $N_FE classifier lncRNA transcripts are missing from the split"
    status=1
  fi
fi

[ "$status" -eq 0 ] && echo "PASS" || echo "FAILED"
exit $status
