#!/usr/bin/env bash
# Stage 07 - annotate the reconstructed transcriptome (gffcompare).
#
# Compares the merged StringTie GTF against the reference annotation, assigning
# a class code to every transcript. Produces merged.annotated.gtf and the
# .tmap table that stage 09b and R01 both read.
#
# Usage:  bash scripts/bash/07_gffcompare.sh
# Env:    STRINGTIE_MERGED, GTF_FILE, GFFCOMPARE_OUTPUT, ENSEMBL_RELEASE
#
# Next: 08_gffread_extract_fasta.sh

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# Ensembl release from metadata/reference_release.tsv (see stage 00).
REL_FILE="$REPO_ROOT/metadata/reference_release.tsv"
RELEASE="${ENSEMBL_RELEASE:-}"
[ -z "$RELEASE" ] && [ -f "$REL_FILE" ] && \
  RELEASE="$(awk -F'\t' 'NR>1 && tolower($2) ~ /^(yes|true|1)$/ {print $1; exit}' "$REL_FILE")"
RELEASE="${RELEASE:-116}"
GTF_FILE="${GTF_FILE:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.${RELEASE}.gtf}"
STRINGTIE_MERGED="${STRINGTIE_MERGED:-results/06_stringtie_assembly/merged.gtf}"
GFFCOMPARE_OUTPUT="${GFFCOMPARE_OUTPUT:-results/07_gffcompare}"

# Ensure output directory exists
mkdir -p "$GFFCOMPARE_OUTPUT"

echo "[INFO] gffcompare Stage 07: Annotate Reconstructed Transcriptome"
echo "[INFO] Reference GTF: $GTF_FILE"
echo "[INFO] Merged StringTie GTF: $STRINGTIE_MERGED"
echo "[INFO] Output directory: $GFFCOMPARE_OUTPUT"
echo ""

# Check if input files exist
if [ ! -f "$GTF_FILE" ]; then
  echo "[ERROR] Reference GTF not found: $GTF_FILE" >&2
  exit 1
fi

if [ ! -f "$STRINGTIE_MERGED" ]; then
  echo "[ERROR] Merged StringTie GTF not found: $STRINGTIE_MERGED" >&2
  echo "[ERROR] Did stage 06 (StringTie) complete successfully?" >&2
  exit 1
fi

# Run gffcompare
echo "[INFO] Running gffcompare..."

gffcompare \
  -r "$GTF_FILE" \
  -o "$GFFCOMPARE_OUTPUT/merged" \
  "$STRINGTIE_MERGED"

# The main output is merged.annotated.gtf which has class_code assignments
ANNOTATED_GTF="$GFFCOMPARE_OUTPUT/merged.annotated.gtf"

if [ -f "$ANNOTATED_GTF" ]; then
  echo "[INFO] Annotated GTF created: $ANNOTATED_GTF"
  echo "[INFO]"
  echo "[INFO] Transcript class codes:"
  echo "[INFO]   = : exact match to reference"
  echo "[INFO]   c : contained in reference"
  echo "[INFO]   j : novel isoform"
  echo "[INFO]   u : intergenic (novel)"
  echo "[INFO]   x : overlapping (antisense/other strand)"
  echo "[INFO]   i : intronic"
  echo "[INFO]   o : other"
  echo "[INFO]   p : pseudo-isoform (dup fragments)"
  echo "[INFO]"

  # Summary: count novel transcripts (class code 'u')
  NOVEL_COUNT=$(grep 'class_code "u"' "$ANNOTATED_GTF" | wc -l)
  echo "[INFO] Novel intergenic transcripts (class u): $NOVEL_COUNT"
else
  echo "[WARNING] Annotated GTF not found at expected location: $ANNOTATED_GTF"
fi

echo "[INFO] gffcompare Stage 07 complete"
echo "[INFO] All outputs saved to: $GFFCOMPARE_OUTPUT"
