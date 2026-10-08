#!/usr/bin/env bash
# Stage 06 - transcript assembly and merge (StringTie).
#
# Phase 1 assembles each BAM in reference-guided mode (-G); phase 2 merges all
# per-sample GTFs into one reconstructed transcriptome.
#
# STRINGTIE_EXTRA passes additional flags to phase 1 and defaults to empty,
# which is plain guided assembly as the manuscript Methods describe it. The
# invocation actually used for the GRCr8 run was never recorded; the recorded
# runs used --conservative, and two early ones also -f 0.05 -m 200 -t -s 3.0.
# See docs/OPEN_ISSUES.md before choosing.
#
# Usage:  bash scripts/bash/06_stringtie_assemble_merge.sh
#         STRINGTIE_EXTRA="--conservative" bash scripts/bash/06_stringtie_assemble_merge.sh
# Env:    BAM_OUTPUT, GTF_FILE, STRINGTIE_OUTPUT, STRINGTIE_MERGED,
#         STRINGTIE_EXTRA, THREADS, ENSEMBL_RELEASE
#
# Next: 07_gffcompare.sh

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
RELEASE="${RELEASE:-115}"
GTF_FILE="${GTF_FILE:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.${RELEASE}.gtf}"
BAM_OUTPUT="${BAM_OUTPUT:-results/05_hisat2_bam}"
THREADS="${THREADS:-8}"
STRINGTIE_OUTPUT="${STRINGTIE_OUTPUT:-results/06_stringtie_assembly}"
STRINGTIE_MERGED="${STRINGTIE_MERGED:-results/06_stringtie_assembly/merged.gtf}"
STRINGTIE_EXTRA="${STRINGTIE_EXTRA:-}"
read -r -a STRINGTIE_EXTRA_ARR <<< "$STRINGTIE_EXTRA"

# Ensure output directory exists
mkdir -p "$STRINGTIE_OUTPUT"

echo "[INFO] StringTie Stage 06: Transcript Assembly and Merge"
echo "[INFO] Reference GTF: $GTF_FILE"
echo "[INFO] BAM input: $BAM_OUTPUT"
echo "[INFO] Output directory: $STRINGTIE_OUTPUT"
echo "[INFO] Threads: $THREADS"
echo "[INFO] Extra assembly flags: ${STRINGTIE_EXTRA:-<none>}"
echo ""

# Check if GTF exists
if [ ! -f "$GTF_FILE" ]; then
  echo "[ERROR] GTF file not found: $GTF_FILE" >&2
  echo "[ERROR] Run 00_download_references.sh first or provide --gtf flag" >&2
  exit 1
fi

# Find all BAM files
BAM_FILES=$(find "$BAM_OUTPUT" -maxdepth 1 -name "*.bam" 2>/dev/null | sort)

if [ -z "$BAM_FILES" ]; then
  echo "[ERROR] No BAM files found in $BAM_OUTPUT" >&2
  echo "[ERROR] Did stage 05 (samtools) complete successfully?" >&2
  exit 1
fi

FILE_COUNT=$(echo "$BAM_FILES" | wc -l)
echo "[INFO] Found $FILE_COUNT BAM file(s) to process"
echo ""

# Phase 1: Assemble transcripts per sample (guided mode)
echo "[INFO] Phase 1: Per-sample transcript assembly..."
declare -a GTF_LIST
for bam_file in $BAM_FILES; do
  sample_name=$(basename "$bam_file" .bam)
  output_gtf="$STRINGTIE_OUTPUT/${sample_name}.gtf"

  echo "[INFO] Assembling $sample_name..."

  stringtie \
    -p "$THREADS" \
    -G "$GTF_FILE" \
    ${STRINGTIE_EXTRA_ARR[@]+"${STRINGTIE_EXTRA_ARR[@]}"} \
    -o "$output_gtf" \
    "$bam_file"

  GTF_LIST+=("$output_gtf")
done

echo "[INFO] Per-sample assembly complete"
echo ""

# Phase 2: Merge all sample GTFs into a single reference transcriptome
echo "[INFO] Phase 2: Merging all sample GTFs..."

# Create a temporary file with list of GTF files
GTF_MERGE_LIST=$(mktemp)
trap "rm -f $GTF_MERGE_LIST" EXIT

printf '%s\n' "${GTF_LIST[@]}" > "$GTF_MERGE_LIST"

stringtie \
  --merge \
  -G "$GTF_FILE" \
  -p "$THREADS" \
  -o "$STRINGTIE_MERGED" \
  "$GTF_MERGE_LIST"

echo "[INFO] Merged GTF saved to: $STRINGTIE_MERGED"
echo "[INFO] StringTie Stage 06 complete"
