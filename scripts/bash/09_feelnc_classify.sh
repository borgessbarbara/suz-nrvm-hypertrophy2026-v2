#!/usr/bin/env bash
# Stage 09: lncRNA identification and classification (FEELnc)
#
#   1. FEELnc_filter.pl     -> candidate transcript set (exon lines only)
#   2. FEELnc_codpot.pl     -> .lncRNA.gtf / .mRNA.gtf / .noORF.gtf
#   3. FEELnc_classifier.pl -> lncRNA_classes_<prefix>_1kb_100kb.txt
#
# Settings are taken verbatim from archive/raw_scripts_original/bash_scripts/bash_scripts.txt,
# which records four successive FEELnc runs, all with:
#     -b transcript_biotype=protein_coding -b transcript_status=KNOWN --mode=intergenic
# and a classifier window of 1 kb / 100 kb (hence the historical filename
# lncRNA_classes_..._1kb_100kb.txt).
#
# WHICH candidate GTF was final is an OPEN item: the notebook shows the input
# evolving through candidate_lncRNA_NE6_NE24.gtf -> candidate_lncRNAs_NE6_NE24_oct.gtf
# -> cand_lncRNAs_NE6_NE24_annotated.gtf -> cand_lncRNAs_ALLPE_annotated.gtf.
# See docs/OPEN_ITEMS.md.
#
# Downstream: stage 09b normalises these outputs into the exons/ layout that
# the Babi preprocessing chain (09c-09e) and R01 expect.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Ensembl release from metadata/reference_release.tsv (see stage 00).
REL_FILE="$REPO_ROOT/metadata/reference_release.tsv"
RELEASE="${ENSEMBL_RELEASE:-}"
[ -z "$RELEASE" ] && [ -f "$REL_FILE" ] && \
  RELEASE="$(awk -F'\t' 'NR>1 && tolower($2) ~ /^(yes|true|1)$/ {print $1; exit}' "$REL_FILE")"
RELEASE="${RELEASE:-115}"
GTF_FILE="${GTF_FILE:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.${RELEASE}.gtf}"
GENOME_FASTA="${GENOME_FASTA:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.primary_assembly.fa}"
GFFCOMPARE_OUTPUT="${GFFCOMPARE_OUTPUT:-results/07_gffcompare}"
FEELNC_OUTPUT="${FEELNC_OUTPUT:-results/09_feelnc_classify}"
# RESOLVED against the raw command notebook. The Methods are wrong.
#
# archive/raw_scripts_original/bash_scripts/bash_scripts.txt records seven
# FEELnc_classifier runs. The FINAL one - the 6+24 h rn8 run, under the
# heading "Para 6 y 24 hrs", whose output filename is exactly what this
# pipeline consumes - is at line 300:
#
#   FEELnc_classifier.pl -i cand_lncRNAs_CNE624_rn8_8_annotated.gtf.lncRNA.gtf \
#     -a .../Rattus_norvegicus.GRCr8.114.gtf \
#     --window=1000 --maxwindow=100000 \
#     > lncRNA_classes_CNE624_rn8_8_annotated_1kb_100kb.txt
#
# --window=200 appears only in the six EARLIER exploratory runs (lines 199,
# 207, 209, 219, 227), every one of which used maxwindow 10000 or 50000 and
# produced a differently named file that nothing consumes.
#
# The Methods' "200 bp to 100 kb window (--window=200 --maxwindow=100000)"
# pairs the window from the discarded runs with the maxwindow from the final
# one. 1000/100000 is correct. See docs/MANUSCRIPT_ALIGNMENT.md.
FEELNC_WINDOW="${FEELNC_WINDOW:-1000}"
FEELNC_MAXWINDOW="${FEELNC_MAXWINDOW:-100000}"
# Dataset prefix. The historical analysis used CNE_6_24_HISAT2 (6 h + 24 h,
# HISAT2). Only the HISAT2 branch is carried forward; see docs/OPEN_ITEMS.md.
PREFIX="${PREFIX:-CNE_6_24_HISAT2}"

mkdir -p "$FEELNC_OUTPUT"

echo "[INFO] FEELnc stage 09"
echo "[INFO]   reference GTF : $GTF_FILE"
echo "[INFO]   genome FASTA  : $GENOME_FASTA"
echo "[INFO]   merged GTF    : $GFFCOMPARE_OUTPUT/merged.annotated.gtf"
echo "[INFO]   output        : $FEELNC_OUTPUT"
echo "[INFO]   prefix        : $PREFIX"
echo "[INFO]   window        : $FEELNC_WINDOW / $FEELNC_MAXWINDOW"

MERGED_GTF="$GFFCOMPARE_OUTPUT/merged.annotated.gtf"
for f in "$GTF_FILE" "$GENOME_FASTA" "$MERGED_GTF"; do
  if [ ! -f "$f" ]; then
    echo "[ERROR] required input not found: $f" >&2
    exit 1
  fi
done

# --- Step 1: filter ---------------------------------------------------------
# FEELnc_filter emits exon lines only, with a minimal attribute set
# (gene_id / transcript_id / exon_number). That is exactly the format of the
# historical exons/*.gtf files.
CANDIDATES="$FEELNC_OUTPUT/${PREFIX}.annotated.gtf"
echo "[INFO] Step 1/3  FEELnc_filter.pl"
FEELnc_filter.pl \
  -i "$MERGED_GTF" \
  -a "$GTF_FILE" \
  -b transcript_biotype=protein_coding \
  -b transcript_status=KNOWN \
  --monoex=1 \
  --outlog "$FEELNC_OUTPUT/feelnc_filter.log" \
  > "$CANDIDATES"
echo "[INFO]   candidates: $CANDIDATES"

# --- Step 2: coding potential ----------------------------------------------
# Produces <outname>.lncRNA.gtf, <outname>.mRNA.gtf and, where ORF detection
# is inconclusive, <outname>.noORF.gtf.
CODPOT_DIR="$FEELNC_OUTPUT/codpot"
mkdir -p "$CODPOT_DIR"
echo "[INFO] Step 2/3  FEELnc_codpot.pl"
FEELnc_codpot.pl \
  -i "$CANDIDATES" \
  -a "$GTF_FILE" \
  -g "$GENOME_FASTA" \
  -b transcript_biotype=protein_coding \
  -b transcript_status=KNOWN \
  --mode=intergenic \
  --outdir "$CODPOT_DIR" \
  --outname "$PREFIX" \
  2>&1 | tee "$FEELNC_OUTPUT/feelnc_codpot.log"

CODPOT_LNC="$CODPOT_DIR/${PREFIX}.lncRNA.gtf"
if [ ! -f "$CODPOT_LNC" ]; then
  echo "[ERROR] FEELnc_codpot did not produce $CODPOT_LNC" >&2
  echo "[ERROR] contents of $CODPOT_DIR:" >&2
  ls -la "$CODPOT_DIR" >&2
  exit 1
fi
echo "[INFO]   lncRNA calls: $CODPOT_LNC"

# --- Step 3: classifier -----------------------------------------------------
CLASSES="$FEELNC_OUTPUT/lncRNA_classes_${PREFIX}_1kb_100kb.txt"
echo "[INFO] Step 3/3  FEELnc_classifier.pl"
FEELnc_classifier.pl \
  -i "$CODPOT_LNC" \
  -a "$GTF_FILE" \
  --window="$FEELNC_WINDOW" \
  --maxwindow="$FEELNC_MAXWINDOW" \
  > "$CLASSES" \
  2> "$FEELNC_OUTPUT/feelnc_classifier.log"

echo "[INFO]   classes: $CLASSES"
echo "[INFO] FEELnc stage 09 complete."
echo "[INFO] Next: scripts/bash/09b_feelnc_split_exon_gtfs.sh"
