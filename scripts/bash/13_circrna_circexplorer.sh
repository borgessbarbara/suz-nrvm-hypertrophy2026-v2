#!/usr/bin/env bash
# Stage 13: circRNA detection with the CIRCexplorer / CLEAR branch.
#
# Consolidates the three batch-specific scripts received from the circRNA
# co-author into one sample-sheet-driven loop:
#   provenance/circRNA_upstream/CIRCexp3_run_all_2023_SebaNECardio.sh    H1-H4
#   provenance/circRNA_upstream/CIRCexp3_run_all_2023_SebaNECardio_2.sh  SQ1-SQ8
#   provenance/circRNA_upstream/CIRCexp3_run_all_2023_SebaNECardio_3.sh  SQ13-SQ16
#
# TOOL NAMING — READ BEFORE EDITING THE METHODS TEXT
# --------------------------------------------------
# The methods prose says "CIRCexplorer3 v1.0.0". The recovered scripts invoke
# `CIRCexplorer2 parse` and `CIRCexplorer2 annotate` together with CLEAR's
# run3.py and circ_quant.py. The string "circexplorer3" appears only as a
# directory name. The actual commands are preserved here unchanged; see
# docs/METHODS_NOTES.md.
#
# WHAT CHANGED FROM THE ORIGINALS
# -------------------------------
#  * One loop over the sample sheet instead of three hard-coded ranges.
#  * Per-sample scratch directory. The originals all decompressed into one
#    shared directory and then globbed `*.fq`, so two concurrent runs would
#    pick up each other's reads. They were only safe to run strictly one at
#    a time, and nothing said so.
#  * Sample name derived from the sample sheet instead of the fragile
#    `${read1: -21:-19}` substring arithmetic.
#  * set -euo pipefail, and `rm -rf` targets are now validated before use.
#  * All absolute /media/storage2/jpablo/... paths are CLI arguments.
#
# fastp settings are unchanged:
#   --detect_adapter_for_pe -f 10 -F 10 -3 --cut_tail_mean_quality 28 -q 30 -w 8

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_circrna

SAMPLE_SHEET="${SAMPLE_SHEET:-$REPO_ROOT/metadata/longRNA_samples.tsv}"
FASTQ_MANIFEST="${FASTQ_MANIFEST:-$REPO_ROOT/metadata/fastq_manifest.tsv}"
OUT_DIR="${OUT_DIR:-results/13_circexplorer}"
SCRATCH="${SCRATCH:-results/13_circexplorer/scratch}"
THREADS="${THREADS:-20}"
HISAT_THREADS="${HISAT_THREADS:-30}"

GENOME_FASTA="${CIRC_GENOME_FASTA:-}"
REFFLAT="${CIRC_REFFLAT:-}"
STAR_INDEX="${STAR_INDEX:-}"
HISAT_INDEX="${HISAT_INDEX:-}"
CLEAR_SRC="${CLEAR_SRC:-}"

usage() {
  cat <<'EOF'
Usage: 13_circrna_circexplorer.sh --genome-fasta F --refflat F --star-index D \
                                  --hisat-index P --clear-src D [options]

Required:
  --genome-fasta PATH   rat genome FASTA (historically rat_primary_assembly_orderr8.fa)
  --refflat PATH        refFlat annotation (historically rn8.refflat.nochr.txt)
  --star-index DIR      STAR index directory (historically STAR_rn8/)
  --hisat-index PREFIX  HISAT2 index prefix (historically .../Ensembl/Rn8)
  --clear-src DIR       directory containing CLEAR's run3.py and circ_quant.py

Options:
  --samples PATH        sample sheet   (default metadata/longRNA_samples.tsv)
  --manifest PATH       FASTQ manifest (default metadata/fastq_manifest.tsv)
  --out-dir DIR         default results/13_circexplorer
  --threads N           STAR/BWA threads (default 20)
  --hisat-threads N     CLEAR run3.py threads (default 30)
  --only SAMPLE[,...]   restrict to these sample_ids

Reference note: rn8 == GRCr8. The refFlat and indexes are rn8-era resources
that are NOT supplied in any package received; see
metadata/circRNA_input_manifest.tsv and docs/OPEN_ITEMS.md.
EOF
}

ONLY=""
while [ $# -gt 0 ]; do
  case "$1" in
    --genome-fasta) GENOME_FASTA="$2"; shift 2 ;;
    --refflat)      REFFLAT="$2";      shift 2 ;;
    --star-index)   STAR_INDEX="$2";   shift 2 ;;
    --hisat-index)  HISAT_INDEX="$2";  shift 2 ;;
    --clear-src)    CLEAR_SRC="$2";    shift 2 ;;
    --samples)      SAMPLE_SHEET="$2"; shift 2 ;;
    --manifest)     FASTQ_MANIFEST="$2"; shift 2 ;;
    --out-dir)      OUT_DIR="$2";      shift 2 ;;
    --threads)      THREADS="$2";      shift 2 ;;
    --hisat-threads) HISAT_THREADS="$2"; shift 2 ;;
    --only)         ONLY="$2";         shift 2 ;;
    --help|-h)      usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

for v in GENOME_FASTA REFFLAT STAR_INDEX HISAT_INDEX CLEAR_SRC; do
  if [ -z "${!v}" ]; then
    echo "[ERROR] --${v,,} is required (see --help)." | tr '_' '-' >&2
    exit 1
  fi
done
[ -f "$GENOME_FASTA" ] || { echo "[ERROR] genome FASTA not found: $GENOME_FASTA" >&2; exit 1; }
[ -f "$REFFLAT" ]      || { echo "[ERROR] refFlat not found: $REFFLAT" >&2; exit 1; }
[ -d "$STAR_INDEX" ]   || { echo "[ERROR] STAR index dir not found: $STAR_INDEX" >&2; exit 1; }
[ -f "$CLEAR_SRC/run3.py" ] || { echo "[ERROR] run3.py not found in $CLEAR_SRC" >&2; exit 1; }
[ -f "$CLEAR_SRC/circ_quant.py" ] || { echo "[ERROR] circ_quant.py not found in $CLEAR_SRC" >&2; exit 1; }
[ -f "$SAMPLE_SHEET" ] || { echo "[ERROR] sample sheet not found: $SAMPLE_SHEET" >&2; exit 1; }
[ -f "$FASTQ_MANIFEST" ] || { echo "[ERROR] FASTQ manifest not found: $FASTQ_MANIFEST" >&2; exit 1; }

mkdir -p "$OUT_DIR" "$SCRATCH"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
SCRATCH="$(cd "$SCRATCH" && pwd)"

# fastq path for sample_id + read, from the manifest
fastq_for() {
  awk -F'\t' -v s="$1" -v r="$2" '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    $h["sample_id"]==s && $h["read"]==r && $h["library_type"]=="longRNA" {
      print ($h["server_path"] != "" && $(h["server_path"]) != "TO_FILL") \
            ? $(h["server_path"]) : $(h["exact_fastq_filename"]); exit }
  ' "$FASTQ_MANIFEST"
}

samples=$(awk -F'\t' 'NR>1 && $1!="" {print $1}' "$SAMPLE_SHEET")
if [ -n "$ONLY" ]; then
  samples=$(echo "$samples" | tr ' ' '\n' | grep -Fx -f <(echo "$ONLY" | tr ',' '\n') || true)
  [ -n "$samples" ] || { echo "[ERROR] --only matched no sample_id" >&2; exit 1; }
fi

echo "[INFO] CIRCexplorer/CLEAR stage 13"
echo "[INFO]   samples : $(echo "$samples" | wc -w)"
echo "[INFO]   out     : $OUT_DIR"

for sample in $samples; do
  echo "=============================================================="
  echo "[INFO] sample: $sample"

  r1_gz=$(fastq_for "$sample" R1)
  r2_gz=$(fastq_for "$sample" R2)
  if [ -z "$r1_gz" ] || [ -z "$r2_gz" ]; then
    echo "[ERROR] no R1/R2 entry for $sample in $FASTQ_MANIFEST" >&2; exit 1
  fi
  for f in "$r1_gz" "$r2_gz"; do
    [ -f "$f" ] || { echo "[ERROR] FASTQ not found: $f" >&2; exit 1; }
  done

  # per-sample scratch: the originals shared one directory and globbed *.fq
  work="$SCRATCH/$sample"
  rm -rf "$work"; mkdir -p "$work"

  echo "[INFO] decompressing"
  gunzip -k -c "$r1_gz" > "$work/${sample}_1.fq"
  gunzip -k -c "$r2_gz" > "$work/${sample}_2.fq"

  echo "[INFO] fastp"
  fastp -i "$work/${sample}_1.fq" -I "$work/${sample}_2.fq" \
        -o "$work/${sample}_1.trimmed.fq" -O "$work/${sample}_2.trimmed.fq" \
        --detect_adapter_for_pe -f 10 -F 10 -3 \
        --cut_tail_mean_quality 28 -q 30 -w 8 \
        --json "$OUT_DIR/${sample}_fastp.json" \
        --html "$OUT_DIR/${sample}_fastp.html"
  rm -f "$work/${sample}_1.fq" "$work/${sample}_2.fq"

  read1="$work/${sample}_1.trimmed.fq"
  read2="$work/${sample}_2.trimmed.fq"
  sdir="$OUT_DIR/${sample}_circexplorer_results"
  mkdir -p "$sdir"
  hisat_out="$work/clear_hisat"

  echo "[INFO] CLEAR run3.py (HISAT2 mapping)"
  python "$CLEAR_SRC/run3.py" -1 "$read1" -2 "$read2" \
         -g "$GENOME_FASTA" -i "$HISAT_INDEX" -G "$REFFLAT" \
         -o "$hisat_out" -p "$HISAT_THREADS"

  unmapped="$hisat_out/unmapped.fq"
  align_bam="$hisat_out/align.bam"
  [ -f "$unmapped" ]  || { echo "[ERROR] CLEAR produced no $unmapped" >&2; exit 1; }
  [ -f "$align_bam" ] || { echo "[ERROR] CLEAR produced no $align_bam" >&2; exit 1; }

  prefix="$work/${sample}_CIRCRNA_"
  echo "[INFO] STAR (chimeric detection on unmapped reads)"
  STAR --chimSegmentMin 10 --runThreadN "$THREADS" \
       --genomeDir "$STAR_INDEX" --readFilesIn "$unmapped" \
       --outFileNamePrefix "$prefix"

  echo "[INFO] CIRCexplorer2 parse"
  CIRCexplorer2 parse -t STAR "${prefix}Chimeric.out.junction" \
    > "$sdir/${sample}_CIRCexplorer2_parse.log"

  # CIRCexplorer2 parse writes back_spliced_junction.bed into the CWD
  bsj="back_spliced_junction.bed"
  [ -f "$bsj" ] || bsj="$(dirname "${prefix}Chimeric.out.junction")/back_spliced_junction.bed"
  [ -f "$bsj" ] || { echo "[ERROR] back_spliced_junction.bed not produced" >&2; exit 1; }

  echo "[INFO] CIRCexplorer2 annotate"
  CIRCexplorer2 annotate -r "$REFFLAT" -g "$GENOME_FASTA" -b "$bsj" \
    -o "$sdir/${sample}_CIRCRNA_circularRNA_known.txt" \
    > "$sdir/${sample}_CIRCRNA_CIRCexplorer2_annotate.log"

  echo "[INFO] CLEAR circ_quant.py"
  python "$CLEAR_SRC/circ_quant.py" \
         -c "$sdir/${sample}_CIRCRNA_circularRNA_known.txt" \
         -b "$align_bam" -t -r "$REFFLAT" \
         -o "$sdir/${sample}_quant.txt"

  cp -f "$bsj" "$sdir/${sample}_back_spliced_junction.bed"

  # the originals removed the 32 GB Aligned.out.sam; keep that behaviour
  rm -rf "$work"
  echo "[INFO] $sample done -> $sdir"
done

echo "[INFO] stage 13 complete -> $OUT_DIR"
