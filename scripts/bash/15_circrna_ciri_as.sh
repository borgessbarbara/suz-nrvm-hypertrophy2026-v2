#!/usr/bin/env bash
# Stage 15: BWA-MEM -> CIRI2 -> CIRI-AS circRNA sequence reconstruction.
#
# Ported from provenance/circRNA_upstream/CIRI_2_run_all_seba.sh.
#
# SAMPLE COVERAGE — DELIBERATELY NOT WIDENED
# ------------------------------------------
# The original loop is `for m in {13..16}`, i.e. SQ13-SQ16 ONLY. There is no
# CIRI2/CIRI-AS loop for H1-H4 or SQ1-SQ8 anywhere in the material received.
# The audit flags this (remaining_open_items.tsv, category CIRI_AS) with the
# instruction: parameterize the generic command for the full sample sheet
# "only if that reflects the actual analysis; otherwise recover the missing
# batch commands."
#
# So --only defaults to SQ13,SQ14,SQ15,SQ16. Running the full sheet is
# possible but must be an explicit choice: pass --all and you will be asked to
# confirm, because doing so silently would change what the methods describe.
#
# NOTE ON fastp: this branch uses DIFFERENT fastp settings from stages 13/14 —
# it adds `-b 140 -B 140 -l 140`, hard-trimming reads to 140 bp. That is
# required by CIRI-AS's fixed-length assumptions and is preserved verbatim.
#
# Reference note: this script uses Rattus_norvegicus.GRCr8.114.gtf, while the
# project lock is GRCr8.115. Both are GRCr8 (releases 114 and 115 share the
# assembly; 113 and earlier are mRatBN7.2), so this is an annotation-release
# gap, not an assembly conflict. See docs/REFERENCE_AUDIT.md.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_circrna

SAMPLE_SHEET="${SAMPLE_SHEET:-$REPO_ROOT/metadata/longRNA_samples.tsv}"
FASTQ_MANIFEST="${FASTQ_MANIFEST:-$REPO_ROOT/metadata/fastq_manifest.tsv}"
OUT_DIR="${OUT_DIR:-results/15_ciri_as}"
SCRATCH="${SCRATCH:-results/15_ciri_as/scratch}"
THREADS="${THREADS:-20}"
GENOME_FASTA="${CIRC_GENOME_FASTA:-}"
GTF="${CIRC_GTF:-}"
CIRI2_PL="${CIRI2_PL:-}"
CIRI_AS_PL="${CIRI_AS_PL:-}"
BWA_MIN_SCORE="${BWA_MIN_SCORE:-19}"   # the original's -T 19
ONLY="SQ13_ln,SQ14_ln,SQ15_ln,SQ16_ln"
RUN_ALL=0
CONFIRM_ALL=0

usage() {
  cat <<'EOF'
Usage: 15_circrna_ciri_as.sh --genome-fasta F --gtf F --ciri2 P --ciri-as P [options]

Required:
  --genome-fasta PATH   historically rat_primary_assembly_orderr8.fa
  --gtf PATH            historically Rattus_norvegicus.GRCr8.114.gtf
  --ciri2 PATH          CIRI2.pl        (v2.0.6; not a conda package)
  --ciri-as PATH        CIRI_AS_v1.2.pl (v1.2;   not a conda package)

Options:
  --samples PATH        default metadata/longRNA_samples.tsv
  --manifest PATH       default metadata/fastq_manifest.tsv
  --out-dir DIR         default results/15_ciri_as
  --threads N           default 20
  --bwa-min-score N     bwa mem -T (default 19, as in the original)
  --only SAMPLE[,...]   default SQ13_ln,SQ14_ln,SQ15_ln,SQ16_ln
  --all                 run every sample in the sheet (requires --yes-widen)
  --yes-widen           acknowledge that --all goes beyond what was originally run
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --genome-fasta)  GENOME_FASTA="$2";  shift 2 ;;
    --gtf)           GTF="$2";           shift 2 ;;
    --ciri2)         CIRI2_PL="$2";      shift 2 ;;
    --ciri-as)       CIRI_AS_PL="$2";    shift 2 ;;
    --samples)       SAMPLE_SHEET="$2";  shift 2 ;;
    --manifest)      FASTQ_MANIFEST="$2"; shift 2 ;;
    --out-dir)       OUT_DIR="$2";       shift 2 ;;
    --threads)       THREADS="$2";       shift 2 ;;
    --bwa-min-score) BWA_MIN_SCORE="$2"; shift 2 ;;
    --only)          ONLY="$2";          shift 2 ;;
    --all)           RUN_ALL=1;          shift ;;
    --yes-widen)     CONFIRM_ALL=1;      shift ;;
    --help|-h)       usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

for v in GENOME_FASTA GTF CIRI2_PL CIRI_AS_PL; do
  [ -n "${!v}" ] || { echo "[ERROR] missing required argument for $v (see --help)" >&2; exit 1; }
done
for f in "$GENOME_FASTA" "$GTF" "$CIRI2_PL" "$CIRI_AS_PL" "$SAMPLE_SHEET" "$FASTQ_MANIFEST"; do
  [ -f "$f" ] || { echo "[ERROR] not found: $f" >&2; exit 1; }
done

if [ "$RUN_ALL" -eq 1 ] && [ "$CONFIRM_ALL" -ne 1 ]; then
  cat >&2 <<'EOF'
[ERROR] --all was given without --yes-widen.

The original CIRI2/CIRI-AS loop covered SQ13-SQ16 only. Running every sample
produces reconstruction output for samples the published analysis never had,
which changes what the methods section describes. If that is intended, pass
--yes-widen as well and record the decision in docs/OPEN_ITEMS.md.
EOF
  exit 1
fi

all_samples=$(awk -F'\t' 'NR>1 && $1!="" {print $1}' "$SAMPLE_SHEET")
if [ "$RUN_ALL" -eq 1 ]; then
  samples="$all_samples"
  echo "[WARN] running ALL samples - this is wider than the original SQ13-SQ16 loop."
else
  samples=$(echo "$all_samples" | tr ' ' '\n' | grep -Fx -f <(echo "$ONLY" | tr ',' '\n') || true)
  [ -n "$samples" ] || {
    echo "[ERROR] --only matched no sample_id in $SAMPLE_SHEET" >&2
    echo "[ERROR] sheet contains: $(echo "$all_samples" | tr '\n' ' ')" >&2
    exit 1; }
fi

mkdir -p "$OUT_DIR/CIRI2_results" "$OUT_DIR/AS_results" "$OUT_DIR/sam_bwa_alignment" "$SCRATCH"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"; SCRATCH="$(cd "$SCRATCH" && pwd)"

fastq_for() {
  awk -F'\t' -v s="$1" -v r="$2" '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    $h["sample_id"]==s && $h["read"]==r && $h["library_type"]=="longRNA" {
      print ($(h["server_path"]) != "" && $(h["server_path"]) != "TO_FILL") \
            ? $(h["server_path"]) : $(h["exact_fastq_filename"]); exit }
  ' "$FASTQ_MANIFEST"
}

echo "[INFO] CIRI2/CIRI-AS stage 15"
echo "[INFO]   samples : $(echo "$samples" | tr '\n' ' ')"

for sample in $samples; do
  echo "=============================================================="
  echo "[INFO] sample: $sample"
  r1_gz=$(fastq_for "$sample" R1); r2_gz=$(fastq_for "$sample" R2)
  [ -n "$r1_gz" ] && [ -n "$r2_gz" ] || { echo "[ERROR] no R1/R2 for $sample" >&2; exit 1; }

  work="$SCRATCH/$sample"; rm -rf "$work"; mkdir -p "$work"
  gunzip -k -c "$r1_gz" > "$work/${sample}_1.fq"
  gunzip -k -c "$r2_gz" > "$work/${sample}_2.fq"

  # NOTE the extra -b/-B/-l 140 relative to stages 13 and 14.
  fastp -i "$work/${sample}_1.fq" -I "$work/${sample}_2.fq" \
        -o "$work/${sample}_1.trimmed.fq" -O "$work/${sample}_2.trimmed.fq" \
        --detect_adapter_for_pe -f 10 -F 10 -3 \
        --cut_tail_mean_quality 28 -q 30 -w 8 -b 140 -B 140 -l 140 \
        --json "$OUT_DIR/${sample}_fastp.json" \
        --html "$OUT_DIR/${sample}_fastp.html"
  rm -f "$work/${sample}_1.fq" "$work/${sample}_2.fq"

  sam="$OUT_DIR/sam_bwa_alignment/${sample}.sam"
  ciri="$OUT_DIR/CIRI2_results/${sample}.ciri"

  echo "[INFO] bwa mem -T $BWA_MIN_SCORE"
  bwa mem -t "$THREADS" -T "$BWA_MIN_SCORE" "$GENOME_FASTA" \
      "$work/${sample}_1.trimmed.fq" "$work/${sample}_2.trimmed.fq" > "$sam"

  echo "[INFO] CIRI2"
  perl "$CIRI2_PL" -I "$sam" -O "$ciri" -F "$GENOME_FASTA" -A "$GTF" -T "$THREADS"

  echo "[INFO] CIRI-AS"
  perl "$CIRI_AS_PL" -S "$sam" -C "$ciri" -F "$GENOME_FASTA" -A "$GTF" \
       -T "$THREADS" -O "$OUT_DIR/AS_results/${sample}.ASout" -D yes

  rm -rf "$work"
  echo "[INFO] $sample done"
done

echo "[INFO] stage 15 complete -> $OUT_DIR"
echo "[INFO] Next: Rscript scripts/R/04b_circRNA_de_temporal.R"
