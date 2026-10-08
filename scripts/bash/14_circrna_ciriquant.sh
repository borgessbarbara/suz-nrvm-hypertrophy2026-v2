#!/usr/bin/env bash
# Stage 14: circRNA detection and quantification with CIRIquant.
#
# Consolidates:
#   archive/final_package_v3_original/02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/CIRIquant_all_run_2023_Seba_NE_Cardio_1.sh  H1-H4
#   archive/final_package_v3_original/02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/CIRIquant_all_run_2023_Seba_NE_Cardio_2.sh  SQ1-SQ8
#   archive/final_package_v3_original/02_CIRCRNA_UPSTREAM_RAW_PROVENANCE/CIRIquant_all_run_2023_Seba_NE_Cardio_3.sh  SQ13-SQ16
#
# The CIRIquant invocation is preserved exactly:
#     CIRIquant -t 30 -1 $read1 -2 $read2 --config rn8.yml -o $sample -p $sample
# and so are the fastp settings. The commented-out `--circ ... --tool
# CIRCexplorer2` variant in the originals is NOT enabled here; if the deposited
# results used it, say so in docs/OPEN_ITEMS.md and switch --use-circexplorer-prior on.
#
# BLOCKING INPUT: rn8.yml
# -----------------------
# CIRIquant needs a YAML config naming the genome FASTA, the GTF, and the
# bwa/hisat2/stringtie/samtools binaries. The original points at
# /media/storage2/jpablo/CIRIquant/rn8.yml, which is not in any package
# received. This stage fails loudly rather than guessing; see
# metadata/circRNA_input_manifest.tsv.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_circrna

SAMPLE_SHEET="${SAMPLE_SHEET:-$REPO_ROOT/metadata/longRNA_samples.tsv}"
FASTQ_MANIFEST="${FASTQ_MANIFEST:-$REPO_ROOT/metadata/fastq_manifest.tsv}"
OUT_DIR="${OUT_DIR:-results/14_ciriquant}"
SCRATCH="${SCRATCH:-results/14_ciriquant/scratch}"
THREADS="${THREADS:-30}"
CIRIQUANT_CONFIG="${CIRIQUANT_CONFIG:-}"
CIRC_PRIOR_DIR="${CIRC_PRIOR_DIR:-}"
ONLY=""

usage() {
  cat <<'EOF'
Usage: 14_circrna_ciriquant.sh --config rn8.yml [options]

Required:
  --config PATH             CIRIquant YAML config (historically rn8.yml).
                            NOT SUPPLIED in any package received - must be
                            recovered or recreated and committed.

Options:
  --samples PATH            default metadata/longRNA_samples.tsv
  --manifest PATH           default metadata/fastq_manifest.tsv
  --out-dir DIR             default results/14_ciriquant
  --threads N               default 30
  --only SAMPLE[,...]       restrict to these sample_ids
  --use-circexplorer-prior DIR
                            enable the commented-out variant from the
                            originals: pass each sample's CIRCexplorer2
                            circularRNA_known.txt via --circ/--tool. Only use
                            this if confirmed as what actually produced the
                            deposited tables.
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --config)    CIRIQUANT_CONFIG="$2"; shift 2 ;;
    --samples)   SAMPLE_SHEET="$2";     shift 2 ;;
    --manifest)  FASTQ_MANIFEST="$2";   shift 2 ;;
    --out-dir)   OUT_DIR="$2";          shift 2 ;;
    --threads)   THREADS="$2";          shift 2 ;;
    --only)      ONLY="$2";             shift 2 ;;
    --use-circexplorer-prior) CIRC_PRIOR_DIR="$2"; shift 2 ;;
    --help|-h)   usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [ -z "$CIRIQUANT_CONFIG" ] || [ ! -f "$CIRIQUANT_CONFIG" ]; then
  cat >&2 <<EOF
[ERROR] CIRIquant config not provided or not found: ${CIRIQUANT_CONFIG:-<unset>}

The historical scripts reference
    /media/storage2/jpablo/CIRIquant/rn8.yml
which was never included in any package received from the circRNA co-author.
It is recorded as a blocking open item:
    metadata/circRNA_input_manifest.tsv   (row: CIRIquant config)
    docs/OPEN_ITEMS.md

Recreating it requires knowing exactly which genome FASTA and GTF the circRNA
branch used - which is itself an open item, because the circRNA scripts carry
GRCr8.114 while the project reference lock is GRCr8.115.

Pass --config once that is resolved.
EOF
  exit 1
fi

[ -f "$SAMPLE_SHEET" ]   || { echo "[ERROR] sample sheet not found: $SAMPLE_SHEET" >&2; exit 1; }
[ -f "$FASTQ_MANIFEST" ] || { echo "[ERROR] FASTQ manifest not found: $FASTQ_MANIFEST" >&2; exit 1; }

mkdir -p "$OUT_DIR" "$SCRATCH"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
SCRATCH="$(cd "$SCRATCH" && pwd)"

fastq_for() {
  awk -F'\t' -v s="$1" -v r="$2" '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    $h["sample_id"]==s && $h["read"]==r && $h["library_type"]=="longRNA" {
      print ($(h["server_path"]) != "" && $(h["server_path"]) != "TO_FILL") \
            ? $(h["server_path"]) : $(h["exact_fastq_filename"]); exit }
  ' "$FASTQ_MANIFEST"
}

samples=$(awk -F'\t' 'NR>1 && $1!="" {print $1}' "$SAMPLE_SHEET")
if [ -n "$ONLY" ]; then
  samples=$(echo "$samples" | tr ' ' '\n' | grep -Fx -f <(echo "$ONLY" | tr ',' '\n') || true)
  [ -n "$samples" ] || { echo "[ERROR] --only matched no sample_id" >&2; exit 1; }
fi

echo "[INFO] CIRIquant stage 14"
echo "[INFO]   config  : $CIRIQUANT_CONFIG"
echo "[INFO]   samples : $(echo "$samples" | wc -w)"

for sample in $samples; do
  echo "=============================================================="
  echo "[INFO] sample: $sample"
  r1_gz=$(fastq_for "$sample" R1); r2_gz=$(fastq_for "$sample" R2)
  [ -n "$r1_gz" ] && [ -n "$r2_gz" ] || { echo "[ERROR] no R1/R2 for $sample" >&2; exit 1; }
  for f in "$r1_gz" "$r2_gz"; do
    [ -f "$f" ] || { echo "[ERROR] FASTQ not found: $f" >&2; exit 1; }
  done

  work="$SCRATCH/$sample"; rm -rf "$work"; mkdir -p "$work"
  gunzip -k -c "$r1_gz" > "$work/${sample}_1.fq"
  gunzip -k -c "$r2_gz" > "$work/${sample}_2.fq"

  fastp -i "$work/${sample}_1.fq" -I "$work/${sample}_2.fq" \
        -o "$work/${sample}_1.trimmed.fq" -O "$work/${sample}_2.trimmed.fq" \
        --detect_adapter_for_pe -f 10 -F 10 -3 \
        --cut_tail_mean_quality 28 -q 30 -w 8 \
        --json "$OUT_DIR/${sample}_fastp.json" \
        --html "$OUT_DIR/${sample}_fastp.html"
  rm -f "$work/${sample}_1.fq" "$work/${sample}_2.fq"

  extra=()
  if [ -n "$CIRC_PRIOR_DIR" ]; then
    prior="$CIRC_PRIOR_DIR/${sample}_circexplorer_results/${sample}_CIRCRNA_circularRNA_known.txt"
    [ -f "$prior" ] || { echo "[ERROR] --use-circexplorer-prior given but $prior not found" >&2; exit 1; }
    extra=(--circ "$prior" --tool CIRCexplorer2)
    echo "[INFO] using CIRCexplorer2 prior: $prior"
  fi

  echo "[INFO] CIRIquant"
  ( cd "$OUT_DIR" && CIRIquant -t "$THREADS" \
      -1 "$work/${sample}_1.trimmed.fq" \
      -2 "$work/${sample}_2.trimmed.fq" \
      --config "$CIRIQUANT_CONFIG" \
      -o "$sample" -p "$sample" "${extra[@]}" )

  # the originals deleted the two large intermediate BAMs
  rm -f "$OUT_DIR/$sample/align/${sample}.sorted.bam" \
        "$OUT_DIR/$sample/circ/${sample}_denovo.sorted.bam"
  rm -rf "$work"
  echo "[INFO] $sample done -> $OUT_DIR/$sample"
done

echo "[INFO] stage 14 complete -> $OUT_DIR"
