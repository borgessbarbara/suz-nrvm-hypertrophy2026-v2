#!/usr/bin/env bash
# Stage 15 - circRNA sequence reconstruction (BWA-MEM, CIRI2, CIRI-AS).
#
# Note the fastp settings differ from stages 13 and 14: this branch adds
# -b 140 -B 140 -l 140, hard-trimming to 140 bp for CIRI-AS.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_circrna

SAMPLE_SHEET="${SAMPLE_SHEET:-$REPO_ROOT/metadata/longRNA_samples.tsv}"
FASTQ_MANIFEST="${FASTQ_MANIFEST:-$REPO_ROOT/metadata/fastq_manifest.tsv}"
OUT_DIR="${OUT_DIR:-results/15_ciri_as}"
SCRATCH="${SCRATCH:-results/15_ciri_as/scratch}"
THREADS="${THREADS:-20}"
REL_FILE="$REPO_ROOT/metadata/reference_release.tsv"
RELEASE="${ENSEMBL_RELEASE:-}"
[ -z "$RELEASE" ] && [ -f "$REL_FILE" ] && \
  RELEASE="$(awk -F'\t' 'NR>1 && tolower($2) ~ /^(yes|true|1)$/ {print $1; exit}' "$REL_FILE")"
RELEASE="${RELEASE:-116}"

GENOME_FASTA="${CIRC_GENOME_FASTA:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.primary_assembly.fa}"
GTF="${CIRC_GTF:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.${RELEASE}.gtf}"
BWA_INDEX="${BWA_INDEX:-$REPO_ROOT/references/bwa_GRCr8/GRCr8}"
CIRI2_PL="${CIRI2_PL:-}"
CIRI_AS_PL="${CIRI_AS_PL:-}"
BWA_MIN_SCORE="${BWA_MIN_SCORE:-19}"   # the original's -T 19
ONLY=""
RUN_ALL=1

usage() {
  cat <<'EOF'
Usage: 15_circrna_ciri_as.sh --genome-fasta F --gtf F --ciri2 P --ciri-as P [options]

Everything has a default; no argument is required.

  --genome-fasta PATH   default: the stage 00 GRCr8 primary assembly
  --gtf PATH            default: the active-release Ensembl GTF
  --bwa-index PREFIX    default: references/bwa_GRCr8/GRCr8 (stage 12b)
  --ciri2 PATH          CIRI2.pl        (default: found on PATH from conda)
  --ciri-as PATH        CIRI_AS.pl      (default: found on PATH from conda)
  --samples PATH        default metadata/longRNA_samples.tsv
  --manifest PATH       default metadata/fastq_manifest.tsv
  --out-dir DIR         default results/15_ciri_as
  --threads N           default 20
  --bwa-min-score N     bwa mem -T (default 19, as in the original)
  --only SAMPLE[,...]   restrict to these samples. Default is every sample in
                        the sheet; pass SQ13_ln,SQ14_ln,SQ15_ln,SQ16_ln to
                        reproduce the original four-sample loop.
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
    --bwa-index)     BWA_INDEX="$2";     shift 2 ;;
    --only)          ONLY="$2";          shift 2 ;;
    --all)           ONLY="";            shift ;;
    --yes-widen)     shift ;;
    --help|-h)       usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# CIRI2 and CIRI-AS come from the bioconda packages ciri2 and ciri-as, at the
# exact versions the Methods state (2.0.6 and 1.2). Their wrappers are on PATH
# once env_circrna is active, so neither path has to be supplied; the names
# have varied across releases, so try the known spellings in turn.
find_tool() {
  local out="" c
  for c in "$@"; do
    out="$(command -v "$c" 2>/dev/null)" && { echo "$out"; return 0; }
  done
  return 1
}
if [ -z "$CIRI2_PL" ]; then
  CIRI2_PL="$(find_tool CIRI2.pl CIRI2 ciri2)" || {
    echo "[ERROR] CIRI2 not found on PATH. It is the bioconda package 'ciri2'," >&2
    echo "[ERROR] pinned in envs/environment_circrna.yml. Activate env_circrna," >&2
    echo "[ERROR] or pass --ciri2 /path/to/CIRI2.pl." >&2; exit 1; }
fi
if [ -z "$CIRI_AS_PL" ]; then
  CIRI_AS_PL="$(find_tool CIRI_AS.pl CIRI_AS_v1.2.pl CIRI_AS ciri-as)" || {
    echo "[ERROR] CIRI-AS not found on PATH. It is the bioconda package" >&2
    echo "[ERROR] 'ciri-as', pinned in envs/environment_circrna.yml. Activate" >&2
    echo "[ERROR] env_circrna, or pass --ciri-as /path/to/CIRI_AS.pl." >&2; exit 1; }
fi
echo "[INFO] CIRI2   : $CIRI2_PL"
echo "[INFO] CIRI-AS : $CIRI_AS_PL"

# A bioconda install may expose a wrapper rather than the raw .pl, and
# `perl <wrapper>` would fail. Run a .pl through perl, anything else directly.
run_perl_tool() {
  local tool="$1"; shift
  case "$tool" in
    *.pl) perl "$tool" "$@" ;;
    *)    "$tool" "$@" ;;
  esac
}

[ -s "${BWA_INDEX}.bwt" ] || {
  echo "[ERROR] BWA index not found at prefix $BWA_INDEX" >&2
  echo "[ERROR] Run scripts/bash/12b_circrna_references.sh, or --bwa-index." >&2
  exit 1; }

for f in "$GENOME_FASTA" "$GTF" "$SAMPLE_SHEET" "$FASTQ_MANIFEST"; do
  [ -f "$f" ] || { echo "[ERROR] not found: $f" >&2; exit 1; }
done

all_samples=$(awk -F'\t' 'NR>1 && $1!="" {print $1}' "$SAMPLE_SHEET")
if [ -z "$ONLY" ]; then
  samples="$all_samples"
  echo "[INFO] running every sample in the sheet."
  echo "[INFO] The original loop was 'for m in {13..16}', i.e. SQ13-SQ16 only."
  echo "[INFO] Widened by decision (B. Borges): restricting reconstruction to"
  echo "[INFO] four of sixteen libraries had no methodological justification."
  echo "[INFO] Use --only to reproduce the original four."
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

  # bwa mem takes the index PREFIX, not the FASTA. The original passed the
  # FASTA because its index had been built in place; stage 12b builds with
  # -p "$BWA_INDEX", so the prefix is what matches here.
  echo "[INFO] bwa mem -T $BWA_MIN_SCORE"
  bwa mem -t "$THREADS" -T "$BWA_MIN_SCORE" "$BWA_INDEX" \
      "$work/${sample}_1.trimmed.fq" "$work/${sample}_2.trimmed.fq" > "$sam"

  # CIRI2 is a prerequisite, not an alternative: CIRI-AS takes its .ciri call
  # set through -C and reconstructs structure only for circles CIRI2 found.
  echo "[INFO] CIRI2"
  run_perl_tool "$CIRI2_PL" -I "$sam" -O "$ciri" -F "$GENOME_FASTA" -A "$GTF" -T "$THREADS"

  echo "[INFO] CIRI-AS"
  run_perl_tool "$CIRI_AS_PL" -S "$sam" -C "$ciri" -F "$GENOME_FASTA" -A "$GTF" \
       -T "$THREADS" -O "$OUT_DIR/AS_results/${sample}.ASout" -D yes

  rm -rf "$work"
  echo "[INFO] $sample done"
done

echo "[INFO] stage 15 complete -> $OUT_DIR"
echo "[INFO] Next: Rscript scripts/R/04b_circRNA_de_temporal.R"
