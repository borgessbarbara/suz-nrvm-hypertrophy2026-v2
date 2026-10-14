#!/usr/bin/env bash
# =============================================================================
# suz_pipeline.bash — top-level orchestrator.
#
# NE vs Control in neonatal rat ventricular myocytes, 6 h and 24 h.
# Long RNA (mRNA + novel lncRNA), small RNA (miRNA), circRNA, and the ceRNA
# networks that tie them together.
#
# Reference lock: Rattus_norvegicus.GRCr8.dna.primary_assembly.fa
#                 Rattus_norvegicus.GRCr8.<REL>.gtf  (release from metadata/reference_release.tsv)
#                 biomaRt archive follows the same release. See docs/CHANGES_FROM_ORIGINALS.md.
#
# Stage groups:
#   00-12   long RNA and small RNA     (env qc_align, mirge, intarna)
#   09b-09e transcript metadata chain  (env qc_align)   <- 09b is RECONSTRUCTED
#   13-16   circRNA module             (env circrna, intarna)
#   R01-R26 analysis and figures       (env r)
#
# See README.md for the full stage table and docs/OPEN_ISSUES.md for what still
# needs confirmation before the repository is frozen.
# =============================================================================

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASH_DIR="$REPO_ROOT/scripts/bash"

LONGRNA_SAMPLES="$REPO_ROOT/metadata/longRNA_samples.tsv"
SMALLRNA_SAMPLES="$REPO_ROOT/metadata/smallRNA_samples.tsv"
FASTQ_MANIFEST="$REPO_ROOT/metadata/fastq_manifest.tsv"
GENOME_FASTA=""
GTF=""
THREADS=16
OUT="results"
DPI=600
STAGES=""
RUN_R=0
SKIP_DOWNLOAD=0
SKIP_CIRC=0
DRY=0

usage() {
  cat <<'EOF'
Usage: ./suz_pipeline.bash [options]

Inputs
  --longrna-samples FILE    default metadata/longRNA_samples.tsv
  --smallrna-samples FILE   default metadata/smallRNA_samples.tsv
  --fastq-manifest FILE     default metadata/fastq_manifest.tsv
  --genome-fasta FILE       default: built by stage 00 (GRCr8 primary assembly)
  --gtf FILE                default: downloaded by stage 00 (GRCr8.115)
  --skip-download-refs      do not run stage 00

Execution
  --threads N               default 16
  --out DIR                 default results
  --dpi N                   figure resolution, default 600
  --stages LIST             e.g. "1-10", "4,7", "13-"  (bash stages only)
  --run-r-analysis          also run the R stages after the bash stages
  --skip-circrna            omit the circRNA module (bash 13-16 and R04/13/20/22)
  --dry-run                 print what would run, run nothing

Notes
  * The circRNA module needs references and tools that are NOT in this
    repository (rn8.refflat, rn8.yml, CIRI2.pl, CIRI_AS.pl, a STAR index).
    See metadata/external_inputs.tsv. Those stages take their own
    arguments; run them directly rather than through this wrapper.
  * Stage 09b is a RECONSTRUCTION of a step whose original command was never
    recovered. Verify it with scripts/checks/check_feelnc_split_counts.sh.
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --longrna-samples) LONGRNA_SAMPLES="$2"; shift 2 ;;
    --smallrna-samples) SMALLRNA_SAMPLES="$2"; shift 2 ;;
    --fastq-manifest) FASTQ_MANIFEST="$2"; shift 2 ;;
    --genome-fasta) GENOME_FASTA="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --skip-download-refs) SKIP_DOWNLOAD=1; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --dpi) DPI="$2"; shift 2 ;;
    --stages) STAGES="$2"; shift 2 ;;
    --run-r-analysis) RUN_R=1; shift ;;
    --skip-circrna) SKIP_CIRC=1; shift ;;
    --dry-run) DRY=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

GENOME_FASTA="${GENOME_FASTA:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.primary_assembly.fa}"
# Ensembl release from metadata/reference_release.tsv (see stage 00).
REL_FILE="$REPO_ROOT/metadata/reference_release.tsv"
RELEASE="${ENSEMBL_RELEASE:-}"
[ -z "$RELEASE" ] && [ -f "$REL_FILE" ] && \
  RELEASE="$(awk -F'\t' 'NR>1 && tolower($2) ~ /^(yes|true|1)$/ {print $1; exit}' "$REL_FILE")"
RELEASE="${RELEASE:-115}"
GTF="${GTF:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.${RELEASE}.gtf}"

for f in "$LONGRNA_SAMPLES" "$SMALLRNA_SAMPLES" "$FASTQ_MANIFEST"; do
  [ -f "$f" ] || { echo "[ERROR] not found: $f" >&2; exit 1; }
done
mkdir -p "$OUT" "$OUT/logs"

export SAMPLE_SHEET="$LONGRNA_SAMPLES"
export SMALLRNA_SAMPLE_SHEET="$SMALLRNA_SAMPLES"
export FASTQ_MANIFEST GENOME_FASTA GTF_FILE="$GTF" THREADS
export RESULTS_DIR="$OUT"

# stage number -> script
declare -a ORDER=(00 01 02 03 04 05 06 07 08 09 09b 09c 09d 09e 10 11 12 12b 13 14 15 16)
declare -A SCRIPT=(
  [00]=00_download_references.sh      [01]=01_qc_raw.sh
  [02]=02_trim_fastp.sh               [03]=03_qc_trimmed.sh
  [04]=04_align_hisat2.sh             [05]=05_sort_index_bam.sh
  [06]=06_stringtie_assemble_merge.sh [07]=07_gffcompare.sh
  [08]=08_gffread_extract_fasta.sh    [09]=09_feelnc_classify.sh
  [09b]=09b_feelnc_split_exon_gtfs.sh [09c]=09c_extract_transcript_ids.sh
  [09d]=09d_subset_merged_gtfs.sh     [09e]=09e_filter_out_lncRNAs.sh
  [10]=10_salmon_index_quant.sh       [11]=11_mirge3_run.sh
  [12]=12_intarna_run.sh              [12b]=12b_circrna_star_index.sh
  [13]=13_circrna_circexplorer.sh
  [14]=14_circrna_ciriquant.sh        [15]=15_circrna_ciri_as.sh
  [16]=16_circrna_intarna.sh
)
CIRC=" 12b 13 14 15 16 "
NEEDS_ARGS=" 13 14 15 "   # require resources this repo cannot supply

num_of() { echo "${1%%b*}" | sed 's/^0*//;s/^$/0/'; }

in_selection() {
  local n; n="$(num_of "$1")"
  [ -z "$STAGES" ] && return 0
  local IFS=','
  for part in $STAGES; do
    case "$part" in
      *-) [ "$n" -ge "$(num_of "${part%-}")" ] && return 0 ;;
      -*) [ "$n" -le "$(num_of "${part#-}")" ] && return 0 ;;
      *-*) [ "$n" -ge "$(num_of "${part%%-*}")" ] && [ "$n" -le "$(num_of "${part##*-}")" ] && return 0 ;;
      *) [ "$n" -eq "$(num_of "$part")" ] && return 0 ;;
    esac
  done
  return 1
}

echo "================================================================"
echo " SUZ NE-NRVM pipeline"
echo "   genome  : $GENOME_FASTA"
echo "   gtf     : $GTF"
echo "   threads : $THREADS"
echo "   out     : $OUT"
echo "   stages  : ${STAGES:-all}"
[ "$SKIP_CIRC" = 1 ] && echo "   circRNA : SKIPPED"
echo "================================================================"

for s in "${ORDER[@]}"; do
  scr="${SCRIPT[$s]}"
  in_selection "$s" || continue
  [ "$s" = "00" ] && [ "$SKIP_DOWNLOAD" = 1 ] && { echo "-- stage 00 skipped"; continue; }
  if [ "$SKIP_CIRC" = 1 ] && [[ "$CIRC" == *" $s "* ]]; then
    echo "-- stage $s $scr  (skipped: circRNA)"; continue
  fi
  if [[ "$NEEDS_ARGS" == *" $s "* ]]; then
    cat <<EOF
-- stage $s $scr  (NOT run by the wrapper)
   This stage needs resources that are not in the repository (refFlat,
   rn8.yml, CIRI2.pl, CIRI_AS.pl, STAR index). Run it directly:
     bash scripts/bash/$scr --help
   See metadata/external_inputs.tsv and docs/OPEN_ISSUES.md.
EOF
    continue
  fi

  echo
  echo "---------------------------------------------------------------"
  echo "-- stage $s : $scr"
  echo "---------------------------------------------------------------"
  if [ "$DRY" = 1 ]; then echo "   [dry-run] bash $BASH_DIR/$scr"; continue; fi
  if ! bash "$BASH_DIR/$scr" 2>&1 | tee "$OUT/logs/stage_${s}.log"; then
    echo "[ERROR] stage $s failed. See $OUT/logs/stage_${s}.log" >&2
    exit 1
  fi
done

if [ "$RUN_R" = 1 ]; then
  echo
  echo "================================================================"
  echo " R analysis"
  echo "================================================================"
  rargs=(--results "$OUT" --dpi "$DPI")
  [ "$SKIP_CIRC" = 1 ] && rargs+=(--skip-circrna)
  [ "$DRY" = 1 ] && rargs+=(--dry-run)
  bash "$BASH_DIR/_run_r_analysis.sh" "${rargs[@]}" || exit 1
fi

echo
echo "[INFO] pipeline finished. Outputs under $OUT"
echo "[INFO] Before treating results as final, read docs/OPEN_ISSUES.md."
