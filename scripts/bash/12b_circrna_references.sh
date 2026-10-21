#!/usr/bin/env bash
# Stage 12b - build every reference the circRNA branch needs.
#
# Builds every circRNA-branch reference that no received script ever produced:
#   <star-index>/                    STAR genome index  (the historical STAR_rn8/)
#   <refflat>                        refFlat for CIRCexplorer2 annotate (stage 13)
#   <bwa-index>.{amb,ann,bwt,pac,sa} BWA index for CIRI2 (15) and CIRIquant (14)
#   <genome>.fai                     samtools index, required by CIRIquant
#
# All derive from the same GRCr8 primary assembly and GTF that stage 00
# produces, so the circRNA branch and the long-RNA branch share one reference.
#
# THE INDEX IS BUILT WITHOUT SPLICE JUNCTIONS BY DEFAULT.
# CIRCexplorer2 detects back-splice junctions from STAR's chimeric output, and
# its own documentation builds the index with no --sjdbGTFfile. The recorded
# run passed no GTF either. Pass --sjdb-gtf to add splice junctions anyway; it
# changes chimeric detection, so do not do it to match the long-RNA branch.
#
# Usage:  bash scripts/bash/12b_circrna_references.sh [options]
#         bash scripts/bash/12b_circrna_references.sh --help
# Env:    GENOME_FASTA, GTF_FILE, STAR_INDEX, REFFLAT_OUT, BWA_INDEX, THREADS,
#         ENSEMBL_RELEASE, STAR_LIMIT_RAM
#
# MEMORY: STAR genomeGenerate on a ~2.6 Gb genome needs roughly 32 GB of RAM.
# It will fail on a smaller machine; --limit-ram lowers STAR's budget but
# cannot remove the requirement.
#
# Next: 13_circrna_circexplorer.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

REL_FILE="$REPO_ROOT/metadata/reference_release.tsv"
RELEASE="${ENSEMBL_RELEASE:-}"
[ -z "$RELEASE" ] && [ -f "$REL_FILE" ] && \
  RELEASE="$(awk -F'\t' 'NR>1 && tolower($2) ~ /^(yes|true|1)$/ {print $1; exit}' "$REL_FILE")"
RELEASE="${RELEASE:-116}"

GENOME_FASTA="${GENOME_FASTA:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.primary_assembly.fa}"
GTF_FILE="${GTF_FILE:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.${RELEASE}.gtf}"
STAR_INDEX="${STAR_INDEX:-$REPO_ROOT/references/STAR_GRCr8_${RELEASE}}"
REFFLAT_OUT="${REFFLAT_OUT:-$REPO_ROOT/references/GRCr8.${RELEASE}.refflat.txt}"
BWA_INDEX="${BWA_INDEX:-$REPO_ROOT/references/bwa_GRCr8/GRCr8}"
THREADS="${THREADS:-8}"
STAR_LIMIT_RAM="${STAR_LIMIT_RAM:-}"
SJDB_GTF=0
SJDB_OVERHANG=149
FORCE=0
SKIP_STAR=0
SKIP_REFFLAT=0
SKIP_BWA=0

usage() {
  cat <<'EOF'
Usage: 12b_circrna_references.sh [options]

  --genome-fasta FILE   GRCr8 primary assembly   (env GENOME_FASTA)
  --gtf FILE            Ensembl GTF              (env GTF_FILE)
  --star-index DIR      STAR index directory     (env STAR_INDEX)
  --refflat FILE        refFlat output           (env REFFLAT_OUT)
  --bwa-index PREFIX    BWA index PREFIX, not a directory (env BWA_INDEX)
  --threads N           default 8                (env THREADS)
  --sjdb-gtf            build a splice-aware index (NOT the default; see header)
  --sjdb-overhang N     only with --sjdb-gtf, default 149 (read length - 1)
  --limit-ram BYTES     STAR --limitGenomeGenerateRAM
  --skip-star           do not build the STAR index
  --skip-refflat        do not build the refFlat
  --skip-bwa            do not build the BWA index or the .fai
  --force               rebuild even if the output already exists
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --genome-fasta)  GENOME_FASTA="$2"; shift 2 ;;
    --gtf)           GTF_FILE="$2"; shift 2 ;;
    --star-index)    STAR_INDEX="$2"; shift 2 ;;
    --refflat)       REFFLAT_OUT="$2"; shift 2 ;;
    --threads)       THREADS="$2"; shift 2 ;;
    --sjdb-gtf)      SJDB_GTF=1; shift ;;
    --sjdb-overhang) SJDB_OVERHANG="$2"; shift 2 ;;
    --limit-ram)     STAR_LIMIT_RAM="$2"; shift 2 ;;
    --bwa-index)     BWA_INDEX="$2"; shift 2 ;;
    --skip-star)     SKIP_STAR=1; shift ;;
    --skip-refflat)  SKIP_REFFLAT=1; shift ;;
    --skip-bwa)      SKIP_BWA=1; shift ;;
    --force)         FORCE=1; shift ;;
    --help|-h)       usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

source "$SCRIPT_DIR/_activate_env.sh" env_circrna

echo "[INFO] Stage 12b: STAR index and refFlat for the circRNA branch"
echo "[INFO] genome      : $GENOME_FASTA"
echo "[INFO] GTF         : $GTF_FILE"
echo "[INFO] STAR index  : $STAR_INDEX"
echo "[INFO] refFlat     : $REFFLAT_OUT"
echo "[INFO] BWA index   : $BWA_INDEX"
echo "[INFO] threads     : $THREADS"
echo ""

for f in "$GENOME_FASTA" "$GTF_FILE"; do
  if [ ! -f "$f" ]; then
    echo "[ERROR] reference not found: $f" >&2
    echo "[ERROR] run scripts/bash/00_download_references.sh first, or pass" >&2
    echo "[ERROR] --genome-fasta / --gtf." >&2
    exit 1
  fi
done

# --- STAR genome index -------------------------------------------------------
if [ "$SKIP_STAR" -eq 0 ]; then
  if [ -s "$STAR_INDEX/SA" ] && [ "$FORCE" -eq 0 ]; then
    echo "[INFO] STAR index already present, skipping: $STAR_INDEX"
    echo "[INFO] use --force to rebuild"
  else
    echo "[INFO] building the STAR index (this needs ~32 GB of RAM)..."
    mkdir -p "$STAR_INDEX"

    declare -a STAR_ARGS=(
      --runMode genomeGenerate
      --runThreadN "$THREADS"
      --genomeDir "$STAR_INDEX"
      --genomeFastaFiles "$GENOME_FASTA"
    )
    if [ "$SJDB_GTF" -eq 1 ]; then
      echo "[WARN] --sjdb-gtf given: building a SPLICE-AWARE index."
      echo "[WARN] CIRCexplorer2 expects an index built without splice"
      echo "[WARN] junctions; this changes chimeric detection."
      STAR_ARGS+=(--sjdbGTFfile "$GTF_FILE" --sjdbOverhang "$SJDB_OVERHANG")
    fi
    [ -n "$STAR_LIMIT_RAM" ] && STAR_ARGS+=(--limitGenomeGenerateRAM "$STAR_LIMIT_RAM")

    STAR "${STAR_ARGS[@]}"
    echo "[INFO] STAR index written to $STAR_INDEX"
  fi
  echo ""
fi

# --- refFlat annotation ------------------------------------------------------
# CIRCexplorer2 annotate -r expects refFlat:
#   geneName  name  chrom  strand  txStart txEnd cdsStart cdsEnd
#   exonCount exonStarts exonEnds
# gtfToGenePred -genePredExt emits genePredExt, whose columns are the same set
# in a different order with the gene name last; the awk below reorders them.
#
# Ensembl chromosome names carry no "chr" prefix, which is what the historical
# filename rn8.refflat.nochr.txt was recording. No renaming is needed as long
# as the GTF and the genome FASTA come from the same Ensembl release.
if [ "$SKIP_REFFLAT" -eq 0 ]; then
  if [ -s "$REFFLAT_OUT" ] && [ "$FORCE" -eq 0 ]; then
    echo "[INFO] refFlat already present, skipping: $REFFLAT_OUT"
    echo "[INFO] use --force to rebuild"
  else
    command -v gtfToGenePred >/dev/null 2>&1 || {
      echo "[ERROR] gtfToGenePred not found. It comes from the bioconda" >&2
      echo "[ERROR] package ucsc-gtftogenepred, pinned in" >&2
      echo "[ERROR] envs/environment_circrna.yml." >&2
      exit 1; }

    mkdir -p "$(dirname "$REFFLAT_OUT")"
    TMP_GP="$(mktemp)"
    trap 'rm -f "$TMP_GP"' EXIT

    echo "[INFO] converting the GTF to refFlat..."
    gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons \
      "$GTF_FILE" "$TMP_GP"

    awk 'BEGIN{OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' \
      "$TMP_GP" > "$REFFLAT_OUT"

    n_rows=$(wc -l < "$REFFLAT_OUT")
    n_chr=$(cut -f3 "$REFFLAT_OUT" | sort -u | wc -l)
    echo "[INFO] refFlat written: $REFFLAT_OUT"
    echo "[INFO]   $n_rows transcripts across $n_chr sequences"

    if [ "$n_rows" -eq 0 ]; then
      echo "[ERROR] the refFlat is empty; the GTF is probably not Ensembl-format" >&2
      exit 1
    fi
    if cut -f3 "$REFFLAT_OUT" | head -1 | grep -q '^chr'; then
      echo "[WARN] chromosome names start with 'chr'. The circRNA genome FASTA"
      echo "[WARN] must use the same convention or CIRCexplorer2 annotate finds"
      echo "[WARN] nothing."
    fi
  fi
fi

# --- BWA index and samtools faidx -------------------------------------------
# CIRI2 (stage 15) aligns with bwa mem, and CIRIquant (stage 14) needs both a
# BWA index prefix and a .fai beside the FASTA.
if [ "$SKIP_BWA" -eq 0 ]; then
  if [ -s "${BWA_INDEX}.bwt" ] && [ "$FORCE" -eq 0 ]; then
    echo "[INFO] BWA index already present, skipping: $BWA_INDEX"
    echo "[INFO] use --force to rebuild"
  else
    echo "[INFO] building the BWA index (slow; single-threaded by design)..."
    mkdir -p "$(dirname "$BWA_INDEX")"
    bwa index -p "$BWA_INDEX" "$GENOME_FASTA"
    echo "[INFO] BWA index written with prefix $BWA_INDEX"
  fi

  if [ -s "${GENOME_FASTA}.fai" ] && [ "$FORCE" -eq 0 ]; then
    echo "[INFO] ${GENOME_FASTA##*/}.fai already present"
  else
    echo "[INFO] indexing the FASTA with samtools faidx..."
    samtools faidx "$GENOME_FASTA"
    echo "[INFO] wrote ${GENOME_FASTA}.fai"
  fi
  echo ""
fi

echo ""
echo "[INFO] stage 12b complete"
echo "[INFO] Next: scripts/bash/13_circrna_circexplorer.sh \\"
echo "[INFO]         --star-index $STAR_INDEX --refflat $REFFLAT_OUT"
echo "[INFO] Stages 14 and 15 pick up the BWA index and .fai automatically."
