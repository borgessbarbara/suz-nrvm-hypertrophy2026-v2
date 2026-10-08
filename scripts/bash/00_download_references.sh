#!/usr/bin/env bash
# Stage 00 - download the Ensembl rat reference.
#
# Builds the GRCr8 PRIMARY ASSEMBLY (not toplevel): prefers a single
# dna.primary_assembly.fa.gz, otherwise downloads the per-chromosome parts and
# concatenates them. Also downloads the matching GTF.
#
# The release comes from metadata/reference_release.tsv, the same file
# scripts/R/lib/ensembl.R reads. Override with ENSEMBL_RELEASE=114.
#
# Usage:  bash scripts/bash/00_download_references.sh [GENOME_FASTA] [GTF_FILE]
# Env:    ENSEMBL_RELEASE, PRIMARY_PARTS
#
# Next: 01_qc_raw.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# The Ensembl release comes from metadata/reference_release.tsv - the same file
# scripts/R/lib/ensembl.R reads - so the GTF and the biomaRt archive can never
# drift apart. Override for one run with ENSEMBL_RELEASE=114.
REL_FILE="$REPO_ROOT/metadata/reference_release.tsv"
if [ -n "${ENSEMBL_RELEASE:-}" ]; then
  RELEASE="$ENSEMBL_RELEASE"
elif [ -f "$REL_FILE" ]; then
  RELEASE="$(awk -F'\t' 'NR>1 && tolower($2) ~ /^(yes|true|1)$/ {print $1; exit}' "$REL_FILE")"
fi
RELEASE="${RELEASE:-115}"
echo "[INFO] Ensembl release: $RELEASE  (metadata/reference_release.tsv)"
GENOME_FASTA="${1:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.primary_assembly.fa}"
GTF_FILE="${2:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.${RELEASE}.gtf}"

DNA_BASE="https://ftp.ensembl.org/pub/release-${RELEASE}/fasta/rattus_norvegicus/dna"
ENSEMBL_FASTA_URL="$DNA_BASE/Rattus_norvegicus.GRCr8.dna.primary_assembly.fa.gz"
ENSEMBL_GTF_URL="https://ftp.ensembl.org/pub/release-${RELEASE}/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.${RELEASE}.gtf.gz"

# Concatenation order when assembling from per-chromosome parts: autosomes in
# numeric order, then X, Y, MT, then any unplaced scaffolds. Primary assembly
# EXCLUDES nonchromosomal by explicit decision (B. Borges, 8 Oct 2026): the
# analysis is chromosome-scoped. Add "nonchromosomal" to PRIMARY_PARTS to include it.
PRIMARY_PARTS="${PRIMARY_PARTS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X Y MT}"

url_exists() {
  curl -sS -o /dev/null -I -w '%{http_code}' --max-time 60 "$1" 2>/dev/null | grep -qE '^(200|226)$'
}

echo "[INFO] Reference download script starting..."
echo "[INFO] Target genome FASTA: $GENOME_FASTA"
echo "[INFO] Target GTF: $GTF_FILE"

# Function to download with curl or wget fallback
download_file() {
  local url="$1"
  local output="$2"

  if command -v curl &>/dev/null; then
    echo "[INFO] Downloading with curl: $url"
    curl -L -o "$output" "$url"
  elif command -v wget &>/dev/null; then
    echo "[INFO] Downloading with wget: $url"
    wget -O "$output" "$url"
  else
    echo "[ERROR] Neither curl nor wget found. Install one and try again." >&2
    return 1
  fi
}

# Download genome FASTA if not present
if [ ! -f "$GENOME_FASTA" ]; then
  mkdir -p "$(dirname "$GENOME_FASTA")"

  if url_exists "$ENSEMBL_FASTA_URL"; then
    echo "[INFO] Concatenated primary assembly is published; downloading it."
    TEMP_FASTA="$GENOME_FASTA.gz"
    download_file "$ENSEMBL_FASTA_URL" "$TEMP_FASTA"
    echo "[INFO] Decompressing..."
    gunzip -v "$TEMP_FASTA"
  else
    echo "[INFO] No single primary_assembly.fa.gz; assembling from parts."
    PARTS_DIR="$(dirname "$GENOME_FASTA")/primary_assembly_parts"
    mkdir -p "$PARTS_DIR"
    GOT=0
    : > "$GENOME_FASTA.tmp"
    for p in $PRIMARY_PARTS; do
      u="$DNA_BASE/Rattus_norvegicus.GRCr8.dna.primary_assembly.${p}.fa.gz"
      if url_exists "$u"; then
        f="$PARTS_DIR/primary_assembly.${p}.fa.gz"
        [ -f "$f" ] || download_file "$u" "$f"
        echo "[INFO]   + ${p}"
        gunzip -c "$f" >> "$GENOME_FASTA.tmp"
        GOT=$((GOT+1))
      fi
    done

    if [ "$GOT" -eq 0 ]; then
      rm -f "$GENOME_FASTA.tmp"
      cat >&2 <<EOF
[ERROR] No primary assembly found for release $RELEASE at
[ERROR]   $DNA_BASE
[ERROR]
[ERROR] Per the Ensembl README, that would mean GRCr8 has no haplotype/patch
[ERROR] regions and toplevel is equivalent. Verify before falling back:
[ERROR]   GENOME_FASTA=<path to toplevel> bash scripts/bash/00_download_references.sh
[ERROR] and record the decision in docs/OPEN_ISSUES.md.
EOF
      exit 1
    fi

    mv "$GENOME_FASTA.tmp" "$GENOME_FASTA"
    echo "[INFO] Concatenated $GOT primary-assembly parts."
    echo "[INFO] Order: $PRIMARY_PARTS"
    echo "[INFO] Parts kept under $PARTS_DIR (delete to reclaim space)."
  fi

  echo "[INFO] Sequences in genome FASTA: $(grep -c '^>' "$GENOME_FASTA")"
  echo "[INFO] Genome FASTA ready: $GENOME_FASTA"
else
  echo "[INFO] Genome FASTA already exists, skipping download."
fi

# Download GTF if not present
if [ ! -f "$GTF_FILE" ]; then
  echo "[INFO] GTF not found, downloading..."
  TEMP_GTF="$GTF_FILE.gz"
  download_file "$ENSEMBL_GTF_URL" "$TEMP_GTF"
  echo "[INFO] Decompressing GTF..."
  gunzip -v "$TEMP_GTF"
  echo "[INFO] GTF ready: $GTF_FILE"
else
  echo "[INFO] GTF already exists, skipping download."
fi

echo "[INFO] Reference download complete."
