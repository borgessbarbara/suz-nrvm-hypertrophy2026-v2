#!/usr/bin/env bash
# Stage 16 - circRNA-miRNA duplex prediction (IntaRNA).
#
# Four runs, pairing UP circRNAs against DOWN miRNAs and vice versa, because a
# sponging relationship predicts opposite directions. All use
# --mode=M --outMode=C. R21 later retains duplexes with MFE <= -20 kcal/mol.
#
# The output names must not change: R21 globs Interact_*.csv and parses the
# state and direction out of them.
#
# Usage:  bash scripts/bash/16_circrna_intarna.sh [--circ-fasta-dir DIR]
#                [--mirna-fasta-dir DIR] [--out-dir DIR] [--threads N]
#
# Targets come from R04b, queries from R14.
#
# Next: Rscript scripts/R/20_circRNA_volcanoes.R

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_intarna

CIRC_FASTA_DIR="${CIRC_FASTA_DIR:-results/R04b_circRNA_de_temporal}"
MIRNA_FASTA_DIR="${MIRNA_FASTA_DIR:-results/13_mirna_mrna_integration/miRNA_FASTA_by_condition}"
OUT_DIR="${OUT_DIR:-results/16_intarna_circrna}"
THREADS="${THREADS:-1}"
MANIFEST="${MANIFEST:-$REPO_ROOT/archive/final_package_v3_original/03_INTARNA_COMMANDS/IntaRNA_manifest.tsv}"

usage() {
  cat <<'EOF'
Usage: 16_circrna_intarna.sh [options]

  --circ-fasta-dir DIR   circRNA target FASTAs  (default results/R04b_circRNA_de_temporal)
  --mirna-fasta-dir DIR  miRNA query FASTAs     (default results/13_.../miRNA_FASTA_by_condition)
  --out-dir DIR          default results/16_intarna_circrna
  --threads N            IntaRNA --threads (default 1; the originals did not set it)

Runs the four circRNA-miRNA IntaRNA commands recorded in
archive/final_package_v3_original/03_INTARNA_COMMANDS/IntaRNA_manifest.tsv.
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --circ-fasta-dir)  CIRC_FASTA_DIR="$2";  shift 2 ;;
    --mirna-fasta-dir) MIRNA_FASTA_DIR="$2"; shift 2 ;;
    --out-dir)         OUT_DIR="$2";         shift 2 ;;
    --threads)         THREADS="$2";         shift 2 ;;
    --help|-h)         usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

mkdir -p "$OUT_DIR"

# state  circ_direction  mirna_direction
runs=(
  "6h   UP   DOWN"
  "6h   DOWN UP"
  "24h  UP   DOWN"
  "24h  DOWN UP"
)

echo "[INFO] IntaRNA stage 16 (circRNA-miRNA)"
echo "[INFO]   circRNA FASTAs : $CIRC_FASTA_DIR"
echo "[INFO]   miRNA FASTAs   : $MIRNA_FASTA_DIR"
echo "[INFO]   output         : $OUT_DIR"

missing=0
for spec in "${runs[@]}"; do
  read -r state cdir mdir <<<"$spec"
  hr="${state%h}"
  target="$CIRC_FASTA_DIR/circRNA_${cdir}_${state}_all_records.fa"
  query="$MIRNA_FASTA_DIR/miRNAs_${mdir}_${hr}hr.fa"
  [ -f "$target" ] || { echo "[ERROR] missing target FASTA: $target" >&2; missing=1; }
  [ -f "$query" ]  || { echo "[ERROR] missing query FASTA : $query"  >&2; missing=1; }
done
if [ "$missing" -ne 0 ]; then
  cat >&2 <<'EOF'

[ERROR] One or more IntaRNA inputs are missing.

The circRNA_{UP,DOWN}_{6,24}h_all_records.fa targets are produced by
scripts/R/04b_circRNA_de_temporal.R, which is itself a RECONSTRUCTION - no
script in any package received generates these four files. The miRNA queries
come from scripts/R/14_mirna_mRNA_integration.R.

See metadata/external_inputs.tsv and docs/OPEN_ISSUES.md.
EOF
  exit 1
fi

for spec in "${runs[@]}"; do
  read -r state cdir mdir <<<"$spec"
  hr="${state%h}"
  target="$CIRC_FASTA_DIR/circRNA_${cdir}_${state}_all_records.fa"
  query="$MIRNA_FASTA_DIR/miRNAs_${mdir}_${hr}hr.fa"
  mi_tag=$([ "$mdir" = "UP" ] && echo "miUP" || echo "miDOWN")
  out="$OUT_DIR/Interact_circ${cdir}_${mi_tag}_${state}_all_modeM.csv"

  echo "[INFO] IntaRNA circ${cdir} x ${mi_tag} @ ${state}"
  echo "[INFO]   -t $(basename "$target")  -q $(basename "$query")"
  IntaRNA --mode=M -t "$target" -q "$query" --outMode=C --out "$out" \
          ${THREADS:+--threads "$THREADS"}

  n=$(( $(wc -l < "$out") - 1 ))
  echo "[INFO]   -> $(basename "$out")  ($n duplexes)"
done

echo "[INFO] stage 16 complete -> $OUT_DIR"
echo "[INFO] Downstream: R20 retains duplexes with MFE <= -20 kcal/mol."
