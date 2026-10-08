#!/usr/bin/env bash
# Stage 16: circRNA-miRNA duplex prediction with IntaRNA.
#
# The four circRNA-miRNA commands, exactly as recorded in
# provenance/intarna/IntaRNA_commands_as_used.sh and
# provenance/intarna/IntaRNA_manifest.tsv:
#
#   IntaRNA --mode=M -t circRNA_UP_6h_all_records.fa    -q miRNAs_DOWN_6hr.fa  --outMode=C --out Interact_circUP_miDOWN_6h_all_modeM.csv
#   IntaRNA --mode=M -t circRNA_DOWN_6h_all_records.fa  -q miRNAs_UP_6hr.fa    --outMode=C --out Interact_circDOWN_miUP_6h_all_modeM.csv
#   IntaRNA --mode=M -t circRNA_UP_24h_all_records.fa   -q miRNAs_DOWN_24hr.fa --outMode=C --out Interact_circUP_miDOWN_24h_all_modeM.csv
#   IntaRNA --mode=M -t circRNA_DOWN_24h_all_records.fa -q miRNAs_UP_24hr.fa   --outMode=C --out Interact_circDOWN_miUP_24h_all_modeM.csv
#
# The pairing is deliberately anti-correlated: UP circRNAs are tested against
# DOWN miRNAs and vice versa, because a sponging ceRNA relationship predicts
# opposite directions. Stage R20 then retains duplexes with MFE <= -20 kcal/mol.
#
# The four target FASTAs come from R18b (the reconstructed circRNA merge/export
# stage) and the four query FASTAs from R13. No script in any package received
# produced the circRNA_*_all_records.fa files; see docs/OPEN_ITEMS.md.
#
# The output names MUST NOT change: R20 globs Interact_*.csv and parses the
# state and direction out of these exact filenames.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_intarna

CIRC_FASTA_DIR="${CIRC_FASTA_DIR:-results/R04b_circRNA_de_temporal}"
MIRNA_FASTA_DIR="${MIRNA_FASTA_DIR:-results/13_mirna_mrna_integration/miRNA_FASTA_by_condition}"
OUT_DIR="${OUT_DIR:-results/16_intarna_circrna}"
THREADS="${THREADS:-1}"
MANIFEST="${MANIFEST:-$REPO_ROOT/provenance/intarna/IntaRNA_manifest.tsv}"

usage() {
  cat <<'EOF'
Usage: 16_circrna_intarna.sh [options]

  --circ-fasta-dir DIR   circRNA target FASTAs  (default results/R04b_circRNA_de_temporal)
  --mirna-fasta-dir DIR  miRNA query FASTAs     (default results/13_.../miRNA_FASTA_by_condition)
  --out-dir DIR          default results/16_intarna_circrna
  --threads N            IntaRNA --threads (default 1; the originals did not set it)

Runs the four circRNA-miRNA IntaRNA commands recorded in
provenance/intarna/IntaRNA_manifest.tsv.
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
come from scripts/R/13_mirna_mRNA_integration.R.

See metadata/circRNA_input_manifest.tsv and docs/OPEN_ITEMS.md.
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
