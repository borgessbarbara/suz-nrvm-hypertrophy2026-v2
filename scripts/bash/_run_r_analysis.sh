#!/usr/bin/env bash
# Run the R analysis stages in dependency order.
#
# The numbering IS the dependency order - scripts/checks/check_io_graph.R
# enforces that no stage reads a file produced by a later one. If you add a
# stage, add it to metadata/stage_io.tsv as well or the check will fail.
#
# Usage:
#   _run_r_analysis.sh [--results DIR] [--dpi N] [--stages LIST] [--dry-run]
#                      [--continue-on-error] [--skip-circrna]
#
#   --stages accepts "1-9", "4,7,12" or "14-" (from 14 to the end).

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
R_DIR="$REPO_ROOT/scripts/R"

RESULTS="results"
DPI=600
STAGES=""
DRY=0
KEEPGOING=0
SKIP_CIRC=0

usage() { sed -n '2,14p' "$0"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --results) RESULTS="$2"; shift 2 ;;
    --dpi) DPI="$2"; shift 2 ;;
    --stages) STAGES="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    --continue-on-error) KEEPGOING=1; shift ;;
    --skip-circrna) SKIP_CIRC=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# circRNA stages, skipped by --skip-circrna
CIRC_STAGES=" 04 04b 13 20 22 "

in_selection() {
  local n="$1"
  [ -z "$STAGES" ] && return 0
  local IFS=','
  for part in $STAGES; do
    case "$part" in
      *-) [ "$((10#$n))" -ge "$((10#${part%-}))" ] && return 0 ;;
      -*) [ "$((10#$n))" -le "$((10#${part#-}))" ] && return 0 ;;
      *-*) [ "$((10#$n))" -ge "$((10#${part%%-*}))" ] &&
           [ "$((10#$n))" -le "$((10#${part##*-}))" ] && return 0 ;;
      *) [ "$((10#$n))" -eq "$((10#$part))" ] && return 0 ;;
    esac
  done
  return 1
}

mapfile -t SCRIPTS < <(find "$R_DIR" -maxdepth 1 -name '[0-9][0-9]*_*.R' | LC_ALL=C sort)
[ ${#SCRIPTS[@]} -gt 0 ] || { echo "[ERROR] no R stages found in $R_DIR" >&2; exit 1; }

echo "================================================================"
echo " R analysis stages"
echo "   results : $RESULTS"
echo "   dpi     : $DPI"
echo "   stages  : ${STAGES:-all}"
[ "$SKIP_CIRC" = 1 ] && echo "   circRNA : SKIPPED"
echo "================================================================"

ran=0; skipped=0; failed=0
FAILED_LIST=""

for s in "${SCRIPTS[@]}"; do
  base="$(basename "$s")"
  n="${base%%_*}"
  n_num="${n%%[a-z]*}"

  if ! in_selection "$n_num"; then skipped=$((skipped+1)); continue; fi
  if [ "$SKIP_CIRC" = 1 ] && [[ "$CIRC_STAGES" == *" $n "* ]]; then
    echo "-- R$n $base  (skipped: circRNA)"
    skipped=$((skipped+1)); continue
  fi

  echo
  echo "---------------------------------------------------------------"
  echo "-- R$n  $base"
  echo "---------------------------------------------------------------"

  args=(--in-dir "$RESULTS" --dpi "$DPI")
  # stages whose CLI has no --dpi
  case "$n" in
    01|04|04b) args=(--in-dir "$RESULTS") ;;
  esac

  if [ "$DRY" = 1 ]; then
    echo "   [dry-run] Rscript $s ${args[*]}"
    ran=$((ran+1)); continue
  fi

  if Rscript "$s" "${args[@]}"; then
    ran=$((ran+1))
  else
    failed=$((failed+1)); FAILED_LIST="$FAILED_LIST R$n"
    echo "[ERROR] R$n ($base) failed." >&2
    if [ "$KEEPGOING" != 1 ]; then
      echo "[ERROR] stopping. Pass --continue-on-error to run the rest anyway." >&2
      exit 1
    fi
  fi
done

echo
echo "================================================================"
echo " ran $ran, skipped $skipped, failed $failed"
[ -n "$FAILED_LIST" ] && echo " failed stages:$FAILED_LIST"
echo "================================================================"
[ "$failed" -eq 0 ] || exit 1
