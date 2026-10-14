#!/usr/bin/env bash
# Create the pipeline's conda environments, then verify and lock them.
#
# Creates env_qc_align, env_mirge, env_intarna, env_circrna and env_r from
# envs/*.yml, checks that the binaries each stage calls are actually on PATH
# inside the environment, and writes envs/*.lock.yml recording what resolved.
#
# The environment files pin no versions, so this installs current releases.
# The lock files are the version record; update
# metadata/software_tool_inventory.tsv from them.
#
# Usage:  bash scripts/commands/create_envs.sh [options]
#
#   --only NAME[,NAME]  create only these (qc_align mirge intarna circrna r)
#   --solver SOLVER     conda solver: libmamba (default) or classic
#   --dry-run           solve without installing, then stop
#   --no-lock           skip writing the .lock.yml files
#   --force             remove an existing environment first
#
# Requires conda. mamba is used automatically when present.
#
# CLEAR, CIRI2 and CIRI-AS are not conda packages. After env_circrna exists
# this script prints the CLEAR install command; CIRI2 and CIRI-AS must be
# supplied to stage 15 with --ciri2-path and --ciri-as-path.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
ENV_DIR="$REPO_ROOT/envs"

ONLY=""
SOLVER="libmamba"
DRY_RUN=0
WRITE_LOCK=1
FORCE=0

usage() { sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
  case "$1" in
    --only)    ONLY="$2"; shift 2 ;;
    --solver)  SOLVER="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --no-lock) WRITE_LOCK=0; shift ;;
    --force)   FORCE=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

command -v conda >/dev/null 2>&1 || {
  echo "[ERROR] conda is not on PATH." >&2; exit 1; }

CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

CREATE="conda"
command -v mamba >/dev/null 2>&1 && CREATE="mamba"

echo "== creating pipeline environments =="
echo "   conda      : $(conda --version)"
echo "   create via : $CREATE"
echo "   solver     : $SOLVER"
echo "   channels   : conda-forge, bioconda, nodefaults (set per file)"
echo ""
echo "Strict channel priority is required by bioconda. Setting it for this run."
conda config --set channel_priority strict

# name -> the binaries that must be callable after activation
binaries_for() {
  case "$1" in
    qc_align) echo "fastqc fastp multiqc hisat2 hisat2-build samtools stringtie gffcompare gffread FEELnc_filter.pl FEELnc_codpot.pl FEELnc_classifier.pl salmon gawk pigz" ;;
    mirge)    echo "miRge3.0 bowtie samtools Rscript" ;;
    intarna)  echo "IntaRNA" ;;
    circrna)  echo "hisat2 STAR bwa samtools fastp CIRCexplorer2 CIRIquant stringtie bedtools perl gawk" ;;
    r)        echo "Rscript R" ;;
  esac
}

ALL="qc_align mirge intarna circrna r"
TARGETS="$ALL"
if [ -n "$ONLY" ]; then TARGETS="$(echo "$ONLY" | tr ',' ' ')"; fi

status=0
declare -a MADE=()

for short in $TARGETS; do
  case " $ALL " in *" $short "*) ;; *)
    echo "[ERROR] unknown environment: $short (expected one of: $ALL)" >&2
    status=1; continue ;;
  esac

  env_name="env_$short"
  yml="$ENV_DIR/environment_${short}.yml"

  echo ""
  echo "---------------------------------------------------------------"
  echo "-- $env_name   <-  envs/environment_${short}.yml"
  echo "---------------------------------------------------------------"

  if [ ! -f "$yml" ]; then
    echo "[ERROR] missing $yml" >&2; status=1; continue
  fi

  # the name inside the file must match what _activate_env.sh activates
  declared="$(awk '/^name:/{print $2; exit}' "$yml")"
  if [ "$declared" != "$env_name" ]; then
    echo "[ERROR] $yml declares 'name: $declared' but the stages activate" >&2
    echo "[ERROR] '$env_name'. Fix the name: field or nothing will activate." >&2
    status=1; continue
  fi

  if conda env list | awk '{print $1}' | grep -qx "$env_name"; then
    if [ "$FORCE" -eq 1 ]; then
      echo "[INFO] --force: removing the existing $env_name"
      conda env remove -n "$env_name" -y >/dev/null || true
    else
      echo "[INFO] $env_name already exists; skipping creation (use --force to rebuild)"
      MADE+=("$short"); continue
    fi
  fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry-run] $CREATE env create -f $yml --solver=$SOLVER --dry-run"
    $CREATE env create -f "$yml" --solver="$SOLVER" --dry-run || {
      echo "[FAIL] $env_name does not solve" >&2; status=1; }
    continue
  fi

  if $CREATE env create -f "$yml" --solver="$SOLVER"; then
    echo "[ok] created $env_name"
    MADE+=("$short")
  else
    echo "[FAIL] could not create $env_name" >&2
    status=1
    continue
  fi
done

[ "$DRY_RUN" -eq 1 ] && { echo ""; echo "dry run complete"; exit $status; }

# ---------------------------------------------------------------------------
# verify the binaries each stage calls, and lock
# ---------------------------------------------------------------------------
for short in ${MADE[@]+"${MADE[@]}"}; do
  env_name="env_$short"
  echo ""
  echo "-- verifying $env_name"
  conda activate "$env_name" 2>/dev/null || {
    echo "  [FAIL] cannot activate $env_name" >&2; status=1; continue; }

  missing=0
  for b in $(binaries_for "$short"); do
    if command -v "$b" >/dev/null 2>&1; then
      printf "  ok      %s\n" "$b"
    else
      printf "  MISSING %s\n" "$b"
      missing=$((missing + 1))
    fi
  done
  [ "$missing" -gt 0 ] && { echo "  [FAIL] $missing binary/binaries missing from $env_name" >&2; status=1; }

  # gawk must really be GNU awk: the metadata chain needs 3-argument match()
  if [ "$short" = "qc_align" ] || [ "$short" = "circrna" ]; then
    if echo | gawk 'BEGIN{ if (match("ab","(a)(b)",m)==0) exit 1 }' 2>/dev/null; then
      echo "  ok      gawk supports 3-argument match()"
    else
      echo "  [FAIL] gawk does not support 3-argument match(); stage 09e would be wrong" >&2
      status=1
    fi
  fi

  if [ "$WRITE_LOCK" -eq 1 ]; then
    lock="$ENV_DIR/environment_${short}.lock.yml"
    conda env export --no-builds -n "$env_name" > "$lock" && \
      echo "  [lock] $(basename "$lock")"
  fi

  conda deactivate
done

# ---------------------------------------------------------------------------
echo ""
echo "================================================================"
case " ${MADE[*]-} " in *" circrna "*)
  cat <<'EOF'
CLEAR is not a conda package and is still missing. Install it into
env_circrna now:

    conda activate env_circrna
    git clone https://github.com/YangLab/CLEAR
    cd CLEAR && python ./setup.py install

CIRI2 (CIRI2.pl) and CIRI-AS (CIRI_AS_v1.2.pl) are standalone Perl scripts.
Pass them to stage 15 with --ciri2-path and --ciri-as-path.
See docs/OPEN_ISSUES.md section 7.4 for the other circRNA resources needed.
EOF
;; esac

if [ "$status" -eq 0 ]; then
  echo "PASS - every requested environment created and verified."
  echo "Next: update metadata/software_tool_inventory.tsv from the .lock.yml files."
else
  echo "FAILED - see the errors above."
fi
echo "================================================================"
exit $status
