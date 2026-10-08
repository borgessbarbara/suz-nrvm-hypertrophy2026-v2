#!/usr/bin/env bash
# Fail if any Ensembl connection in the runnable pipeline is unpinned.
#
# The whole point of lib/ensembl.R is that there is exactly ONE place where a
# biomaRt connection is made, and it carries an explicit archive host. The
# historical scripts did not: they used useEnsembl(..., version = 114), whose
# host is resolved at runtime against www.ensembl.org, or bare
# useMart("ensembl"), which is not pinned at all.
#
# provenance/ and archive/ are exempt by design - they are preserved verbatim.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1

# the expected host is derived from metadata/reference_release.tsv, not fixed here
status=0

echo "== biomaRt pinning check =="

# 1. no connection call outside lib/ensembl.R
offenders=$(grep -rnE '\b(useMart|useEnsembl|useDataset)\s*\(' scripts/R \
            --include='*.R' 2>/dev/null \
            | grep -v '^scripts/R/lib/ensembl.R:' \
            | grep -vE ':[0-9]+:\s*#')
if [ -n "$offenders" ]; then
  echo "FAIL: biomaRt connection opened outside scripts/R/lib/ensembl.R:"
  echo "$offenders" | sed 's/^/  /'
  echo "  -> use connect_ensembl() instead."
  status=1
else
  echo "  ok: all connections go through lib/ensembl.R"
fi

# 2. the active release resolves to a known archive, and the GTF follows it
if [ ! -f metadata/reference_release.tsv ]; then
  echo "FAIL: metadata/reference_release.tsv is missing; the release is unset."
  status=1
else
  ACTIVE=$(awk -F'\t' 'NR>1 && tolower($2) ~ /^(yes|true|1)$/ {print $1; exit}' \
           metadata/reference_release.tsv)
  N_ACTIVE=$(awk -F'\t' 'NR>1 && tolower($2) ~ /^(yes|true|1)$/' \
             metadata/reference_release.tsv | wc -l)
  if [ "$N_ACTIVE" -ne 1 ]; then
    echo "FAIL: metadata/reference_release.tsv has $N_ACTIVE active releases; exactly 1 is required."
    status=1
  elif ! grep -qE "\"$ACTIVE\"[[:space:]]*=[[:space:]]*\"https://[a-z0-9]+\.archive\.ensembl\.org\"" \
          scripts/R/lib/ensembl.R; then
    echo "FAIL: release $ACTIVE has no archive host in SUZ_ENSEMBL_ARCHIVES."
    grep -n 'SUZ_ENSEMBL_ARCHIVES' -A4 scripts/R/lib/ensembl.R | sed 's/^/  /'
    status=1
  else
    HOST=$(awk -F'\t' -v r="$ACTIVE" 'NR>1 && $1==r {print $3; exit}' metadata/reference_release.tsv)
    echo "  ok: active release $ACTIVE -> $HOST"
    # the bash stages must read the same file rather than hard-coding a release
    # Only an ASSIGNMENT that bakes in a release is a problem. Help text that
    # mentions one is fine and often necessary - stage 15 documents that the
    # circRNA branch historically used GRCr8.114.
    BAKED_RX='^[^:]+:[0-9]+:[^#]*=[^=]*GRCr8\.[0-9]+\.gtf'
    if grep -rnE 'GRCr8\.[0-9]+\.gtf' scripts/bash suz_pipeline.bash 2>/dev/null \
         | grep -vE '\$\{RELEASE\}' | grep -qE "$BAKED_RX"; then
      echo "FAIL: a bash stage hard-codes an Ensembl release instead of reading"
      echo "      metadata/reference_release.tsv:"
      grep -rnE 'GRCr8\.[0-9]+\.gtf' scripts/bash suz_pipeline.bash \
        | grep -vE '\$\{RELEASE\}' | grep -E "$BAKED_RX" | sed 's/^/  /'
      status=1
    else
      echo "  ok: every bash stage takes the release from the same file"
    fi
  fi
fi

# 3. nobody reintroduced version= as a substitute for host=
if grep -rnE 'useEnsembl\s*\([^)]*version[[:space:]]*=' scripts/R --include='*.R' \
     | grep -vE ':[0-9]+:\s*#' >/dev/null 2>&1; then
  echo "FAIL: useEnsembl(..., version = ...) found. With an archive host, version="
  echo "      is redundant and its resolution has changed between biomaRt releases."
  grep -rnE 'useEnsembl\s*\([^)]*version[[:space:]]*=' scripts/R --include='*.R' \
    | grep -vE ':[0-9]+:\s*#' | sed 's/^/  /'
  status=1
else
  echo "  ok: no version= pin competing with the host"
fi

# 4. getBM with values= but no filters= (the 00_annotation_biomart.Rmd bug)
bad=$(grep -rn 'getBM' scripts/R --include='*.R' | grep -vE ':[0-9]+:\s*#' \
      | grep 'values' | grep -v 'filters' | grep -v 'suz_getBM')
if [ -n "$bad" ]; then
  echo "FAIL: getBM() with values= but no filters= (biomaRt silently ignores values):"
  echo "$bad" | sed 's/^/  /'
  status=1
else
  echo "  ok: no getBM(values=) without filters="
fi

[ "$status" -eq 0 ] && echo "PASS" || echo "FAILED"
exit $status
