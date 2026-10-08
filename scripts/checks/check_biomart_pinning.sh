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

EXPECTED_HOST="sep2025.archive.ensembl.org"
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

# 2. the pinned default is the expected archive
if grep -qE "SUZ_ENSEMBL_HOST_DEFAULT[[:space:]]*<-[[:space:]]*\"https://${EXPECTED_HOST}\"" \
     scripts/R/lib/ensembl.R; then
  echo "  ok: default host is https://${EXPECTED_HOST} (Ensembl release 115)"
else
  echo "FAIL: the default host in scripts/R/lib/ensembl.R is not https://${EXPECTED_HOST}"
  grep -n 'SUZ_ENSEMBL_HOST_DEFAULT' scripts/R/lib/ensembl.R | sed 's/^/  /'
  status=1
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
