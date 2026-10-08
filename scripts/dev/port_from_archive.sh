#!/usr/bin/env bash
# =============================================================================
# port_from_archive.sh — turn an archived original script into a pipeline stage.
#
# The large analysis scripts (up to 2,564 lines) are ported by SURGICAL
# TRANSFORMATION rather than by retyping, so the analysis itself stays
# byte-for-byte identical to what produced the manuscript figures and the diff
# against the archive is small enough to read.
#
# Exactly five classes of change are applied:
#
#   1. A standard preamble is prepended: shebang, provenance header, optparse
#      CLI, lib/ sourcing, and base_dir / out_dir / dpi taken from the CLI.
#   2. Lines assigning base_dir or out_dir to a hard-coded absolute path are
#      commented out (the preamble has already set them).
#   3. Default function arguments of the form base = "C:/Users/sebau/..." are
#      rewritten to base = base_dir.
#   4. Browser-download duplicate filenames -- "foo(1).csv", "bar(2).txt" --
#      are deleted from candidate vectors. They were listed as FIRST
#      preference in six scripts and would silently select a stale copy.
#   5. The "/mnt/data" container search path is dropped.
#
# Everything else is passed through untouched. Run with --check to see the
# resulting diff without writing.
#
# Usage:
#   port_from_archive.sh --src <archived.R> --dst <scripts/R/NN_name.R> \
#                        --stage NN --desc "one line" [--out-default DIR]
#                        [--check]
# =============================================================================

set -euo pipefail

SRC=""; DST=""; STAGE=""; DESC=""; OUT_DEFAULT=""; CHECK=0
while [ $# -gt 0 ]; do
  case "$1" in
    --src) SRC="$2"; shift 2 ;;
    --dst) DST="$2"; shift 2 ;;
    --stage) STAGE="$2"; shift 2 ;;
    --desc) DESC="$2"; shift 2 ;;
    --out-default) OUT_DEFAULT="$2"; shift 2 ;;
    --check) CHECK=1; shift ;;
    --help|-h) sed -n '2,32p' "$0"; exit 0 ;;
    *) echo "[ERROR] unknown argument: $1" >&2; exit 1 ;;
  esac
done
[ -n "$SRC" ] && [ -n "$DST" ] && [ -n "$STAGE" ] || { echo "[ERROR] --src, --dst and --stage are required" >&2; exit 1; }
[ -f "$SRC" ] || { echo "[ERROR] source not found: $SRC" >&2; exit 1; }
OUT_DEFAULT="${OUT_DEFAULT:-results/R${STAGE}_out}"

SRC_BASE="$(basename "$SRC")"
SRC_LINES="$(wc -l < "$SRC")"
TMP="$(mktemp)"; BODY="$(mktemp)"; SRC_LF="$(mktemp)"
trap 'rm -f "$TMP" "$BODY" "$SRC_LF"' EXIT

# Every archived script has CRLF line terminators (they came off a Windows
# machine). Normalise to LF first, otherwise the trailing \r defeats the "$"
# anchors below and makes the audit diff show every line as changed.
tr -d '\r' < "$SRC" > "$SRC_LF"

# ---- body transformation (rules 2-5) ---------------------------------------
sed -E \
  -e 's@^([[:space:]]*)(base_dir|out_dir)([[:space:]]*)<-[[:space:]]*"?[A-Za-z]:/Users/sebau[^\n]*$@\1## [ported] \2 is set from the CLI in the preamble. Original line:\n\1# &@' \
  -e 's@^([[:space:]]*)out_dir[[:space:]]*<-[[:space:]]*file\.path\(base_dir,[^)]*\)[[:space:]]*$@\1## [ported] out_dir is set from the CLI in the preamble. Original line:\n\1# &@' \
  -e 's@^([[:space:]]*)(out_dir|outdir|output_dir)([[:space:]]*)<-([[:space:]]*)"([^"/][^"]*)"[[:space:]]*$@\1\2\3<-\4ensure_out_dir(file.path(.suz_out_root, "\5"))  ## [ported] was a CWD-relative path@' \
  -e "s@^([[:space:]]*)(out_dir|outdir|output_dir)([[:space:]]*)<-([[:space:]]*)file\.path\(getwd\(\),[[:space:]]*'([^']*)'\)[[:space:]]*\$@\1\2\3<-\4ensure_out_dir(file.path(.suz_out_root, '\5'))  ## [ported] was getwd()-relative@" \
  -e 's@base[[:space:]]*=[[:space:]]*"[A-Za-z]:/Users/sebau[^"]*"@base = base_dir@g' \
  -e 's@"[^"]*\([0-9]\)\.(csv|txt|tsv)",[[:space:]]*@@g' \
  -e 's@,[[:space:]]*"[^"]*\([0-9]\)\.(csv|txt|tsv)"@@g' \
  -e 's@c\("\.",[[:space:]]*"/mnt/data"\)@c(".")@g' \
  -e 's@"/mnt/data",[[:space:]]*@@g' \
  -e 's@,[[:space:]]*"/mnt/data"@@g' \
  "$SRC_LF" > "$BODY"

# NOTE: rule 4 deliberately only removes an artifact filename TOGETHER with an
# adjacent comma, i.e. only where it is one element of a candidate vector.
# An earlier version also rewrote a lone artifact string to "" and then
# cleaned up empty elements -- that was too blunt and silently ate legitimate
# empty-string arguments such as str_replace(x, pattern, "") and na_if(out, "").
# If any artifact filename survives as a sole element, the post-check below
# reports it rather than guessing.
if grep -qE '\([0-9]\)\.(csv|txt|tsv)' "$BODY"; then
  echo "[WARN] browser-download duplicate filename(s) survived the transform in $SRC_BASE:" >&2
  grep -nE '\([0-9]\)\.(csv|txt|tsv)' "$BODY" | sed 's/^/[WARN]   /' >&2
  echo "[WARN] remove them by hand before using this stage." >&2
fi

# strip a leading shebang from the body (the preamble supplies one)
sed -i -E '1{/^#!/d}' "$BODY"

# ---- preamble (rule 1) ------------------------------------------------------
cat > "$TMP" <<PREAMBLE
#!/usr/bin/env Rscript
# =============================================================================
# R${STAGE} — ${DESC}
#
# PORTED, NOT REWRITTEN.
# Source: archive/.../${SRC_BASE} (${SRC_LINES} lines)
# Produced by scripts/dev/port_from_archive.sh, which applies exactly five
# classes of change and passes everything else through unmodified:
#   1. this preamble (optparse CLI, lib/ sourcing, base_dir/out_dir/dpi)
#   2. hard-coded C:/Users/sebau/... base_dir and out_dir commented out
#   3. base = "C:/Users/sebau/..." default arguments -> base = base_dir
#   4. browser-download duplicates "name(1).csv" removed from candidate lists
#   5. the /mnt/data container search path removed
#
# Diff this file against the archived original to audit the port:
#   diff archive/.../${SRC_BASE} scripts/R/$(basename "$DST")
# =============================================================================

suppressPackageStartupMessages(library(optparse))

.here <- dirname(normalizePath(sub("^--file=", "",
          commandArgs(trailingOnly = FALSE)[grepl("^--file=",
          commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
source(file.path(.here, "lib", "io.R"))
source(file.path(.here, "lib", "plotting.R"))

.opt <- parse_args(OptionParser(
  option_list = list(
    make_option("--in-dir", type = "character", default = "results",
                help = "root directory searched for inputs [default %default]"),
    make_option("--out-dir", type = "character", default = "${OUT_DEFAULT}",
                help = "output directory [default %default]"),
    make_option("--dpi", type = "integer", default = 600,
                help = "raster figure resolution [default %default]")),
  description = "\nR${STAGE}: ${DESC}\nPorted from ${SRC_BASE}."))

base_dir <- normalizePath(.opt\$\`in-dir\`, winslash = "/", mustWork = FALSE)
out_dir  <- ensure_out_dir(.opt\$\`out-dir\`)
.suz_out_root <- out_dir   ## bodies that rebind out_dir are nested under this
png_dpi  <- .opt\$dpi
dpi_png  <- .opt\$dpi
message("[R${STAGE}] in-dir : ", base_dir)
message("[R${STAGE}] out-dir: ", out_dir)

# ----------------------------- ported body ----------------------------------
PREAMBLE

cat "$BODY" >> "$TMP"

if [ "$CHECK" = "1" ]; then
  echo "=== would write $DST ($(wc -l < "$TMP") lines, from $SRC_LINES) ==="
  diff -u "$SRC_LF" "$BODY" || true
  exit 0
fi

mkdir -p "$(dirname "$DST")"
cp "$TMP" "$DST"
echo "[port] $SRC_BASE ($SRC_LINES lines) -> $DST ($(wc -l < "$DST") lines)"
