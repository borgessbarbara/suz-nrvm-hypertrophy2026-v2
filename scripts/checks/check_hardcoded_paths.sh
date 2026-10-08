#!/usr/bin/env bash
# Fail if a machine-specific path or a browser-download duplicate filename has
# crept back into the runnable pipeline.
#
# Between them the original scripts carried five different machine roots:
#   C:/Users/sebau/OneDrive/Documentos/...        (11 raw + 6 final scripts)
#   C:/Users/sebau/Downloads                      (18_curated_Fig4A_Fig6C)
#   C:/Users/Barbara/OneDrive/.../seba-short-reads (the three .qmd execute-dirs)
#   /home/vmaracaja/shared-projects/SUZ_sr_2025   (two metadata-chain .sh)
#   /home/vmaracaja/shared-projects/sebauz_sr     (filter_out_lncRNAs.sh)
#   D:/ and E:/ drive letters                     (the circRNA .Rmd notebooks)
#   /mnt/data                                     (36b_, a container path)
#
# and six final scripts listed "name(1).csv" style browser-download duplicates
# as FIRST-preference input candidates.
#
# archive/ is exempt - it is every received package, preserved verbatim on purpose.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1

status=0
SCOPE=(scripts/R scripts/bash scripts/commands scripts/checks metadata suz_pipeline.bash)

echo "== hard-coded path check =="
echo "   scope: ${SCOPE[*]}"
echo "   exempt: archive/ (every received package, verbatim by design)"

# Lines that are comments are allowed: the porter deliberately preserves the
# original absolute paths as commented provenance.
scan() {
  local label="$1" pattern="$2"
  local hits
  hits=$(grep -rnE "$pattern" "${SCOPE[@]}" 2>/dev/null \
         | grep -vE ':[0-9]+:[[:space:]]*(#|##)' \
         | grep -vE '^scripts/checks/check_hardcoded_paths\.sh:')
  if [ -n "$hits" ]; then
    echo "FAIL: $label"
    echo "$hits" | sed 's/^/  /'
    status=1
  else
    echo "  ok: no $label"
  fi
}

scan "Windows user paths (C:/Users/...)"            'C:/Users/[A-Za-z]'
scan "other Windows drive letters (D:/ E:/)"        '[[:space:]"(][DE]:/[A-Za-z]'
scan "vmaracaja server roots"                       '/home/vmaracaja'
scan "the /mnt/data container path"                 '/mnt/data'
scan "browser-download duplicates, name(N).ext"     '\([0-9]\)\.(csv|txt|tsv|xlsx)'

# The circRNA stages must never bake in a /media/storage path, but they are
# allowed to NAME one inside a heredoc or an echo, because telling the user
# where the historical resource lived is exactly what those messages are for
# (e.g. 14_circrna_ciriquant.sh points at the missing rn8.yml). So strip
# heredoc bodies before scanning.
strip_heredocs() {
  awk '
    # entering a heredoc: <<EOF, <<-EOF, <<"EOF", <<'"'"'EOF'"'"'
    !inhd && match($0, /<<-?[[:space:]]*['"'"'"]?[A-Za-z_][A-Za-z0-9_]*['"'"'"]?/) {
      tag = substr($0, RSTART, RLENGTH)
      gsub(/^<<-?[[:space:]]*/, "", tag); gsub(/['"'"'"]/, "", tag)
      inhd = 1; hdtag = tag; print ""; next
    }
    inhd {
      t = $0; gsub(/^[[:space:]\t]+/, "", t)
      if (t == hdtag) { inhd = 0 }
      print ""; next
    }
    { print }
  ' "$1"
}

hits=""
while IFS= read -r file; do
  found=$(strip_heredocs "$file" | grep -nE '/media/storage[0-9]' \
          | grep -vE '^[0-9]+:[[:space:]]*(#|##)' || true)
  [ -n "$found" ] && hits="${hits}${file}:${found}"$'\n'
done < <(find scripts -type f \( -name '*.sh' -o -name '*.R' \) ! -path 'scripts/checks/*')
hits=$(printf '%s' "$hits" | sed '/^$/d')
if [ -n "$hits" ]; then
  echo "FAIL: /media/storage... path outside a comment"
  echo "$hits" | sed 's/^/  /'
  status=1
else
  echo "  ok: no literal /media/storage paths"
fi

[ "$status" -eq 0 ] && echo "PASS" || echo "FAILED"
exit $status
