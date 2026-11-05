#!/bin/bash
# Extract windowed observables (v1) from every dataset zip of figshare 10.6084/m9.figshare.5969335.
# Usage: empirical/run_all.sh <directory containing the 16 dataset zips>
set -e
SRC="${1:?give the directory that holds the dataset zips}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/observables"
mkdir -p "$OUT"
for z in "$SRC"/*.zip; do
  d="$(basename "$z" .zip)"
  [ -f "$OUT/$d/meta.json" ] && continue
  echo "$(date +%T) $d"
  python3 "$HERE/extract.py" "$z" "$OUT/$d"
done
