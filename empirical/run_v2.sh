#!/bin/bash
# Continuous tracking with per-cell fluorescence readings for the wild-type transfer datasets.
# Usage: empirical/run_v2.sh <directory containing the dataset zips> [parallel jobs, default 2]
set -e
SRC="${1:?give the directory that holds the dataset zips}"
JOBS="${2:-2}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/observables"
ls "$SRC"/"WT YFP-CFP"*.zip | while read -r z; do
  d="$(basename "$z" .zip)"
  [ -f "$OUT/$d/meta_full.json" ] && continue
  printf '%s\0' "$z"
done | xargs -0 -P "$JOBS" -I{} sh -c 'd=$(basename "{}" .zip); python3 "'"$HERE"'/extract_v2.py" "{}" "'"$OUT"'/$d" > /dev/null && echo "done $d"'
