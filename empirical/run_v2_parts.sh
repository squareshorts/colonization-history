#!/bin/bash
# Run every part of the extended empirical analysis (resumable), then merge the results.
# Usage: empirical/run_v2_parts.sh [bootstrap resamples, default 200]
set -e
B="${1:-200}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/.."
( for p in time_matched P1_first3h S2 P2_P3; do python3 empirical/run_v2_final.py $p "$B"; done ) &
( for p in occupancy_adjusted P1_extended; do python3 empirical/run_v2_final.py $p "$B"; done ) &
wait
python3 empirical/run_v2_final.py combine
