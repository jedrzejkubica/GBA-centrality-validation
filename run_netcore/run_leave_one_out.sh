#!/usr/bin/env bash

# Usage:
# ./run_leave_one_out.sh <interactome> <seeds_dir> <permutation_dir> <out_dir>

set -euo pipefail

INTERACTOME="$1"
SEEDS_DIR="$2"
PERM_DIR="$3"
OUT_DIR="$4"

for seedfile in "$SEEDS_DIR"/seeds_*.txt ;
    do
    seed=$(basename "$seedfile" .txt | sed 's/^seeds_//')

    OUT_LEFT_OUT="${OUT_DIR}/output_${seed}/"
    mkdir -p "$OUT_LEFT_OUT"

    python run_netcore/NetCore/netcore/netcore.py \
        -e "$INTERACTOME" \
        -s "$seedfile" \
        -pd "$PERM_DIR" \
        -o "$OUT_LEFT_OUT" \
        1>"$OUT_LEFT_OUT"/log.txt \
        2>&1 ;
    done