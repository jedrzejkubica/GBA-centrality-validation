############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2024-2026
#
# This file was written by Jędrzej Kubica and Nicolas Thierry-Mieg
# (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
############################################################################################

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