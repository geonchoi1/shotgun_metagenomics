#!/bin/bash
# === Step 1: verify MGE annotations exist in all 3 tracks ===
# Checks ISEScan + IntegronFinder + ICEberg3 outputs in plasmid/MAG/UB tracks.
# Output: $PROJECT/04_cross_track/mobile_arg/step1/availability.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

OUT=$PROJECT/04_cross_track/mobile_arg/step1
mkdir -p "$OUT"

# (source_label, track_root, isescan, integronfinder, iceberg3, amrfinder)
declare -a SRC=(plasmid mag unbinned)

declare -A IS=(
    [plasmid]="$PROJECT/01_plasmid_track/20_isescan/all"
    [mag]="$PROJECT/02_mag_track/19_isescan/all"
    [unbinned]="$PROJECT/03_unbinned_track/19_isescan/all"
)
declare -A IF=(
    [plasmid]="$PROJECT/01_plasmid_track/21_integronfinder/all"
    [mag]="$PROJECT/02_mag_track/20_integronfinder/all"
    [unbinned]="$PROJECT/03_unbinned_track/20_integronfinder/all"
)
declare -A IC=(
    [plasmid]="$PROJECT/01_plasmid_track/22_iceberg3/all/iceberg3.tsv"
    [mag]="$PROJECT/02_mag_track/21_iceberg3/all/iceberg3.tsv"
    [unbinned]="$PROJECT/03_unbinned_track/21_iceberg3/all/iceberg3.tsv"
)
declare -A AMR=(
    [plasmid]="$PROJECT/01_plasmid_track/09_amrfinder/all/amrfinder.tsv"
    [mag]="$PROJECT/02_mag_track/08_amrfinder/all/amrfinder.tsv"
    [unbinned]="$PROJECT/03_unbinned_track/08_amrfinder/all/amrfinder.tsv"
)

avail="$OUT/availability.tsv"
{
    printf "source\tisescan_ok\tintegronfinder_ok\ticeberg3_ok\tamrfinder_ok\n"
    fail=0
    for s in "${SRC[@]}"; do
        is_ok=$(ls "${IS[$s]}"/*.tsv 2>/dev/null | head -1)
        if_ok=$(find "${IF[$s]}" -mindepth 1 -name "Results_Integron_Finder_*" 2>/dev/null | head -1)
        ic_ok=$([ -s "${IC[$s]}" ] && echo yes || echo no)
        ar_ok=$([ -s "${AMR[$s]}" ] && echo yes || echo no)
        [ -n "$is_ok" ] && is_v=yes || is_v=no
        [ -n "$if_ok" ] && if_v=yes || if_v=no
        printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$is_v" "$if_v" "$ic_ok" "$ar_ok"
        for v in "$is_v" "$if_v" "$ic_ok" "$ar_ok"; do
            [ "$v" = "no" ] && fail=1
        done
    done
    if [ "$fail" = "1" ]; then
        echo "WARNING: some tracks missing MGE/ARG outputs — Step 1 partial" >&2
    fi
} > "$avail"

cat "$avail"
echo "[$(date '+%F %T')] Step 1 DONE — $avail"
