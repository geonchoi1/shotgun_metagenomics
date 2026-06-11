#!/bin/bash
# === 22 CRISPRCasTyper per-MAG + aggregated NR100 spacers ===
# Input:  $PROJECT/02_mag_track/01_raw_fasta/{circ,frag}/*.fa
# Output: $PROJECT/02_mag_track/22_cctyper/per_mag/<MAG>/  (cctyper output)
#         $PROJECT/02_mag_track/22_cctyper/all_mag_spacers.fa  (NR100 via cd-hit-est)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN_BASE=$PROJECT/02_mag_track/01_raw_fasta
OUT=$PROJECT/02_mag_track/22_cctyper
PER_DIR=$OUT/per_mag
mkdir -p "$PER_DIR"

activate_env "$ENV_CCTYPER"

echo "[$(date '+%F %T')] cctyper per MAG"
for topo in circ frag; do
    for fa in "$IN_BASE/$topo"/*.fa; do
        [ -f "$fa" ] || continue
        mag=$(basename "$fa" .fa)
        out="$PER_DIR/$mag"
        if [ -s "$out/cctyper.log" ] || [ -s "$out/CRISPR_Cas.tab" ] || [ -d "$out/spacers" ]; then
            echo "  $mag already done"
            continue
        fi
        mkdir -p "$(dirname "$out")"
        echo "  $topo/$mag"
        # cctyper requires the output dir to NOT exist
        rm -rf "$out"
        cctyper "$fa" "$out" --threads "$THREADS" --prodigal meta \
            > "$out.log" 2>&1 || { echo "  FAIL $mag — see ${out}.log"; continue; }
        [ -f "$out.log" ] && mv "$out.log" "$out/cctyper.log" || true
    done
done

# Aggregate spacers across MAGs
echo "[$(date '+%F %T')] aggregate spacers"
AGG_RAW="$OUT/all_mag_spacers.raw.fa"
: > "$AGG_RAW"
for sp_dir in "$PER_DIR"/*/spacers; do
    [ -d "$sp_dir" ] || continue
    mag=$(basename "$(dirname "$sp_dir")")
    for f in "$sp_dir"/*.fa; do
        [ -f "$f" ] || continue
        awk -v m="$mag" '/^>/{sub(/^>/,">"m"|"); print; next}{print}' "$f" >> "$AGG_RAW"
    done
done

n_raw=$(grep -c '^>' "$AGG_RAW" 2>/dev/null || echo 0)
if [ "$n_raw" -gt 0 ]; then
    echo "  NR100 via cd-hit-est ($n_raw raw spacers)"
    cd-hit-est -i "$AGG_RAW" -o "$OUT/all_mag_spacers.fa" -c 1.0 -n 8 -M 0 -T "$THREADS" -d 0 > "$OUT/cdhit.log" 2>&1
else
    echo "  no spacers found"
    : > "$OUT/all_mag_spacers.fa"
fi
n_nr=$(grep -c '^>' "$OUT/all_mag_spacers.fa" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] DONE — spacers raw=$n_raw NR100=$n_nr"
