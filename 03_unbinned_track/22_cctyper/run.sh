#!/bin/bash
# === 22 CRISPRCasTyper on UB contigs ===
# Detect Cas operons + extract CRISPR spacers from unbinned chromosomal contigs.
# Contig-based (not ORF-based) — independent of Bakta re-runs.
#
# Input:  $PROJECT/unbinned/01_raw_fasta/{circ,frag}/unbinned.fna
# Output: $PROJECT/unbinned/22_cctyper/{circ,frag}/  (cctyper output)
#         $PROJECT/unbinned/22_cctyper/all_ub_spacers.fa  (NR100 via cd-hit-est)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN_BASE=$PROJECT/unbinned/01_raw_fasta
OUT=$PROJECT/unbinned/22_cctyper
mkdir -p "$OUT"

activate_env "$ENV_CCTYPER"

echo "[$(date '+%F %T')] cctyper on UB (circ + frag)"
for topo in circ frag; do
    fa="$IN_BASE/$topo/unbinned.fna"
    out="$OUT/$topo"
    [ -s "$fa" ] || { echo "  $topo: no input"; continue; }
    if [ -s "$out/CRISPR_Cas.tab" ] || [ -d "$out/spacers" ]; then
        echo "  $topo already done"
        continue
    fi
    rm -rf "$out"
    echo "  running cctyper on $topo"
    cctyper "$fa" "$out" --threads "$THREADS" --prodigal meta \
        > "$out.log" 2>&1 || { echo "  FAIL $topo — see ${out}.log"; continue; }
    [ -f "$out.log" ] && mv "$out.log" "$out/cctyper.log" || true
done

# Aggregate spacers across topologies
echo "[$(date '+%F %T')] aggregate UB spacers"
AGG_RAW="$OUT/all_ub_spacers.raw.fa"
: > "$AGG_RAW"
for topo in circ frag; do
    sp_dir="$OUT/$topo/spacers"
    [ -d "$sp_dir" ] || continue
    for f in "$sp_dir"/*.fa; do
        [ -f "$f" ] || continue
        awk -v t="$topo" '/^>/{sub(/^>/,">UB_"t"|"); print; next}{print}' "$f" >> "$AGG_RAW"
    done
done

n_raw=$(grep -c '^>' "$AGG_RAW" 2>/dev/null || echo 0)
if [ "$n_raw" -gt 0 ]; then
    echo "  NR100 via cd-hit-est ($n_raw raw spacers)"
    cd-hit-est -i "$AGG_RAW" -o "$OUT/all_ub_spacers.fa" -c 1.0 -n 8 -M 0 -T "$THREADS" -d 0 > "$OUT/cdhit.log" 2>&1
else
    echo "  no spacers found"
    : > "$OUT/all_ub_spacers.fa"
fi
n_nr=$(grep -c '^>' "$OUT/all_ub_spacers.fa" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] DONE — spacers raw=$n_raw NR100=$n_nr"
