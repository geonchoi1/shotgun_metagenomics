#!/bin/bash
# === 05 Topology split using USER-provided circ/frag map ===
# Input:  $PROJECT/00_input/circ_frag_map.tsv   (2 cols: contig_id<TAB>{circ|frag})
#                                                contig_id must match <SAMPLE>|<seq> in all_putative.fna
#         $PROJECT/04_plasmid/all_putative.fna
# Output: $PROJECT/05_plasmid_split/circ/all.fna
#         $PROJECT/05_plasmid_split/frag/all.fna

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

MAP=$PROJECT/00_input/circ_frag_map.tsv
IN_FA=$PROJECT/04_plasmid/all_putative.fna
OUT_BASE=$PROJECT/05_plasmid_split

[ -f "$MAP" ]   || { echo "ERROR: USER map missing: $MAP" >&2; exit 2; }
[ -f "$IN_FA" ] || { echo "ERROR: input fasta missing: $IN_FA" >&2; exit 2; }

mkdir -p "$OUT_BASE/circ" "$OUT_BASE/frag"
CIRC_FA=$OUT_BASE/circ/all.fna
FRAG_FA=$OUT_BASE/frag/all.fna

if [ -s "$CIRC_FA" ] && [ -s "$FRAG_FA" ]; then
    echo "[$(date '+%F %T')] 05 split already done"
    exit 0
fi

echo "[$(date '+%F %T')] split all_putative.fna by user circ_frag_map.tsv"
awk -F'\t' -v circ_out="$CIRC_FA" -v frag_out="$FRAG_FA" '
    FNR==NR { topo[$1]=$2; next }
    /^>/ {
        id=$1; sub(/^>/,"",id)
        t = (id in topo) ? topo[id] : ""
        if (t=="circ")      { write="circ" }
        else if (t=="frag") { write="frag" }
        else                { write="" }
        if (write=="circ") print > circ_out
        else if (write=="frag") print > frag_out
        next
    }
    { if (write=="circ") print > circ_out
      else if (write=="frag") print > frag_out }
' "$MAP" "$IN_FA"

nc=$(grep -c '^>' "$CIRC_FA" 2>/dev/null || echo 0)
nf=$(grep -c '^>' "$FRAG_FA" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] DONE — circ:$nc, frag:$nf"
