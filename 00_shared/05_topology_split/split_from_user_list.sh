#!/bin/bash
# === 05 Topology split: split plasmids (all_putative.fna) into circ / frag ===
# Topology is read from the MASTER contig_topology.tsv, auto-generated from metaFlye
# assembly_info (col 4 "circ." Y/N). No user-provided map needed — the same master is
# the single source of truth for plasmid (here), MAG and unbinned topology.
#
# Input:  $PROJECT/00_shared/02_assembly/contig_topology.tsv  (auto; <sample>|<contig> TAB circ|frag)
#         $PROJECT/00_shared/04_genomad_plasmid/all_putative.fna
# Output: $PROJECT/00_shared/05_topology_split/{circ,frag}/all.fna

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

MASTER=$PROJECT/00_shared/02_assembly/contig_topology.tsv
IN_FA=$PROJECT/00_shared/04_genomad_plasmid/all_putative.fna
OUT_BASE=$PROJECT/00_shared/05_topology_split

# Auto-generate the master topology from metaFlye assembly_info if missing.
if [ ! -s "$MASTER" ]; then
    echo "[$(date '+%F %T')] generating master contig_topology.tsv from assembly_info"
    : > "$MASTER"
    for info in "$PROJECT"/00_shared/02_assembly/metaflye/*/assembly_info.txt; do
        [ -s "$info" ] || continue
        s=$(basename "$(dirname "$info")")
        awk -F'\t' -v s="$s" 'NR>1{print s"|"$1"\t"($4=="Y"?"circ":"frag")}' "$info" >> "$MASTER"
    done
fi

[ -s "$MASTER" ] || { echo "ERROR: topology map empty: $MASTER (no assembly_info?)" >&2; exit 2; }
[ -f "$IN_FA" ]  || { echo "ERROR: input fasta missing: $IN_FA" >&2; exit 2; }

mkdir -p "$OUT_BASE/circ" "$OUT_BASE/frag"
CIRC_FA=$OUT_BASE/circ/all.fna
FRAG_FA=$OUT_BASE/frag/all.fna

if [ -s "$CIRC_FA" ] && [ -s "$FRAG_FA" ]; then
    echo "[$(date '+%F %T')] 05 split already done"
    exit 0
fi

echo "[$(date '+%F %T')] split all_putative.fna by master contig_topology.tsv"
awk -F'\t' -v circ_out="$CIRC_FA" -v frag_out="$FRAG_FA" '
    FNR==NR { topo[$1]=$2; next }
    /^>/ {
        id=$1; sub(/^>/,"",id)
        t = (id in topo) ? topo[id] : "frag"   # default frag if absent
        write = (t=="circ") ? "circ" : "frag"
        print > (write=="circ" ? circ_out : frag_out)
        next
    }
    { print > (write=="circ" ? circ_out : frag_out) }
' "$MASTER" "$IN_FA"

nc=$(grep -c '^>' "$CIRC_FA" 2>/dev/null || echo 0)
nf=$(grep -c '^>' "$FRAG_FA" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] DONE — circ:$nc, frag:$nf"
