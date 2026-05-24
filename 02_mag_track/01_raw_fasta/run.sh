#!/bin/bash
# === 01 Raw FASTA: symlink dereplicated MAGs into circ/frag folders ===
# Input:  $PROJECT/07_mag/08_drep_species/dereplicated_genomes/*.fa
#         $PROJECT/07_mag/circ_mag.tsv  (optional: MAG<TAB>circ|frag; default frag)
# Output: $PROJECT/mag/01_raw_fasta/{circ,frag}/<MAG>.fa  (symlinks)
#         $PROJECT/mag/01_raw_fasta/mag_topology.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

SRC=$PROJECT/07_mag/08_drep_species/dereplicated_genomes
MAP=$PROJECT/07_mag/circ_mag.tsv
OUT=$PROJECT/mag/01_raw_fasta
mkdir -p "$OUT/circ" "$OUT/frag"

[ -d "$SRC" ] || { echo "ERROR: no dereplicated_genomes at $SRC" >&2; exit 1; }

echo "[$(date '+%F %T')] symlink MAGs into circ/frag"
topo="$OUT/mag_topology.tsv"
: > "$topo"

# Build a lookup: MAG -> circ|frag (default frag)
declare -A TOPO
if [ -s "$MAP" ]; then
    while IFS=$'\t' read -r m t; do
        [ -z "$m" ] && continue
        TOPO["$m"]="$t"
    done < "$MAP"
    echo "  loaded topology map ($MAP)"
else
    echo "  no $MAP — defaulting all to frag"
fi

n_circ=0; n_frag=0
for fa in "$SRC"/*.fa "$SRC"/*.fasta "$SRC"/*.fna; do
    [ -f "$fa" ] || continue
    mag=$(basename "$fa")
    mag="${mag%.fa}"; mag="${mag%.fasta}"; mag="${mag%.fna}"
    topo_label="${TOPO[$mag]:-frag}"
    [ "$topo_label" = "circ" ] || topo_label="frag"
    dest="$OUT/$topo_label/${mag}.fa"
    [ -L "$dest" ] || ln -sf "$fa" "$dest"
    echo -e "$mag\t$topo_label" >> "$topo"
    if [ "$topo_label" = "circ" ]; then n_circ=$((n_circ+1)); else n_frag=$((n_frag+1)); fi
done

echo "[$(date '+%F %T')] DONE — circ=$n_circ frag=$n_frag"
