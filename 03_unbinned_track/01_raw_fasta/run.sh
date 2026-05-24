#!/bin/bash
# === 01 Raw FASTA (UB): symlink unbinned chromosomal all.fna ===
# Input:  $PROJECT/08_chromosomal_split/unbinned/all.fna  (output from 00_shared/08)
#         $PROJECT/08_chromosomal_split/unbinned/circ.ids (optional list of circular contigs)
# Output: $PROJECT/unbinned/01_raw_fasta/all/unbinned.fna   (symlink)
#         optionally: $PROJECT/unbinned/01_raw_fasta/{circ,frag}/*.fna (when circ.ids present)
#         $PROJECT/unbinned/01_raw_fasta/contig_topology.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

SRC=$PROJECT/08_chromosomal_split/unbinned/all.fna
CIRC_IDS=$PROJECT/08_chromosomal_split/unbinned/circ.ids
OUT=$PROJECT/unbinned/01_raw_fasta
mkdir -p "$OUT/all"

[ -s "$SRC" ] || { echo "ERROR: missing $SRC" >&2; exit 1; }

echo "[$(date '+%F %T')] symlink unbinned all.fna -> $OUT/all/unbinned.fna"
ln -sf "$SRC" "$OUT/all/unbinned.fna"

topo="$OUT/contig_topology.tsv"
: > "$topo"

if [ -s "$CIRC_IDS" ]; then
    echo "  circular flag list found -> splitting circ/frag"
    mkdir -p "$OUT/circ" "$OUT/frag"
    awk 'BEGIN{while((getline l < ARGV[1])>0){ids[l]=1} ARGV[1]=""}
         /^>/{name=substr($1,2); is_circ=(name in ids)?1:0;
              outf=(is_circ?"'"$OUT"'/circ/"name".fna":"'"$OUT"'/frag/"name".fna");
              print name"\t"(is_circ?"circ":"frag") >> "'"$topo"'";
              print > outf; next}
         {print > outf}' "$CIRC_IDS" "$SRC"
    n_circ=$(awk -F'\t' '$2=="circ"' "$topo" | wc -l)
    n_frag=$(awk -F'\t' '$2=="frag"' "$topo" | wc -l)
    echo "  circ=$n_circ frag=$n_frag"
else
    echo "  no $CIRC_IDS — treating all contigs as frag (default)"
    grep '^>' "$SRC" | awk '{name=substr($1,2); print name"\tfrag"}' > "$topo"
fi

n_total=$(grep -c '^>' "$SRC" || true)
echo "[$(date '+%F %T')] DONE — $n_total contigs symlinked"
