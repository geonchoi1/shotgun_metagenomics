#!/bin/bash
# === 01 Raw FASTA (UB): split into circ/+frag/ multi-FASTA ===
# Input:  $PROJECT/00_shared/08_split_binned_unbinned/unbinned/all.fna
#         $PROJECT/00_shared/08_split_binned_unbinned/unbinned/circ.ids (optional one-per-line)
# Output: $PROJECT/03_unbinned_track/01_raw_fasta/circ/unbinned.fna   (multi-FASTA)
#         $PROJECT/03_unbinned_track/01_raw_fasta/frag/unbinned.fna   (multi-FASTA)
#         $PROJECT/03_unbinned_track/01_raw_fasta/contig_topology.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

SRC=$PROJECT/00_shared/08_split_binned_unbinned/unbinned/all.fna
CIRC_IDS=$PROJECT/00_shared/08_split_binned_unbinned/unbinned/circ.ids
OUT=$PROJECT/03_unbinned_track/01_raw_fasta
mkdir -p "$OUT"/{circ,frag}

[ -s "$SRC" ] || { echo "ERROR: missing $SRC" >&2; exit 1; }

topo="$OUT/contig_topology.tsv"
: > "$topo"

if [ -s "$CIRC_IDS" ]; then
    echo "[$(date '+%F %T')] Splitting circ + frag by $CIRC_IDS"
    python3 - <<PYEOF
import os
SRC = "$SRC"; CIRC_IDS = "$CIRC_IDS"; OUT = "$OUT"
circ = set(line.strip() for line in open(CIRC_IDS) if line.strip())
circ_fh = open(f"{OUT}/circ/unbinned.fna", "w")
frag_fh = open(f"{OUT}/frag/unbinned.fna", "w")
topo = open(f"{OUT}/contig_topology.tsv", "w")
cur, name, is_circ = [], None, False
def flush():
    if cur and name:
        (circ_fh if is_circ else frag_fh).writelines(cur)
        topo.write(f"{name}\t{'circ' if is_circ else 'frag'}\n")
for line in open(SRC):
    if line.startswith(">"):
        flush(); cur=[]
        name = line[1:].split()[0]
        is_circ = name in circ
    cur.append(line)
flush()
circ_fh.close(); frag_fh.close(); topo.close()
nc = sum(1 for _ in open(f"{OUT}/contig_topology.tsv") if "\tcirc\n" in open(f"{OUT}/contig_topology.tsv").read()[:0]+"\tcirc\n" or True)
import subprocess
print(f"  circ: {subprocess.run(['grep','-c','^>',f'{OUT}/circ/unbinned.fna'],capture_output=True,text=True).stdout.strip()}")
print(f"  frag: {subprocess.run(['grep','-c','^>',f'{OUT}/frag/unbinned.fna'],capture_output=True,text=True).stdout.strip()}")
PYEOF
else
    echo "[$(date '+%F %T')] No $CIRC_IDS — all → frag"
    : > "$OUT/circ/unbinned.fna"
    ln -sf "$SRC" "$OUT/frag/unbinned.fna"
    grep '^>' "$SRC" | awk '{name=substr($1,2); print name"\tfrag"}' > "$topo"
fi

n_circ=$(grep -c '^>' "$OUT/circ/unbinned.fna" 2>/dev/null || echo 0)
n_frag=$(grep -c '^>' "$OUT/frag/unbinned.fna" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] DONE — circ=$n_circ, frag=$n_frag"
