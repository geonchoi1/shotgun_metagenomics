#!/bin/bash
# === ALT for mob (relaxase) ===
# Pfam PF03432 (MOB Relaxase) HMM single-profile search.
# Catches MOB-like proteins that broader MOBscan family HMMs might split.

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

THREADS=${THREADS:-16}
INPUT_FAA=${INPUT_FAA:-../inputs/plasmidome.master.faa}
ORF2CONTIG=${ORF2CONTIG:-../inputs/orf2contig.tsv}
PFAM_DIR=${PFAM_DIR:-/mnt/nas/DB/geon/pfam_hallmark}
OUT_DIR=${OUT_DIR:-../outputs/02_mob_pfam}
mkdir -p $OUT_DIR

if [ ! -f $PFAM_DIR/PF03432.hmm.h3i ]; then
  hmmpress $PFAM_DIR/PF03432.hmm 2>/dev/null || true
fi

echo "[$(date '+%F %T')] PF03432 (MOB Relaxase) hmmsearch"
hmmsearch --cut_ga --cpu $THREADS \
  --tblout $OUT_DIR/PF03432.tblout \
  $PFAM_DIR/PF03432.hmm $INPUT_FAA > $OUT_DIR/PF03432.stdout

python3 <<PYEOF
import re
from collections import defaultdict
orf2c = {}
with open("$ORF2CONTIG") as f:
    for line in f: o, c = line.rstrip("\n").split("\t")[:2]; orf2c[o] = c
contig_mob = defaultdict(int)
with open("$OUT_DIR/PF03432.tblout") as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        orf = re.split(r"\s+", line.rstrip())[0]
        c = orf2c.get(orf)
        if c: contig_mob[c] += 1
with open("$OUT_DIR/contig_mob.tsv", "w") as fout:
    for c in contig_mob:
        fout.write(f"{c}\tPF03432_MOB\n")
print(f"  Plasmid with PF03432 hit: {len(contig_mob)}")
PYEOF
