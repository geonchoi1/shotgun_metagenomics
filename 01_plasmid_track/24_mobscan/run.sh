#!/bin/bash
# === DEFAULT for mob (relaxase) ===
# MOBscan: 9 MOB family HMM (MOBF, MOBP1, MOBP3, MOBQ, MOBV, MOBH, MOBC, MOBT, MOBB, MOBM).
# Coluzzi 2022 framework; superior for environmental plasmid (vs mob_typer's DIAMOND DB which is clinical-biased).

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

THREADS=${THREADS:-16}
INPUT_FAA=${INPUT_FAA:-../inputs/plasmidome.master.faa}
ORF2CONTIG=${ORF2CONTIG:-../inputs/orf2contig.tsv}
MOBSCAN_HMM=${MOBSCAN_HMM:-/mnt/nas/DB/geon/mobscan_db/MOBfamDB}
OUT_DIR=${OUT_DIR:-../outputs/02_mob_mobscan}
mkdir -p $OUT_DIR

if [ ! -f ${MOBSCAN_HMM}.h3i ]; then
  hmmpress $MOBSCAN_HMM
fi

echo "[$(date '+%F %T')] MOBscan HMM search"
hmmsearch -E 0.01 --domE 0.01 --cpu $THREADS \
  --tblout $OUT_DIR/mobscan.tblout \
  $MOBSCAN_HMM $INPUT_FAA > $OUT_DIR/mobscan.stdout

# ORF -> contig, family filter (HMM coverage > 60%, E-value < 0.01)
python3 <<PYEOF
import re
from collections import defaultdict
orf2c = {}
with open("$ORF2CONTIG") as f:
    for line in f:
        o, c = line.rstrip("\n").split("\t")[:2]
        orf2c[o] = c

contig_mob = defaultdict(set)
with open("$OUT_DIR/mobscan.tblout") as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        parts = re.split(r"\s+", line.rstrip())
        orf, family = parts[0], parts[2]
        c = orf2c.get(orf)
        if c: contig_mob[c].add(family)

with open("$OUT_DIR/contig_mob.tsv", "w") as fout:
    for c, fams in contig_mob.items():
        for f in fams:
            fout.write(f"{c}\t{f}\n")

print(f"  Plasmid with relaxase: {len(contig_mob)}")
print(f"  Family breakdown:")
from collections import Counter
all_fam = Counter()
for fams in contig_mob.values():
    for f in fams: all_fam[f] += 1
for f, n in all_fam.most_common():
    print(f"    {f}: {n}")
PYEOF
