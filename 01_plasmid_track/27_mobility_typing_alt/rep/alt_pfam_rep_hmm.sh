#!/bin/bash
# === ALT for rep ===
# Pfam Rep HMM consensus: PF03090 (RepA), PF07042 (TrfA).
# Plus T4CP (PF12696) which is mobility but related to rep operon.
# Broad-coverage rep-domain detection — catches plasmid that PlasmidFinder misses.

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

THREADS=${THREADS:-16}
INPUT_FAA=${INPUT_FAA:-../inputs/plasmidome.master.faa}
ORF2CONTIG=${ORF2CONTIG:-../inputs/orf2contig.tsv}
PFAM_DIR=${PFAM_DIR:-/mnt/nas/DB/geon/pfam_hallmark}
OUT_DIR=${OUT_DIR:-../outputs/01_rep_pfam_hmm}
mkdir -p $OUT_DIR

# Rep-specific Pfams (PF03090 RepA, PF07042 TrfA)
for hmm in PF03090 PF07042; do
  if [ ! -f $PFAM_DIR/${hmm}.hmm.h3i ]; then
    hmmpress $PFAM_DIR/${hmm}.hmm 2>/dev/null || true
  fi
  echo "[$(date '+%F %T')] hmmsearch $hmm"
  hmmsearch --cut_ga --cpu $THREADS --tblout $OUT_DIR/${hmm}.tblout \
    $PFAM_DIR/${hmm}.hmm $INPUT_FAA > $OUT_DIR/${hmm}.stdout
done

# Join ORF -> contig
python3 <<PYEOF
import re
from collections import defaultdict
orf2c = {}
with open("$ORF2CONTIG") as f:
    for line in f:
        o, c = line.rstrip("\n").split("\t")[:2]
        orf2c[o] = c

rep_per_contig = defaultdict(set)
for hmm in ["PF03090", "PF07042"]:
    with open(f"$OUT_DIR/{hmm}.tblout") as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            orf = re.split(r"\s+", line.rstrip())[0]
            c = orf2c.get(orf)
            if c: rep_per_contig[c].add(hmm)

with open("$OUT_DIR/contig_rep_pfam.tsv", "w") as fout:
    for c, hmms in rep_per_contig.items():
        for h in hmms:
            fout.write(f"{c}\t{h}\n")

print(f"  Unique plasmid with Rep Pfam: {len(rep_per_contig)}")
print(f"  Total contig-rep rows: {sum(len(h) for h in rep_per_contig.values())}")
PYEOF
