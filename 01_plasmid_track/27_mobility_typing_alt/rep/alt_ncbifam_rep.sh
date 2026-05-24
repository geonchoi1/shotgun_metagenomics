#!/bin/bash
# === ALT for rep ===
# NCBIfam plasmid Rep-initiator curated catalog (20 initiator profiles, 24 with partition).
# Very precise — but small library; mostly catches well-known incompatibility groups.

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

THREADS=${THREADS:-16}
INPUT_FAA=${INPUT_FAA:-../inputs/plasmidome.master.faa}
ORF2CONTIG=${ORF2CONTIG:-../inputs/orf2contig.tsv}
NCBIFAM_HMM=${NCBIFAM_HMM:-/mnt/nas/DB/geon/ncbifam/hmm_PGAP.LIB}
OUT_DIR=${OUT_DIR:-../outputs/01_rep_ncbifam}
mkdir -p $OUT_DIR

echo "[$(date '+%F %T')] hmmsearch NCBIfam (full library)"
hmmsearch --cut_ga --cpu $THREADS \
  --tblout $OUT_DIR/ncbifam.tblout \
  $NCBIFAM_HMM $INPUT_FAA > $OUT_DIR/ncbifam.stdout

# Filter for plasmid-rep-initiator profiles (NF00000-NF99999 range, but specifically Plasmid_Rep_Initiator family)
# Known curated NF accession prefixes for plasmid replication initiators
python3 <<PYEOF
import re
from collections import defaultdict
orf2c = {}
with open("$ORF2CONTIG") as f:
    for line in f:
        o, c = line.rstrip("\n").split("\t")[:2]
        orf2c[o] = c

# Filter HMM description for "rep" / "initiator" / "Inc" keywords
rep_per_contig = defaultdict(set)
with open("$OUT_DIR/ncbifam.tblout") as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        parts = re.split(r"\s+", line.rstrip(), maxsplit=18)
        if len(parts) < 4: continue
        orf, hmm_name = parts[0], parts[2]
        desc = parts[-1] if len(parts) > 18 else ""
        # Keep only plasmid-replication-related profiles
        if not re.search(r"(plasmid.*rep|initiator|RepA|RepB|RepC|Inc[A-Z])", hmm_name + " " + desc, re.I):
            continue
        c = orf2c.get(orf)
        if c: rep_per_contig[c].add(hmm_name)

with open("$OUT_DIR/contig_rep_ncbifam.tsv", "w") as fout:
    for c, hmms in rep_per_contig.items():
        for h in hmms:
            fout.write(f"{c}\t{h}\n")

print(f"  Plasmid with NCBIfam rep initiator: {len(rep_per_contig)}")
PYEOF
