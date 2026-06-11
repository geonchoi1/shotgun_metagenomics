#!/bin/bash
# === PLSDB lookup: per-plasmid Mash + NUCmer + ecosystem positioning ===
# - Mash dist vs PLSDB sketch (D≤0.1, P≤0.1)
# - NUCmer (--delta-filter -i 90 -l 500 -1) validation
# - Join PLSDB metadata (nuccore + taxonomy + biosample)
# - Simka Bray-Curtis k-mer distance (our plasmids vs PLSDB grouped by ecosystem) for NMDS
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

PLASMID=$PROJECT/01_plasmid_track/02_drep/dereplicated.fna
OUT=$PROJECT/01_plasmid_track/60_plsdb_lookup
mkdir -p $OUT/mash $OUT/nucmer $OUT/simka

############ Mash ############
if [ ! -s $OUT/mash/mash_filtered.tsv ]; then
  activate_env "$ENV_MASH"
  echo "[$(date '+%F %T')] Mash vs PLSDB"
  mash dist -p $THREADS -d 0.1 -v 0.1 $PLSDB_MASH_SKETCH $PLASMID > $OUT/mash/mash_all.tsv
  awk -F'\t' '$3<=0.1 && $4<=0.1' $OUT/mash/mash_all.tsv > $OUT/mash/mash_filtered.tsv
  echo "  Mash hits: $(wc -l < $OUT/mash/mash_filtered.tsv)"
fi

############ NUCmer validation ############
if [ ! -s $OUT/nucmer/nucmer_validated.tsv ]; then
  echo "[$(date '+%F %T')] NUCmer validation"
  activate_env "$ENV_NUCMER"
  awk -F'\t' '{print $1}' $OUT/mash/mash_filtered.tsv | sort -u > $OUT/nucmer/plsdb_candidates.txt
  python3 - <<PYEOF
ids=set(open("$OUT/nucmer/plsdb_candidates.txt").read().split())
out=open("$OUT/nucmer/plsdb_subset.fna","w"); keep=False
with open("$PLSDB_FASTA") as f:
    for line in f:
        if line.startswith('>'):
            keep = line[1:].split()[0] in ids
        if keep: out.write(line)
out.close()
PYEOF
  cd $OUT/nucmer
  nucmer --threads $THREADS -p nucmer plsdb_subset.fna $PLASMID > nucmer.log 2>&1 || true
  delta-filter -i 90 -l 500 -1 nucmer.delta > nucmer_filtered.delta
  show-coords -rclT nucmer_filtered.delta > nucmer_validated.tsv
  echo "  NUCmer validated rows: $(wc -l < $OUT/nucmer/nucmer_validated.tsv)"
fi

############ Metadata join ############
echo "[$(date '+%F %T')] PLSDB metadata join"
python3 - <<PYEOF
import os, csv
from collections import defaultdict
OUT="$OUT"
PLSDB_DIR="$PLSDB_DIR"
# Locate the metadata file
meta_candidates=[f"{PLSDB_DIR}/metadata/plsdb.csv",
                 f"{PLSDB_DIR}/plsdb.csv",
                 f"{PLSDB_DIR}/metadata/nuccore.csv"]
META=None
for m in meta_candidates:
    if os.path.exists(m): META=m; break
if not META:
    print("PLSDB metadata not found"); raise SystemExit(0)

# Build acc → row
rows={}
with open(META) as fh:
    r=csv.DictReader(fh)
    for row in r:
        acc=row.get('NUCCORE_ACC') or row.get('ACC_NUCCORE') or row.get('accession') or row.get('nuccore.accession')
        if acc: rows[acc]=row

# Parse mash + nucmer validated hits
plasmid_to_refs=defaultdict(list)
with open(f"{OUT}/mash/mash_filtered.tsv") as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)<5: continue
        ref, qry, dist, pval = p[0], p[1], p[2], p[3]
        plasmid_to_refs[qry].append((ref, float(dist), float(pval)))

validated=set()
with open(f"{OUT}/nucmer/nucmer_validated.tsv") as f:
    for line in f:
        if line.startswith('/') or not line.strip() or line.startswith('NUCMER'): continue
        p=line.split('\t')
        if len(p)<13: continue
        validated.add((p[-1].strip(), p[-2]))

with open(f"{OUT}/per_plasmid_plsdb_lookup.tsv","w") as o:
    o.write("plasmid\tref_acc\tmash_dist\tmash_pval\tnucmer_validated\thost_genus\thost_species\tecosystem\tbiosample\n")
    for q, hits in plasmid_to_refs.items():
        for ref, d, pv in sorted(hits, key=lambda x: x[1]):
            meta=rows.get(ref, {})
            o.write("\t".join([
                q, ref, f"{d:.4f}", f"{pv:.2g}",
                "yes" if (q,ref) in validated else "no",
                meta.get('TAXONOMY_genus',''),
                meta.get('TAXONOMY_species',''),
                meta.get('BIOSAMPLE_ECOSYSTEM','') or meta.get('biosample.ecosystem',''),
                meta.get('BIOSAMPLE_UID','') or meta.get('biosample.uid',''),
            ])+"\n")
print(f"plasmids with ≥1 lookup: {len(plasmid_to_refs)}")
PYEOF

############ Simka ecosystem positioning ############
SK=$OUT/simka
if [ ! -s $SK/mat_abundance_braycurtis.csv ]; then
  echo "[$(date '+%F %T')] Simka ecosystem Bray-Curtis"
  activate_env "$ENV_SIMKA"
  # Build inputs.txt: per-plasmid (our) + per-ecosystem-group PLSDB subset
  mkdir -p $SK/fasta_ours
  python3 - <<PYEOF
import os, csv
ours_seq={}; cur=None
with open("$PLASMID") as f:
    for line in f:
        if line.startswith('>'):
            cur=line[1:].split()[0]; ours_seq[cur]=[]
        else: ours_seq[cur].append(line)
SK="$SK"
os.makedirs(f"{SK}/fasta_ours", exist_ok=True)
inp=open(f"{SK}/inputs.txt","w")
for c,lines in ours_seq.items():
    safe=c.replace('|','_').replace('/','_')
    with open(f"{SK}/fasta_ours/{safe}.fna","w") as o:
        o.write(f">{c}\n"); o.write("".join(lines))
    inp.write(f"ours_{safe}: {SK}/fasta_ours/{safe}.fna\n")
inp.close()
PYEOF
  # Add PLSDB ecosystem groups if metadata available
  python3 - <<PYEOF
import os, csv
from collections import defaultdict
SK="$SK"
PLSDB_DIR="$PLSDB_DIR"
META=None
for m in [f"{PLSDB_DIR}/metadata/plsdb.csv", f"{PLSDB_DIR}/plsdb.csv"]:
    if os.path.exists(m): META=m; break
if not META:
    print("no PLSDB meta — Simka ours-only"); raise SystemExit(0)
acc_eco={}
with open(META) as fh:
    r=csv.DictReader(fh)
    for row in r:
        acc=row.get('NUCCORE_ACC') or row.get('ACC_NUCCORE') or row.get('accession')
        eco=row.get('BIOSAMPLE_ECOSYSTEM') or row.get('biosample.ecosystem','')
        if acc and eco: acc_eco[acc]=eco.split(';')[0].strip().replace(' ','_')
# Bin PLSDB seqs by ecosystem
eco_out={}
import gzip
def open_any(p):
    return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
cur_acc=None; cur_eco=None; cur_buf=[]
with open_any("$PLSDB_FASTA") as f:
    for line in f:
        if line.startswith('>'):
            if cur_acc and cur_eco:
                eco_out.setdefault(cur_eco, open(f"{SK}/plsdb_{cur_eco}.fna","a")).write(">"+cur_acc+"\n"+"".join(cur_buf))
            cur_acc=line[1:].split()[0]
            cur_eco=acc_eco.get(cur_acc)
            cur_buf=[]
        else:
            cur_buf.append(line)
    if cur_acc and cur_eco:
        eco_out.setdefault(cur_eco, open(f"{SK}/plsdb_{cur_eco}.fna","a")).write(">"+cur_acc+"\n"+"".join(cur_buf))
for fh in eco_out.values(): fh.close()
with open(f"{SK}/inputs.txt","a") as inp:
    for eco in eco_out:
        inp.write(f"plsdb_{eco}: {SK}/plsdb_{eco}.fna\n")
print(f"PLSDB ecosystem groups: {len(eco_out)}")
PYEOF
  simka -in $SK/inputs.txt -out $SK -out-tmp $SK/tmp -nb-cores $THREADS -kmer-size 21 -abundance-min 1 || true
  echo "  Simka done"
fi

echo "[$(date '+%F %T')] DONE"
