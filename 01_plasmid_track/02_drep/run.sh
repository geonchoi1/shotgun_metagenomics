#!/bin/bash
# === RBH 100/100 dereplication on putative plasmids (Fiamenghi 2025) ===
# - all-vs-all megablast on all_putative.fna
# - filter: pident == 100, aln_length/qlen == 1.0, RBH (both directions)
# - one rep per connected component (longest)
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

CIRC_FNA=$PROJECT/05_plasmid_split/circ/all.fna
FRAG_FNA=$PROJECT/05_plasmid_split/frag/all.fna
OUT=$PROJECT/plasmid/02_drep
mkdir -p $OUT

ALL_FNA=$OUT/all_putative.fna
if [ ! -s $ALL_FNA ]; then
  cat $CIRC_FNA $FRAG_FNA > $ALL_FNA
fi

# circ_frag_map.tsv: contig_id -> topology(circ|frag)
if [ ! -s $OUT/circ_frag_map.tsv ]; then
  grep '^>' $CIRC_FNA | awk '{print substr($1,2)"\tcirc"}' >  $OUT/circ_frag_map.tsv
  grep '^>' $FRAG_FNA | awk '{print substr($1,2)"\tfrag"}' >> $OUT/circ_frag_map.tsv
fi

N_IN=$(grep -c '^>' $ALL_FNA)
echo "[$(date '+%F %T')] drep input: $N_IN putative plasmids"

# ---- BLAST all-vs-all ----
activate_env "$ENV_DIAMOND"   # base env has BLAST+
if [ ! -s $OUT/blastn_full.tsv ]; then
  makeblastdb -in $ALL_FNA -dbtype nucl -out $OUT/plasmids 2>&1 | tail -3
  blastn -task megablast \
    -query $ALL_FNA -db $OUT/plasmids \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
    -evalue 1e-3 -num_threads $THREADS_BLAST \
    -max_target_seqs 25000 -perc_identity 0.0 \
    -out $OUT/blastn_full.tsv
fi
echo "[$(date '+%F %T')] BLAST hits: $(wc -l < $OUT/blastn_full.tsv)"

# ---- RBH 100/100 filter + connected components ----
if [ ! -s $OUT/dereplicated.fna ]; then
python3 - <<PYEOF
import os
from collections import defaultdict
import networkx as nx

OUT = "$OUT"
ALL_FNA = "$ALL_FNA"

# Read seq lengths
plen, seq = {}, {}
cur = None
with open(ALL_FNA) as f:
    for line in f:
        if line.startswith('>'):
            cur = line[1:].split()[0]
            seq[cur] = []
        else:
            seq[cur].append(line.strip())
for k,v in seq.items():
    plen[k] = sum(len(x) for x in v)
print(f"Sequences: {len(plen)}")

# Aggregate BLAST per (q,s): total aln, aln-weighted pident
pair = defaultdict(lambda: {'aln':0, 'pid_sum':0.0})
with open(f"{OUT}/blastn_full.tsv") as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        q, s, pid, aln, ql = p[0], p[1], float(p[2]), int(p[3]), int(p[12])
        if q == s: continue
        pair[(q,s)]['aln'] += aln
        pair[(q,s)]['pid_sum'] += aln * pid

# Compute metrics: tANI%, AF% (per query)
metrics = {}
for (q,s),d in pair.items():
    tANI = d['pid_sum']/d['aln'] if d['aln']>0 else 0
    af = 100.0 * d['aln'] / plen[q]
    metrics[(q,s)] = (tANI, af)

# RBH 100/100: both (q,s) and (s,q) with tANI>=99.99 AND AF>=99.99
RBH = []
for (q,s),(tANI,af) in metrics.items():
    if q >= s: continue
    if (s,q) not in metrics: continue
    tANI2, af2 = metrics[(s,q)]
    if tANI>=99.99 and af>=99.99 and tANI2>=99.99 and af2>=99.99:
        RBH.append((q,s))
print(f"RBH 100/100 pairs: {len(RBH)}")

# Connected components → keep longest as rep
G = nx.Graph()
for p in plen: G.add_node(p)
for a,b in RBH: G.add_edge(a,b)

reps = []
members = []
for ci, comp in enumerate(nx.connected_components(G)):
    rep = max(comp, key=lambda x: plen[x])
    reps.append(rep)
    for m in comp:
        members.append((rep, m, plen[m]))
print(f"Clusters / reps: {len(reps)}")

# Write outputs
reps_set = set(reps)
with open(f"{OUT}/dereplicated.fna","w") as fout:
    for r in sorted(reps_set):
        fout.write(f">{r}\n")
        for chunk in seq[r]:
            fout.write(chunk+"\n")
with open(f"{OUT}/drep_members.tsv","w") as fout:
    fout.write("rep\tmember\tlength\n")
    for r,m,l in members:
        fout.write(f"{r}\t{m}\t{l}\n")
print("Wrote dereplicated.fna + drep_members.tsv")
PYEOF
fi

# ---- Split dereplicated.fna by topology (circ / frag) ----
python3 - <<PYEOF
OUT = "$OUT"
topo = {}
with open(f"{OUT}/circ_frag_map.tsv") as f:
    for line in f:
        c,t = line.rstrip('\n').split('\t')
        topo[c] = t
with open(f"{OUT}/dereplicated.fna") as f, \
     open(f"{OUT}/circ.fna","w") as oc, \
     open(f"{OUT}/frag.fna","w") as ofg:
    out=None
    for line in f:
        if line.startswith('>'):
            cid = line[1:].split()[0]
            out = oc if topo.get(cid)=='circ' else ofg
        out.write(line)
PYEOF

echo "[$(date '+%F %T')] DONE — dereplicated: $(grep -c '^>' $OUT/dereplicated.fna) (circ=$(grep -c '^>' $OUT/circ.fna), frag=$(grep -c '^>' $OUT/frag.fna))"
