#!/bin/bash
# === RBH 100/100 dereplication on putative plasmids (Fiamenghi 2025 + Camargo formula) ===
# - all-vs-all megablast on all_putative.fna
# - ANI/AF calculation: Camargo anicalc.py (HSP union AF + pruned tANI, more accurate
#   than simple HSP sum which over-counts internal repeats)
# - filter: ani ≥ 99.99% AND qcov ≥ 99.99% AND tcov ≥ 99.99% (= RBH 100/100 both directions)
# - one rep per connected component (longest)
#
# Refs:
#   - Fiamenghi 2025 Nat Comm (10.1038/s41467-025-65102-6): RBH 100/100 dereplication
#   - Camargo 2024 Nat Microbiol: anicalc.py formula
#     (github.com/apcamargo/bioinformatics-snakemake-pipelines/contig-ani-leiden-clustering-pipeline)
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

# ---- BLAST all-vs-all (Camargo-compatible 11-column outfmt) ----
activate_env "$ENV_DIAMOND"   # base env has BLAST+
if [ ! -s $OUT/blastn_full.tsv ]; then
  makeblastdb -in $ALL_FNA -dbtype nucl -out $OUT/plasmids 2>&1 | tail -3
  blastn -task megablast \
    -query $ALL_FNA -db $OUT/plasmids \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue qlen slen' \
    -evalue 1e-3 -num_threads $THREADS_BLAST \
    -max_target_seqs 25000 -perc_identity 0.0 \
    -out $OUT/blastn_full.tsv
fi
echo "[$(date '+%F %T')] BLAST hits: $(wc -l < $OUT/blastn_full.tsv)"

# ---- ANI/AF via Camargo anicalc.py (HSP union + prune) ----
CAMARGO_PIPELINE=${CAMARGO_PIPELINE:-$HOME/tools/bioinformatics-snakemake-pipelines/contig-ani-leiden-clustering-pipeline}
if [ ! -s $OUT/ani.tsv ]; then
  echo "[$(date '+%F %T')] Computing ANI/AF via Camargo anicalc.py logic"
  python3 - <<PYEOF
from collections import namedtuple
Hsp = namedtuple("Hsp", ["qname","tname","pid","len","qcoords","tcoords","evalue","qlen","tlen"])

def parse_blast(handle):
    for line in handle:
        p = line.split()
        yield Hsp(p[0], p[1], float(p[2])/100, float(p[3]),
                  sorted([int(p[4]), int(p[5])]),
                  sorted([int(p[6]), int(p[7])]),
                  float(p[8]), float(p[9]), float(p[10]))

def yield_alignment_blocks(handle):
    key, alns = None, None
    for aln in parse_blast(handle):
        key = (aln.qname, aln.tname); alns = [aln]; break
    for aln in parse_blast(handle):
        if (aln.qname, aln.tname) == key:
            alns.append(aln)
        else:
            yield alns
            key = (aln.qname, aln.tname); alns = [aln]
    yield alns

def prune_alns(alns, max_evalue=1e-3):
    keep, cur_aln = [], 0
    qry_len = alns[0].qlen
    for aln in alns:
        qcoords = aln.qcoords
        aln_len = max(qcoords) - min(qcoords) + 1
        if aln.evalue > max_evalue: continue
        if cur_aln >= qry_len or aln_len + cur_aln >= 1.10 * qry_len: break
        keep.append(aln); cur_aln += aln_len
    return keep

def compute_ani(alns):
    return round(sum(a.len*a.pid for a in alns)/sum(a.len for a in alns), 4)

def compute_cov(alns):
    def union(coords_list, total_len):
        coords = sorted(coords_list)
        nr = [coords[0]]
        for s,e in coords[1:]:
            if s <= nr[-1][1]+1: nr[-1][1] = max(nr[-1][1], e)
            else: nr.append([s,e])
        return round(sum(e-s+1 for s,e in nr)/total_len, 4)
    qcov = union([a.qcoords for a in alns], alns[0].qlen)
    tcov = union([a.tcoords for a in alns], alns[0].tlen)
    return qcov, tcov

with open(f"$OUT/blastn_full.tsv") as fin, open(f"$OUT/ani.tsv","w") as fo:
    fo.write("contig_1\tcontig_2\tnum_alns\tani\tqcov\ttcov\n")
    for alns in yield_alignment_blocks(fin):
        alns = prune_alns(alns)
        if not alns: continue
        ani = compute_ani(alns); qcov, tcov = compute_cov(alns)
        fo.write(f"{alns[0].qname}\t{alns[0].tname}\t{len(alns)}\t{ani}\t{qcov}\t{tcov}\n")
print("ani.tsv done")
PYEOF
fi

# ---- RBH 100/100 filter + connected components ----
if [ ! -s $OUT/dereplicated.fna ]; then
python3 - <<PYEOF
import os
import networkx as nx
from collections import defaultdict

# Read sequences + lengths
plen, seq = {}, {}
cur = None
with open("$ALL_FNA") as f:
    for line in f:
        if line.startswith('>'):
            cur = line[1:].split()[0]; seq[cur] = []
        else: seq[cur].append(line.strip())
for k,v in seq.items(): plen[k] = sum(len(x) for x in v)
print(f"Sequences: {len(plen)}")

# Read ani.tsv → RBH 100/100 pairs
# Both qcov and tcov ≥ 99.99% (= AND both directions, not OR)
RBH_pairs = []
with open("$OUT/ani.tsv") as f:
    next(f)  # header
    for line in f:
        c1, c2, _, ani, qcov, tcov = line.rstrip().split('\t')
        if c1 == c2: continue
        if c1 >= c2: continue  # unique unordered pair
        ani, qcov, tcov = float(ani), float(qcov), float(tcov)
        if ani >= 0.9999 and qcov >= 0.9999 and tcov >= 0.9999:
            RBH_pairs.append((c1, c2))
print(f"RBH 100/100 pairs: {len(RBH_pairs)}")

# Connected components → keep longest as rep
G = nx.Graph()
for p in plen: G.add_node(p)
for a,b in RBH_pairs: G.add_edge(a,b)

reps, members = [], []
for comp in nx.connected_components(G):
    rep = max(comp, key=lambda x: plen[x])
    reps.append(rep)
    for m in comp: members.append((rep, m, plen[m]))
print(f"Clusters / reps: {len(reps)}")

# Write outputs
reps_set = set(reps)
with open("$OUT/dereplicated.fna","w") as fout:
    for r in sorted(reps_set):
        fout.write(f">{r}\n")
        for chunk in seq[r]: fout.write(chunk+"\n")
with open("$OUT/drep_members.tsv","w") as fout:
    fout.write("rep\tmember\tlength\n")
    for r,m,l in members: fout.write(f"{r}\t{m}\t{l}\n")
print("Wrote dereplicated.fna + drep_members.tsv")
PYEOF
fi

# ---- Split dereplicated.fna by topology (circ / frag) ----
python3 - <<PYEOF
topo = {}
with open("$OUT/circ_frag_map.tsv") as f:
    for line in f:
        c,t = line.rstrip('\n').split('\t'); topo[c] = t
with open("$OUT/dereplicated.fna") as f, \
     open("$OUT/circ.fna","w") as oc, \
     open("$OUT/frag.fna","w") as ofg:
    out=None
    for line in f:
        if line.startswith('>'):
            cid = line[1:].split()[0]
            out = oc if topo.get(cid)=='circ' else ofg
        out.write(line)
PYEOF

echo "[$(date '+%F %T')] DONE — dereplicated: $(grep -c '^>' $OUT/dereplicated.fna) (circ=$(grep -c '^>' $OUT/circ.fna), frag=$(grep -c '^>' $OUT/frag.fna))"
