#!/bin/bash
# === Track 1 + 2 + 3 clustering + PTU co-membership validation ===
# Track 1: BLAST-Leiden pOTU on our 1873 (Camargo / Fiamenghi)
# Track 2: COPLA on circular only → PTU
# Track 3: combined our + PLSDB ref → indirect PTU labels for our plasmids
# Validation: classify Track 3 clusters as Pure/Mixed/Novel
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

OURS=$PROJECT/plasmid/02_drep/dereplicated.fna
CIRC=$PROJECT/plasmid/02_drep/circ.fna
OUT=$PROJECT/plasmid/30_clustering
mkdir -p $OUT/track1 $OUT/track2_copla $OUT/track3 $OUT/validation

############ TRACK 1: BLAST + Leiden on ours ############
T1=$OUT/track1
if [ ! -s $T1/pOTU_membership.tsv ]; then
  echo "[$(date '+%F %T')] Track1 — all-vs-all megablast"
  activate_env "$ENV_DIAMOND"
  makeblastdb -in $OURS -dbtype nucl -out $T1/plasmids
  blastn -task megablast -query $OURS -db $T1/plasmids \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
    -evalue 1e-3 -num_threads $THREADS_BLAST \
    -max_target_seqs 25000 -perc_identity 0.0 \
    -out $T1/blastn.tsv

  python3 - <<PYEOF
import networkx as nx, igraph as ig, leidenalg
from collections import defaultdict
T1 = "$T1"; OURS="$OURS"
plen = {}; cur=None; ln=0
with open(OURS) as f:
    for line in f:
        if line.startswith('>'):
            if cur: plen[cur]=ln
            cur=line[1:].split()[0]; ln=0
        else: ln += len(line.strip())
plen[cur]=ln
pair = defaultdict(lambda: {'aln':0,'pid_sum':0.0})
with open(f"{T1}/blastn.tsv") as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        q,s=p[0],p[1]
        if q==s: continue
        pid=float(p[2]); aln=int(p[3])
        pair[(q,s)]['aln']+=aln; pair[(q,s)]['pid_sum']+=aln*pid
G=nx.Graph()
for k in plen: G.add_node(k)
for (q,s),d in pair.items():
    tANI=d['pid_sum']/d['aln'] if d['aln']>0 else 0
    af=100.0*d['aln']/plen[q]
    if af>=70.0 and tANI>=70.0:
        w = tANI*af/10000.0
        if G.has_edge(q,s):
            G[q][s]['weight']=max(G[q][s]['weight'],w)
        else:
            G.add_edge(q,s,weight=w)
nodes=list(G.nodes()); ni={n:i for i,n in enumerate(nodes)}
g=ig.Graph(n=len(nodes), edges=[(ni[a],ni[b]) for a,b in G.edges()],
           edge_attrs={'weight':[G[a][b]['weight'] for a,b in G.edges()]})
part=leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                              weights='weight', resolution_parameter=1.0, seed=42)
with open(f"{T1}/pOTU_membership.tsv","w") as o:
    o.write("contig\tpOTU\n")
    for i,c in enumerate(part):
        for nid in c:
            o.write(f"{nodes[nid]}\tpOTU_{i+1:05d}\n")
print(f"Track1 pOTU: {len(part)} clusters")
PYEOF
fi

############ TRACK 2: COPLA on circular only ############
T2=$OUT/track2_copla
if [ -s $CIRC ] && [ ! -s $T2/copla_summary.tsv ]; then
  echo "[$(date '+%F %T')] Track2 — COPLA on circular"
  activate_env "$ENV_COPLA"
  cd $COPLA_DIR
  python3 bin/copla.py $CIRC $T2 \
    --threads $THREADS || echo "WARN: COPLA failed — check inputs"
  cd $SCRIPT_DIR
fi

############ TRACK 3: combined (ours ∪ PLSDB ref) BLAST + Leiden ############
T3=$OUT/track3
if [ ! -s $T3/combined_clusters.tsv ]; then
  echo "[$(date '+%F %T')] Track3 — combined dereplicated + PLSDB ref"
  mkdir -p $T3
  cat $OURS $PLSDB_FASTA > $T3/combined.fna
  activate_env "$ENV_DIAMOND"
  makeblastdb -in $T3/combined.fna -dbtype nucl -out $T3/combined
  blastn -task megablast -query $T3/combined.fna -db $T3/combined \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
    -evalue 1e-3 -num_threads $THREADS_BLAST \
    -max_target_seqs 25000 -perc_identity 0.0 \
    -out $T3/combined_blast.tsv

  python3 - <<PYEOF
import networkx as nx, igraph as ig, leidenalg
from collections import defaultdict
T3="$T3"
plen={}; cur=None; ln=0
with open(f"{T3}/combined.fna") as f:
    for line in f:
        if line.startswith('>'):
            if cur: plen[cur]=ln
            cur=line[1:].split()[0]; ln=0
        else: ln+=len(line.strip())
plen[cur]=ln
pair=defaultdict(lambda:{'aln':0,'pid_sum':0.0})
with open(f"{T3}/combined_blast.tsv") as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        q,s=p[0],p[1]
        if q==s: continue
        pair[(q,s)]['aln']+=int(p[3]); pair[(q,s)]['pid_sum']+=int(p[3])*float(p[2])
G=nx.Graph()
for k in plen: G.add_node(k)
for (q,s),d in pair.items():
    tANI=d['pid_sum']/d['aln'] if d['aln']>0 else 0
    af=100.0*d['aln']/plen[q]
    if af>=70.0 and tANI>=70.0:
        w=tANI*af/10000.0
        if G.has_edge(q,s): G[q][s]['weight']=max(G[q][s]['weight'],w)
        else: G.add_edge(q,s,weight=w)
nodes=list(G.nodes()); ni={n:i for i,n in enumerate(nodes)}
g=ig.Graph(n=len(nodes), edges=[(ni[a],ni[b]) for a,b in G.edges()],
           edge_attrs={'weight':[G[a][b]['weight'] for a,b in G.edges()]})
part=leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                              weights='weight', resolution_parameter=1.0, seed=42)
with open(f"{T3}/combined_clusters.tsv","w") as o:
    o.write("contig\tpOTU_T3\n")
    for i,c in enumerate(part):
        for nid in c:
            o.write(f"{nodes[nid]}\tpOTU_T3_{i+1:05d}\n")
print(f"Track3 clusters: {len(part)}")
PYEOF
fi

############ VALIDATION — Pure / Mixed / Novel ############
echo "[$(date '+%F %T')] Validation — Pure/Mixed/Novel"
python3 - <<PYEOF
from collections import defaultdict
OURS_IDS=set()
with open("$OURS") as f:
    for line in f:
        if line.startswith('>'): OURS_IDS.add(line[1:].split()[0])
cl=defaultdict(list)
with open("$OUT/track3/combined_clusters.tsv") as f:
    next(f)
    for line in f:
        c,p=line.rstrip('\n').split('\t')
        cl[p].append(c)
out=open("$OUT/validation/cluster_classification.tsv","w")
out.write("pOTU_T3\tn_total\tn_ours\tn_ref\tclass\tindirect_ptu_ref_members\n")
ours_to_ref={}
for p,ms in cl.items():
    ours=[m for m in ms if m in OURS_IDS]
    ref=[m for m in ms if m not in OURS_IDS]
    if ours and ref: cls="Mixed"
    elif ours and not ref: cls="Novel"
    else: cls="PureRef"
    out.write(f"{p}\t{len(ms)}\t{len(ours)}\t{len(ref)}\t{cls}\t{';'.join(ref[:5])}\n")
    if cls=="Mixed":
        # Assign indirect-PTU = mode-ref-member name
        for o in ours:
            ours_to_ref[o]=p
out.close()
with open("$OUT/validation/our_plasmid_indirect_ptu.tsv","w") as o:
    o.write("contig\tpOTU_T3\n")
    for c,p in ours_to_ref.items():
        o.write(f"{c}\t{p}\n")
print("validation done")
PYEOF
echo "[$(date '+%F %T')] DONE"
