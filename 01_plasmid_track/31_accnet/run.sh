#!/bin/bash
# === AcCNET protein-community via kClust + bipartite + Louvain ===
# Internal (our ORFs only) and External (our ∪ PLSDB ref ORFs)
# NMI vs Track 1 (internal) and Track 3 (external)
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/plasmid/04_master_orf/all/master.faa
ORF2CONTIG=$PROJECT/plasmid/04_master_orf/orf2contig.tsv
T1=$PROJECT/plasmid/30_clustering/track1/pOTU_membership.tsv
T3=$PROJECT/plasmid/30_clustering/track3/combined_clusters.tsv
OUT=$PROJECT/plasmid/31_accnet
mkdir -p $OUT/internal $OUT/external

KCLUST=${KCLUST_BIN:-$HOME/cfp/tools/kClust/kClust}
[ -x "$KCLUST" ] || KCLUST=$(command -v kClust || true)
[ -n "$KCLUST" ] || { echo "ERROR: kClust binary not found"; exit 1; }

############ INTERNAL ############
IN=$OUT/internal
if [ ! -s $IN/clusters.dmp ]; then
  echo "[$(date '+%F %T')] kClust internal"
  $KCLUST -i $FAA -d $IN -s 5.03 -c 0.8 -e 1e-4 -M $MEM_KCLUST 2>&1 | tail -10
fi

############ EXTERNAL: prodigal PLSDB ref + combine ############
EX=$OUT/external
if [ ! -s $EX/combined_master.faa ]; then
  echo "[$(date '+%F %T')] Prodigal on PLSDB ref"
  activate_env "$ENV_DIAMOND"  # base has prodigal usually
  command -v prodigal >/dev/null || { echo "ERROR: prodigal missing"; exit 1; }
  prodigal -i $PLSDB_FASTA -a $EX/reference_orfs.faa -o $EX/reference_orfs.gff \
           -p meta -q -f gff > $EX/prodigal.log 2>&1
  cat $FAA $EX/reference_orfs.faa > $EX/combined_master.faa
fi
if [ ! -s $EX/clusters.dmp ]; then
  echo "[$(date '+%F %T')] kClust external"
  $KCLUST -i $EX/combined_master.faa -d $EX -s 5.03 -c 0.8 -e 1e-4 -M $MEM_KCLUST 2>&1 | tail -10
fi

############ Bipartite + Louvain + NMI ############
python3 - <<PYEOF
import os, re
from collections import defaultdict
import networkx as nx
try:
    import community as community_louvain  # python-louvain
except ImportError:
    import networkx.algorithms.community as nxc
    community_louvain = None
from sklearn.metrics import normalized_mutual_info_score as NMI

def parse_kclust(out_dir):
    id2orf={}
    with open(f"{out_dir}/headers.dmp") as f:
        for line in f:
            m=re.match(r'^(\d+)\s+>(\S+)', line)
            if m: id2orf[int(m.group(1))]=m.group(2)
    orf2c={}
    with open(f"{out_dir}/clusters.dmp") as f:
        for line in f:
            if line.startswith('#'): continue
            p=line.split()
            if len(p)<2: continue
            nid=int(p[0]); cid=int(p[1])
            if nid in id2orf: orf2c[id2orf[nid]]=cid
    return orf2c

def orf2contig_map(path):
    o2c={}
    with open(path) as f:
        for line in f:
            p=line.rstrip('\n').split('\t')
            if len(p)>=2: o2c[p[0]]=p[1]
    return o2c

def ref_orf2contig(gff):
    o2c={}
    cur=None
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                m=re.search(r'seqhdr="([^"]+)"', line)
                if m: cur=m.group(1).split()[0]
                continue
            p=line.rstrip('\n').split('\t')
            if len(p)<9: continue
            attrs=dict(kv.split('=',1) for kv in p[8].split(';') if '=' in kv)
            oid=attrs.get('ID','')
            if cur and oid:
                # Prodigal -p meta IDs are like "1_1","1_2" — match the actual FAA header which is "<contig>_<n>"
                # Build candidate: prodigal FAA header is "<seqhdr>_<n>"
                num=oid.split('_')[-1]
                o2c[f"{cur}_{num}"]=cur
    return o2c

def louvain_clusters(G):
    if community_louvain is not None:
        part=community_louvain.best_partition(G, weight='weight', random_state=42)
        cl=defaultdict(list)
        for n,c in part.items(): cl[c].append(n)
        return list(cl.values())
    else:
        # Fallback: greedy modularity
        return list(nxc.greedy_modularity_communities(G, weight='weight'))

def build_bipartite_and_communities(orf2cluster, orf2contig):
    p2cl=defaultdict(set)
    for orf,cl in orf2cluster.items():
        c=orf2contig.get(orf)
        if c: p2cl[c].add(cl)
    # Plasmid graph: edge if share ≥1 cluster, weight = Jaccard
    plasmids=list(p2cl.keys())
    G=nx.Graph()
    for p in plasmids: G.add_node(p)
    # Inverted index: cluster→plasmids
    cl2p=defaultdict(set)
    for p,cls in p2cl.items():
        for c in cls: cl2p[c].add(p)
    pair_share=defaultdict(int)
    for c,ps in cl2p.items():
        ps=list(ps)
        for i in range(len(ps)):
            for j in range(i+1,len(ps)):
                a,b=sorted((ps[i],ps[j]))
                pair_share[(a,b)]+=1
    for (a,b),s in pair_share.items():
        union=len(p2cl[a] | p2cl[b])
        w=s/union if union else 0
        if w>0:
            G.add_edge(a,b,weight=w)
    return G

def nmi_vs(comm, ref_tsv, label_col=1):
    ref={}
    with open(ref_tsv) as f:
        h=next(f)
        for line in f:
            p=line.rstrip('\n').split('\t')
            if len(p)>label_col: ref[p[0]]=p[label_col]
    # Build comm dict
    cmap={}
    for ci,memb in enumerate(comm):
        for m in memb: cmap[m]=str(ci)
    common=sorted(set(cmap) & set(ref))
    if len(common)<2: return float('nan'), 0
    y1=[cmap[m] for m in common]; y2=[ref[m] for m in common]
    return NMI(y1,y2), len(common)

# ---- Internal ----
IN="$IN"
orf2c_int=parse_kclust(IN)
orf2contig=orf2contig_map("$ORF2CONTIG")
G_int=build_bipartite_and_communities(orf2c_int, orf2contig)
print(f"Internal graph: {G_int.number_of_nodes()} nodes, {G_int.number_of_edges()} edges")
comm_int=louvain_clusters(G_int)
with open(f"{IN}/plasmid_communities.tsv","w") as o:
    o.write("contig\tcommunity\n")
    for ci,m in enumerate(comm_int):
        for x in m: o.write(f"{x}\tcomm_{ci+1:05d}\n")
nmi1,n1 = nmi_vs(comm_int, "$T1")
with open(f"{IN}/nmi_vs_track1.txt","w") as o:
    o.write(f"NMI\t{nmi1:.4f}\nN_common\t{n1}\n")
print(f"NMI(internal vs Track1) = {nmi1:.4f} on n={n1}")

# ---- External ----
EX="$EX"
orf2c_ext=parse_kclust(EX)
orf2contig_ext=dict(orf2contig)
orf2contig_ext.update(ref_orf2contig(f"{EX}/reference_orfs.gff"))
G_ext=build_bipartite_and_communities(orf2c_ext, orf2contig_ext)
print(f"External graph: {G_ext.number_of_nodes()} nodes, {G_ext.number_of_edges()} edges")
comm_ext=louvain_clusters(G_ext)
with open(f"{EX}/plasmid_communities.tsv","w") as o:
    o.write("contig\tcommunity\n")
    for ci,m in enumerate(comm_ext):
        for x in m: o.write(f"{x}\tcomm_{ci+1:05d}\n")
nmi3,n3 = nmi_vs(comm_ext, "$T3")
with open(f"{EX}/nmi_vs_track3.txt","w") as o:
    o.write(f"NMI\t{nmi3:.4f}\nN_common\t{n3}\n")
print(f"NMI(external vs Track3) = {nmi3:.4f} on n={n3}")
PYEOF

echo "[$(date '+%F %T')] DONE"
