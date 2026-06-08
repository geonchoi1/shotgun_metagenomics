#!/usr/bin/env python
"""Fig2e/f networks with global force-directed (igraph) layout instead of phyllotaxis packing.
Official Leiden clusters; intra-cluster AF70 edges. Track1=compartment, Track3=ours black+ref phylum."""
import sys, pandas as pd, numpy as np, json
import igraph as ig
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
A="/home/gchoi/wwtp_plasmidome/analysis/plasmid_catalog"
C="/home/gchoi/wwtp_plasmidome/01_plasmid_track/30_clustering"
df=pd.read_csv(f"{A}/plasmid_master.tsv",sep="\t")
ZONES=["IN","Anaerobic","Anoxic","Oxic","RAS","EF"]
DISP={"IN":"Inf","Anaerobic":"Ana","Anoxic":"Anx","Oxic":"Oxi","RAS":"RAS","EF":"Eff"}
ZCOL={"IN":"#4C72B0","Anaerobic":"#55A868","Anoxic":"#C44E52","Oxic":"#8172B3","RAS":"#CCB974","EF":"#64B5CD"}
samp=dict(zip(df.contig,df["sample"]))

def load_clusters(path):
    cl={}; members={}
    for i,ln in enumerate(open(path)):
        mem=ln.rstrip("\n").split("\t")[1].split(","); members[i]=mem
        for m in mem: cl[m]=i
    return cl,members
def intra_edges(ani_path, cl, mnodes):
    a=pd.read_csv(ani_path,sep="\t",usecols=["contig_1","contig_2","ani","qcov","tcov"])
    a=a[(a.contig_1!=a.contig_2)&((a.qcov>=0.7)|(a.tcov>=0.7))]
    a=a[a.contig_1.isin(mnodes)&a.contig_2.isin(mnodes)]
    a=a[a.contig_1.map(cl)==a.contig_2.map(cl)]
    return list(zip(a.contig_1,a.contig_2))

def render(track):
    if track=="e":
        cl,members=load_clusters(f"{C}/track1/plasmid_1869_clusters.tsv")
        ani=f"{C}/track1/plasmid_1869_ani.tsv"
    else:
        cl,members=load_clusters(f"{C}/track3/combined_11763_clusters.tsv")
        ani=f"{C}/track3/combined_11763_ani.tsv"
    multi={c:m for c,m in members.items() if len(m)>1}
    mnodes=[n for m in multi.values() for n in m]
    edges=intra_edges(ani,cl,set(mnodes))
    nodes=sorted(set(mnodes))
    idx={n:i for i,n in enumerate(nodes)}
    g=ig.Graph(n=len(nodes),edges=[(idx[a],idx[b]) for a,b in edges],directed=False); g.simplify()
    lay=g.layout_fruchterman_reingold(niter=400)
    xy=np.array(lay.coords); xy=(xy-xy.mean(0))
    sc=np.percentile(np.abs(xy),99) or 1.0; xy=xy/sc
    fig,ax=plt.subplots(figsize=(9,9) if track=="e" else (11,11))
    ei=np.array([(e.source,e.target) for e in g.es])
    seg=np.stack([xy[ei[:,0]],xy[ei[:,1]]],axis=1)
    ax.add_collection(LineCollection(seg,colors="#9a9a9a",linewidths=0.5,alpha=0.7,zorder=1))
    is_our=lambda n:"|" in n
    if track=="e":
        for z in ZONES:
            m=[i for i,n in enumerate(nodes) if samp.get(n)==z]
            ax.scatter(xy[m,0],xy[m,1],s=22,color=ZCOL[z],edgecolor="white",linewidths=0.3,alpha=0.95,zorder=3,label=DISP[z])
        singl=[c for c,m in members.items() if len(m)==1]
        scoll=[members[c][0] for c in singl]
        np.random.seed(0); ang=np.random.uniform(0,2*np.pi,len(scoll)); R=1.12+np.random.uniform(0,0.18,len(scoll))
        ax.scatter(R*np.cos(ang),R*np.sin(ang),s=6,c=[ZCOL[samp.get(s,"IN")] for s in scoll],alpha=0.6,linewidths=0,zorder=2)
        hand=[Line2D([0],[0],marker="o",color="w",markerfacecolor=ZCOL[z],markersize=9,label=DISP[z]) for z in ZONES]
        ax.legend(handles=hand,fontsize=10,frameon=False,ncol=6,loc="lower center",bbox_to_anchor=(0.5,-0.02))
        out="Fig2e_potu_network.png"
    else:
        acc2phy=json.load(open(f"{A}/plsdb_acc2phylum.json"))
        TOP=["Pseudomonadota","Bacillota","Actinomycetota","Bacteroidota","Spirochaetota","Campylobacterota"]
        PCOL={p:plt.cm.tab10(i) for i,p in enumerate(TOP)}
        def phy(n):
            p=acc2phy.get(n) or acc2phy.get(n.split(".")[0]); return p if p in TOP else "Other"
        ref=[i for i,n in enumerate(nodes) if not is_our(n)]
        for ph in TOP+["Other"]:
            m=[i for i in ref if phy(nodes[i])==ph]
            ax.scatter(xy[m,0],xy[m,1],s=8,color=PCOL.get(ph,"#bcbcbc"),linewidths=0,alpha=0.8,zorder=2,label=ph)
        m=[i for i,n in enumerate(nodes) if is_our(n)]
        ax.scatter(xy[m,0],xy[m,1],s=14,color="black",linewidths=0,alpha=0.95,zorder=3)
        hand=[Line2D([0],[0],marker="o",color="w",markerfacecolor="black",markersize=10,label="our plasmids")]+[Line2D([0],[0],marker="o",color="w",markerfacecolor=PCOL.get(p,"#bcbcbc"),markersize=9,label=p) for p in TOP+["Other"]]
        ax.legend(handles=hand,fontsize=9,frameon=False,ncol=4,loc="lower center",bbox_to_anchor=(0.5,-0.04))
        out="Fig2f_track3_phylum.png"
    for t in ax.get_legend().get_texts(): t.set_fontweight("bold")
    ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(f"{A}/{out}",dpi=500,bbox_inches="tight"); plt.close(fig)
    print(f"saved {out} | nodes={len(nodes)} edges={g.ecount()}")
    from PIL import Image; im=Image.open(f"{A}/{out}"); w,h=im.size; im.resize((1000,int(1000*h/w))).save(f"/tmp/fig2{track}_ig.png")

render(sys.argv[1] if len(sys.argv)>1 else "e")
