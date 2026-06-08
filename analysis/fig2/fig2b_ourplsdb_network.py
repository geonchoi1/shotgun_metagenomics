#!/usr/bin/env python
"""our-plasmid -- PLSDB network (our-PLSDB edges only), GTDB phylum colour, ours=black.
Node overlap removed (KDTree push-apart). Legend on the right (shorter figure, larger text)."""
import networkx as nx, numpy as np, json, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from scipy.spatial import cKDTree
A="/home/gchoi/wwtp_plasmidome/analysis/plasmid_catalog"; N=f"{A}/plsdb_network"
phy=json.load(open(f"{A}/plsdb_acc2gtdbphylum.json"))
def P(a): return phy.get(a) or phy.get(a.split(".")[0]) or "Unknown"
isour=lambda s:"|" in s
op=pd.read_csv(f"{N}/ours_plsdb_edges.tsv",sep="\t",header=None)
G=nx.Graph()
for a,b in zip(op[0],op[1]): G.add_edge(a,b,phylum=P(b))
nodes=list(G)
PHY=["Pseudomonadota","Bacteroidota","Bacillota","Actinomycetota","Campylobacterota","Desulfobacterota","Unknown"]
PCOL=dict(zip(PHY,list(plt.cm.tab10(np.linspace(0,1,6)))+[(0.7,0.7,0.7,1)]))
def pcol(p): return PCOL.get(p,(0.7,0.7,0.7,1))
pos=nx.spring_layout(G,seed=1,k=2.5/np.sqrt(len(nodes)),iterations=60)
Pp=np.array([pos[n] for n in nodes]); Pp-=Pp.mean(0); Pp/=np.abs(Pp).max()
R=0.010
for it in range(600):
    pr=cKDTree(Pp).query_pairs(2*R,output_type='ndarray')
    if len(pr)==0: break
    d=Pp[pr[:,1]]-Pp[pr[:,0]]; dist=np.hypot(d[:,0],d[:,1]); dist[dist<1e-9]=1e-9
    push=d/dist[:,None]*((2*R-dist)/2)[:,None]
    disp=np.zeros_like(Pp); np.add.at(disp,pr[:,0],-push); np.add.at(disp,pr[:,1],push); Pp+=disp
pos={n:Pp[i] for i,n in enumerate(nodes)}
xr=Pp[:,0].max()-Pp[:,0].min(); FW=16
ppu=FW*72/(xr*1.05); d_pts=2*R*ppu*0.42; s=np.pi*(d_pts/2)**2
fig,ax=plt.subplots(figsize=(FW,FW))
segs=[[pos[u],pos[v]] for u,v in G.edges()]; cols=[pcol(d["phylum"]) for _,_,d in G.edges(data=True)]
ax.add_collection(LineCollection(segs,colors=cols,linewidths=0.5,alpha=0.55,zorder=1))
ref=[n for n in nodes if not isour(n)]
for ph in PHY:
    ns=[n for n in ref if P(n)==ph]
    if ns: ax.scatter([pos[n][0] for n in ns],[pos[n][1] for n in ns],s=s,color=pcol(ph),linewidths=0,zorder=2)
ours=[n for n in nodes if isour(n)]
ax.scatter([pos[n][0] for n in ours],[pos[n][1] for n in ours],s=s*1.8,color="black",linewidths=0,zorder=3)
ax.axis("off"); ax.set_aspect("equal")
hand=[Line2D([0],[0],marker="o",color="w",markerfacecolor="black",markersize=15,label="This study")]+\
     [Line2D([0],[0],marker="o",color="w",markerfacecolor=pcol(p),markersize=14,label=p) for p in PHY if any(P(n)==p for n in ref)]
lg=ax.legend(handles=hand,fontsize=18,frameon=False,loc="center left",bbox_to_anchor=(1.0,0.5),
             handletextpad=0.4,labelspacing=0.7); plt.setp(lg.get_texts(),fontweight="bold")
fig.savefig(f"{A}/Fig2b_network.png",dpi=300,bbox_inches="tight"); plt.close()
from PIL import Image; im=Image.open(f"{A}/Fig2b_network.png"); w,h=im.size; im.resize((1250,int(1250*h/w))).save("/tmp/nov4.png")
print(f"saved | {len(nodes)} nodes ({len(ours)} ours), aspect w/h={w/h:.2f}")
