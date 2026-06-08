#!/usr/bin/env python
"""Fig4a (alt) — ARG/MGE co-occurrence NETWORK across 1,873 plasmids.
Edge = significant positive co-occurrence (CooccurrenceAffinity / Mainali alpha_mle > 0, p<0.05);
edge width ~ affinity. Nodes coloured by element type (ARG/Metal/Integron/IS/ICE/IME). Modules emerge
as cliques: mer (Hg) operon, class-1 integron (qacEdelta1-sul1), sul2-IS91, ICE-IME."""
import pandas as pd, numpy as np, networkx as nx
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
A="/home/gchoi/wwtp_plasmidome/analysis/04_arg_mge"
S="/home/gchoi/wwtp_plasmidome/04_cross_track/01_mobile_arg/cooccur_test/affinity_prev5_all.tsv"
d=pd.read_csv(S,sep="\t")
d["a"]=pd.to_numeric(d["alpha_mle"],errors="coerce")
d=d[(d["p_value"]<0.05)&(d["a"]>0)].copy()
prev={}   # node prevalence = # plasmids carrying that element
for _,r in d.iterrows(): prev[r["entity_1"]]=r["entity_1_count_mA"]; prev[r["entity_2"]]=r["entity_2_count_mB"]
METAL=set("merA merB merD merE merF merP merR merT arsA arsB arsC copA pcoA silA cadA".split())
INTG={"qacEdelta1","sul1"}
def typ(e):
    p,n=e.split(":",1)
    if p=="IS": return "IS"
    if e=="MGE:ICE": return "ICE"
    if e=="MGE:IME": return "IME"
    if p=="INT": return "Integron"
    if n in INTG: return "Integron"
    if n in METAL: return "Metal"
    return "ARG"
TCOL={"ARG":"#d62728","Metal":"#8c564b","Integron":"#ff7f0e","IS":"#2166ac","ICE":"#1b7837","IME":"#5ab4ac"}
from scipy.spatial import cKDTree
G=nx.Graph()
for _,r in d.iterrows():
    if typ(r["entity_1"])=="Metal" or typ(r["entity_2"])=="Metal": continue   # drop metal/biocide
    G.add_edge(r["entity_1"],r["entity_2"],w=min(r["a"],10))
nodes=list(G.nodes())
pos=nx.spring_layout(G,seed=7,k=2.6/np.sqrt(len(G)),iterations=400)
P=np.array([pos[n] for n in nodes]); P-=P.mean(0); P/=np.abs(P).max()
R=0.090                                                    # min node separation -> de-overlap
for _ in range(500):
    pr=cKDTree(P).query_pairs(2*R,output_type='ndarray')
    if len(pr)==0: break
    dd=P[pr[:,1]]-P[pr[:,0]]; dist=np.hypot(dd[:,0],dd[:,1]); dist[dist<1e-9]=1e-9
    push=dd/dist[:,None]*((2*R-dist)/2)[:,None]
    disp=np.zeros_like(P); np.add.at(disp,pr[:,0],-push); np.add.at(disp,pr[:,1],push); P+=disp
pos={n:P[i] for i,n in enumerate(nodes)}
fig,ax=plt.subplots(figsize=(8.4,7.0))
for u,v,ed in G.edges(data=True):
    ax.plot([pos[u][0],pos[v][0]],[pos[u][1],pos[v][1]],color="0.62",lw=0.5+ed["w"]/10*2.8,alpha=0.5,zorder=1,solid_capstyle="round")
def nsize(n): return 140+min(prev.get(n,5),90)*9   # size ~ prevalence (# plasmids), capped
for n in nodes:
    t=typ(n)
    ax.scatter(*pos[n],s=nsize(n),color=TCOL[t],edgecolor="black",lw=1.0,zorder=3)
    lab=n.split(":")[1]; fs=6.8 if len(lab)<=4 else (5.6 if len(lab)<=7 else 4.6)   # shrink long labels
    ax.annotate(lab,pos[n],fontsize=fs,fontweight="bold",ha="center",va="center",color="black",zorder=4,
                path_effects=[pe.withStroke(linewidth=1.6,foreground="white")])   # black text + white halo
ax.axis("off"); ax.set_aspect("equal")
LEG=[("ARG","#d62728"),("Integron","#ff7f0e"),("IS","#2166ac"),("ICE","#1b7837"),("IME","#5ab4ac")]
hand=[Line2D([0],[0],marker="o",color="w",markerfacecolor=c,markeredgecolor="black",markersize=11,label=l) for l,c in LEG]
# node-size key (prevalence) + affinity width key
hand+=[Line2D([0],[0],marker="o",color="w",markerfacecolor="0.7",markeredgecolor="black",
       markersize=np.sqrt(140+min(p,90)*9)/2.3,label=f"{p} plasmids") for p in (10,50,90)]
hand+=[Line2D([0],[0],color="0.6",lw=0.4+a/10*2.6,label=f"α={a}") for a in (3,10)]
lg=ax.legend(handles=hand,loc="center left",bbox_to_anchor=(1.0,0.5),frameon=False,fontsize=9.5,handletextpad=0.5,labelspacing=0.9)
plt.setp(lg.get_texts(),fontweight="bold")
fig.savefig(f"{A}/Fig4a_network.png",dpi=600,bbox_inches="tight"); plt.close(fig)
from PIL import Image; im=Image.open(f"{A}/Fig4a_network.png"); w,h=im.size; im.resize((780,int(780*h/w))).save("/tmp/fig4a_net.png")
print(f"nodes={G.number_of_nodes()} edges={G.number_of_edges()}")
