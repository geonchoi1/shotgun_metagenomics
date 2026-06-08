#!/usr/bin/env python
"""Plasmid functional comparison (TPM-weighted abundance view) across the A2O train.
(a) PCA ordination of 6 zones on TPM-weighted KO abundance (StandardScaler + PCA2; coords+variance
    recomputed from the same env matrix that fed GSEA, so internally consistent). Train trajectory drawn.
(b) Zone-specific enriched KEGG pathways (gseapy prerank, env-vs-rest log2FC ranking). Dot = FDR<0.05 &
    NES>0; size=-log10(FDR), colour=NES. Only signal-bearing zones/pathways are named (no full matrix)."""
import pandas as pd, numpy as np, glob, os
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.lines import Line2D

F="/home/gchoi/wwtp_plasmidome/01_plasmid_track/32_functional_comparison"
GMT="/home/gchoi/wwtp_plasmidome/legacy/plasmidome_downstream_legacy/43_functional_comparison_1873_base/kegg_pathways.gmt"
OUT="/home/gchoi/wwtp_plasmidome/analysis/01_functional_pfam_ko"
ZONE=["IN","Anaerobic","Anoxic","Oxic","RAS","EF"]; DISP=["Inf","Ana","Anx","Oxi","RAS","Eff"]
ZCOL=dict(zip(ZONE,plt.cm.viridis(np.linspace(0.05,0.92,6))))   # train gradient Inf->Eff

# ---- (a) PCA on TPM-weighted KO matrix ----
M=pd.read_csv(f"{F}/abundance_weighted_ko_env.tsv",sep="\t",index_col=0).loc[ZONE]
X=StandardScaler().fit_transform(M.values)
pca=PCA(n_components=2); Y=pca.fit_transform(X); ve=pca.explained_variance_ratio_*100

# ---- (b) GSEA: gather per-env significant enriched pathways ----
id2name={}
for ln in open(GMT):
    p=ln.rstrip("\n").split("\t")
    if len(p)>=2: id2name[p[0]]=p[1]
rows=[]
for f in glob.glob(f"{F}/abundance_weighted_ko_gsea_*.tsv"):
    env=os.path.basename(f).replace("abundance_weighted_ko_gsea_","").replace(".tsv","")
    d=pd.read_csv(f,sep="\t")
    for _,r in d.iterrows():
        rows.append((env,r["Term"],float(r["NES"]),float(r["FDR q-val"])))
g=pd.DataFrame(rows,columns=["env","term","NES","FDR"])
sig=g[(g.FDR<0.05)&(g.NES>0)].copy()
sig=sig[sig.term.map(lambda t:t in id2name and not t.startswith("map05"))]   # drop human-disease/unnamed artifacts
terms=sorted(sig.term.unique())
# order pathways by the zone of their strongest enrichment (train order), then NES
best=sig.loc[sig.groupby("term").NES.idxmax()].set_index("term")
order=sorted(terms,key=lambda t:(ZONE.index(best.loc[t,"env"]),-best.loc[t,"NES"]))
SHORT={"Chlorocyclohexane and chlorobenzene degradation":"Chlorocyclohexane/-benzene degr.",
       "Microbial metabolism in diverse environments":"Microbial metab. (diverse env.)",
       "Degradation of aromatic compounds":"Aromatic compound degradation",
       "Ascorbate and aldarate metabolism":"Ascorbate/aldarate metabolism"}
lab={t:SHORT.get(id2name.get(t,t),id2name.get(t,t)) for t in order}

# ==== (a) PCA — separate figure (descriptive ordination, no connecting line, no title) ====
figa,ax=plt.subplots(figsize=(5.0,4.6))
ax.axhline(0,color="0.85",lw=0.8,zorder=0); ax.axvline(0,color="0.85",lw=0.8,zorder=0)
for i,z in enumerate(ZONE):
    ax.scatter(Y[i,0],Y[i,1],s=340,color=ZCOL[z],edgecolor="black",lw=1.4,zorder=3)
    ax.annotate(DISP[i],(Y[i,0],Y[i,1]),fontsize=11,fontweight="bold",ha="center",va="center",
                color="white" if i<3 else "black",zorder=4)
ax.set_xlabel(f"PC1 ({ve[0]:.0f}%)",fontsize=12,fontweight="bold")
ax.set_ylabel(f"PC2 ({ve[1]:.0f}%)",fontsize=12,fontweight="bold")
ax.tick_params(labelsize=9)
for s in ["top","right"]: ax.spines[s].set_visible(False)
figa.savefig(f"{OUT}/Fig_functional_pca.png",dpi=600,bbox_inches="tight"); plt.close(figa)

# ==== (b) GSEA bubble — separate figure, portrait (taller than wide), top/right spines removed, all-bold ====
figb,ax2=plt.subplots(figsize=(5.2,7.2))
nes=sig.set_index(["term","env"]).NES.to_dict(); fdr=sig.set_index(["term","env"]).FDR.to_dict()
norm=mc.Normalize(sig.NES.min(),sig.NES.max()); cmap=plt.cm.Reds
for yi,t in enumerate(order):
    for xi,z in enumerate(ZONE):
        if (t,z) in nes:
            sz=(-np.log10(max(fdr[(t,z)],1e-3)))*80+35
            ax2.scatter(xi,yi,s=sz,color=cmap(norm(nes[(t,z)])),edgecolor="black",lw=0.9,zorder=3)
ax2.set_xticks(range(6)); ax2.set_xticklabels(DISP,fontsize=11,fontweight="bold")
ax2.set_yticks(range(len(order))); ax2.set_yticklabels([lab[t] for t in order],fontsize=10,fontweight="bold")
ax2.set_ylim(-0.6,len(order)-0.4); ax2.set_xlim(-0.6,5.6); ax2.invert_yaxis()
ax2.grid(True,color="0.9",lw=0.6,zorder=0); ax2.set_axisbelow(True)
for s in ["top","right"]: ax2.spines[s].set_visible(False)
cb=figb.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax2,shrink=0.45,pad=0.02,aspect=14)
cb.set_label("NES",fontsize=11,fontweight="bold"); cb.ax.tick_params(labelsize=9)
for t in cb.ax.get_yticklabels(): t.set_fontweight("bold")
sl=[Line2D([0],[0],marker="o",color="w",markerfacecolor="0.6",markeredgecolor="black",
           markersize=np.sqrt((-np.log10(q))*80+35),label=f"{q:g}") for q in (0.05,0.01,0.001)]
lg=ax2.legend(handles=sl,title="FDR",fontsize=9,title_fontsize=10,frameon=False,
              loc="upper center",bbox_to_anchor=(0.5,-0.07),ncol=3,columnspacing=2.6,handletextpad=0.5)
plt.setp(lg.get_title(),fontweight="bold"); plt.setp(lg.get_texts(),fontweight="bold")
figb.savefig(f"{OUT}/Fig_functional_gsea.png",dpi=600,bbox_inches="tight"); plt.close(figb)

from PIL import Image
for n in ("Fig_functional_pca","Fig_functional_gsea"):
    im=Image.open(f"{OUT}/{n}.png"); w,h=im.size; im.resize((620,int(620*h/w))).save(f"/tmp/{n}.png")
print(f"PCA var: PC1={ve[0]:.1f}% PC2={ve[1]:.1f}% | sig pathways={len(order)}")
print("pathways:",[f"{t}:{lab[t]}" for t in order])
