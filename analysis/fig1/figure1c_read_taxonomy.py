#!/usr/bin/env python
"""Read-based community composition (Kraken2/Bracken).
Fig 1d = Genus (wide, horizontal). Fig S2 = Phylum. Unified fonts: axis 12 / ticks 11 / legend 10 (bold)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
SAMP=["IN","Anaerobic","Anoxic","Oxic","RAS","EF"]
DISP={"IN":"Inf","Anaerobic":"Ana","Anoxic":"Anx","Oxic":"Oxi","RAS":"RAS","EF":"Eff"}
XL=[DISP[s] for s in SAMP]
EXCLUDE={"Chordata","Homo"}
def matrix(rank):
    m={}
    for s in SAMP:
        d={}
        for ln in open(f"/home/gchoi/wwtp_plasmidome/00_shared/09_kraken2_community/{s}/bracken.G.report"):
            p=ln.rstrip("\n").split("\t")
            if len(p)>=6 and p[3]==rank:
                name=p[5].strip()
                if name in EXCLUDE: continue
                d[name]=float(p[0])
        m[s]=d
    M=pd.DataFrame(m).fillna(0.0); return M/M.sum()*100
OUT="/home/gchoi/wwtp_plasmidome/analysis/taxonomy"
def stackbar(ax,M,ntop):
    top=M.mean(axis=1).sort_values(ascending=False).head(ntop).index.tolist()
    plot=M.loc[top].copy(); plot.loc["Other"]=M.drop(index=top).sum(); plot=plot[SAMP]
    cols=list(plt.cm.tab20(np.linspace(0,1,20)))[:len(plot)-1]+[(0.72,0.72,0.72,1)]
    base=np.zeros(len(SAMP))
    for (nm,row),c in zip(plot.iterrows(),cols):
        ax.barh(np.arange(len(SAMP)),row.values,left=base,label=nm,color=c,height=0.8,edgecolor="black",linewidth=0.5)
        base+=row.values
    ax.set_yticks(np.arange(len(SAMP))); ax.set_yticklabels(XL); ax.invert_yaxis()
    ax.set_xlabel("Read-based relative abundance (%)",fontsize=12,fontweight="bold"); ax.set_xlim(0,100)
    ax.tick_params(axis="both",labelsize=11)
    for lab in ax.get_xticklabels()+ax.get_yticklabels(): lab.set_fontweight("bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(bbox_to_anchor=(1.005,1),loc="upper left",fontsize=10,handlelength=1.1,labelspacing=0.3,frameon=False)
MP,MG=matrix("P"),matrix("G")
# Fig 1d — Genus (wide)
fig,ax=plt.subplots(figsize=(10,4.2)); stackbar(ax,MG,15)
fig.tight_layout(); fig.savefig(f"{OUT}/Fig1d_read_taxonomy.png",dpi=600,bbox_inches="tight"); plt.close(fig)
MG.to_csv(f"{OUT}/read_taxonomy_genus_matrix.tsv",sep="\t")
# Fig S2 — Phylum
fig,ax=plt.subplots(figsize=(10,3.6)); stackbar(ax,MP,12)
fig.tight_layout(); fig.savefig(f"{OUT}/FigureS2_phylum_taxonomy.png",dpi=600,bbox_inches="tight"); plt.close(fig)
MP.to_csv(f"{OUT}/read_taxonomy_phylum_matrix.tsv",sep="\t")
print("saved Fig1d_read_taxonomy.png (Genus) + FigureS2_phylum_taxonomy.png (Phylum)")
