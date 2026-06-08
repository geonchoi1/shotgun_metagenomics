#!/usr/bin/env python
"""Fig 1c — ARG donut (MAGs carrying ARG) + phylum x antibiotic drug-class stacked bar.
AMRFinder Type=AMR only. Fonts: axis 12 / ticks 11 / legend 10 (bold); donut center 16, %% 14, caption 12."""
import pandas as pd, numpy as np, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
A="/home/gchoi/wwtp_plasmidome/analysis/mag"
pm=pd.read_csv(f"{A}/inputs/arg_vfg_per_mag.tsv",sep="\t")
nA=int((pm.ARG>0).sum()); N=len(pm)
amr=pd.read_csv("/home/gchoi/wwtp_plasmidome/02_mag_track/08_amrfinder/amr.tsv",sep="\t")
cls=amr[amr["Type"]=="AMR"]["Class"].value_counts()
labels=[c.title() for c in cls.index]; pal=plt.cm.tab20(np.linspace(0,1,len(cls)))
clpal={labels[i]:pal[i] for i in range(len(labels))}

def arg_donut(ax):
    sizes=[nA,N-nA]; labels=["ARG-carrying MAGs","No ARG"]; colors=["#e31a1c","#969696"]
    w,_=ax.pie(sizes,colors=colors,startangle=90,counterclock=False,
               wedgeprops=dict(width=0.42,edgecolor="black",lw=1.1))
    for we,sz in zip(w,sizes):
        ang=np.deg2rad((we.theta1+we.theta2)/2)
        ax.text(0.79*np.cos(ang),0.79*np.sin(ang),f"{100*sz/N:.0f}%",ha="center",va="center",
                fontsize=11,fontweight="bold",color="white")
    # labels placed top / bottom (vertical) so the donut stays large in the narrow subplot
    place=[(0.55,1.55,"bottom"),(-0.55,-1.55,"top")]
    for (we,lab,sz),(tx,ty,va) in zip(zip(w,labels,sizes),place):
        ang=np.deg2rad((we.theta1+we.theta2)/2); x,y=np.cos(ang),np.sin(ang)
        ax.annotate(f"{lab} ({sz})",xy=(x*0.98,y*0.98),xytext=(tx,ty),
                    ha="center",va=va,fontsize=10.5,fontweight="bold",
                    arrowprops=dict(arrowstyle="-",color="black",lw=0.8))
    ax.text(0,0,f"n={N}",ha="center",va="center",fontsize=13,fontweight="bold")
    ax.set_aspect("equal"); ax.set_xlim(-1.25,1.25); ax.set_ylim(-1.9,1.9)

# phylum x drug-class (counts)
meta=pd.read_csv(f"{A}/inputs/tree_meta.tsv",sep="\t")
meta["phylum"]=meta["phylum"].map(lambda p:re.sub(r"_[A-Z]+$","",str(p)))
mag2phy=dict(zip(meta["label"],meta["phylum"]))
orf2mag={}
for ln in open("/home/gchoi/wwtp_plasmidome/02_mag_track/03_master_orf/mag.master.gff"):
    if "\tCDS\t" not in ln: continue
    pp=ln.split("\t"); m=re.search(r"ID=([^;]+)",pp[8])
    if m: orf2mag[m.group(1)]=pp[0].split("|")[0]
am=amr[amr["Type"]=="AMR"].copy()
am["phylum"]=am["Protein id"].map(orf2mag).map(mag2phy)
am["Class_t"]=am["Class"].str.title()
piv=am.dropna(subset=["phylum"]).pivot_table(index="phylum",columns="Class_t",aggfunc="size",fill_value=0)
piv=piv.loc[piv.sum(1).sort_values(ascending=True).index]
piv=piv.reindex(columns=[c for c in labels if c in piv.columns],fill_value=0)

fig,(axd,axb)=plt.subplots(1,2,figsize=(12,5),gridspec_kw=dict(width_ratios=[1,2]))
arg_donut(axd)
base=np.zeros(len(piv))
for c in piv.columns:
    axb.barh(range(len(piv)),piv[c].values,left=base,color=clpal[c],label=c,edgecolor="black",lw=0.5,height=0.8)
    base+=piv[c].values
axb.set_yticks(range(len(piv))); axb.set_yticklabels(piv.index,fontsize=11,fontweight="bold")
axb.set_xlabel("ARG genes (n)",fontsize=12,fontweight="bold"); axb.tick_params(labelsize=11)
for lab in axb.get_xticklabels(): lab.set_fontweight("bold")
axb.spines["top"].set_visible(False); axb.spines["right"].set_visible(False)
lg=axb.legend(bbox_to_anchor=(1.0,1.0),loc="upper left",fontsize=10,frameon=False,title="Drug class")
plt.setp(lg.get_texts(),fontweight="bold"); plt.setp(lg.get_title(),fontweight="bold")
fig.tight_layout(); fig.savefig(f"{A}/Fig1c.png",dpi=600,bbox_inches="tight"); plt.close(fig)
print(f"saved Fig1c | ARG {nA}/{N} ({100*nA/N:.0f}%), {len(cls)} classes, {int(cls.sum())} AMR genes")
