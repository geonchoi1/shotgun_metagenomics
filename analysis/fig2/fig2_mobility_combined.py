#!/usr/bin/env python
"""Fig2 mobility combined: (left) 5-tier stacked bar Putative vs Complete circular;
(right) tier x plasmid-length boxplot with pairwise MWU+BH brackets. Shared 5-tier colour key
(boxplot x-axis names every tier) -> no separate legend."""
import pandas as pd, numpy as np
import scipy.stats as st
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A="/home/gchoi/wwtp_plasmidome/analysis/plasmid_catalog"
df=pd.read_csv(f"{A}/plasmid_master.tsv",sep="\t"); circ=df[df.topology=="circular"]
TIERS=["pCONJ","pdCONJ","pMOB","pOriT","pNT"]
TCOL={"pCONJ":"#08519c","pdCONJ":"#3182bd","pMOB":"#6baed6","pOriT":"#bdd7e7","pNT":"#d9d9d9"}
B=dict(fontweight="bold")
fig,(axL,axR)=plt.subplots(1,2,figsize=(10.5,5.0),gridspec_kw=dict(width_ratios=[1,1.85],wspace=0.28))

# ---- LEFT: stacked bar ----
def tier_n(sub): c=sub.tier.value_counts(); return [c.get(t,0) for t in TIERS]
groups=[(f"Putative\n(n={len(df):,})",df),(f"Complete\ncircular (n={len(circ)})",circ)]
x=np.arange(len(groups)); base=np.zeros(len(groups))
counts=np.array([tier_n(g[1]) for g in groups],dtype=float); frac=counts/counts.sum(1,keepdims=True)*100
for j,t in enumerate(TIERS):
    axL.bar(x,frac[:,j],bottom=base,width=0.62,color=TCOL[t],edgecolor="black",lw=0.6)
    for k in range(len(groups)):
        if frac[k,j]>=4.5:
            axL.text(k,base[k]+frac[k,j]/2,f"{frac[k,j]:.0f}%",ha="center",va="center",
                     fontsize=9,color=("white" if t in("pCONJ","pdCONJ") else "black"),**B)
    base+=frac[:,j]
for k,(lab,sub) in enumerate(groups):
    axL.text(k,103,f"{100*sub.mobile.mean():.1f}%\nmobile",ha="center",fontsize=9.5,color="#08519c",**B)
axL.set_xticks(x); axL.set_xticklabels([g[0] for g in groups],fontsize=10,**B)
axL.set_ylim(0,116); axL.set_ylabel("Plasmids (%)",fontsize=12,**B); axL.tick_params(labelsize=10)
for l in axL.get_yticklabels(): l.set_fontweight("bold")
axL.spines["top"].set_visible(False); axL.spines["right"].set_visible(False)

# ---- RIGHT: boxplot tier x length ----
data=[df[df.tier==t].length.values/1000 for t in TIERS]            # kb
bp=axR.boxplot(data,positions=range(len(TIERS)),widths=0.6,patch_artist=True,showfliers=False,
               medianprops=dict(color="black",lw=1.4),whiskerprops=dict(color="black"),
               capprops=dict(color="black"),boxprops=dict(edgecolor="black",lw=0.8))
for patch,t in zip(bp["boxes"],TIERS): patch.set_facecolor(TCOL[t])
np.random.seed(0)
for i,d in enumerate(data):
    axR.scatter(np.random.normal(i,0.07,len(d)),d,s=4,color="black",alpha=0.18,zorder=1,linewidths=0)
axR.set_yscale("log"); axR.set_ylabel("Plasmid length (kb)",fontsize=12,**B)
axR.set_xticks(range(len(TIERS)))
axR.set_xticklabels([f"{t}\n(n={len(data[i])})" for i,t in enumerate(TIERS)],fontsize=10,**B)
axR.tick_params(labelsize=10)
for l in axR.get_yticklabels(): l.set_fontweight("bold")
axR.spines["top"].set_visible(False); axR.spines["right"].set_visible(False)
# pairwise MWU + BH
pairs=[(i,j) for i in range(5) for j in range(i+1,5)]
pv=[st.mannwhitneyu(data[i],data[j],alternative="two-sided").pvalue for i,j in pairs]
rej,padj,_,_=multipletests(pv,method="fdr_bh")
def star(p): return "***" if p<1e-3 else ("**" if p<1e-2 else ("*" if p<5e-2 else None))
sig=[(i,j,star(padj[k])) for k,(i,j) in enumerate(pairs) if rej[k] and star(padj[k])]
# greedy non-overlapping levels
placed=[]
for i,j,s in sorted(sig,key=lambda p:abs(p[1]-p[0])):
    lvl=0
    while any(l==lvl and not(j<=a or i>=b) for a,b,l in placed): lvl+=1
    placed.append((i,j,lvl))
starmap={(i,j):s for i,j,s in sig}
ymax=max(d.max() for d in data); STEP=1.28
for i,j,lvl in placed:
    s=starmap[(i,j)]
    yy=ymax*(STEP**(lvl+1))
    axR.plot([i,i,j,j],[yy/1.05,yy,yy,yy/1.05],lw=1.0,color="black")
    axR.text((i+j)/2,yy*1.01,s,ha="center",va="bottom",fontsize=9,**B)
axR.set_ylim(top=ymax*(STEP**(max(l for *_,l in placed)+2.4)))

# ---- shared (unified) 5-tier legend ----
from matplotlib.patches import Patch
hand=[Patch(facecolor=TCOL[t],edgecolor="black",lw=0.6,label=t) for t in TIERS]
lg=fig.legend(handles=hand,ncol=5,loc="lower center",bbox_to_anchor=(0.5,-0.04),frameon=False,
              fontsize=11,columnspacing=1.4,handletextpad=0.5)
plt.setp(lg.get_texts(),fontweight="bold")
fig.savefig(f"{A}/Fig2d_mobility.png",dpi=600,bbox_inches="tight"); plt.close(fig)
from PIL import Image; im=Image.open(f"{A}/Fig2d_mobility.png"); w,h=im.size; im.resize((1100,int(1100*h/w))).save("/tmp/mobcomb.png")
print(f"saved | sig pairs {len(sig)}")
