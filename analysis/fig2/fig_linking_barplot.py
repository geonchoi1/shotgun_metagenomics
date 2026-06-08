#!/usr/bin/env python
"""Number of our-plasmid -- PLSDB links per GTDB phylum (He-2024 Fig2b style), vertical bars.
One y-axis break (Pseudomonadota above); double-slash cut marks; n on bars."""
import pandas as pd, numpy as np, json
from collections import Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A="/home/gchoi/wwtp_plasmidome/analysis/plasmid_catalog"
phy=json.load(open(f"{A}/plsdb_acc2gtdbphylum.json"))
def P(a): return phy.get(a) or phy.get(a.split(".")[0]) or "Unknown"
op=pd.read_csv(f"{A}/plsdb_network/ours_plsdb_edges.tsv",sep="\t",header=None)
c=Counter(P(b) for b in op[1])
ORDER=["Pseudomonadota","Bacteroidota","Bacillota","Actinomycetota","Campylobacterota","Desulfobacterota"]
vals=[c[p] for p in ORDER]
PCOL=dict(zip(ORDER,plt.cm.tab10(np.linspace(0,1,6))))
x=np.arange(len(ORDER)); colors=[PCOL[p] for p in ORDER]
TOP=(7120,7500); BOT=(0,1120)

fig,(axT,axB)=plt.subplots(2,1,sharex=True,figsize=(4.8,4.3),
                           gridspec_kw=dict(height_ratios=[0.7,3],hspace=0.10))
for ax in (axT,axB):
    ax.bar(x,vals,width=0.66,color=colors,edgecolor="black",lw=0.7)
    ax.spines["right"].set_visible(False)
axT.set_ylim(*TOP); axB.set_ylim(*BOT)
axT.spines["bottom"].set_visible(False); axT.spines["top"].set_visible(False)
axB.spines["top"].set_visible(False)
axT.tick_params(bottom=False)
axT.set_yticks([7300]); axB.set_yticks([0,500,1000])
# double-slash cut marks at left of the break
d=.6
kw=dict(marker=[(-1,-d),(1,d)],markersize=8,linestyle="none",color="k",mec="k",mew=1.1,clip_on=False)
for ox in (-0.013,0.013):
    axT.plot([ox],[0],transform=axT.transAxes,**kw)
    axB.plot([ox],[1],transform=axB.transAxes,**kw)
# n labels above each bar in its panel
for xi,(p,v) in enumerate(zip(ORDER,vals)):
    ax=axT if v>1120 else axB; rng=ax.get_ylim()[1]-ax.get_ylim()[0]
    ax.text(xi,v+rng*0.04,f"{v:,}",ha="center",va="bottom",fontsize=9.5,fontweight="bold")
axB.set_xticks(x); axB.set_xticklabels(ORDER,rotation=35,ha="right",fontsize=10,fontweight="bold")
for ax in (axT,axB):
    for lab in ax.get_yticklabels(): lab.set_fontweight("bold"); lab.set_fontsize(9)
axB.set_ylabel("Number of links",fontsize=11,fontweight="bold")
axB.yaxis.set_label_coords(-0.17,0.72)
fig.savefig(f"{A}/Fig2c_links.png",dpi=600,bbox_inches="tight"); plt.close(fig)
from PIL import Image; im=Image.open(f"{A}/Fig2c_links.png"); w,h=im.size; im.resize((720,int(720*h/w))).save("/tmp/barlink.png")
print("saved")
