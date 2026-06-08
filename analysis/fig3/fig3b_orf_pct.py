#!/usr/bin/env python
"""Fig3b — cargo composition (ARG / metal-biocide / virulence) of ORFs, plasmid vs chromosomal.
Stacked bar normalized to 100% (relative composition); % labels on segments; monochrome blue."""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A="/home/gchoi/wwtp_plasmidome/analysis/plasmid_cargo"
CATS=["ARG","Metal / biocide resistance","Virulence (VF)"]
COL=["#08519c","#4292c6","#c6dbef"]                       # monochrome blue (dark->light)
# ORF counts: plasmid, chromosomal(MAG+unbinned)
PL=np.array([186,89,123]); CH=np.array([866,323,7113])
plp=PL/PL.sum()*100; chp=CH/CH.sum()*100
groups=["Plasmid","Chromosomal"]; data=np.array([plp,chp]); y=np.arange(2)
fig,ax=plt.subplots(figsize=(5.6,3.0))
base=np.zeros(2)
for j,c in enumerate(CATS):
    ax.barh(y,data[:,j],left=base,height=0.6,color=COL[j],edgecolor="black",lw=0.8,label=c)
    for k in range(2):
        if data[k,j]>=4:
            ax.text(base[k]+data[k,j]/2,k,f"{data[k,j]:.0f}%",ha="center",va="center",
                    fontsize=10,fontweight="bold",color=("white" if j==0 else "black"))
    base+=data[:,j]
ax.set_yticks(y); ax.set_yticklabels([f"{g}\n(n={n:,})" for g,n in zip(groups,[PL.sum(),CH.sum()])],fontsize=11,fontweight="bold")
ax.set_xlabel("Cargo ORF composition (%)",fontsize=12,fontweight="bold"); ax.set_xlim(0,100); ax.tick_params(labelsize=10)
for l in ax.get_xticklabels(): l.set_fontweight("bold")
ax.invert_yaxis()
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
lg=ax.legend(fontsize=9,frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.22),ncol=3,
             columnspacing=1.4,handletextpad=0.5)
plt.setp(lg.get_texts(),fontweight="bold")
fig.savefig(f"{A}/Fig3b_orf_pct.png",dpi=600,bbox_inches="tight"); plt.close(fig)
from PIL import Image; im=Image.open(f"{A}/Fig3b_orf_pct.png"); w,h=im.size; im.resize((520,int(520*h/w))).save("/tmp/fig3b_v2.png")
print("plasmid comp%:",plp.round(1),"chrom comp%:",chp.round(1))
