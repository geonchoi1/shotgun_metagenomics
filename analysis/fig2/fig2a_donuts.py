#!/usr/bin/env python
"""Fig2a — two donuts: (left) Complete circular vs Putative; (right) PLSDB hit vs No PLSDB hit.
% in white inside ring, n=1,869 in centre, capitalised leader labels. PLSDB hit = red."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A="/home/gchoi/wwtp_plasmidome/analysis/plasmid_catalog"
df=pd.read_csv(f"{A}/plasmid_master.tsv",sep="\t"); N=len(df)
nCirc=int((df.topology=="circular").sum()); nPut=N-nCirc
nHit=int(df.plsdb_match.sum());            nNo=N-nHit
BLUE="#2166ac"; GREY="#969696"; RED="#e31a1c"

def donut(ax,sizes,labels,colors):
    w,_=ax.pie(sizes,colors=colors,startangle=90,counterclock=False,
               wedgeprops=dict(width=0.42,edgecolor="black",lw=1.1))
    for we,sz in zip(w,sizes):
        ang=np.deg2rad((we.theta1+we.theta2)/2)
        ax.text(0.79*np.cos(ang),0.79*np.sin(ang),f"{100*sz/N:.1f}%",ha="center",va="center",
                fontsize=11,fontweight="bold",color="white")
    place=[(0.0,1.55,"bottom"),(0.0,-1.55,"top")]               # top / bottom leader labels
    for (we,lab,sz),(tx,ty,va) in zip(zip(w,labels,sizes),place):
        ang=np.deg2rad((we.theta1+we.theta2)/2); x,y=np.cos(ang),np.sin(ang)
        ax.annotate(f"{lab} ({sz})",xy=(x*0.98,y*0.98),xytext=(tx,ty),ha="center",va=va,
                    fontsize=10.5,fontweight="bold",
                    arrowprops=dict(arrowstyle="-",color="black",lw=0.8))
    ax.text(0,0,f"n={N:,}",ha="center",va="center",fontsize=13,fontweight="bold")
    ax.set_aspect("equal"); ax.set_xlim(-1.25,1.25); ax.set_ylim(-1.9,1.9)

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7,4.2))
donut(ax1,[nCirc,nPut],["Complete circular","Putative"],[BLUE,GREY])
donut(ax2,[nHit,nNo],["PLSDB hit","No PLSDB hit"],[RED,GREY])
fig.savefig(f"{A}/Fig2a_donut.png",dpi=600,bbox_inches="tight"); plt.close(fig)
print(f"saved Fig2a_donut.png | circular {nCirc}/{nPut} | PLSDB hit {nHit}(red)/{nNo}")
