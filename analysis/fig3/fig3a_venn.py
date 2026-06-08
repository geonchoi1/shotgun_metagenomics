#!/usr/bin/env python
"""Fig3a — shared ORFs (diamond blastp >=90% id, >=90% q&s cov), plasmid vs chromosomal.
Totals inside circles; shared pulled out with a leader line. Blue (chromosomal) larger, shifted right."""
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
A="/home/gchoi/wwtp_plasmidome/analysis/plasmid_cargo"
PLT,CHT,SH=60798,4383003,24368        # plasmid total, chromosomal total, shared (plasmid ORFs w/ >=90% homolog)
RED="#e41a1c"; BLUE="#377eb8"
fig,ax=plt.subplots(figsize=(6.6,4.8))
ax.add_patch(Circle((0.34,0.50),0.24,facecolor=RED,alpha=0.45,edgecolor="black",lw=1.6,zorder=2))
ax.add_patch(Circle((0.62,0.50),0.34,facecolor=BLUE,alpha=0.45,edgecolor="black",lw=1.6,zorder=1))
# totals inside circles
ax.text(0.205,0.50,f"Plasmid\n{PLT:,}",ha="center",va="center",fontsize=12,fontweight="bold",color="black",zorder=4)
ax.text(0.74,0.50,f"Chromosomal\n{CHT:,}",ha="center",va="center",fontsize=12,fontweight="bold",color="black",zorder=4)
# shared pulled out with leader line to the overlap region
ax.annotate(f"Shared\n{SH:,}",xy=(0.44,0.50),xytext=(0.44,0.97),ha="center",va="center",
            fontsize=12,fontweight="bold",color="black",zorder=5,
            arrowprops=dict(arrowstyle="-",color="black",lw=1.1))
ax.set_xlim(0,1); ax.set_ylim(0.10,1.04); ax.set_aspect("equal"); ax.axis("off")
fig.savefig(f"{A}/Fig3a_venn.png",dpi=600,bbox_inches="tight"); plt.close(fig)
from PIL import Image; im=Image.open(f"{A}/Fig3a_venn.png"); w,h=im.size; im.resize((600,int(600*h/w))).save("/tmp/fig3a_v5.png")
print("done")
