#!/usr/bin/env python
"""Fig3 (combined) — three plasmid heatmaps in ONE row, all '% of plasmid pool', transposed
(6 samples Inf..Eff on y, function on x):
  (1) Plasmid abundance : mobility tier, CONTIG TPM (per-plasmid molecule abundance)
  (2) Cargo             : ARG/Metal/VF, ORF TPM (per-gene abundance; backbone not over-counted)
  (3) Drug class        : AMRFinder Class, ORF TPM
Panel 1 normalized by total plasmid CONTIG TPM; panels 2-3 by total plasmid ORF mean-depth (D_s cancels).
in-cell %, black borders, per-panel LogNorm (Reds), no colorbars."""
import pandas as pd, numpy as np, re
from collections import Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mc
W="/home/gchoi/wwtp_plasmidome"; A=f"{W}/analysis/plasmid_cargo"; PT=f"{W}/01_plasmid_track"; BC=f"{W}/analysis/mag/cov/bedcov"
SAMP=["IN","Anaerobic","Anoxic","Oxic","RAS","EF"]; DISP=["Inf","Ana","Anx","Oxi","RAS","Eff"]; TPM=[f"{z}_TPM" for z in SAMP]
TIERS=["pCONJ","pdCONJ","pMOB","pOriT","pNT"]
cg=pd.read_csv(f"{A}/plasmid_cargo.tsv",sep="\t")

# ---- (1) mobility: sample x tier, CONTIG TPM, % of plasmid contig pool ----
MM=np.array([[cg[cg.tier==t][s].sum() for t in TIERS] for s in TPM])
totc=np.array([cg[s].sum() for s in TPM])            # total plasmid contig TPM per sample
MM=MM/totc[:,None]*100

# ---- plasmid ORF set + cargo ORF sets (plasmid only) ----
porf=set()
for ln in open(f"{PT}/04_master_orf/plasmidome.master.gff"):
    if "\tCDS\t" not in ln: continue
    m=re.search(r"ID=([^;]+)",ln.split("\t")[8])
    if m: porf.add(m.group(1))
amr=pd.read_csv(f"{PT}/09_amrfinder/amr.tsv",sep="\t")
arg=amr[amr.Type=="AMR"]
arg_orfs=set(arg["Protein id"])&porf
stress=set(amr[amr.Type=="STRESS"]["Protein id"])
bm=set(pd.read_csv(f"{PT}/10_bacmet/bacmet.tsv",sep="\t",header=None)[0])
metal_orfs=(bm|stress)&porf                                    # metal = BacMet U STRESS (matches has_metal)
vf_orfs=set(pd.read_csv(f"{PT}/11_vfdb/vfdb.tsv",sep="\t",header=None)[0])&porf
# ARG ORF -> drug class(es)
o2cls={}
for o,c in zip(arg["Protein id"],arg["Class"]):
    if o in porf and isinstance(c,str): o2cls.setdefault(o,set()).add(c.title())
cnt=Counter(c for cs in o2cls.values() for c in cs)
CL=[c for c,_ in cnt.most_common(12)]
allcargo=arg_orfs|metal_orfs|vf_orfs|set(o2cls)

# ---- per-sample bedcov: plasmid ORF pool + cargo/class ORF mean-depth ----
MC=[]; MD=[]
for s in SAMP:
    pool=0.0; dep={}
    for ln in open(f"{BC}/{s}.bedcov"):
        p=ln.split("\t"); orf=p[3]
        if orf not in porf: continue
        L=int(p[2])-int(p[1])
        if L<=0: continue
        md=int(p[4])/L; pool+=md
        if orf in allcargo: dep[orf]=md
    MC.append([sum(dep.get(o,0) for o in S)/pool*100 for S in (arg_orfs,metal_orfs,vf_orfs)])
    MD.append([sum(dep.get(o,0) for o,cs in o2cls.items() if cl in cs)/pool*100 for cl in CL])
MC=np.array(MC); MD=np.array(MD)
CATS=["ARG","Metal","VF"]

panels=[("Plasmid abundance",MM,TIERS),("Cargo",MC,CATS),("Drug class",MD,CL)]
fig,axes=plt.subplots(1,3,figsize=(9.6,3.6),gridspec_kw=dict(width_ratios=[5,3,12],wspace=0.07))
for k,(ax,(title,M,xl)) in enumerate(zip(axes,panels)):
    im=ax.imshow(M,aspect="auto",cmap="Reds",norm=mc.LogNorm(vmin=max(M[M>0].min(),1e-3),vmax=M.max()))
    ax.set_xticks(range(len(xl))); ax.set_xticklabels(xl,fontsize=9,fontweight="bold",rotation=40,ha="right")
    ax.set_yticks(range(6))
    ax.set_yticklabels(DISP if k==0 else [],fontsize=10,fontweight="bold")
    ax.set_xticks(np.arange(-.5,len(xl),1),minor=True); ax.set_yticks(np.arange(-.5,6,1),minor=True)
    ax.grid(which="minor",color="black",lw=1.0); ax.tick_params(which="minor",length=0)
    for s in ax.spines.values(): s.set_edgecolor("black"); s.set_linewidth(1.2)
    thr=np.percentile(M[M>0],68)
    for i in range(6):
        for j in range(len(xl)):
            v=M[i,j]
            if v>0: ax.text(j,i,(f"{v:.0f}" if v>=10 else (f"{v:.1f}" if v>=1 else f"{v:.2f}")),ha="center",va="center",
                            fontsize=7.5,fontweight="bold",color="white" if v>thr else "black")
    ax.set_title(title,fontsize=12,fontweight="bold")
fig.savefig(f"{A}/Fig3cd_combined.png",dpi=600,bbox_inches="tight"); plt.close(fig)
from PIL import Image; im=Image.open(f"{A}/Fig3cd_combined.png"); w,h=im.size; im.resize((1100,int(1100*h/w))).save("/tmp/fig3row.png")
print("plasmid ORFs:",len(porf),"| ARG/Metal/VF ORFs:",len(arg_orfs),len(metal_orfs),len(vf_orfs),"| classes:",CL)
