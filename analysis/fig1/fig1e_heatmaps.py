#!/usr/bin/env python
"""Render Fig 1e heatmaps (horizontal: samples rows x 14 tiles cols). No title, no category strip.
Versions: relative (Greys), relative_zscore (diverging), tpm (Greys) — renders whichever TSVs exist."""
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
A="/home/gchoi/wwtp_plasmidome/analysis/mag"
SAMP=["Inf","Ana","Anx","Oxi","RAS","Eff"]
TORDER=["Acetate","Propionate","Butyrate","Ammonia oxidation","Nitrite oxidation",
        "Nitrite reduction","NO reduction","N2O reduction","poly-P kinase","PHA storage",
        "Pi transport (pst)","GH","PL","CE"]
import re as _re
_mk=pd.read_csv(f"{A}/function_ko_gene.tsv",sep="\t")
_mkd=_mk[_mk["marker"]=="Y"].groupby("function")["gene_symbol"].apply(lambda s:"/".join(dict.fromkeys(s)))
def _lab(t):
    if t in ("GH","PL","CE"): return t
    return f"{_re.sub(r'\s*\(.*\)$','',t)} ({_mkd.get(t,'')})"
LABELS=[_lab(t) for t in TORDER]

def render(tsv,out,label,cmap,diverging=False,log=False,zscore=False,power=None):
    if not os.path.exists(tsv): print("skip (no file):",tsv); return
    d=pd.read_csv(tsv,sep="\t",index_col=0).reindex(TORDER)
    M=d[SAMP].values.T.astype(float)          # rows=samples, cols=tiles
    if log: M=np.log10(M+1.0)
    if zscore:                                 # per-tile (column) z-score across samples
        M=(M-np.nanmean(M,axis=0,keepdims=True))/(np.nanstd(M,axis=0,keepdims=True)+1e-12)
        diverging=True
    fig,ax=plt.subplots(figsize=(11,3.0))
    if diverging:
        v=np.nanmax(np.abs(M)); im=ax.imshow(M,aspect="auto",cmap="RdBu_r",vmin=-v,vmax=v)
    elif power is not None:
        im=ax.imshow(M,aspect="auto",cmap=cmap,norm=PowerNorm(gamma=power,vmin=0,vmax=np.nanmax(M)))
    else:
        im=ax.imshow(M,aspect="auto",cmap=cmap,vmin=0,vmax=np.nanmax(M))
    ax.set_xticks(range(len(TORDER))); ax.set_xticklabels(LABELS,rotation=40,ha="right",fontsize=10,fontweight="bold")
    ax.set_yticks(range(len(SAMP))); ax.set_yticklabels(SAMP,fontweight="bold",fontsize=11)
    for lab in ax.get_yticklabels()+ax.get_xticklabels(): lab.set_fontweight("bold")
    cb_bold=True
    # black cell borders
    ax.set_xticks(np.arange(-0.5,len(TORDER),1),minor=True)
    ax.set_yticks(np.arange(-0.5,len(SAMP),1),minor=True)
    ax.grid(which="minor",color="black",linewidth=0.6); ax.tick_params(which="minor",length=0)
    for sp in ax.spines.values(): sp.set_edgecolor("black"); sp.set_linewidth(0.8)
    cb=fig.colorbar(im,ax=ax,shrink=0.85,pad=0.012); cb.set_label(label,fontsize=12,fontweight="bold")
    for t in cb.ax.get_yticklabels(): t.set_fontweight("bold")
    fig.tight_layout(); fig.savefig(out,dpi=600,bbox_inches="tight"); plt.close(fig)
    print("saved",os.path.basename(out))

render(f"{A}/Fig1e_tpm.tsv",f"{A}/Fig1e_tpm.png","TPM (log₁₀)","Greys",log=True)
