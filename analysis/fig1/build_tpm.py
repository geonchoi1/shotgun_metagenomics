#!/usr/bin/env python
"""Fig 1e v1 (TPM): gene-level TPM from per-ORF coverage (bedcov).
mean_depth_orf = bedcov_sum / ORF_length ; TPM_orf = mean_depth / sum(all-ORF mean_depth) * 1e6 (per sample).
tile value = SUM of TPM over the tile's ORFs in that sample (total functional abundance)."""
import pandas as pd, numpy as np
A="/home/gchoi/wwtp_plasmidome/analysis/mag"
SAMP=["IN","Anaerobic","Anoxic","Oxic","RAS","EF"]
DISP=dict(zip(SAMP,["Inf","Ana","Anx","Oxi","RAS","Eff"]))
TILES=[("Acetate","C: VFA"),("Propionate","C: VFA"),("Butyrate","C: VFA"),
       ("Ammonia oxidation","N: Nitrification"),("Nitrite oxidation","N: Nitrification"),
       ("Nitrite reduction","N: Denitrification"),("NO reduction","N: Denitrification"),("N2O reduction","N: Denitrification"),
       ("poly-P kinase","P: PAO"),("PHA storage","P: PAO"),("Pi transport (pst)","P: PAO"),
       ("GH","C: Carb. degradation"),("PL","C: Carb. degradation"),("CE","C: Carb. degradation")]
TORDER=[t for t,_ in TILES]; TCAT=dict(TILES)

# orf -> set(tiles), per sample (from build/orf_tile_sample.tsv)
L=pd.read_csv(f"{A}/build/orf_tile_sample.tsv",sep="\t")
orf2tiles={}
for orf,t in zip(L["orf"],L["tile"]): orf2tiles.setdefault(orf,set()).add(t)
tile_orfs=set(orf2tiles)

# per sample: D_s = sum of all-ORF mean depth ; store mean depth for tile ORFs only
Ds={}; mdepth={}                       # mdepth[(s,orf)] = mean depth
for s in SAMP:
    tot=0.0
    for ln in open(f"{A}/cov/bedcov/{s}.bedcov"):
        p=ln.rstrip("\n").split("\t")
        Ln=int(p[2])-int(p[1]); md=float(p[4])/Ln if Ln>0 else 0.0
        tot+=md
        if p[3] in tile_orfs: mdepth[(s,p[3])]=md
    Ds[s]=tot
print("D_s (sum all-ORF mean depth):",{DISP[s]:round(Ds[s],1) for s in SAMP})

# tile x sample = sum TPM over tile ORFs assigned to that sample
M=np.zeros((len(TORDER),len(SAMP)))
ti={t:i for i,t in enumerate(TORDER)}
for orf,sp,t in zip(L["orf"],L["sample"],L["tile"]):
    md=mdepth.get((sp,orf))
    if md is None: continue
    M[ti[t],SAMP.index(sp)] += md/Ds[sp]*1e6
out=pd.DataFrame(M.round(2),index=TORDER,columns=[DISP[s] for s in SAMP])
out.insert(0,"category",[TCAT[t] for t in TORDER])
out.to_csv(f"{A}/Fig1e_tpm.tsv",sep="\t")
print("\n=== Fig1e_tpm.tsv (sum TPM per tile) ===\n",out.to_string())
