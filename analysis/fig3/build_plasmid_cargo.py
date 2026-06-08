#!/usr/bin/env python
"""Per-plasmid cargo table: map ORF(locus tag)->contig via plasmidome.master.gff, then attach
AMRFinder ARG/stress, BacMet metal, VFDB virulence + master(sample,tier,length) + TPM(6 zones)."""
import pandas as pd, re, json
P="/home/gchoi/wwtp_plasmidome/01_plasmid_track"; A="/home/gchoi/wwtp_plasmidome/analysis/plasmid_cargo"; CAT="/home/gchoi/wwtp_plasmidome/analysis/plasmid_catalog"
# ORF -> contig
orf2c={}
for ln in open(f"{P}/04_master_orf/plasmidome.master.gff"):
    if "\tCDS\t" not in ln: continue
    f=ln.split("\t"); m=re.search(r"ID=([^;]+)",f[8])
    if m: orf2c[m.group(1)]=f[0]
print("ORF->contig:",len(orf2c))
def orfcol(df): return df.iloc[:,0].map(orf2c)
# AMRFinder (ARG + stress)
amr=pd.read_csv(f"{P}/09_amrfinder/amr.tsv",sep="\t")
amr["contig"]=amr["Protein id"].map(orf2c)
arg=amr[amr.Type=="AMR"].dropna(subset=["contig"])
stress=amr[amr.Type=="STRESS"].dropna(subset=["contig"])
arg_by=arg.groupby("contig")["Class"].apply(lambda s:";".join(sorted(set(s))))
# BacMet metal, VFDB
bm=pd.read_csv(f"{P}/10_bacmet/bacmet.tsv",sep="\t",header=None); bm["contig"]=bm[0].map(orf2c)
vf=pd.read_csv(f"{P}/11_vfdb/vfdb.tsv",sep="\t",header=None); vf["contig"]=vf[0].map(orf2c)
# master + TPM
m=pd.read_csv(f"{CAT}/plasmid_master.tsv",sep="\t")
tpm=pd.read_csv(f"{P}/40_quantification/tpm_matrix.tsv",sep="\t").rename(columns={"Contig":"contig"})
m=m.merge(tpm,on="contig",how="left")
m["n_arg"]=m.contig.map(arg.groupby("contig").size()).fillna(0).astype(int)
m["n_metal"]=m.contig.map(bm.dropna(subset=["contig"]).groupby("contig").size()).fillna(0).astype(int)
m["n_metal_stress"]=m.contig.map(stress.groupby("contig").size()).fillna(0).astype(int)
m["n_vf"]=m.contig.map(vf.dropna(subset=["contig"]).groupby("contig").size()).fillna(0).astype(int)
m["arg_classes"]=m.contig.map(arg_by).fillna("")
m["has_arg"]=m.n_arg>0; m["has_metal"]=(m.n_metal>0)|(m.n_metal_stress>0); m["has_vf"]=m.n_vf>0
m.to_csv(f"{A}/plasmid_cargo.tsv",sep="\t",index=False)
print(f"plasmids {len(m)} | ARG-carrying {m.has_arg.sum()} | metal {m.has_metal.sum()} | VF {m.has_vf.sum()}")
print("by zone (sample) ARG-carriers:",m[m.has_arg].groupby('sample').size().to_dict())
