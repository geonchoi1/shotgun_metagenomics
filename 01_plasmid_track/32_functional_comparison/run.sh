#!/bin/bash
# === Functional comparison across environments ===
# Counts Pfam/KO per environment (sample-of-origin parsed from contig prefix)
# Filter ≥100 total / ≥2 env
# 3 tests x 2 weighting modes = PCA + GSEA + Fisher (richness + TPM-weighted)
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

PFAM_TBLOUT=$PROJECT/plasmid/05_pfam/pfam.tblout
KOFAM_MAPPER=$PROJECT/plasmid/07_kofamscan/kofam.mapper.tsv
ORF2CONTIG=$PROJECT/plasmid/04_master_orf/orf2contig.tsv
COVERM_DIR=$PROJECT/plasmid/40_quantification
SAMPLE_ENV=${SAMPLE_ENV:-$PROJECT/sample_to_env.tsv}   # optional user-supplied: sample_id<TAB>env
OUT=$PROJECT/plasmid/32_functional_comparison
mkdir -p $OUT

KEGG_GMT=${KEGG_GMT:-$DB_ROOT/gsea/kegg_pathway.gmt}

# Auto-build KEGG GMT from REST API if missing (one-time, ~30s)
if [ ! -s "$KEGG_GMT" ]; then
    echo "[$(date +%F\ %T)] KEGG GMT not found — building from KEGG REST API"
    mkdir -p "$(dirname $KEGG_GMT)"
    python3 "$SCRIPT_DIR/build_kegg_gmt.py" --out "$KEGG_GMT"
fi

activate_env "$ENV_DIAMOND"  # base — has python sklearn/pandas/scipy/gseapy

python3 - <<PYEOF
import os, re, sys
from collections import defaultdict, Counter
import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from scipy.stats import fisher_exact

OUT="$OUT"
os.makedirs(OUT, exist_ok=True)

# ---- 1. ORF → contig → sample → env ----
orf2c={}
with open("$ORF2CONTIG") as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)>=2: orf2c[p[0]]=p[1]

# sample id parsed from contig prefix before "|" or "_" (try both)
def sample_of(contig):
    if '|' in contig: return contig.split('|')[0]
    return contig.split('_')[0]

sample_env={}
if os.path.exists("$SAMPLE_ENV"):
    with open("$SAMPLE_ENV") as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                p=line.rstrip('\n').split('\t')
                if len(p)>=2: sample_env[p[0]]=p[1]
else:
    # fallback: env = sample (each sample is its own env)
    pass

def env_of(orf):
    c=orf2c.get(orf)
    if not c: return None
    s=sample_of(c)
    return sample_env.get(s, s)

# ---- 2. Pfam annotation per ORF ----
orf2pfam=defaultdict(set)
with open("$PFAM_TBLOUT") as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        parts=re.split(r'\s+', line.rstrip())
        orf=parts[0]; pfam=parts[3]  # accession column
        orf2pfam[orf].add(pfam)

# KO
orf2ko={}
with open("$KOFAM_MAPPER") as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)>=2 and p[1]: orf2ko[p[0]]=p[1]

# ---- 3. TPM table per contig (CoverM output) ----
contig_tpm=defaultdict(dict)
import glob
for tf in glob.glob("$COVERM_DIR/coverm_*.tsv"):
    sample=os.path.basename(tf).replace('coverm_','').replace('.tsv','')
    with open(tf) as f:
        header=f.readline().split('\t')
        # Locate TPM column
        tpm_idx=None
        for i,h in enumerate(header):
            if 'TPM' in h: tpm_idx=i; break
        if tpm_idx is None: continue
        for line in f:
            p=line.rstrip('\n').split('\t')
            try: contig_tpm[p[0]][sample]=float(p[tpm_idx])
            except: pass

# ---- 4. Build env × feature matrices (richness + TPM-weighted) ----
def build_matrix(orf2feat):
    envs=set()
    feat_env_richness=defaultdict(lambda: defaultdict(int))
    feat_env_tpm=defaultdict(lambda: defaultdict(float))
    for orf,feats in orf2feat.items() if isinstance(orf2feat, dict) else []:
        e=env_of(orf)
        if not e: continue
        envs.add(e)
        if isinstance(feats, set):
            for f_ in feats:
                feat_env_richness[f_][e]+=1
        else:
            feat_env_richness[feats][e]+=1
        # TPM-weighted: sum TPM across samples in this env mapped to contig
        c=orf2c.get(orf)
        if c and c in contig_tpm:
            for s,t in contig_tpm[c].items():
                if sample_env.get(s,s)==e:
                    if isinstance(feats, set):
                        for f_ in feats: feat_env_tpm[f_][e]+=t
                    else:
                        feat_env_tpm[feats][e]+=t
    envs=sorted(envs)
    return envs, feat_env_richness, feat_env_tpm

def write_pca(mat_dict, envs, prefix, label):
    feats=sorted(mat_dict.keys())
    M=np.array([[mat_dict[f].get(e,0) for e in envs] for f in feats], dtype=float)
    # Filter: feature ≥100 total AND present in ≥2 envs
    keep=[i for i,row in enumerate(M) if row.sum()>=100 and (row>0).sum()>=2]
    if len(keep)<3 or len(envs)<2:
        with open(f"{OUT}/{prefix}.pca.tsv","w") as o: o.write("(skipped — too few)\n")
        return None,None
    Mk=M[keep,:]
    fk=[feats[i] for i in keep]
    Z=(Mk - Mk.mean(axis=1, keepdims=True))/(Mk.std(axis=1, keepdims=True)+1e-9)
    pca=PCA(n_components=min(3,Mk.shape[1],Mk.shape[0]))
    coords=pca.fit_transform(Z.T)
    pd.DataFrame(coords, index=envs, columns=[f"PC{i+1}" for i in range(coords.shape[1])])\
      .to_csv(f"{OUT}/{prefix}.pca.tsv", sep='\t')
    pd.DataFrame(Mk, index=fk, columns=envs).to_csv(f"{OUT}/{prefix}.matrix.tsv", sep='\t')
    return Mk, fk

def fisher_envs(mat_dict, envs, prefix, label):
    feats=sorted(mat_dict.keys())
    rows=[]
    M=np.array([[mat_dict[f].get(e,0) for e in envs] for f in feats], dtype=float)
    keep=[i for i,row in enumerate(M) if row.sum()>=100 and (row>0).sum()>=2]
    M=M[keep,:]; feats=[feats[i] for i in keep]
    total=M.sum()
    for fi,f_ in enumerate(feats):
        for ei,e in enumerate(envs):
            a=M[fi,ei]
            b=M[fi,:].sum()-a
            c=M[:,ei].sum()-a
            d=total-a-b-c
            if a+b<1 or c+d<1: continue
            try:
                OR,p=fisher_exact([[a,b],[c,d]], alternative='greater')
                rows.append((f_, e, int(a), int(b), int(c), int(d), OR, p))
            except: pass
    pd.DataFrame(rows, columns=['feature','env','a','b','c','d','OR','p'])\
      .to_csv(f"{OUT}/{prefix}.fisher.tsv", sep='\t', index=False)

def gsea_run(mat_rich, mat_tpm, envs, prefix):
    try:
        import gseapy
    except ImportError:
        with open(f"{OUT}/{prefix}.gsea.tsv","w") as o: o.write("gseapy missing\n")
        return
    if not os.path.exists("$KEGG_GMT"):
        with open(f"{OUT}/{prefix}.gsea.tsv","w") as o: o.write("KEGG_GMT missing\n")
        return
    # Build KO rank by env diff: env-of-interest mean vs others
    # For each env, run prerank with rank = z(env) - mean(z(other envs))
    feats=sorted(mat_rich.keys())
    M=np.array([[mat_rich[f].get(e,0) for e in envs] for f in feats], dtype=float)
    out_all=[]
    for ei,e in enumerate(envs):
        rank=M[:,ei] - M[:,[j for j in range(len(envs)) if j!=ei]].mean(axis=1)
        df=pd.DataFrame({'feat':feats,'rank':rank}).sort_values('rank', ascending=False)
        try:
            r=gseapy.prerank(rnk=df, gene_sets="$KEGG_GMT", outdir=None, min_size=5, max_size=2000, seed=42)
            res=r.res2d.copy(); res['env']=e
            out_all.append(res)
        except Exception as ex:
            print(f"gsea {e} fail: {ex}")
    if out_all:
        pd.concat(out_all).to_csv(f"{OUT}/{prefix}.gsea.tsv", sep='\t', index=False)

# ---- 5. Run all combinations ----
for label, src in [("pfam", orf2pfam), ("ko", orf2ko)]:
    envs, mat_rich, mat_tpm = build_matrix(src)
    if not envs: continue
    for mode, mat in [("richness", mat_rich), ("tpm", mat_tpm)]:
        prefix=f"{label}_{mode}"
        write_pca(mat, envs, prefix, label)
        fisher_envs(mat, envs, prefix, label)
        if label=="ko":
            gsea_run(mat_rich, mat_tpm, envs, prefix)
print("functional comparison done")
PYEOF
echo "[$(date '+%F %T')] DONE"
