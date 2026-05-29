#!/bin/bash
# === Functional comparison across environments ===
# Counts Pfam/KO per environment (sample-of-origin parsed from contig prefix)
# Filter ≥MIN_CARRIERS / ≥2 env (statsmodels power analysis: OR=4, BH-FDR α≈0.001, power 0.80)
# 3 tests x 2 weighting modes = PCA + Fisher + GSEA
#   (richness MAIN + TPM-weighted SUPPL ; all plasmids)
#
# Methods aligned with Fiamenghi 2025 (Nat Commun) plasmidome framework:
#   * Normalize raw counts by env plasmid count (proportion per env)
#   * sklearn StandardScaler row-wise (per function) before PCA
#   * scipy.stats.fisher_exact (two-sided) + scipy.stats.odds_ratio (conditional MLE)
#   * log_OR = natural log (matches paper); log2OR also reported for convenience
#   * scipy.stats.false_discovery_control (Benjamini-Hochberg FDR, q<0.05)
#   * GSEA: rank input = Fisher log_OR (DEG-like pipeline, equivalent to RNA-seq log2FC→GSEA)
#   * gseapy.prerank (Subramanian 2005 GSEA; equivalent to R fgsea/hypeR)
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

PFAM_TBLOUT=$PROJECT/plasmid/05_pfam/pfam.tblout
KOFAM_MAPPER=$PROJECT/plasmid/07_kofamscan/kofam.mapper.tsv
ORF2CONTIG=$PROJECT/plasmid/04_master_orf/orf2contig.tsv
COVERM_DIR=$PROJECT/plasmid/40_quantification
SAMPLE_ENV=${SAMPLE_ENV:-$PROJECT/sample_to_env.tsv}   # sample_id<TAB>env
OUT=$PROJECT/plasmid/32_functional_comparison
mkdir -p $OUT

KEGG_GMT=${KEGG_GMT:-$DB_ROOT/gsea/kegg_pathway.gmt}

# === User MUST set the filter threshold ===
# MIN_CARRIERS = minimum number of plasmids carrying a function to include
#   it in the analysis. Power analysis (Fleiss approximation; n1=188 IN, n2=300 EF;
#   two-sided; BH-FDR α≈0.001; power=0.80; OR=4) suggests ≥50; Fiamenghi 2025 used 100.
#   Choose based on (a) target effect size, (b) multiple-testing burden,
#   (c) trade-off between sensitivity and false-positive control.
#   Example: export MIN_CARRIERS=50
: ${MIN_CARRIERS:?ERROR: please export MIN_CARRIERS — e.g. export MIN_CARRIERS=50 (see header comment for guidance)}

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
from sklearn.preprocessing import StandardScaler
from scipy.stats import fisher_exact, false_discovery_control
try:
    from scipy.stats import odds_ratio
    HAS_OR_FN = True
except ImportError:
    HAS_OR_FN = False  # scipy < 1.10 fallback to sample OR

MIN_CARRIERS = $MIN_CARRIERS  # user-supplied via export MIN_CARRIERS=<int>

OUT="$OUT"
os.makedirs(OUT, exist_ok=True)

# ---- 1. ORF → contig → sample → env ----
orf2c={}
with open("$ORF2CONTIG") as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)>=2: orf2c[p[0]]=p[1]

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
        orf=parts[0]; pfam=parts[3]
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
    for orf,feats in (orf2feat.items() if isinstance(orf2feat, dict) else []):
        c=orf2c.get(orf)
        e=env_of(orf)
        if not e: continue
        envs.add(e)
        if isinstance(feats, set):
            for f_ in feats: feat_env_richness[f_][e]+=1
        else:
            feat_env_richness[feats][e]+=1
        if c and c in contig_tpm:
            for s,t in contig_tpm[c].items():
                if sample_env.get(s,s)==e:
                    if isinstance(feats, set):
                        for f_ in feats: feat_env_tpm[f_][e]+=t
                    else:
                        feat_env_tpm[feats][e]+=t
    envs=sorted(envs)
    return envs, feat_env_richness, feat_env_tpm

# ---- 4.5. Plasmid count per env (for proportion normalization) ----
def env_plasmid_counts():
    env_p=defaultdict(set)
    for orf,c in orf2c.items():
        s=sample_of(c)
        e=sample_env.get(s,s)
        env_p[e].add(c)
    return {e:len(s) for e,s in env_p.items()}

def write_pca(mat_dict, envs, prefix, label, n_per_env):
    """PCA after (1) divide by env plasmid count, (2) row-wise StandardScaler.
    Filtering is applied upstream (richness-count basis) before this call.
    """
    feats=sorted(mat_dict.keys())
    if len(feats)<3 or len(envs)<2:
        with open(f"{OUT}/{prefix}.pca.tsv","w") as o: o.write("(skipped — too few)\n")
        return None,None
    Mk=np.array([[mat_dict[f].get(e,0) for e in envs] for f in feats], dtype=float)
    fk=feats
    # Step 1: normalize by env plasmid count → proportion per env (Fiamenghi 2025)
    n_vec=np.array([max(n_per_env.get(e,1),1) for e in envs], dtype=float)
    Pk=Mk / n_vec[np.newaxis, :]
    # Step 2: row-wise StandardScaler (per function); sklearn scales columns → transpose
    Z=StandardScaler().fit_transform(Pk.T).T
    pca=PCA(n_components=min(3,Mk.shape[1],Mk.shape[0]))
    coords=pca.fit_transform(Z.T)
    pd.DataFrame(coords, index=envs, columns=[f"PC{i+1}" for i in range(coords.shape[1])])\
      .to_csv(f"{OUT}/{prefix}.pca.tsv", sep='\t')
    # Loadings (Fiamenghi 2025 Suppl Data 9 equivalent)
    pd.DataFrame(pca.components_.T, index=fk, columns=[f"PC{i+1}" for i in range(pca.n_components_)])\
      .to_csv(f"{OUT}/{prefix}.pca_loadings.tsv", sep='\t')
    pd.DataFrame({'PC':[f"PC{i+1}" for i in range(pca.n_components_)],
                  'var_explained':pca.explained_variance_ratio_})\
      .to_csv(f"{OUT}/{prefix}.pca_variance.tsv", sep='\t', index=False)
    pd.DataFrame(Mk, index=fk, columns=envs).to_csv(f"{OUT}/{prefix}.matrix.tsv", sep='\t')
    pd.DataFrame(Pk, index=fk, columns=envs).to_csv(f"{OUT}/{prefix}.proportion.tsv", sep='\t')
    return Mk, fk

def fisher_envs(mat_dict, envs, prefix, label):
    """Per-env vs rest Fisher exact (two-sided) + conditional MLE OR + BH-FDR.
    Reports log_OR (natural log, matches Fiamenghi 2025 paper text) AND log2OR (convenience).
    Filtering is applied upstream (richness-count basis) before this call.
    """
    feats=sorted(mat_dict.keys())
    rows=[]
    M=np.array([[mat_dict[f].get(e,0) for e in envs] for f in feats], dtype=float)
    total=M.sum()
    for fi,f_ in enumerate(feats):
        for ei,e in enumerate(envs):
            a=M[fi,ei]; b=M[fi,:].sum()-a
            c=M[:,ei].sum()-a; d=total-a-b-c
            if a+b<1 or c+d<1: continue
            try:
                tab=[[a,b],[c,d]]
                _, p = fisher_exact(tab, alternative='two-sided')
                if HAS_OR_FN:
                    OR = odds_ratio(tab).statistic  # conditional MLE (paper-grade)
                else:
                    OR = (a*d) / max(b*c, 1e-9)
                if OR and OR > 0 and np.isfinite(OR):
                    log_OR  = np.log(OR)            # natural log — matches Fiamenghi text
                    log2OR  = np.log2(OR)           # convenience
                else:
                    log_OR  = np.nan
                    log2OR  = np.nan
                rows.append((f_, e, int(a), int(b), int(c), int(d), OR, log_OR, log2OR, p))
            except: pass
    df = pd.DataFrame(rows, columns=['feature','env','a','b','c','d','OR','log_OR','log2OR','p'])
    # BH-FDR per env
    df['FDR'] = np.nan
    for e in envs:
        mask = (df['env'] == e) & df['p'].notna()
        if mask.sum() > 0:
            df.loc[mask, 'FDR'] = false_discovery_control(df.loc[mask, 'p'].values)
    df.to_csv(f"{OUT}/{prefix}.fisher.tsv", sep='\t', index=False)
    return df  # used by gsea_run() for log_OR rank metric

def gsea_run(fisher_df, envs, prefix):
    """GSEA per env via gseapy.prerank (equivalent to fgsea/hypeR).
    Rank = Fisher log_OR (DEG-like pipeline; matches Fiamenghi 2025 method).
    """
    try:
        import gseapy
    except ImportError:
        with open(f"{OUT}/{prefix}.gsea.tsv","w") as o: o.write("gseapy missing\n")
        return
    if not os.path.exists("$KEGG_GMT"):
        with open(f"{OUT}/{prefix}.gsea.tsv","w") as o: o.write("KEGG_GMT missing\n")
        return
    if fisher_df is None or len(fisher_df) == 0:
        with open(f"{OUT}/{prefix}.gsea.tsv","w") as o: o.write("(skipped — fisher empty)\n")
        return
    out_all=[]
    for e in envs:
        sub = fisher_df[fisher_df['env'] == e].copy()
        # Drop NaN/Inf, then clip extreme log_OR to ±10 for GSEA stability
        sub = sub[sub['log_OR'].notna() & np.isfinite(sub['log_OR'])]
        if len(sub) < 5: continue
        sub['log_OR'] = sub['log_OR'].clip(-10, 10)
        rnk = (sub[['feature','log_OR']]
               .sort_values('log_OR', ascending=False)
               .rename(columns={'feature':'gene','log_OR':'rank'}))
        try:
            r = gseapy.prerank(rnk=rnk, gene_sets="$KEGG_GMT", outdir=None,
                               min_size=5, max_size=2000, seed=42)
            res = r.res2d.copy(); res['env'] = e
            out_all.append(res)
        except Exception as ex:
            print(f"gsea {e} fail: {ex}")
    if out_all:
        pd.concat(out_all).to_csv(f"{OUT}/{prefix}.gsea.tsv", sep='\t', index=False)

# ---- 5. Run all combinations: 2 annotations × 2 modes ----
n_per_env = env_plasmid_counts()
print(f"\nplasmids per env: {dict(sorted(n_per_env.items()))}", flush=True)

for label, src in [("pfam", orf2pfam), ("ko", orf2ko)]:
    envs, mat_rich, mat_tpm = build_matrix(src)
    if not envs: continue
    # === Filter on RICHNESS-COUNT basis (richness ≥ MIN_CARRIERS AND ≥2 envs) ===
    # Apply the SAME feature set to both richness and TPM modes so they are directly comparable.
    feats_all=sorted(mat_rich.keys())
    M_rich=np.array([[mat_rich[f].get(e,0) for e in envs] for f in feats_all], dtype=float)
    keep_feats=set(feats_all[i] for i,row in enumerate(M_rich)
                   if row.sum()>=MIN_CARRIERS and (row>0).sum()>=2)
    mat_rich_f={f:v for f,v in mat_rich.items() if f in keep_feats}
    mat_tpm_f ={f:v for f,v in mat_tpm.items()  if f in keep_feats}
    print(f"  {label}: filtered {len(keep_feats)} / {len(feats_all)} features (MIN_CARRIERS={MIN_CARRIERS})", flush=True)
    for mode, mat in [("richness", mat_rich_f), ("tpm", mat_tpm_f)]:
        prefix=f"{label}_{mode}"
        write_pca(mat, envs, prefix, label, n_per_env)
        fisher_df = fisher_envs(mat, envs, prefix, label)   # returns df with log_OR
        if label=="ko":
            gsea_run(fisher_df, envs, prefix)               # GSEA ranked by log_OR
print("functional comparison done")
PYEOF
echo "[$(date '+%F %T')] DONE"
