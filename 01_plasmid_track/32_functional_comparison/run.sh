#!/bin/bash
# === Functional comparison across environments ===
# Counts Pfam/KO per environment (sample-of-origin parsed from contig prefix)
# Filter ≥MIN_CARRIERS / ≥2 env (statsmodels power analysis: OR=4, BH-FDR α≈0.001, power 0.80)
# Two modes share the SAME filter (richness-count basis) but use mode-appropriate tests:
#   * richness MAIN — Fisher exact (binary count) + log_OR + GSEA(log_OR rank)
#                     Fiamenghi 2025-style: divide by env plasmid count -> proportion -> StandardScaler -> PCA
#   * TPM-weighted SUPPL — Mann-Whitney U (continuous) + log2FC + GSEA(log2FC rank)
#                     TPM is already sample-level normalized; no further normalization before PCA/test
#
# Packages:
#   * statsmodels.stats.power.NormalIndPower  — MIN_CARRIERS helper (compute_min_carriers.py)
#   * sklearn StandardScaler + decomposition.PCA
#   * scipy.stats.fisher_exact + scipy.stats.odds_ratio (conditional MLE)
#   * scipy.stats.mannwhitneyu                 — TPM test
#   * scipy.stats.false_discovery_control       — BH-FDR
#   * gseapy.prerank (Subramanian 2005; equivalent to R fgsea/hypeR)
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
#   it in the analysis. Power analysis (statsmodels NormalIndPower; n1=188 IN, n2=300 EF;
#   two-sided; BH-FDR α≈0.001; power=0.80; OR=4) suggests ≥50; Fiamenghi 2025 used 100.
#   Use compute_min_carriers.py to derive for your data.
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
from scipy.stats import fisher_exact, mannwhitneyu, false_discovery_control
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

# ---- 3. TPM table per contig (CoverM output; already sample-level normalized to 1e6) ----
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

# ---- 4. Build env × feature structures (richness count, TPM sum, TPM distribution) ----
# mat_rich[f][e]      : ORF count (richness)
# mat_tpm_sum[f][e]   : sum of contig TPM across plasmids carrying f in env e (for PCA)
# mat_tpm_dist[f][e]  : list of per-ORF contig TPMs (for Mann-Whitney U distribution comparison)
def build_matrix(orf2feat):
    envs=set()
    mat_rich=defaultdict(lambda: defaultdict(int))
    mat_tpm_sum=defaultdict(lambda: defaultdict(float))
    mat_tpm_dist=defaultdict(lambda: defaultdict(list))
    for orf,feats in (orf2feat.items() if isinstance(orf2feat, dict) else []):
        c=orf2c.get(orf)
        e=env_of(orf)
        if not e: continue
        envs.add(e)
        iter_feats = feats if isinstance(feats, set) else [feats]
        # ORF richness
        for f_ in iter_feats: mat_rich[f_][e] += 1
        # TPM aggregation per ORF (only env-matched samples)
        if c and c in contig_tpm:
            tpm_in_env = sum(t for s,t in contig_tpm[c].items() if sample_env.get(s,s)==e)
            for f_ in iter_feats:
                mat_tpm_sum[f_][e] += tpm_in_env
                if tpm_in_env > 0:
                    mat_tpm_dist[f_][e].append(tpm_in_env)
    envs=sorted(envs)
    return envs, mat_rich, mat_tpm_sum, mat_tpm_dist

# ---- 4.5. Plasmid count per env (for richness proportion only) ----
def env_plasmid_counts():
    env_p=defaultdict(set)
    for orf,c in orf2c.items():
        s=sample_of(c)
        e=sample_env.get(s,s)
        env_p[e].add(c)
    return {e:len(s) for e,s in env_p.items()}

def write_pca(mat_dict, envs, prefix, label, n_per_env=None):
    """PCA with row-wise StandardScaler.
    If n_per_env given (richness mode), divide by env plasmid count first -> proportion.
    If None (TPM mode), use values as-is (TPM is already sample-level normalized).
    """
    feats=sorted(mat_dict.keys())
    if len(feats)<3 or len(envs)<2:
        with open(f"{OUT}/{prefix}.pca.tsv","w") as o: o.write("(skipped — too few)\n")
        return None,None
    Mk=np.array([[mat_dict[f].get(e,0) for e in envs] for f in feats], dtype=float)
    fk=feats
    if n_per_env is not None:
        # Richness: divide by env plasmid count -> proportion (Fiamenghi 2025)
        n_vec=np.array([max(n_per_env.get(e,1),1) for e in envs], dtype=float)
        Pk=Mk / n_vec[np.newaxis, :]
        pd.DataFrame(Pk, index=fk, columns=envs).to_csv(f"{OUT}/{prefix}.proportion.tsv", sep='\t')
    else:
        # TPM: use raw aggregated TPM directly
        Pk=Mk
    # Row-wise StandardScaler then PCA
    Z=StandardScaler().fit_transform(Pk.T).T
    pca=PCA(n_components=min(3,Mk.shape[1],Mk.shape[0]))
    coords=pca.fit_transform(Z.T)
    pd.DataFrame(coords, index=envs, columns=[f"PC{i+1}" for i in range(coords.shape[1])])\
      .to_csv(f"{OUT}/{prefix}.pca.tsv", sep='\t')
    pd.DataFrame(pca.components_.T, index=fk, columns=[f"PC{i+1}" for i in range(pca.n_components_)])\
      .to_csv(f"{OUT}/{prefix}.pca_loadings.tsv", sep='\t')
    pd.DataFrame({'PC':[f"PC{i+1}" for i in range(pca.n_components_)],
                  'var_explained':pca.explained_variance_ratio_})\
      .to_csv(f"{OUT}/{prefix}.pca_variance.tsv", sep='\t', index=False)
    pd.DataFrame(Mk, index=fk, columns=envs).to_csv(f"{OUT}/{prefix}.matrix.tsv", sep='\t')
    return Mk, fk

def fisher_envs(mat_dict, envs, prefix, label):
    """Richness mode: Fisher exact (two-sided) + conditional MLE OR + BH-FDR.
    Reports log_OR (natural log, matches Fiamenghi 2025) AND log2OR (convenience).
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
                    OR = odds_ratio(tab).statistic
                else:
                    OR = (a*d) / max(b*c, 1e-9)
                if OR and OR > 0 and np.isfinite(OR):
                    log_OR  = np.log(OR)
                    log2OR  = np.log2(OR)
                else:
                    log_OR  = np.nan
                    log2OR  = np.nan
                rows.append((f_, e, int(a), int(b), int(c), int(d), OR, log_OR, log2OR, p))
            except: pass
    df = pd.DataFrame(rows, columns=['feature','env','a','b','c','d','OR','log_OR','log2OR','p'])
    df['FDR'] = np.nan
    for e in envs:
        mask = (df['env'] == e) & df['p'].notna()
        if mask.sum() > 0:
            df.loc[mask, 'FDR'] = false_discovery_control(df.loc[mask, 'p'].values)
    df.to_csv(f"{OUT}/{prefix}.fisher.tsv", sep='\t', index=False)
    return df  # consumed by gsea_run (rank by log_OR)

def mwu_envs(mat_dist, envs, prefix, label):
    """TPM mode: Mann-Whitney U (env vs rest) + log2FC + rank-biserial + BH-FDR.
    Input mat_dist[f][e] is a list of per-ORF contig TPMs.
    """
    feats=sorted(mat_dist.keys())
    rows=[]
    for f_ in feats:
        for e in envs:
            X1 = mat_dist[f_].get(e, [])
            X2 = []
            for ef in envs:
                if ef != e:
                    X2.extend(mat_dist[f_].get(ef, []))
            n1, n2 = len(X1), len(X2)
            if n1 < 3 or n2 < 3:
                continue
            try:
                U, p = mannwhitneyu(X1, X2, alternative='two-sided')
                rb = 1.0 - 2.0 * U / (n1 * n2)   # rank-biserial r
                m1 = float(np.mean(X1)) if X1 else 0.0
                m2 = float(np.mean(X2)) if X2 else 0.0
                # log2 FC with pseudocount to handle zero means
                log2FC = float(np.log2((m1 + 1e-9) / (m2 + 1e-9)))
                rows.append((f_, e, n1, n2, float(U), m1, m2, log2FC, rb, p))
            except Exception:
                continue
    df = pd.DataFrame(rows, columns=['feature','env','n1','n2','U','mean1','mean2','log2FC','rank_biserial','p'])
    df['FDR'] = np.nan
    for e in envs:
        mask = (df['env'] == e) & df['p'].notna()
        if mask.sum() > 0:
            df.loc[mask, 'FDR'] = false_discovery_control(df.loc[mask, 'p'].values)
    df.to_csv(f"{OUT}/{prefix}.mwu.tsv", sep='\t', index=False)
    return df  # consumed by gsea_run (rank by log2FC)

def gsea_run(test_df, envs, prefix, rank_col):
    """GSEA per env via gseapy.prerank (equivalent to fgsea/hypeR).
    rank_col = 'log_OR' (richness/Fisher) or 'log2FC' (TPM/MW-U).
    Clipped to ±10 for stability.
    """
    try:
        import gseapy
    except ImportError:
        with open(f"{OUT}/{prefix}.gsea.tsv","w") as o: o.write("gseapy missing\n")
        return
    if not os.path.exists("$KEGG_GMT"):
        with open(f"{OUT}/{prefix}.gsea.tsv","w") as o: o.write("KEGG_GMT missing\n")
        return
    if test_df is None or len(test_df) == 0 or rank_col not in test_df.columns:
        with open(f"{OUT}/{prefix}.gsea.tsv","w") as o: o.write("(skipped — test result empty)\n")
        return
    out_all=[]
    for e in envs:
        sub = test_df[test_df['env'] == e].copy()
        sub = sub[sub[rank_col].notna() & np.isfinite(sub[rank_col])]
        if len(sub) < 5: continue
        sub[rank_col] = sub[rank_col].clip(-10, 10)
        rnk = (sub[['feature', rank_col]]
               .sort_values(rank_col, ascending=False)
               .rename(columns={'feature':'gene', rank_col:'rank'}))
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
    envs, mat_rich, mat_tpm_sum, mat_tpm_dist = build_matrix(src)
    if not envs: continue
    # === Filter on RICHNESS-COUNT basis (same feature set for both modes) ===
    feats_all=sorted(mat_rich.keys())
    M_rich=np.array([[mat_rich[f].get(e,0) for e in envs] for f in feats_all], dtype=float)
    keep_feats=set(feats_all[i] for i,row in enumerate(M_rich)
                   if row.sum()>=MIN_CARRIERS and (row>0).sum()>=2)
    mat_rich_f      = {f:v for f,v in mat_rich.items()     if f in keep_feats}
    mat_tpm_sum_f   = {f:v for f,v in mat_tpm_sum.items()  if f in keep_feats}
    mat_tpm_dist_f  = {f:v for f,v in mat_tpm_dist.items() if f in keep_feats}
    print(f"  {label}: filtered {len(keep_feats)} / {len(feats_all)} features (MIN_CARRIERS={MIN_CARRIERS})", flush=True)

    # Richness branch — Fisher + log_OR + GSEA(log_OR rank)
    prefix_rich=f"{label}_richness"
    write_pca(mat_rich_f, envs, prefix_rich, label, n_per_env=n_per_env)
    fisher_df = fisher_envs(mat_rich_f, envs, prefix_rich, label)
    if label=="ko":
        gsea_run(fisher_df, envs, prefix_rich, rank_col='log_OR')

    # TPM branch — MW-U + log2FC + GSEA(log2FC rank)
    prefix_tpm=f"{label}_tpm"
    write_pca(mat_tpm_sum_f, envs, prefix_tpm, label, n_per_env=None)  # no normalize
    mwu_df = mwu_envs(mat_tpm_dist_f, envs, prefix_tpm, label)
    if label=="ko":
        gsea_run(mwu_df, envs, prefix_tpm, rank_col='log2FC')

print("functional comparison done")
PYEOF
echo "[$(date '+%F %T')] DONE"
