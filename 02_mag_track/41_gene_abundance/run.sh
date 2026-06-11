#!/bin/bash
# === 41 Gene-level abundance (Method 2 of doi:10.1186/s40793-026-00892-w) ===
# Input:  $PROJECT/02_mag_track/03_master_orf/all/master.ffn   (gene nucleotides)
#         $PROJECT/02_mag_track/03_master_orf/all/master.gff3
#         $PROJECT/00_shared/01_read_qc/dehuman/*_clean.fastq.gz
#         $PROJECT/02_mag_track/04_pfam/pfam.tblout
#         $PROJECT/02_mag_track/06_kofamscan/kofam_mapper.tsv
# Output: $PROJECT/02_mag_track/41_gene_abundance/{bam,counts,tpm,pfam,kofam,pathway}/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FFN=$PROJECT/02_mag_track/03_master_orf/all/master.ffn
GFF=$PROJECT/02_mag_track/03_master_orf/all/master.gff3
READS_DIR=$PROJECT/00_shared/01_read_qc/dehuman
PFAM_TBL=$PROJECT/02_mag_track/04_pfam/pfam.tblout
KOFAM_TSV=$PROJECT/02_mag_track/06_kofamscan/kofam_mapper.tsv
OUT=$PROJECT/02_mag_track/41_gene_abundance
mkdir -p "$OUT"/{bam,counts,tpm,pfam,kofam,pathway}

[ -s "$FFN" ] || { echo "ERROR: master.ffn missing" >&2; exit 1; }
[ -s "$GFF" ] || { echo "ERROR: master.gff3 missing" >&2; exit 1; }

# --- Build SAF from GFF3 (CDS rows: chr,start,end,strand,gene_id=locus_tag) ---
SAF="$OUT/master.saf"
if [ ! -s "$SAF" ]; then
    echo "[$(date '+%F %T')] GFF3 -> SAF"
    python3 - <<PYEOF
import re
out = open("$SAF","w")
out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
with open("$GFF") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip(): continue
        if line.startswith(">"): break
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        attrs = dict(re.findall(r"([^=;]+)=([^;]+)", p[8]))
        lt = attrs.get("locus_tag") or attrs.get("ID","")
        if not lt: continue
        out.write(f"{lt}\t{lt}\t1\t{int(p[4])-int(p[3])+1}\t{p[6]}\n")
out.close()
PYEOF
fi

# --- minimap2 each sample to master.ffn (HiFi); samtools sort BAM ---
activate_env "$ENV_MINIMAP2"  # coverm env has minimap2 + samtools
echo "[$(date '+%F %T')] minimap2 map-hifi to master.ffn"
BAMS=()
for fq in "$READS_DIR"/*_clean.fastq.gz; do
    [ -f "$fq" ] || continue
    sample=$(basename "$fq" _clean.fastq.gz)
    bam="$OUT/bam/${sample}.bam"
    if [ -s "$bam" ]; then
        echo "  $sample already mapped"
        BAMS+=("$bam"); continue
    fi
    echo "  $sample"
    minimap2 -ax map-hifi -t "$THREADS_MINIMAP2" --secondary=no -I 50G "$FFN" "$fq" 2> "$OUT/bam/${sample}.mm2.log" \
        | samtools sort -@ 8 -o "$bam" -
    samtools index "$bam"
    BAMS+=("$bam")
done

# --- featureCounts (long-read mode -L, fractional counting) ---
activate_env "$ENV_HMMER"  # base env: featureCounts (subread) usually available; fallback below
if ! command -v featureCounts >/dev/null 2>&1; then
    # featureCounts often lives in coverm env
    activate_env "$ENV_COVERM"
fi
command -v featureCounts >/dev/null 2>&1 || { echo "ERROR: featureCounts not found; install subread" >&2; exit 1; }

echo "[$(date '+%F %T')] featureCounts (long-read, fractional)"
featureCounts \
    -a "$SAF" -F SAF \
    -L -T "$THREADS" -O --fraction --minOverlap 20 \
    -o "$OUT/counts/gene_counts.tsv" \
    "${BAMS[@]}"

# --- TPM, Pfam/KO aggregation, pathway TPM matrix, row-centered ---
echo "[$(date '+%F %T')] TPM + Pfam/KO aggregation + pathway"
python3 - <<PYEOF
import os
OUT = "$OUT"
counts_f = os.path.join(OUT, "counts/gene_counts.tsv")
# parse featureCounts (header + skip first '#' line)
with open(counts_f) as fh:
    line = fh.readline()
    while line.startswith("#"): line = fh.readline()
    header = line.rstrip("\n").split("\t")
    samples = [os.path.basename(h).replace(".bam","") for h in header[6:]]
    gene_len = {}; counts = {}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        gid, length = p[0], int(p[5])
        vals = [float(x) for x in p[6:]]
        gene_len[gid] = length
        counts[gid] = vals

# TPM = (count / length_kb) * 1e6 / sum
import math
n = len(samples)
tpm = {g: [0.0]*n for g in counts}
sums = [0.0]*n
rpk = {}
for g, vals in counts.items():
    lkb = max(gene_len[g],1) / 1000.0
    rpk[g] = [v/lkb for v in vals]
    for i,v in enumerate(rpk[g]): sums[i] += v
for g, vals in rpk.items():
    for i,v in enumerate(vals):
        tpm[g][i] = (v/sums[i]*1e6) if sums[i]>0 else 0.0

with open(os.path.join(OUT,"tpm/gene_tpm.tsv"),"w") as fh:
    fh.write("gene\t" + "\t".join(samples) + "\n")
    for g in sorted(tpm):
        fh.write(g + "\t" + "\t".join(f"{v:.4f}" for v in tpm[g]) + "\n")

# Pfam aggregation: gene -> pfam family (from pfam.tblout col1=gene, col3=pfam name, col2=acc)
pfam_map = {}
ptbl = "$PFAM_TBL"
if os.path.isfile(ptbl):
    seen = set()
    with open(ptbl) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p = line.split()
            gene, fam = p[0], p[2]
            key = (gene, fam)
            if key in seen: continue
            seen.add(key)
            pfam_map.setdefault(gene, []).append(fam)

from collections import defaultdict
def aggregate(map_g2f, outpath):
    agg = defaultdict(lambda: [0.0]*n)
    for g, fams in map_g2f.items():
        if g not in tpm: continue
        for f in fams:
            for i,v in enumerate(tpm[g]):
                agg[f][i] += v
    with open(outpath,"w") as fh:
        fh.write("feature\t" + "\t".join(samples) + "\n")
        for k in sorted(agg):
            fh.write(k + "\t" + "\t".join(f"{v:.4f}" for v in agg[k]) + "\n")
    return agg

pfam_agg = aggregate(pfam_map, os.path.join(OUT,"pfam/pfam_tpm.tsv"))

# KOfam aggregation: kofam_mapper.tsv = gene<TAB>KO (sparse)
ko_map = {}
ktsv = "$KOFAM_TSV"
if os.path.isfile(ktsv):
    with open(ktsv) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[1]: continue
            ko_map.setdefault(p[0], []).append(p[1])

ko_agg = aggregate(ko_map, os.path.join(OUT,"kofam/ko_tpm.tsv"))

# Pathway-level: KO -> pathway via $KOFAM_KO_LIST (best-effort: tab-delim, col 'definition' or 'pathway')
# Use the KO list file from config if present, else skip pathway step
ko_list = os.environ.get("KOFAM_KO_LIST","$KOFAM_KO_LIST")
ko2path = defaultdict(set)
# Try the conventional KEGG link file: ko_pathway.list (KO\tpath:map#####)
pathlinks = [os.path.join(os.path.dirname(ko_list),"ko_pathway.list"),
             os.path.join(os.path.dirname(ko_list),"links/ko2pathway.tsv")]
found_link = next((p for p in pathlinks if os.path.isfile(p)), None)
if found_link:
    with open(found_link) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2: continue
            ko = p[0].replace("ko:","")
            pth = p[1].replace("path:","")
            if pth.startswith("map"):
                ko2path[ko].add(pth)

if ko2path:
    # mean TPM per pathway across KOs
    path_tpm = {}
    for path in {p for s in ko2path.values() for p in s}:
        kos = [k for k,ps in ko2path.items() if path in ps and k in ko_agg]
        if not kos: continue
        means = [sum(ko_agg[k][i] for k in kos)/len(kos) for i in range(n)]
        path_tpm[path] = means
    with open(os.path.join(OUT,"pathway/pathway_tpm.tsv"),"w") as fh:
        fh.write("pathway\t" + "\t".join(samples) + "\n")
        for k in sorted(path_tpm):
            fh.write(k + "\t" + "\t".join(f"{v:.4f}" for v in path_tpm[k]) + "\n")
    # Row-centered heatmap matrix
    with open(os.path.join(OUT,"pathway/pathway_tpm_rowcentered.tsv"),"w") as fh:
        fh.write("pathway\t" + "\t".join(samples) + "\n")
        for k in sorted(path_tpm):
            row = path_tpm[k]; mu = sum(row)/len(row)
            fh.write(k + "\t" + "\t".join(f"{v-mu:.4f}" for v in row) + "\n")
    print(f"  pathway matrix: {len(path_tpm)} rows")
else:
    print("  WARN no ko_pathway.list found — skipping pathway aggregation")

print(f"  Genes={len(tpm)} Pfam={len(pfam_agg)} KO={len(ko_agg)} Samples={n}")
PYEOF

echo "[$(date '+%F %T')] DONE"
