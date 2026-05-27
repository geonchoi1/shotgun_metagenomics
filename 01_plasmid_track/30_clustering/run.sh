#!/bin/bash
# === Track 1 + 2 + 3 clustering + PTU label propagation ===
#
# Track 1: Camargo snakemake on dereplicated own set → pOTU (own only)
# Track 2: COPLA on circular only → direct PTU
# Track 3: Camargo snakemake on combined (ours ∪ PLSDB PTU members) → indirect PTU labels
#          via co-clustering with reference plasmids
#
# Algorithm: Camargo's contig-ani-leiden-clustering-pipeline
#   - blastn megablast (-max_target_seqs 25000 -perc_identity 0)
#   - anicalc.py: HSP union AF + pruned tANI
#   - aniclust.py: edge if (qcov OR tcov ≥ min_cov), weight = ani × max(qcov, tcov)
#   - Leiden community detection (resolution = 1)
#
# Refs:
#   - Fiamenghi 2025 Nat Comm (10.1038/s41467-025-65102-6)
#   - Camargo 2024 (github.com/apcamargo/bioinformatics-snakemake-pipelines)
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

OURS=$PROJECT/plasmid/02_drep/dereplicated.fna
CIRC=$PROJECT/plasmid/02_drep/circ.fna
OUT=$PROJECT/plasmid/30_clustering
mkdir -p $OUT/track1 $OUT/track2_copla $OUT/track3 $OUT/validation

CAMARGO_PIPELINE=${CAMARGO_PIPELINE:-$HOME/tools/bioinformatics-snakemake-pipelines/contig-ani-leiden-clustering-pipeline}
SMK=$CAMARGO_PIPELINE/contig-ani-leiden-clustering-pipeline.smk

############ TRACK 1: Camargo snakemake on OUR set ############
T1=$OUT/track1
if [ ! -s $T1/pOTU_membership.tsv ]; then
  echo "[$(date '+%F %T')] Track1 — Camargo snakemake on dereplicated own set"
  activate_env "$ENV_SNAKEMAKE_CAMARGO"
  cp $OURS $T1/own.fna
  cd $T1
  snakemake --snakefile $SMK \
    --config input=own.fna leiden_resolution=1 min_ani=0 min_cov=0.70 blast_threads=$THREADS_BLAST \
    --cores $THREADS_BLAST
  # convert clusters.tsv → pOTU_membership.tsv
  python3 - <<PYEOF
with open("own_clusters.tsv") as f, open("pOTU_membership.tsv","w") as o:
    o.write("contig\tpOTU\n")
    for i, line in enumerate(f, 1):
        rep, members = line.strip().split('\t')
        for m in members.split(','):
            o.write(f"{m}\tpOTU_T1_{i:05d}\n")
PYEOF
  cd $SCRIPT_DIR
  echo "[$(date '+%F %T')] Track1 pOTU: $(awk 'NR>1{print $2}' $T1/pOTU_membership.tsv | sort -u | wc -l)"
fi

############ TRACK 2: COPLA on circular only ############
T2=$OUT/track2_copla
if [ -s $CIRC ] && [ ! -s $T2/copla_summary.tsv ]; then
  echo "[$(date '+%F %T')] Track2 — COPLA on circular"
  activate_env "$ENV_COPLA"
  cd $COPLA_DIR
  python3 bin/copla.py $CIRC $T2 --threads $THREADS || echo "WARN: COPLA failed (often all unassigned due to fragmented data)"
  cd $SCRIPT_DIR
fi

############ TRACK 3: combined (OUR ∪ PLSDB PTU refs) → Camargo + post-hoc PTU labeling ############
T3=$OUT/track3
if [ ! -s $T3/combined_clusters.tsv ]; then
  echo "[$(date '+%F %T')] Track3 — combined Camargo snakemake"
  mkdir -p $T3
  # Build combined FASTA: ours + PLSDB PTU reference plasmids
  # PLSDB plasmid IDs do not contain '|' — easy to distinguish from our IDs (sample|contig_N)
  awk '!/^>/' $OURS > /dev/null   # validate
  cat $OURS $PLSDB_FASTA > $T3/combined.fna
  N_COMBINED=$(grep -c '^>' $T3/combined.fna)
  echo "  combined input: $N_COMBINED plasmids (ours + PLSDB)"

  activate_env "$ENV_SNAKEMAKE_CAMARGO"
  cd $T3
  snakemake --snakefile $SMK \
    --config input=combined.fna leiden_resolution=1 min_ani=0 min_cov=0.70 blast_threads=$THREADS_BLAST \
    --cores $THREADS_BLAST
  cd $SCRIPT_DIR
fi

############ POST-HOC: PTU label propagation (Track 3) ############
echo "[$(date '+%F %T')] Post-hoc PTU label propagation"
python3 - <<PYEOF
from collections import defaultdict, Counter

# Load OUR plasmid IDs
OURS_IDS = set()
with open("$OURS") as f:
    for line in f:
        if line.startswith('>'):
            OURS_IDS.add(line[1:].split()[0])

# Load PLSDB PTU table (id → PTU label) if available
PLSDB_PTU = {}
plsdb_ptu_file = "$PLSDB_PTU_TABLE"
import os
if os.path.exists(plsdb_ptu_file):
    with open(plsdb_ptu_file) as f:
        next(f)
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                PLSDB_PTU[parts[0]] = parts[1]
print(f"  PLSDB PTU labels: {len(PLSDB_PTU)}")

# Load combined cluster
clusters = defaultdict(list)
with open("$T3/combined_clusters.tsv") as f:
    for line in f:
        rep, members = line.strip().split('\t')
        for m in members.split(','):
            clusters[rep].append(m)

# Classify each cluster + propagate PTU
out_cluster = open("$OUT/validation/cluster_classification.tsv", "w")
out_cluster.write("pOTU_T3_rep\tn_total\tn_ours\tn_ref\tclass\tassigned_ptu\tref_members\n")
out_ptu = open("$OUT/validation/our_plasmid_indirect_ptu.tsv", "w")
out_ptu.write("contig\tindirect_PTU\tassignment_type\n")

n_pure_novel = n_mixed_with_ptu = n_mixed_no_ptu = n_pure_ref = 0
for rep, members in clusters.items():
    ours = [m for m in members if m in OURS_IDS]
    ref = [m for m in members if m not in OURS_IDS]
    if ours and ref:
        # Mixed cluster — propagate PTU label
        ref_ptu_labels = [PLSDB_PTU[r] for r in ref if r in PLSDB_PTU]
        if ref_ptu_labels:
            # Use most common PTU; mark Mixed if multiple
            ptu_counter = Counter(ref_ptu_labels)
            if len(ptu_counter) == 1:
                assigned = ptu_counter.most_common(1)[0][0]
                cls = "Mixed_single_PTU"
            else:
                assigned = ptu_counter.most_common(1)[0][0]
                cls = "Mixed_multi_PTU"
            for o in ours:
                out_ptu.write(f"{o}\t{assigned}\t{cls}\n")
            n_mixed_with_ptu += 1
        else:
            assigned = "no_PTU_in_cluster"
            cls = "Mixed_no_PTU_label"
            for o in ours:
                out_ptu.write(f"{o}\t{assigned}\t{cls}\n")
            n_mixed_no_ptu += 1
    elif ours and not ref:
        cls = "PureNovel"; assigned = "Novel"
        for o in ours:
            out_ptu.write(f"{o}\tNovel\tPureNovel\n")
        n_pure_novel += 1
    else:
        cls = "PureRef"; assigned = "-"
        n_pure_ref += 1
    out_cluster.write(f"{rep}\t{len(members)}\t{len(ours)}\t{len(ref)}\t{cls}\t{assigned}\t{';'.join(ref[:3])}\n")

out_cluster.close(); out_ptu.close()
print(f"  PureNovel clusters: {n_pure_novel}")
print(f"  Mixed with PTU label: {n_mixed_with_ptu}")
print(f"  Mixed without PTU: {n_mixed_no_ptu}")
print(f"  PureRef clusters: {n_pure_ref}")

# Count our plasmid PTU assignment
ptu_counts = Counter()
with open("$OUT/validation/our_plasmid_indirect_ptu.tsv") as f:
    next(f)
    for line in f:
        _, ptu, _ = line.rstrip().split('\t')
        ptu_counts[ptu] += 1
print(f"\n  Our plasmid PTU assignment summary:")
for ptu, n in ptu_counts.most_common(20):
    print(f"    {ptu}: {n}")
PYEOF

echo "[$(date '+%F %T')] DONE — see $OUT/validation/"
