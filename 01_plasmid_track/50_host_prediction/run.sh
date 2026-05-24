#!/bin/bash
# === Plasmid host prediction — Track A (per-MAG cctyper) + C (iPHoP CRISPR ref) + D (Mash+NUCmer) ===
# Track A: needs $PROJECT/mag/22_cctyper/per_mag/*/spacers/*.fa (or run here per MAG)
# Track C: BLAST plasmid vs iPHoP CRISPR DB
# Track D: Mash vs PLSDB sketch (D≤0.1, P≤0.1) + NUCmer validation (id≥90, len≥500)
# Integration: union → per-plasmid host range + pOTU aggregation
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

PLASMID=$PROJECT/plasmid/02_drep/dereplicated.fna
OUT=$PROJECT/plasmid/50_host_prediction
mkdir -p $OUT/trackA $OUT/trackC $OUT/trackD

PLSDB_META=${PLSDB_META:-$PLSDB_DIR/metadata/plsdb.csv}

############ TRACK A — per-MAG cctyper spacer → plasmid ############
TA=$OUT/trackA
CCTYPER_DIR=${CCTYPER_DIR:-$PROJECT/mag/22_cctyper/per_mag}
if [ -d "$CCTYPER_DIR" ] && [ ! -s $TA/spacer_plasmid_filtered.tsv ]; then
  echo "[$(date '+%F %T')] TrackA — aggregate MAG spacers"
  > $TA/all_mag_spacers.fa
  for mag_dir in $CCTYPER_DIR/*/; do
    mag=$(basename $mag_dir)
    [ -d $mag_dir/spacers ] || continue
    for sf in $mag_dir/spacers/*.fa; do
      [ -s $sf ] || continue
      sed "s/^>/>${mag}__/" $sf >> $TA/all_mag_spacers.fa
    done
  done
  activate_env "$ENV_DIAMOND"
  # cd-hit-est nr100
  if command -v cd-hit-est >/dev/null; then
    cd-hit-est -i $TA/all_mag_spacers.fa -o $TA/mag_spacers_nr100.fa -c 1.0 -n 8 -M 16000 -T $THREADS > $TA/cdhit.log 2>&1 || cp $TA/all_mag_spacers.fa $TA/mag_spacers_nr100.fa
  else
    cp $TA/all_mag_spacers.fa $TA/mag_spacers_nr100.fa
  fi
  makeblastdb -in $PLASMID -dbtype nucl -out $TA/plasmid_db
  blastn -task blastn-short -dust no -word_size 8 -max_target_seqs 1000 \
    -query $TA/mag_spacers_nr100.fa -db $TA/plasmid_db \
    -evalue 0.01 -num_threads $THREADS_BLAST \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen evalue bitscore' \
    -out $TA/spacer_plasmid.tsv
  # Fiamenghi filter: aln≥25, mm+gap≤1, aln/spacer≥0.95
  awk -F'\t' '$4>=25 && ($5+$6)<=1 && ($4/$7)>=0.95' $TA/spacer_plasmid.tsv > $TA/spacer_plasmid_filtered.tsv
  echo "  TrackA filtered hits: $(wc -l < $TA/spacer_plasmid_filtered.tsv)"
else
  echo "[$(date '+%F %T')] TrackA — skip (no MAG cctyper dir or output exists)"
fi

############ TRACK C — iPHoP CRISPR DB direct BLAST ############
TC=$OUT/trackC
if [ ! -s $TC/iphop_spacer_filtered.tsv ]; then
  echo "[$(date '+%F %T')] TrackC — BLAST plasmid vs iPHoP spacer DB"
  activate_env "$ENV_DIAMOND"
  if [ ! -f ${IPHOP_SPACER_DB}.nhr ]; then
    makeblastdb -in $IPHOP_SPACER_DB -dbtype nucl -out $IPHOP_SPACER_DB || true
  fi
  blastn -task blastn-short -dust no -word_size 8 -max_target_seqs 1000 \
    -query $IPHOP_SPACER_DB -db <(makeblastdb -in $PLASMID -dbtype nucl -out $TC/plasmid_db -parse_seqids 2>/dev/null; echo $TC/plasmid_db) \
    -evalue 0.01 -num_threads $THREADS_BLAST \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen evalue bitscore' \
    -out $TC/iphop_spacer.tsv 2>/dev/null || \
  blastn -task blastn-short -dust no -word_size 8 -max_target_seqs 1000 \
    -query $IPHOP_SPACER_DB -db $TC/plasmid_db \
    -evalue 0.01 -num_threads $THREADS_BLAST \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen evalue bitscore' \
    -out $TC/iphop_spacer.tsv
  awk -F'\t' '$4>=25 && ($5+$6)<=1 && ($4/$7)>=0.95' $TC/iphop_spacer.tsv > $TC/iphop_spacer_filtered.tsv
  echo "  TrackC filtered hits: $(wc -l < $TC/iphop_spacer_filtered.tsv)"
fi

############ TRACK D — Mash + NUCmer vs PLSDB ############
TD=$OUT/trackD
if [ ! -s $TD/mash_filtered.tsv ]; then
  echo "[$(date '+%F %T')] TrackD — Mash vs PLSDB"
  activate_env "$ENV_MASH"
  mash dist -p $THREADS -d 0.1 -v 0.1 $PLSDB_MASH_SKETCH $PLASMID > $TD/mash_all.tsv
  awk -F'\t' '$3<=0.1 && $4<=0.1' $TD/mash_all.tsv > $TD/mash_filtered.tsv
  echo "  Mash candidates: $(wc -l < $TD/mash_filtered.tsv)"
fi
if [ ! -s $TD/nucmer_validated.tsv ]; then
  echo "[$(date '+%F %T')] TrackD — NUCmer validation"
  activate_env "$ENV_NUCMER"
  # Build subset PLSDB FASTA from top Mash hits
  awk -F'\t' '{print $1}' $TD/mash_filtered.tsv | sort -u > $TD/plsdb_candidates.txt
  python3 - <<PYEOF
ids=set(open("$TD/plsdb_candidates.txt").read().split())
out=open("$TD/plsdb_subset.fna","w"); keep=False
with open("$PLSDB_FASTA") as f:
    for line in f:
        if line.startswith('>'):
            keep = line[1:].split()[0] in ids
        if keep: out.write(line)
out.close()
PYEOF
  cd $TD
  nucmer --threads $THREADS -p nucmer $TD/plsdb_subset.fna $PLASMID > nucmer.log 2>&1 || true
  delta-filter -i 90 -l 500 -1 nucmer.delta > nucmer_filtered.delta
  show-coords -rclT nucmer_filtered.delta > nucmer_validated.tsv
  echo "  NUCmer validated lines: $(wc -l < $TD/nucmer_validated.tsv)"
fi

############ INTEGRATION — union per plasmid + pOTU aggregation ############
echo "[$(date '+%F %T')] Integration"
python3 - <<PYEOF
import os, csv
from collections import defaultdict
OUT="$OUT"
PROJ="$PROJECT"
plasmid_hosts=defaultdict(lambda: defaultdict(set))  # plasmid → track → {host strings}

# TrackA: spacer → plasmid; need MAG → GTDB tax map
mag_tax={}
for fp in [f"{PROJ}/mag/08_gtdbtk/gtdbtk.bac120.summary.tsv",
           f"{PROJ}/mag/08_gtdbtk/gtdbtk.ar53.summary.tsv"]:
    if os.path.exists(fp):
        with open(fp) as fh:
            next(fh)
            for line in fh:
                p=line.rstrip('\n').split('\t')
                if len(p)>=2: mag_tax[p[0]]=p[1]
ta_file=f"{OUT}/trackA/spacer_plasmid_filtered.tsv"
if os.path.exists(ta_file):
    with open(ta_file) as f:
        for line in f:
            p=line.rstrip('\n').split('\t')
            spacer=p[0]; plasmid=p[1]
            mag=spacer.split('__')[0]
            tax=mag_tax.get(mag, f"MAG:{mag}")
            plasmid_hosts[plasmid]["A_MAG_cctyper"].add(tax)

# TrackC: iPHoP spacer hits — qid encodes spacer source taxonomy in iPHoP DB headers
tc_file=f"{OUT}/trackC/iphop_spacer_filtered.tsv"
if os.path.exists(tc_file):
    with open(tc_file) as f:
        for line in f:
            p=line.rstrip('\n').split('\t')
            sp=p[0]; plasmid=p[1]
            # iPHoP spacer header: typically <accession>__<taxonomy_tokens>
            tax=sp.split('__',1)[1] if '__' in sp else sp
            plasmid_hosts[plasmid]["C_iPHoP_CRISPR"].add(tax)

# TrackD: Mash+NUCmer-validated PLSDB hits → metadata host
plsdb_host={}
PLSDB_META="$PLSDB_META"
if os.path.exists(PLSDB_META):
    with open(PLSDB_META) as fh:
        r=csv.DictReader(fh)
        for row in r:
            acc=row.get('NUCCORE_ACC') or row.get('ACC_NUCCORE') or row.get('accession')
            host=row.get('TAXONOMY_genus') or row.get('TAXONOMY_species') or row.get('host')
            if acc and host: plsdb_host[acc]=host
validated_pairs=set()
tdv=f"{OUT}/trackD/nucmer_validated.tsv"
if os.path.exists(tdv):
    with open(tdv) as f:
        for line in f:
            if line.startswith('/') or not line.strip() or line.startswith('NUCMER'): continue
            p=line.split('\t')
            if len(p)<13: continue
            # show-coords -T: refid, queryid at last two cols
            ref=p[-2]; qry=p[-1].strip()
            validated_pairs.add((qry,ref))
for plasmid, ref in validated_pairs:
    if ref in plsdb_host:
        plasmid_hosts[plasmid]["D_Mash_NUCmer"].add(plsdb_host[ref])

# Write per-plasmid integration
with open(f"{OUT}/per_plasmid_hosts.tsv","w") as o:
    o.write("plasmid\ttrackA\ttrackC\ttrackD\tunion_hosts\n")
    for p, tracks in plasmid_hosts.items():
        union=set()
        for t in tracks.values(): union |= t
        o.write(f"{p}\t{';'.join(sorted(tracks.get('A_MAG_cctyper',[])))}\t{';'.join(sorted(tracks.get('C_iPHoP_CRISPR',[])))}\t{';'.join(sorted(tracks.get('D_Mash_NUCmer',[])))}\t{';'.join(sorted(union))}\n")

# pOTU host range aggregation
potu_file=f"{PROJ}/plasmid/30_clustering/track1/pOTU_membership.tsv"
if os.path.exists(potu_file):
    potu_hosts=defaultdict(set)
    contig2potu={}
    with open(potu_file) as f:
        next(f)
        for line in f:
            c,p=line.rstrip('\n').split('\t')
            contig2potu[c]=p
    for c, tracks in plasmid_hosts.items():
        potu=contig2potu.get(c)
        if not potu: continue
        for s in tracks.values(): potu_hosts[potu] |= s
    with open(f"{OUT}/pOTU_host_range.tsv","w") as o:
        o.write("pOTU\tn_hosts\thosts\n")
        for p,h in potu_hosts.items():
            o.write(f"{p}\t{len(h)}\t{';'.join(sorted(h))}\n")
print(f"Plasmids with any host: {len(plasmid_hosts)}")
PYEOF
echo "[$(date '+%F %T')] DONE"
