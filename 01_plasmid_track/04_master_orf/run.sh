#!/bin/bash
# === Combine Bakta circ+frag FAA/GFF/FFN into master ORF set ===
# Outputs:
#   $PROJECT/01_plasmid_track/04_master_orf/{circ,frag,all}/master.{faa,gff,ffn}
#   $PROJECT/01_plasmid_track/04_master_orf/orf2contig.tsv      (locus_tag -> contig)
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

BAKTA=$PROJECT/01_plasmid_track/03_bakta
OUT=$PROJECT/01_plasmid_track/04_master_orf
mkdir -p $OUT/circ $OUT/frag $OUT/all

copy_one() {
  local src_pre=$1  # e.g. $BAKTA/circ/plasmid_circ
  local dst=$2      # e.g. $OUT/circ
  for ext in faa gff3 ffn; do
    if [ -s ${src_pre}.${ext} ]; then
      out_ext=$([ "$ext" = "gff3" ] && echo "gff" || echo "$ext")
      cp -f ${src_pre}.${ext} $dst/master.${out_ext}
    fi
  done
}

copy_one $BAKTA/circ/plasmid_circ  $OUT/circ
copy_one $BAKTA/frag/plasmid_frag  $OUT/frag

# Combined all/master.* (frag first, then circ — matches existing convention)
cat $OUT/frag/master.faa  $OUT/circ/master.faa  > $OUT/all/master.faa
cat $OUT/frag/master.ffn  $OUT/circ/master.ffn  > $OUT/all/master.ffn
{ cat $OUT/frag/master.gff; grep -v '^#' $OUT/circ/master.gff; } > $OUT/all/master.gff

# Build orf2contig.tsv from GFF (locus_tag <- ID/locus_tag attr; contig <- column 1)
python3 - <<PYEOF
import os, re
OUT = "$OUT"
def parse(gff, fout):
    with open(gff) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'CDS': continue
            contig = p[0]
            attrs = dict(kv.split('=',1) for kv in p[8].split(';') if '=' in kv)
            lt = attrs.get('locus_tag') or attrs.get('ID','').split('-')[-1]
            if lt:
                fout.write(f"{lt}\t{contig}\n")
with open(f"{OUT}/orf2contig.tsv","w") as fout:
    for sub in ("circ","frag"):
        gff = f"{OUT}/{sub}/master.gff"
        if os.path.exists(gff): parse(gff, fout)
PYEOF

echo "[$(date '+%F %T')] DONE — all ORFs: $(grep -c '^>' $OUT/all/master.faa), orf2contig: $(wc -l < $OUT/orf2contig.tsv)"
