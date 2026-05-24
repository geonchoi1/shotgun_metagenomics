#!/bin/bash
# === DEFAULT for oriT — Step c ===
# Merge per-workspace results + un-sanitize Acc column (Sample_contig_N → Sample|contig_N).

set -e
OUT_DIR=${OUT_DIR:-../outputs/03_orit_oritfinder2}
cd $OUT_DIR
mkdir -p output

for kind in oriT relaxase t4cp t4ss auxiliary_pro ARG VF metal degradation symbiosis anti-CRISPR matrix; do
  out=output/oriTfinder2_${kind}_summary.txt
  first=1
  for ws in workspace/ws_*; do
    f=$ws/oriTfinder2_${kind}_summary.txt
    [ -f $f ] || continue
    if [ $first -eq 1 ]; then cp $f $out; first=0
    else tail -n +2 $f >> $out; fi
  done
  echo "Merged $kind: $(wc -l < $out 2>/dev/null || echo 0)"
done

# Un-sanitize: SampleX_contig_NUMBER → SampleX|contig_NUMBER
for f in output/oriTfinder2_*_summary.txt; do
  python3 - <<PYEOF
import re
data = open("$f").read()
data2 = re.sub(r"^([A-Za-z]+)_(contig_\d+)", r"\1|\2", data, flags=re.MULTILINE)
open("$f", "w").write(data2)
PYEOF
done

awk -F'\t' 'NR>1 && $2!="-"{print $1}' output/oriTfinder2_oriT_summary.txt | sort -u > output/orit_positive_contigs.txt
echo "oriT-positive plasmid: $(wc -l < output/orit_positive_contigs.txt)"
echo "Output: $OUT_DIR/output/"
