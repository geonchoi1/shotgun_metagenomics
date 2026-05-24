#!/bin/bash
# === DEFAULT for oriT — Step b ===
# 16-way parallel oriTfinder2. Each workspace has its own scripts/tools/data symlinks +
# own ./input + own output summary files. Source tool dir is read-only.

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

PARALLEL=${PARALLEL:-16}
OUT_DIR=${OUT_DIR:-../outputs/03_orit_oritfinder2}
TOOL_DIR=${TOOL_DIR:-/home/gchoi/tools/oriTfinder2_linux}

cd $OUT_DIR
mkdir -p workspace chunks logs

rm -f chunks/list_*
split -d -n l/$PARALLEL -a 2 all_plasmid_ids.list chunks/list_
echo "Created $(ls chunks/ | wc -l) chunk files"

for i in $(seq -w 0 $((PARALLEL - 1))); do
  ws=workspace/ws_$i
  rm -rf $ws
  mkdir -p $ws/{input,tmp}
  ln -sf $TOOL_DIR/scripts $ws/scripts
  ln -sf $TOOL_DIR/tools $ws/tools
  ln -sf $TOOL_DIR/data $ws/data
  ln -sf $TOOL_DIR/run_oriTfinder2_local_batch.pl $ws/run_oriTfinder2_local_batch.pl
  chunk=chunks/list_$i
  [ -f $chunk ] || continue
  while read fname; do
    ln -sf $(realpath .)/input_gbk/$fname $ws/input/$fname
  done < $chunk
  cp $chunk $ws/list.txt
done

cat > run_one_workspace.sh <<'WSEOF'
#!/bin/bash
WS=$1
cd "$WS"
START=$(date +%s)
perl run_oriTfinder2_local_batch.pl list.txt > run.log 2>&1
echo "DONE $WS  elapsed=$(($(date +%s)-START))s  $(date '+%F %T')"
WSEOF
chmod +x run_one_workspace.sh

echo "[$(date '+%F %T')] Launching $PARALLEL parallel oriTfinder2"
ls -d workspace/ws_* | parallel -j $PARALLEL --joblog logs/parallel.joblog ./run_one_workspace.sh {}
echo "[$(date '+%F %T')] All workspaces complete"
