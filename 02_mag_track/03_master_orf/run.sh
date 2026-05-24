#!/bin/bash
# === 03 Master ORF aggregation ===
# Input:  $PROJECT/mag/02_bakta/{circ,frag}/<MAG>/<MAG>.{faa,ffn,gff3}
# Output: $PROJECT/mag/03_master_orf/{circ,frag,all}/master.{faa,ffn,gff3}
#         $PROJECT/mag/03_master_orf/per_mag/<MAG>/{mag.faa,mag.ffn,mag.gff3}
#         $PROJECT/mag/03_master_orf/all/orf2contig.tsv  (locus_tag<TAB>contig<TAB>MAG)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

BAKTA_BASE=$PROJECT/mag/02_bakta
OUT=$PROJECT/mag/03_master_orf
mkdir -p "$OUT"/{circ,frag,all,per_mag}

echo "[$(date '+%F %T')] master ORF builder"

python3 - <<PYEOF
import os, glob, re, gzip
BAKTA = "$BAKTA_BASE"
OUT   = "$OUT"

groups = ["circ", "frag"]
per_group = {g: {"faa": [], "ffn": [], "gff": []} for g in groups}
o2c_rows = []

for g in groups:
    base = os.path.join(BAKTA, g)
    if not os.path.isdir(base): continue
    for mag in sorted(os.listdir(base)):
        d = os.path.join(base, mag)
        faa = os.path.join(d, f"{mag}.faa")
        ffn = os.path.join(d, f"{mag}.ffn")
        gff = os.path.join(d, f"{mag}.gff3")
        if not (os.path.isfile(faa) and os.path.isfile(gff)):
            continue
        # per-MAG copies
        pm = os.path.join(OUT, "per_mag", mag); os.makedirs(pm, exist_ok=True)
        for src, dstname in [(faa,"mag.faa"),(ffn,"mag.ffn"),(gff,"mag.gff3")]:
            dst = os.path.join(pm, dstname)
            if os.path.isfile(src) and not os.path.islink(dst):
                if os.path.exists(dst): os.remove(dst)
                os.symlink(src, dst)
        per_group[g]["faa"].append(faa)
        if os.path.isfile(ffn): per_group[g]["ffn"].append(ffn)
        per_group[g]["gff"].append(gff)

        # orf2contig parse from GFF3 (CDS lines)
        with open(gff) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip(): continue
                if line.startswith(">"): break
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9: continue
                if parts[2] != "CDS": continue
                contig = parts[0]
                attrs = dict(re.findall(r"([^=;]+)=([^;]+)", parts[8]))
                lt = attrs.get("locus_tag") or attrs.get("ID","").split("_cds")[0]
                if lt:
                    o2c_rows.append((lt, contig, mag))

def concat(files, dst):
    with open(dst, "wb") as out:
        for f in files:
            with open(f, "rb") as fh:
                out.write(fh.read())

for g in groups:
    concat(per_group[g]["faa"], os.path.join(OUT, g, "master.faa"))
    concat(per_group[g]["ffn"], os.path.join(OUT, g, "master.ffn"))
    concat(per_group[g]["gff"], os.path.join(OUT, g, "master.gff3"))

# all = circ + frag
all_faa = per_group["circ"]["faa"] + per_group["frag"]["faa"]
all_ffn = per_group["circ"]["ffn"] + per_group["frag"]["ffn"]
all_gff = per_group["circ"]["gff"] + per_group["frag"]["gff"]
concat(all_faa, os.path.join(OUT, "all", "master.faa"))
concat(all_ffn, os.path.join(OUT, "all", "master.ffn"))
concat(all_gff, os.path.join(OUT, "all", "master.gff3"))

with open(os.path.join(OUT, "all", "orf2contig.tsv"), "w") as fh:
    for lt, c, m in o2c_rows:
        fh.write(f"{lt}\t{c}\t{m}\n")

print(f"  MAGs circ={len(per_group['circ']['faa'])} frag={len(per_group['frag']['faa'])}")
print(f"  ORFs total={len(o2c_rows)}")
PYEOF

echo "[$(date '+%F %T')] DONE — outputs in $OUT"
