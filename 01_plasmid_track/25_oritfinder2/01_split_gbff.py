#!/usr/bin/env python3
"""
=== DEFAULT for oriT — Step a ===
Split Bakta gbff into per-plasmid GBK files (one per LOCUS).
Sanitizes 'Sample|contig_N' → 'Sample_contig_N.gbk' (oriTfinder2 perl tool dislikes '|').
"""
import os, re, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--circ_gbff", default="../inputs/bakta_circ.gbff",
                help="Bakta --complete output for circular plasmid")
ap.add_argument("--frag_gbff", default="../inputs/bakta_fragmented.gbff",
                help="Bakta default output for fragmented plasmid")
ap.add_argument("--keep_ids", default="../inputs/dereplicated_ids.txt",
                help="Optional: filter to keep only these original IDs ('Sample|contig_N')")
ap.add_argument("--out_dir", default="../outputs/03_orit_oritfinder2")
args = ap.parse_args()

OUT_GBK = os.path.join(args.out_dir, "input_gbk")
os.makedirs(OUT_GBK, exist_ok=True)
sanitize_re = re.compile(r"[^A-Za-z0-9_.-]")

keep = None
if args.keep_ids and os.path.exists(args.keep_ids):
    with open(args.keep_ids) as f:
        keep = set(line.strip() for line in f if line.strip())
    print(f"  ID filter: {len(keep)}")

def split_file(path):
    if not os.path.exists(path):
        print(f"  SKIP (not found): {path}"); return 0
    n = sk = 0
    cur = []; cur_id = None
    with open(path) as f:
        for line in f:
            cur.append(line)
            if line.startswith("LOCUS"):
                cur_id = line.split()[1].strip()
            elif line.startswith("//"):
                if cur_id:
                    if keep is None or cur_id in keep:
                        safe = sanitize_re.sub("_", cur_id)
                        with open(os.path.join(OUT_GBK, safe + ".gbk"), "w") as fout:
                            fout.writelines(cur)
                        n += 1
                    else: sk += 1
                cur, cur_id = [], None
    print(f"  {path}: wrote {n}, skipped {sk}"); return n

total = split_file(args.circ_gbff) + split_file(args.frag_gbff)
print(f"Total GBK files: {total}")

with open(os.path.join(args.out_dir, "all_plasmid_ids.list"), "w") as fout:
    for fn in sorted(os.listdir(OUT_GBK)):
        fout.write(fn + "\n")
print(f"ID list: {args.out_dir}/all_plasmid_ids.list")
