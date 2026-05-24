#!/usr/bin/env python3
"""
=== ALT for mob (relaxase) ===
Extracts MOB system from CONJScan best_solution.tsv (model_fqn endswith "/MOB").
Cross-check vs MOBscan. CONJScan integrates relaxase + MPF context.
"""
import os, argparse, re
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--conjscan_tsv", default="../outputs/04_mpf_conjscan/result/best_solution.tsv")
ap.add_argument("--out", default="../outputs/02_mob_conjscan/contig_mob.tsv")
args = ap.parse_args()
os.makedirs(os.path.dirname(args.out), exist_ok=True)

def sanitize(s): return re.sub(r"[^A-Za-z0-9]", "_", s)

contig_mob = defaultdict(set)
with open(args.conjscan_tsv) as f:
    for line in f:
        if line.startswith("#") or line.startswith("replicon\t") or not line.strip(): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5: continue
        replicon, model = parts[0], parts[4]
        if model.endswith("/MOB"):
            contig_mob[replicon].add(model)

with open(args.out, "w") as fout:
    for c, ms in contig_mob.items():
        for m in ms:
            fout.write(f"{c}\t{m}\n")
print(f"  CONJScan MOB-system plasmid (sanitized contig names): {len(contig_mob)}")
print("  NB: replicon names are CONJScan-sanitized; remap to original via sanitize().")
