#!/usr/bin/env python3
"""
=== ALT for mob (relaxase) ===
Extracts relaxase column from mob_typer report (col index 7 = relaxase_type).
Cross-check vs MOBscan HMM. mob_typer uses DIAMOND so it can miss diverged environmental relaxases.
"""
import os, argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--mobtyper_tsv", default="../outputs/01_rep_mob_typer/mobtyper_report.tsv")
ap.add_argument("--out", default="../outputs/02_mob_mob_typer_relaxase/contig_mob.tsv")
args = ap.parse_args()
os.makedirs(os.path.dirname(args.out), exist_ok=True)

contig_mob = defaultdict(set)
with open(args.mobtyper_tsv) as f:
    for line in f:
        if line.lower().startswith(("sample_id", "id\t")): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8: continue
        contig = parts[0]
        relaxase = parts[7]
        if relaxase and relaxase != "-":
            for r in relaxase.split(","):
                contig_mob[contig].add(r)

with open(args.out, "w") as fout:
    for c, fams in contig_mob.items():
        for f in fams: fout.write(f"{c}\t{f}\n")
print(f"  mob_typer relaxase-positive plasmid: {len(contig_mob)}")
