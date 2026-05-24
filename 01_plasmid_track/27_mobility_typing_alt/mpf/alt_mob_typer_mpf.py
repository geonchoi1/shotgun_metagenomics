#!/usr/bin/env python3
"""
=== ALT for MPF (T4SS) ===
Extracts mpf_type column from mob_typer report (col 9 = mpf_type).
Cross-check vs CONJScan; mob_typer's MPF DB is smaller/clinical-biased.
"""
import os, argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--mobtyper_tsv", default="../outputs/01_rep_mob_typer/mobtyper_report.tsv")
ap.add_argument("--out", default="../outputs/04_mpf_mob_typer/contig_mpf.tsv")
args = ap.parse_args()
os.makedirs(os.path.dirname(args.out), exist_ok=True)

contig_mpf = defaultdict(set)
with open(args.mobtyper_tsv) as f:
    for line in f:
        if line.lower().startswith(("sample_id", "id\t")): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10: continue
        mpf = parts[9]
        if mpf and mpf != "-":
            for m in mpf.split(","):
                contig_mpf[parts[0]].add(m)

with open(args.out, "w") as fout:
    for c, ms in contig_mpf.items():
        for m in ms: fout.write(f"{c}\t{m}\n")
print(f"  mob_typer MPF-positive plasmid: {len(contig_mpf)}")
