#!/usr/bin/env python3
"""
=== ALT for oriT ===
Extracts oriT column from mob_typer report (col 11 = orit_type).
Very few hits expected — mob_typer's oriT DB is small.
"""
import os, argparse
ap = argparse.ArgumentParser()
ap.add_argument("--mobtyper_tsv", default="../outputs/01_rep_mob_typer/mobtyper_report.tsv")
ap.add_argument("--out", default="../outputs/03_orit_mob_typer/contigs_with_orit.txt")
args = ap.parse_args()
os.makedirs(os.path.dirname(args.out), exist_ok=True)

contigs = set()
with open(args.mobtyper_tsv) as f:
    for line in f:
        if line.lower().startswith(("sample_id", "id\t")): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12: continue
        if parts[11] and parts[11] != "-":
            contigs.add(parts[0])

with open(args.out, "w") as fout:
    for c in sorted(contigs): fout.write(c + "\n")
print(f"  mob_typer oriT-positive: {len(contigs)}")
