#!/usr/bin/env python3
"""
=== Step 5: 5-tier mobility classification ===
Consumes the 4 DEFAULT tool outputs:
  rep   = mob_typer report     → outputs/01_rep_mob_typer/mobtyper_report.tsv
  mob   = MOBscan              → outputs/02_mob_mobscan/contig_mob.tsv
  orit  = oriTfinder2          → outputs/03_orit_oritfinder2/output/orit_positive_contigs.txt
  MPF   = CONJScan             → outputs/04_mpf_conjscan/result/best_solution.tsv

5-tier logic (Coluzzi 2022 + Ares-Arroyo 2023):
  pCONJ   = full MPF (CONJScan T4SS_type*) + MOB
  pdCONJ  = degraded MPF (CONJScan dCONJ_type*) + MOB
  pMOB    = MOB only
  pOriT   = oriT only (no MOB, no MPF)
  pNT     = none
"""
import os, re, argparse
from collections import defaultdict, Counter

ap = argparse.ArgumentParser()
ap.add_argument("--input_fna", default="inputs/dereplicated.fna")
ap.add_argument("--circ_fna", default=None, help="Optional: circular subset FASTA")
ap.add_argument("--mobtyper_tsv", default="outputs/01_rep_mob_typer/mobtyper_report.tsv")
ap.add_argument("--mobscan_tsv", default="outputs/02_mob_mobscan/contig_mob.tsv")
ap.add_argument("--orit_positive", default="outputs/03_orit_oritfinder2/output/orit_positive_contigs.txt")
ap.add_argument("--conjscan_tsv", default="outputs/04_mpf_conjscan/result/best_solution.tsv")
ap.add_argument("--out_tsv", default="outputs/05_5tier/mobility_5tier.tsv")
args = ap.parse_args()
os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)

def sanitize(s): return re.sub(r"[^A-Za-z0-9]", "_", s)

all_p = set()
with open(args.input_fna) as f:
    for line in f:
        if line.startswith(">"): all_p.add(line[1:].split()[0])
print(f"[1] Universe: {len(all_p)}")

circ = set()
if args.circ_fna and os.path.exists(args.circ_fna):
    with open(args.circ_fna) as f:
        for line in f:
            if line.startswith(">"): circ.add(line[1:].split()[0])
    print(f"    Circular subset: {len(circ)}  (∩ universe = {len(circ & all_p)})")

mob_set = defaultdict(set)
with open(args.mobscan_tsv) as f:
    for line in f:
        c, m = line.rstrip("\n").split("\t")[:2]
        mob_set[c].add(m.replace("T4SS_", "").replace("profile_", ""))
print(f"[3] MOB (MOBscan): {len(mob_set)}")

rep_set = defaultdict(set)
with open(args.mobtyper_tsv) as f:
    first = f.readline()
    if first.lower().startswith(("sample_id", "id\t")): pass
    else:
        p = first.rstrip("\n").split("\t")
        if len(p) >= 6 and p[5] != "-":
            for r in p[5].split(","): rep_set[p[0]].add(r)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 6 and p[5] != "-":
            for r in p[5].split(","): rep_set[p[0]].add(r)
print(f"[4] Rep (mob_typer): {len(rep_set)}")

orit_set = set()
if os.path.exists(args.orit_positive):
    with open(args.orit_positive) as f:
        for line in f: orit_set.add(line.strip())
print(f"[5] oriT-positive: {len(orit_set)}")

sanitized_to_real = {sanitize(c): c for c in all_p}
full_mpf  = defaultdict(set)
deg_mpf   = defaultdict(set)
conj_mob  = defaultdict(set)
with open(args.conjscan_tsv) as f:
    for line in f:
        if line.startswith("#") or line.startswith("replicon\t") or not line.strip(): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 5: continue
        real = sanitized_to_real.get(p[0]); model = p[4]
        if not real: continue
        if   "/T4SS_type"  in model: full_mpf[real].add(model)
        elif "/dCONJ_type" in model: deg_mpf[real].add(model)
        elif model.endswith("/MOB"): conj_mob[real].add(model)
print(f"[6] CONJScan: full_MPF={len(full_mpf)}  dCONJ={len(deg_mpf)}  MOB_via_CONJ={len(conj_mob)}")

tier, detail = {}, {}
for c in sorted(all_p):
    has_mob  = (c in mob_set) or (c in conj_mob)
    has_full = c in full_mpf
    has_deg  = c in deg_mpf
    has_ori  = c in orit_set
    if   has_full and has_mob: t = "pCONJ"
    elif has_deg  and has_mob: t = "pdCONJ"
    elif has_mob:              t = "pMOB"
    elif has_ori:              t = "pOriT"
    else:                      t = "pNT"
    tier[c] = t
    detail[c] = (has_mob, has_full, has_deg, has_ori)

with open(args.out_tsv, "w") as f:
    f.write("contig\ttopology\ttier\tmob\tfull_mpf\tdeg_mpf\torit\tmob_families\trep_families\torit_check\tfull_mpf_types\tdeg_mpf_types\n")
    for c in sorted(all_p):
        topo = "circular" if c in circ else "linear"
        hm, hf, hd, ho = detail[c]
        mob_f = ",".join(sorted(mob_set.get(c, set()))) or "-"
        rep_f = ",".join(sorted(rep_set.get(c, set()))) or "-"
        full_t = ",".join(sorted(m.split("/")[-1] for m in full_mpf.get(c, set()))) or "-"
        deg_t  = ",".join(sorted(m.split("/")[-1] for m in deg_mpf.get(c, set())))  or "-"
        f.write(f"{c}\t{topo}\t{tier[c]}\t{int(hm)}\t{int(hf)}\t{int(hd)}\t{int(ho)}\t{mob_f}\t{rep_f}\t{'Y' if ho else '-'}\t{full_t}\t{deg_t}\n")

def report(name, refs):
    n = len(refs)
    if n == 0: return
    cnt = Counter(tier[x] for x in refs if x in tier)
    print(f"\n=== {name} (n={n}) ===")
    for t in ["pCONJ","pdCONJ","pMOB","pOriT","pNT"]:
        c = cnt.get(t,0); print(f"  {t:<8}{c:>6}  {100*c/n:>5.1f}%")
    mob = sum(cnt.get(t,0) for t in ["pCONJ","pdCONJ","pMOB","pOriT"])
    print(f"  {'Mobile':<8}{mob:>6}  {100*mob/n:>5.1f}%")

report("ALL putative", all_p)
if circ: report("COMPLETE circular ∩ universe", circ & all_p)
print(f"\nOutput: {args.out_tsv}")
