#!/usr/bin/env python3
"""
Camargo dRep collapse validation (Fig S3):
  (a) Collapsed pair count per dRep cluster (histogram)
  (b) BLAST identity distribution of collapsed pairs (should be ~100%)
  (c) AF (HSP-union) distribution of collapsed pairs (should be ~100% under
      Camargo formula vs simple-sum which dilutes under repeats)
  (d) Length-difference distribution of collapsed pairs (validates that the
      rep is representative)

Output: Fig S3 multi-panel + collapse_summary.tsv per cluster.

This is DIFFERENT from pOTU threshold validation (30_clustering/validation/):
  - 02_drep: validates 100/100 cutoff applied to putative (1881 -> drep set)
  - 30_clustering: validates the 70% AF biological cutoff between species-
    equivalent pOTU members (Fiamenghi-style)

Refs:
  - Fiamenghi 2025 Nat Comm Fig S2 (10.1038/s41467-025-65102-6)
  - Camargo 2024 anicalc (HSP-union AF + prune_alns)
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ani', required=True, help='anicalc output for collapsed pairs')
    ap.add_argument('--clusters', required=True, help='aniclust output (rep,members)')
    ap.add_argument('--fasta', required=True, help='input FASTA (for length info)')
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Parse cluster member multiplicity
    cluster_sizes = []
    rep_members = {}
    with open(args.clusters) as f:
        for line in f:
            rep, members = line.rstrip().split('\t')
            mems = members.split(',')
            cluster_sizes.append(len(mems))
            rep_members[rep] = mems

    # Parse FASTA lengths
    lengths = {}
    with open(args.fasta) as f:
        cur, seq = None, []
        for line in f:
            if line.startswith('>'):
                if cur:
                    lengths[cur] = sum(len(s) for s in seq)
                cur = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if cur:
            lengths[cur] = sum(len(s) for s in seq)

    # Parse ANI of collapsed pairs (only those in same cluster)
    pair_ani, pair_af, pair_lendiff = [], [], []
    rep_of = {}
    for rep, mems in rep_members.items():
        for m in mems:
            rep_of[m] = rep

    with open(args.ani) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) < 5:
                continue
            q, s = parts[0], parts[1]
            try:
                ani = float(parts[3])
                af = float(parts[4])
            except ValueError:
                continue
            if rep_of.get(q) == rep_of.get(s) and q != s:
                pair_ani.append(ani)
                pair_af.append(af)
                if q in lengths and s in lengths:
                    pair_lendiff.append(abs(lengths[q] - lengths[s]))

    # Plot
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))

    # (a) cluster size histogram
    multi = [c for c in cluster_sizes if c > 1]
    if multi:
        axes[0].hist(multi, bins=range(2, max(multi) + 2),
                     color='steelblue', alpha=0.8, rwidth=0.8)
    axes[0].set_xlabel('Cluster size (members)')
    axes[0].set_ylabel('# clusters')
    axes[0].set_title(f'(a) Multi-member clusters (n={len(multi)})')
    axes[0].set_yscale('log')

    # (b) BLAST identity
    if pair_ani:
        axes[1].hist(pair_ani, bins=50, color='coral', alpha=0.8)
        axes[1].axvline(100, color='black', ls='--', label='100% threshold')
        axes[1].set_xlabel('BLAST identity (%)')
        axes[1].set_ylabel('# pairs')
        axes[1].set_title(f'(b) Identity distribution\n(n={len(pair_ani)} collapsed pairs)')
        axes[1].set_xlim(99, 100.1)
        axes[1].legend()

    # (c) AF (HSP-union)
    if pair_af:
        axes[2].hist(pair_af, bins=50, color='seagreen', alpha=0.8)
        axes[2].axvline(100, color='black', ls='--', label='100% threshold')
        axes[2].set_xlabel('AF (HSP-union, %)')
        axes[2].set_ylabel('# pairs')
        axes[2].set_title('(c) AF distribution')
        axes[2].set_xlim(99, 100.1)
        axes[2].legend()

    # (d) length difference
    if pair_lendiff:
        axes[3].hist(pair_lendiff, bins=30, color='mediumpurple', alpha=0.8)
    axes[3].set_xlabel('|len(query) - len(subject)| (bp)')
    axes[3].set_ylabel('# pairs')
    axes[3].set_title('(d) Length difference')

    plt.tight_layout()
    plt.savefig(f'{args.out_dir}/fig_S3_drep_collapse_validation.png', dpi=200)
    plt.close()

    # Summary table
    with open(f'{args.out_dir}/collapse_summary.tsv', 'w') as o:
        o.write('cluster_rep\tn_members\trep_length\n')
        for rep, mems in rep_members.items():
            if len(mems) > 1:
                o.write(f'{rep}\t{len(mems)}\t{lengths.get(rep, "NA")}\n')

    print(f'Saved: {args.out_dir}/fig_S3_drep_collapse_validation.png')
    print(f'Saved: {args.out_dir}/collapse_summary.tsv')
    print(f'Multi-member clusters: {len(multi)}, total collapsed pairs analyzed: {len(pair_ani)}')


if __name__ == '__main__':
    main()
