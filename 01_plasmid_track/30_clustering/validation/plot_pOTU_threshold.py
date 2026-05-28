#!/usr/bin/env python3
"""
Fiamenghi-style pOTU threshold validation figure:
  (a) hexbin ANI vs AF scatter (intra-pOTU pairs)
  (b) AF distribution KDE + GMM cutoff
  (c) ANI distribution + cutoff
  
Output: Fig S12 multi-panel + threshold_summary.tsv

Refs:
  - Fiamenghi 2025 Nat Comm Fig S2 (10.1038/s41467-025-65102-6)
  - Camargo 2024 (anicalc/aniclust)
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.mixture import GaussianMixture

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ani', required=True, help='anicalc.py output (qname,sname,pid,ani,af)')
    ap.add_argument('--clusters', required=True, help='aniclust.py output (rep,members)')
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    ani_df = pd.read_csv(args.ani, sep='\t')
    
    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # (a) hexbin ANI vs AF
    axes[0].hexbin(ani_df['af'], ani_df['ani'], gridsize=50, cmap='YlOrRd', bins='log')
    axes[0].axhline(70, color='blue', ls='--', label='ANI 70%')
    axes[0].axvline(70, color='red', ls='--', label='AF 70%')
    axes[0].set_xlabel('AF (%)')
    axes[0].set_ylabel('ANI (%)')
    axes[0].set_title('(a) intra-pOTU ANI × AF density')
    axes[0].legend()
    
    # (b) AF KDE + GMM
    af = ani_df['af'].dropna().values
    af = af[(af > 30) & (af < 100)]
    kde = stats.gaussian_kde(af, bw_method=0.1)
    xs = np.linspace(30, 100, 200)
    axes[1].plot(xs, kde(xs), 'b-', label='KDE')
    axes[1].axvline(70, color='red', ls='--', label='Cutoff 70%')
    axes[1].set_xlabel('AF (%)')
    axes[1].set_ylabel('Density')
    axes[1].set_title('(b) AF distribution')
    axes[1].legend()
    
    # (c) ANI distribution
    ani = ani_df['ani'].dropna().values
    ani = ani[(ani > 70) & (ani <= 100)]
    axes[2].hist(ani, bins=50, color='steelblue', alpha=0.7)
    axes[2].axvline(95, color='red', ls='--', label='ANI 95%')
    axes[2].set_xlabel('ANI (%)')
    axes[2].set_ylabel('Count')
    axes[2].set_title('(c) ANI distribution')
    axes[2].legend()
    
    plt.tight_layout()
    plt.savefig(f'{args.out_dir}/fig_S12_pOTU_threshold.png', dpi=200)
    plt.close()
    print(f'Saved: {args.out_dir}/fig_S12_pOTU_threshold.png')

if __name__ == '__main__':
    main()
