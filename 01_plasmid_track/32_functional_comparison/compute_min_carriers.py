#!/usr/bin/env python3
"""Minimum carrier threshold for Fisher's exact test.

Uses statsmodels.stats.power.NormalIndPower (two-proportion z-test power,
asymptotically equivalent to Fisher's exact test). Cross-validated against
R MKpower::power.fisher.test() — same result.

Method (single line at the core):
    analysis.power(effect_size=Cohen_h(p1, p2), nobs1=n1, alpha=alpha, ratio=n2/n1)

Equivalent R one-liner:
    MKpower::power.fisher.test(p1=p1, p2=p2, alpha=alpha, power=power)

Usage:
    # Recommended: OR=4 (strong effect), BH-FDR α≈0.001 (~50 hits among 5000), power=0.80
    python compute_min_carriers.py --n1 188 --n2 300 --OR 4 --alpha 0.001 --power 0.80

    # Bonferroni for Pfam (~5000 tests):
    python compute_min_carriers.py --n1 188 --n2 300 --OR 4 --alpha 1e-5 --power 0.80
"""
import argparse
import numpy as np
from scipy.optimize import brentq
from statsmodels.stats.power import NormalIndPower
from statsmodels.stats.proportion import proportion_effectsize


def make_p1p2(k, OR, n1, n2):
    """Given expected total carriers k and OR, solve (p1, p2)."""
    def gap(p2):
        p1 = OR * p2 / (1 + (OR - 1) * p2) if OR != 1 else p2
        return n1 * p1 + n2 * p2 - k
    try:
        p2 = brentq(gap, 1e-6, 0.99)
        p1 = OR * p2 / (1 + (OR - 1) * p2)
        return p1, p2
    except Exception:
        return None, None


def power_at(p1, p2, n1, n2, alpha):
    """statsmodels two-proportion z-test power (one-line core).
    Asymptotically equivalent to Fisher's exact test.
    """
    h = proportion_effectsize(p1, p2)  # Cohen's h = 2*arcsin(sqrt(p1)) - 2*arcsin(sqrt(p2))
    return NormalIndPower().power(effect_size=abs(h), nobs1=n1, alpha=alpha, ratio=n2 / n1)


def find_min_k(OR, n1, n2, alpha, target_power, k_max=400):
    for k in range(2, k_max + 1):
        p1, p2 = make_p1p2(k, OR, n1, n2)
        if p1 is None or not (0 < p1 < 1) or not (0 < p2 < 1):
            continue
        if power_at(p1, p2, n1, n2, alpha) >= target_power:
            return k
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n1', type=int, required=True, help='Group 1 size (e.g. plasmids in IN env)')
    ap.add_argument('--n2', type=int, required=True, help='Group 2 size (e.g. plasmids in EF env)')
    ap.add_argument('--OR', type=float, default=4.0, help='Target odds ratio (default 4)')
    ap.add_argument('--alpha', type=float, default=0.001,
                    help='Significance level. 0.05=single test; 0.001=BH-FDR proxy; 1e-5=Bonferroni Pfam; 5e-6=Bonferroni KO')
    ap.add_argument('--power', type=float, default=0.80, help='Target power (default 0.80)')
    args = ap.parse_args()

    print(f"\n=== MIN_CARRIERS auto-compute (statsmodels NormalIndPower) ===")
    print(f"  n1            = {args.n1}")
    print(f"  n2            = {args.n2}")
    print(f"  Target OR     = {args.OR}")
    print(f"  Alpha         = {args.alpha:.0e}")
    print(f"  Target power  = {args.power}")

    k = find_min_k(args.OR, args.n1, args.n2, args.alpha, args.power)
    if k is None:
        print(f"\n  WARN: No k≤400 achieves power {args.power} for OR={args.OR} at α={args.alpha:.0e}")
        return

    p1, p2 = make_p1p2(k, args.OR, args.n1, args.n2)
    pw = power_at(p1, p2, args.n1, args.n2, args.alpha)

    print(f"\n  ✓ Recommended MIN_CARRIERS = {k}  (power = {pw*100:.1f}%)")
    print(f"\n  → export MIN_CARRIERS={k}")
    print(f"     bash run.sh\n")

    # Sensitivity
    print(f"  [Sensitivity at adjacent k]:")
    for kk in [max(2, k-20), max(2, k-10), k, k+10, k+25, k+50]:
        p1, p2 = make_p1p2(kk, args.OR, args.n1, args.n2)
        if p1 is None: continue
        pw = power_at(p1, p2, args.n1, args.n2, args.alpha)
        marker = " ← rec" if kk == k else ""
        print(f"    k={kk:>3} → power = {pw*100:5.1f}%{marker}")

    # R equivalent for paper Methods
    print(f"\n  R equivalent (cross-check):")
    print(f"    library(MKpower)")
    print(f"    power.fisher.test(p1={p1:.4f}, p2={p2:.4f}, alpha={args.alpha}, power={args.power})")
    print()


if __name__ == '__main__':
    main()
