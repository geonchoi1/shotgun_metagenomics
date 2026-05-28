# pOTU validation

Two-axis validation that the Camargo Leiden pOTU clustering produces biologically meaningful units (not just sequence-similarity noise). Implementation lives in `30_clustering/run.sh` (post-hoc PTU propagation) and `32_functional_comparison/run.sh`; this README documents the rationale.

---

## 1) PTU co-membership analysis (taxonomy-level external validation)

**Input**: `30_clustering/track3/combined_clusters.tsv` — Track 3 output where our plasmids and PLSDB PTU reference plasmids are co-clustered by Camargo aniclust.

**Logic** (per cluster, member-list parsing):
```
For each cluster, count how many PTU-reference plasmids it contains:
  → PureRef         : cluster contains only PLSDB ref plasmids (no ours)
  → PureNovel       : cluster contains only our plasmids (0 PTU ref)
  → Mixed_single_PTU: cluster has ≥1 our + ≥1 ref, all ref share the same PTU
  → Mixed_multi_PTU : cluster has ≥1 our + ≥2 ref from different PTUs (bridging)
  → Mixed_unknown_PTU: cluster has our + ref, but the ref plasmids have no PTU label
```

**Example**:
```
cluster_X = {ours_1, ours_2, ours_3, NC_IncF-001, NC_IncF-002}
  → Mixed_single_PTU (IncF) — propagate IncF label to ours_1..3

cluster_Y = {ours_4, ours_5}
  → PureNovel — our-only family, no known PTU

cluster_Z = {ours_6, ours_7, NC_IncP-ref, NC_IncQ-ref}
  → Mixed_multi_PTU (IncP+IncQ) — bridges two known PTUs; ambiguous label
```

**Two simultaneous outcomes**:
1. **External validation**: a high *Mixed_single_PTU* + *PureRef* signal means Camargo Leiden clustering reproduces paper-validated PTU taxonomy → "not random noise, real plasmid-family units."
2. **Label acquisition**: COPLA-unassigned fragmented plasmids inherit an **indirect PTU label** via shared cluster with a labeled reference — recovers PTU coverage that direct COPLA misses (especially on environmental, non-circular contigs).

**No extra DB/tool**: pure post-processing of `combined_clusters.tsv`. Done in `30_clustering/run.sh` step "POST-HOC PTU label propagation" and emits `30_clustering/track3/our_plasmid_indirect_ptu.tsv` + per-cluster classification.

---

## 2) Functional category enrichment (PCA) — internal validation

Already produced by `32_functional_comparison/run.sh`; no extra step:
- Per-pOTU Pfam / KOfam annotation count
- Filter: ≥100 total count, present in ≥2 environments
- Normalize by # plasmids per sample
- StandardScaler + PCA `fit_transform` → 2D embedding

If pOTU members / environment compartments separate visually in PC1×PC2 → clusters carry coherent gene function (not just nucleotide similarity).

---

## 3) Joint interpretation

| Axis | Source | Meaning |
|------|--------|---------|
| PTU co-membership | Sequence/taxonomy (PLSDB reference) | Cluster aligns with known plasmid family — external |
| Functional PCA | Pfam/KOfam gene content | Cluster has coherent gene repertoire — internal |

Both positive → the pOTU is a biologically meaningful unit at both sequence-similarity AND functional levels.

If only PTU co-membership is positive: clusters are taxonomically valid but functionally heterogeneous → may need finer-grained subclustering.

If only Functional PCA is positive: clusters share function but bridge taxa → suggests horizontal gene capsule or convergent function rather than vertical inheritance — interesting biology but not a "plasmid taxon."

---

## Related Fig S panels

- Fig S3 — see `02_drep/validation/plot_drep_collapse.py` (1881→drep set 100/100 collapse)
- Fig S12 — see `plot_pOTU_threshold.py` (1869→1233 pOTU 70% AF cutoff KDE/hexbin)
- Fig 3b alluvial — PTU co-membership classes as ribbon colors
