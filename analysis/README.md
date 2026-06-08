# analysis — figure-generation scripts (Fig 1–4)

Scripts used to produce the confirmed main figures (Fig 1–4) of the WWTP plasmidome / MAG study.
Grouped by figure; one figure may use several scripts (and a data-build script shared across panels).

## Fig 1 — Treatment train · community context · plasmidome recovery  (`fig1/`)
- `mag_tree.R` — (b) 169-MAG phylogenomic tree
- `arg_panels.py` — (c) MAG tree concentric rings (taxonomy / quality / TPM / AMR)
- `figure1c_read_taxonomy.py` — (d) read-based community phylum composition
- `build_tpm.py`, `fig1e_heatmaps.py` — (e) MAG TPM heatmap across zones
- (a) WWTP A2O schematic = drawn externally (no script)

## Fig 2 — Plasmidome catalog architecture  (`fig2/`)
- `fig2a_donuts.py` — (a) size / topology donut
- `fig2b_ourplsdb_network.py` — (b) ours + PLSDB plasmid network
- `fig_linking_barplot.py` — (c) linking barplot
- `fig2_mobility_combined.py` — (d) 5-tier mobility composition
- `fig2ef_igraph.py` — (e) MOB×MPF / replicon heatmap

## Fig 3 — Cargo landscape & functional comparison  (`fig3/`)
- `build_plasmid_cargo.py` — cargo/abundance data table (used by b, c)
- `fig3a_venn.py` — (a) cargo Venn
- `fig3b_orf_pct.py` — (b) ARG class composition
- `fig3cd_combined.py` — (c) abundance / cargo / drug-class heatmaps
- `fig_functional_tpm.py` — (d) functional GSEA across zones

## Fig 4 — ARG–MGE mobility  (`fig4/`)
- `fig4a_network.py` — (a) ARG–MGE co-occurrence network
- `cassette_sharing_pool.py` — shared-cassette detection/verification (input to b)
- `fig4b_main.py` — (b) plasmid↔chromosome cassette synteny

All figures rendered at 600 dpi. Input paths are referenced inside each script (project tracks `00_shared/`, `01_plasmid_track/`, `02_mag_track/`, `03_unbinned_track/`, `04_cross_track/`).
