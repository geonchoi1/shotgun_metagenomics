# shotgun_metagenomics

End-to-end shotgun metagenomics pipeline producing **MAGs + Unbinned chromosomal gene catalog + Plasmidome** from a single project, with a downstream cross-track Mobile ARG analysis. Works on both Illumina short reads (metaSPAdes) and PacBio HiFi long reads (metaFlye). Generic — runs on any shotgun metagenome FASTA.

---

## Pipeline at a glance

```
Raw reads (Illumina or HiFi)
      │
      ▼ 00_shared/01_read_qc/        fastp + bowtie2 dehuman  (Illumina)
      │                              minimap2 dehuman          (HiFi)
      ▼ 00_shared/02_assembly/       metaSPAdes / metaFlye
      │
      ▼ 00_shared/03_genomad_virus/  geNomad default            → VIRUS
      ▼ 00_shared/04_genomad_plasmid/ -s 4.8 --relaxed + F1-F5  → PLASMID
      │
      ▼ 05_topology_split/           uses user-provided circ_frag_map.tsv
      ▼ 06_chromosomal_extract/      assembly − plasmid − virus = chromosomal
      ▼ 07_mag_production/           mapping → binning → DAS_Tool → CheckM2
      │                              → barrnap + tRNAscan → dRep species (95%)
      │                              → dRep strain (99% supp)
      │                              → MIMAG classify → GTDB-Tk
      ▼ 08_split_binned_unbinned/    chromosomal → MAG-binned + UNBINNED
      ▼ 09_kraken2_community/        read-level Kraken2 + Bracken (covers everything)
      │
      ├───────────┬────────────┬────────────────┐
      ▼           ▼            ▼                │
01_plasmid   02_mag       03_unbinned           │
   track       track          track             │
      │           │            │                │
      └───────────┴────────────┴────────────────┘
                        ▼
              04_cross_track/01_mobile_arg/
        ARG × MGE co-localization + plasmid↔MAG↔UB
              cross-source HGT network
```

## Top-level layout

```
shotgun_metagenomics/
├── run.sh                  # master orchestrator (runs all 5 sections)
├── run_initial.sh          # 00_shared
├── run_plasmid.sh          # 01_plasmid_track
├── run_mag.sh              # 02_mag_track
├── run_unbinned.sh         # 03_unbinned_track
├── run_cross.sh            # 04_cross_track
├── config.sh               # single source of truth: DB paths + conda envs + threads
├── environment.yml         # conda env for analysis Python scripts
├── requirements.txt        # pip deps
├── setup_envs.sh           # create all tool-specific conda envs (idempotent)
├── README.md
├── LICENSE
├── .gitignore
├── 00_shared/              # 18 step scripts (reads → assembly → classify → bin)
├── 01_plasmid_track/       # 27 step scripts (plasmid pipeline)
├── 02_mag_track/           # 25 step scripts (MAG pipeline)
├── 03_unbinned_track/      # 22 step scripts (UB annotation)
├── 04_cross_track/         # 9 step scripts (Mobile ARG)
└── tools/                  # external standalone tools (install scripts only)
    ├── install_oritfinder2.sh
    └── install_ares_arroyo.sh
```

Each numbered subfolder is a self-contained pipeline step with its own `run.sh` (or per-tool script) that:
- `source`s `config.sh` at the repo root
- accepts `$PROJECT` env var pointing to a per-project working dir
- assumes input from the prior step's standard output path under `$PROJECT/`
- writes to its own `$PROJECT/<section>/<step>/` output dir
- is idempotent (skips if expected output exists)
- can be re-run independently if a prior step's input is available

## Prerequisites

### Reference databases (default under `/mnt/nas/DB/geon/`, override via `config/db_paths.sh`)

| DB | Path | Used by |
|---|---|---|
| GRCh38 (human) | `GRCh38.p14/` | dehuman (Illumina bowtie2 + HiFi minimap2) |
| geNomad | `genomad_db/` | virus + plasmid classification |
| Rfam CMs | `rfam/` | 5-filter F5 (rRNA exclusion) |
| Bakta | `bakta_db/db/` | core ORF annotation |
| Pfam-A | `pfam/Pfam-A.hmm` | Pfam-A HMM search |
| Pfam hallmark | `pfam_hallmark/` | plasmid hallmark HMMs (5 markers) |
| NCBIfam | `ncbifam/hmm_PGAP.LIB` | NCBIfam HMM |
| KofamScan | `kofam_db/` | KEGG KO annotation |
| eggNOG | `eggnog_db/` | functional/COG annotation |
| dbAPIS | `dbAPIS/` | anti-defense HMM |
| AMRFinder | `amrfinder_db/latest` | AMR |
| BacMet | `bacmet_db/` | metal/biocide resistance |
| VFDB | `vfdb_db/` | virulence |
| TADB | `tadb_db/` | toxin-antitoxin |
| dbCAN | `dbcan_db/` | CAZyme |
| DefenseFinder models | `defense-finder-models-v3.1/` | bacterial defense systems |
| ISEScan | (conda env DB) | insertion sequences |
| IntegronFinder | (conda env DB) | integrons |
| ICEberg3 | `iceberg3/ICE_combined` | ICE/IME |
| MOB-suite | `mob_suite/` | rep typing |
| MOBscan | `mobscan_db/MOBfamDB` | relaxase HMM |
| PlasmidFinder | `plasmidfinder_db2/` | replicon (alt) |
| CONJScan | `conjscan_models/` | MPF (T4SS) |
| COPLA | `copla_install/COPLA/` | direct PTU |
| iPHoP | `iphop_db/Jun_2025_pub_rw/` | host prediction (CRISPR spacer DB) |
| PLSDB | `plsdb/` | per-plasmid lookup + ecosystem positioning |
| MMseqs2 GTDB | `mmseqs2_db/gtdbAA_DB` | UB contig taxonomy |
| GTDB-Tk | `gtdbtk_db/release232/` | MAG taxonomy |
| CheckM2 | `checkm2_db/` | MAG completeness/contamination |

See `config/db_paths.sh` for full list and override pattern.

### External standalone tools (not committed)

Installed locally via `tools/install_*.sh`:

| Tool | Install script | Default location |
|---|---|---|
| oriTfinder2 | `tools/install_oritfinder2.sh` | `$HOME/tools/oriTfinder2_linux/` |
| Ares-Arroyo 91-oriT BLAST DB | `tools/install_ares_arroyo.sh` | `$HOME/tools/ares_arroyo_oriT/` |

### Conda envs (default names; override via `config/envs.sh`)

`fastp`, `bowtie2`, `coverm` (with minimap2), `spades`, `metaflye`, `bakta`, `infernal`, `kofamscan`, `eggnog-mapper`, `amrfinderplus`, `dbcan`, `macrel`, `defense-finder`, `antismash`, `bigscape`, `phage` (DBSCAN-SWA + barrnap), `isescan`, `integronfinder`, `metabinner`, `metadecoder`, `semibin`, `metabat2`, `das_tool`, `checkm2`, `drep`, `gtdbtk`, `mob_suite`, `plasmidfinder`, `macsyfinder`, `copla`, `mummer4`, `simka`, `iphop`, `cctyper`, `metabolic`, `mmseqs2`, `r-base`.

---

## Inputs the user must provide

Each project has its own working directory (the repo is the **tool set**, not the data location):

```
my_project/
├── 00_input/
│   ├── reads/                        # *.fastq.gz (per-sample)
│   └── circ_frag_map.tsv             # USER-provided: contig→topology
├── 01_qc/                            # auto-created by pipeline
├── 02_assembly/                      # auto-created
└── ...
```

### `circ_frag_map.tsv` format

Two-column TSV (header required):

```
contig_id	topology
IN|contig_11059	circ
IN|contig_11447	circ
IN|contig_42177	linear
Anaerobic|contig_10318	linear
...
```

- `contig_id` must match assembly FASTA headers (and downstream geNomad output)
- `topology` ∈ {`circ`, `linear`}

#### Generation example (PacBio HiFi metaFlye)

`assembly_info.txt` has `circ.` column (Y/N):

```bash
for s in SAMPLES; do
  tail -n +2 assembly/$s/assembly_info.txt | \
    awk -v s=$s -v OFS='\t' '{ t=($4=="Y")?"circ":"linear"; print s"|"$1, t }'
done > circ_frag_map.tsv
```

#### For Illumina (metaSPAdes)

metaSPAdes rarely emits circular flags. Default everything to `linear`:

```bash
grep "^>" assembly/contigs.fasta | sed 's/^>//;s/ .*//' | \
  awk -v OFS='\t' '{print $1, "linear"}' > circ_frag_map.tsv
```

---

## Track structure (consistent across all 3 tracks)

Each track follows the **ORF-first** pattern:

```
<track>/
├── run.sh                            # track-level orchestrator
├── 01_raw_fasta/                     # symlink target: circ/, frag/ subfolders
├── 02_drep/  (plasmid only)          # RBH 100/100 dereplication
├── 03_bakta/                         # ORF prediction first
│   ├── circ_complete.sh              # Bakta --complete on circular
│   └── frag_default.sh               # Bakta default on linear
├── 04_master_orf/                    # combined master FAA + GFF
│   └── circ/, frag/, all/            # 3 subfolders
├── 05–22/                            # annotation tools (all run on master.faa)
└── (track-specific 23+)              # mobility/clustering/host/PLSDB (plasmid), cctyper/METABOLIC-G/CoverM/gene-abundance (MAG), MMseqs2-LCA (UB)
```

**Key design fix**: ORF prediction (Bakta) happens **first** at step 02-03 (after dRep for plasmid), then `master_orf` at 04 with `circ/`, `frag/`, `all/` subfolders. All downstream annotation tools (04–22) read from `04_master_orf/all/` so they see every ORF exactly once.

---

## 00_shared — Shared upstream

| Step | Folder | Tool / Purpose |
|---|---|---|
| 01 | `01_read_qc/` | Illumina: `fastp` trim → `bowtie2` dehuman vs GRCh38<br>HiFi: `minimap2 -ax map-hifi` dehuman |
| 02 | `02_assembly/` | `metaSPAdes` (Illumina) / `metaFlye --meta --pacbio-hifi` (HiFi) |
| 03 | `03_genomad_virus/` | geNomad **default** → virus contigs |
| 04 | `04_genomad_plasmid/` | geNomad **-s 4.8 --relaxed --enable-score-calibration** → 5-filter:<br>F1 score ≥ 0.7, F2 FDR < 0.05, F3 hallmark ≥ 1, F4 USCG ≤ 1, **F5 no rRNA (Infernal cmscan vs Rfam — separate step)** |
| 05 | `05_topology_split/` | Uses user-provided `circ_frag_map.tsv` → plasmid `circ/` + `frag/` |
| 06 | `06_chromosomal_extract/` | assembly − plasmid − virus = chromosomal |
| 07 | `07_mag_production/` | Mapping (minimap2 long / BWA short) → depth (jgi) → binning → DAS_Tool → CheckM2 → Barrnap rRNA → tRNAscan-SE → dRep species (95%, main) + strain (99%, supp) → MIMAG → GTDB-Tk |
| 08 | `08_split_binned_unbinned/` | chromosomal contigs → MAG-binned set + unbinned set |
| 09 | `09_kraken2_community/` | Kraken2 + Bracken on clean reads → community-level taxonomy matrix (sample × species/genus). Complements MAG GTDB-Tk by covering reads that did not bin |

### Binner choice

- **HiFi (long reads)**: MetaBinner + MetaDecoder + SemiBin2 (Han 2025 top-3) → DAS_Tool consensus
- **Illumina (short reads)**: MetaBAT2 default + user can add SemiBin2/MaxBin2/VAMB (project-specific)

See `00_shared/07_mag_production/03_binner/README.md`.

---

## 01_plasmid_track — Plasmid analysis

| Step | Folder | Tool |
|---|---|---|
| 01 | `01_raw_fasta/` | symlink: `circ/`+`frag/` from `00_shared/05` |
| 02 | `02_drep/` | RBH 100/100 ANI + AF dereplication (Fiamenghi 2025 method) |
| 03 | `03_bakta/` | Bakta `--complete` on circ + default on frag |
| 04 | `04_master_orf/` | combined master FAA + GFF + FFN; `circ/`, `frag/`, `all/` |
| 05–08 | broad HMM | Pfam, NCBIfam, KofamScan, eggNOG-mapper |
| 09 | `09_amrfinder/` | AMRFinderPlus `--plus` (AMR + stress + biocide) |
| 10–12 | DIAMOND | BacMet (metal/biocide), VFDB (virulence), TADB (toxin-antitoxin) |
| 14 | `14_dbcan/` | dbCAN CAZyme |
| 15 | `15_macrel/` | Macrel AMP |
| 16 | `16_defensefinder/` | DefenseFinder bacterial defense systems |
| 17 | `17_antismash/` | antiSMASH BGC |
| 18 | `18_bigscape/` | BiG-SCAPE GCF clustering + MIBiG novelty |
| 19 | `19_dbscan_swa/` | DBSCAN-SWA prophage |
| 20–22 | MGE | ISEScan, IntegronFinder v2, ICEberg3 BLAST |
| 23 | `23_mob_typer/` | **rep typing (default)** — MOB-suite mob_typer |
| 24 | `24_mobscan/` | **relaxase / mob (default)** — MOBscan HMM (9 MOB families) |
| 25 | `25_oritfinder2/` | **oriT (default)** — oriTfinder2 16-way parallel (split + parallel + finalize) |
| 26 | `26_conjscan/` | **MPF (default)** — MacSyFinder CONJScan/Plasmids gembase |
| 27 | `27_mobility_typing_alt/` | Alt tools: PlasmidFinder, Pfam Rep HMM, NCBIfam Rep, Pfam MOB (PF03432), CONJScan MOB, Ares-Arroyo oriT BLAST, mob_typer columns |
| 28 | `28_5tier_classification/` | pCONJ / pdCONJ / pMOB / pOriT / pNT |
| 30 | `30_clustering/` | Track 1 Camargo BLAST+Leiden pOTU + Track 2 COPLA + Track 3 indirect PTU (combined) + PTU co-membership validation |
| 31 | `31_accnet/` | Protein community network: internal (own set) + external (combined with reference) → NMI vs Track 1 + Track 3 |
| 32 | `32_functional_comparison/` | PCA + GSEA (gseapy KEGG) + Fisher per environment, both richness + TPM-weighted (Fiamenghi 2025 main + supplementary) |
| 40 | `40_quantification/` | minimap2 + CoverM (TPM/RPKM/mean) |
| 50 | `50_host_prediction/` | Track A (community MAG cctyper) + Track C (iPHoP CRISPR direct BLAST) + Track D (PLSDB Mash + NUCmer); pOTU host range aggregation |
| 60 | `60_plsdb_lookup/` | Per-plasmid Mash + NUCmer validation; Simka ecosystem NMDS |

### 5-tier mobility classification logic

```
pCONJ   = full MPF (CONJScan T4SS_type*) + MOB (MOBscan relaxase)
pdCONJ  = degraded MPF (CONJScan dCONJ_type*) + MOB
pMOB    = MOB only
pOriT   = oriT (oriTfinder2) + no MOB + no MPF
pNT     = none
```

Reference: Coluzzi & Rocha 2022 NAR (10.1093/nar/gkac1079); Ares-Arroyo et al. 2023 NAR (10.1093/nar/gkad084).

---

## 02_mag_track — MAG analysis

| Step | Folder | Tool |
|---|---|---|
| 01 | `01_raw_fasta/` | symlink: dereplicated MAGs from `00_shared/07` (`circ/`+`frag/` if any complete circular) |
| 02 | `02_bakta/` | Bakta on each MAG (`--complete` for circular MAGs, default for linear) |
| 03 | `03_master_orf/` | combined master with `circ/`, `frag/`, `all/`, `per_mag/` subfolders |
| 04–21 | annotation | **Reduced set vs plasmid** — Pfam, NCBIfam, KofamScan, eggNOG, AMRFinder, dbCAN, DBSCAN-SWA, ISEScan, IntegronFinder, ICEberg3 (BacMet/VFDB/TADB/dbAPIS/Macrel/DefenseFinder/antiSMASH/BigSCAPE skipped — plasmid only) |
| 22 | `22_cctyper/` | MAG-specific: CRISPRCasTyper (CRISPR-Cas system + spacer extraction; spacers feed Track A of plasmid host prediction) |
| 30 | `30_metabolic_g/` | MAG-specific: METABOLIC-G per-MAG × pathway matrix + N/S/C/P cycling heatmap + biogeochemical diagram |
| 40 | `40_coverm/` | Per-MAG per-sample TPM matrix (`coverm genome` with `--min-read-aligned-percent 75 --min-read-percent-identity 95`) |
| 41 | `41_gene_abundance/` | Gene-level abundance pipeline (minimap2 → featureCounts → TPM normalization → HMM-profile aggregation → KEGG pathway-level abundance → row-centered heatmap); reference: Method 2 from doi:10.1186/s40793-026-00892-w |

---

## 03_unbinned_track — Unbinned chromosomal

| Step | Folder | Tool |
|---|---|---|
| 01 | `01_raw_fasta/` | symlink: unbinned chromosomal contigs (`circ/`+`frag/`) |
| 02 | `02_bakta/` | Bakta |
| 03 | `03_master_orf/` | combined master |
| 04–21 | annotation | **Same set as MAG track 04-21** (Pfam, NCBIfam, KofamScan, eggNOG, AMRFinder, dbCAN, DBSCAN-SWA, ISEScan, IntegronFinder, ICEberg3) |
| 22 | `22_cctyper/` | CRISPRCasTyper on UB contigs (parallel to MAG cctyper) |
| 22b | `22_mmseqs2_lca_taxonomy/` | UB-specific: MMseqs2 `easy-taxonomy --tax-lineage 1` LCA against GTDB r214 (Priest 2025 method) for contig-level taxonomy |

---

## 04_cross_track — Cross-track analyses

### `01_mobile_arg/` — ARG mobilization network

Spans all 3 tracks. 7-step pipeline (Zheng 2026 framework, doi:10.1186/s40168-025-02297-2):

1. **MGE annotation** — uses each track's `19_isescan/`, `20_integronfinder/`, `21_iceberg3/` outputs
2. **ARG × MGE presence/absence matrix** per source (plasmid / MAG / UB)
3. **CooccurrenceAffinity** (R, MLE α > 2; together ≥ 10; individual ≥ 20)
4. **ARG-MGE coordinate distance** (Bakta GFF) — 3-tier: ±2 kb strict / ±5 kb standard / ±10 kb broad
5. **Module extraction** — ARG + ±5 kb MGE + ±2 kb flanking
6. **Cross-source BLAST** — modules vs (plasmid ∪ MAG ∪ UB) combined DB; 3-tier identity (≥95% recent / ≥99% very recent / 100% clonal); filter aln/qlen ≥ 0.95, length ≥ 1000 bp, self-hit removed
7. **Mobility pathway classification** — plasmid↔plasmid (cross-pOTU) / plasmid↔MAG / plasmid↔UB / MAG↔MAG
8. **ARG mobility network output** — Cytoscape / Gephi compatible

---

## Quick start

```bash
# 1) Configure (edit paths + envs for your environment)
vim config.sh                     # all DB paths + conda env names + threads

# 2) Create conda envs (idempotent — skips existing)
bash setup_envs.sh                # all envs
# or selective:
bash setup_envs.sh bakta mob_suite

# 3) Install Python analysis env
conda env create -f environment.yml      # creates `shotgun-analysis`
conda activate shotgun-analysis
pip install -r requirements.txt

# 4) Install external standalone tools
bash tools/install_oritfinder2.sh
bash tools/install_ares_arroyo.sh

# 5) Stage your project inputs
mkdir -p ~/my_project/00_input/reads
cp /path/to/reads/*.fastq.gz ~/my_project/00_input/reads/
# Provide circ_frag_map.tsv (see format above)
vim ~/my_project/00_input/circ_frag_map.tsv

# 6) Run the full pipeline
export READ_TYPE=hifi             # or "illumina"
bash /path/to/shotgun_metagenomics/run.sh ~/my_project
```

### Run a single section (skip the rest)

```bash
# Initial only (reads → MAG production)
bash /path/to/shotgun_metagenomics/run_initial.sh  ~/my_project

# Plasmid track only
bash /path/to/shotgun_metagenomics/run_plasmid.sh  ~/my_project

# MAG track only
bash /path/to/shotgun_metagenomics/run_mag.sh      ~/my_project

# Unbinned chromosomal only
bash /path/to/shotgun_metagenomics/run_unbinned.sh ~/my_project

# Cross-track (Mobile ARG) — requires all 3 tracks done first
bash /path/to/shotgun_metagenomics/run_cross.sh    ~/my_project
```

### Skip individual sections in master run.sh

```bash
SKIP_INITIAL=1 SKIP_UNBINNED=1 bash run.sh ~/my_project
```

### Run a single step

```bash
export PROJECT=~/my_project
bash /path/to/shotgun_metagenomics/01_plasmid_track/23_mob_typer/run.sh
```

Every step's `run.sh` is self-contained and re-runnable as long as its required upstream outputs are present.

---

## Build status

| Section | Scripts | Status |
|---|---|---|
| `00_shared/` | 18 | ✅ |
| `01_plasmid_track/` | 27 | ✅ |
| `02_mag_track/` | 25 | ✅ |
| `03_unbinned_track/` | 22 | ✅ |
| `04_cross_track/` | 9 | ✅ |
| Root orchestrators + config + setup | 7 | ✅ |
| **Total** | **108** scripts | **all implemented** |

Each is fully implemented, idempotent, parameterizable via `config.sh`. Some tool integrations may need per-environment tweaking (e.g. BiG-SCAPE v1 vs v2, MetaBinner kmer profile path, custom conda env names).

---

## Reference papers (key references for paper writeup)

| Topic | Paper | DOI |
|---|---|---|
| Plasmid identification | geNomad — Camargo 2024 Nat Biotechnol | 10.1038/s41587-023-02053-7 |
| Plasmidome benchmark | Fiamenghi 2025 Nat Commun (GSPR) | 10.1038/s41467-025-65102-6 |
| Plasmid mobility framework | Coluzzi 2022 NAR | 10.1093/nar/gkac1079 |
| Plasmid mobility (modern review) | Coluzzi & Rocha 2025 Sci Adv | 10.1126/sciadv.adk5811 |
| oriT framework | Ares-Arroyo 2023 NAR | 10.1093/nar/gkad084 |
| PTU framework | Redondo-Salvo 2020 Nat Commun | 10.1038/s41467-020-17278-2 |
| Plasmid defense reservoir (WWTP) | Zheng 2026 Microbiome | 10.1186/s40168-025-02297-2 |
| ARG islands | Wang & Dagan 2024 NAR Genom Bioinform | 10.1093/nargab/lqae003 |
| BSI plasmidome (clustering, NMI) | Lipworth 2024 Nat Commun | 10.1038/s41467-024-45761-7 |
| WWTP plasmid communities | Smyth 2025 npj AMR | 10.1038/s44259-025-00151-x |
| HiFi MAG | Kim 2022 Nat Commun | 10.1038/s41467-022-34149-0 |
| HiFi binner benchmark | Han 2025 Nat Commun | 10.1038/s41467-025-57957-6 |
| Lineage-resolved MAGs (HiFi+Hi-C) | Bickhart 2022 Nat Biotechnol | 10.1038/s41587-022-XXX |
| MAG gene abundance | Lemos 2026 Environ Microbiome | 10.1186/s40793-026-00892-w |
| AcCNET protein network | Lanza 2017 Mol Biol Evol | 10.1093/molbev/msac115 |

---

## License

TBD
