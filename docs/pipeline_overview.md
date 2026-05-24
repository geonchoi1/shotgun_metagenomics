# Pipeline overview

End-to-end shotgun metagenomics pipeline. Same upstream + 3 track split + cross-track analyses.

## Top-level layout

```
shotgun_metagenomics/
├── config/                          # DB paths + conda envs + threads
├── 00_shared/                       # reads → assembly → classify → bin
├── 01_plasmid_track/                # plasmid analysis
├── 02_mag_track/                    # MAG analysis
├── 03_unbinned_track/               # unbinned chromosomal analysis
├── 04_cross_track/                  # Mobile ARG flow + other cross-track
├── tools/                           # external standalone (oriTfinder2 etc.)
└── docs/                            # this dir
```

## 00_shared — pre-track-separation

| Step | Folder | Purpose |
|---|---|---|
| 01 | `01_read_qc/` | Illumina fastp + bowtie2 dehuman / HiFi minimap2 dehuman |
| 02 | `02_assembly/` | metaSPAdes / metaFlye |
| 03 | `03_genomad_virus/` | geNomad **default** for virus extraction |
| 04 | `04_genomad_plasmid/` | geNomad **-s 4.8 --relaxed** + 5-filter (F1-F4 + F5 rRNA) |
| 05 | `05_topology_split/` | uses user-provided `circ_frag_map.tsv` → splits plasmid into `circ/`+`frag/` |
| 06 | `06_chromosomal_extract/` | assembly − plasmid − virus = chromosomal |
| 07 | `07_mag_production/` | mapping + binning + DAS_Tool + CheckM2 + Barrnap + tRNAscan + dRep + MIMAG + GTDB-Tk |
| 08 | `08_split_binned_unbinned/` | chromosomal → MAG-binned + unbinned |

## 01_plasmid_track — plasmid annotation + mobility + clustering

| Step | Folder | Tool / purpose |
|---|---|---|
| 01 | `01_raw_fasta/` | symlink: `circ/` + `frag/` from `00_shared/05` |
| 02 | `02_drep/` | RBH 100/100 ANI + AF dereplication |
| 03 | `03_bakta/` | Bakta `--complete` on circ + default on frag |
| 04 | `04_master_orf/` | combined master FAA + GFF (`circ/`, `frag/`, `all/`) |
| 05–08 | functional broad HMM | Pfam, NCBIfam, KofamScan, eggNOG |
| 09–15 | specific annotation | AMRFinder, BacMet, VFDB, TADB, dbAPIS, dbCAN, Macrel |
| 16 | defense | DefenseFinder |
| 17–18 | BGC | antiSMASH, BiG-SCAPE |
| 19 | prophage | DBSCAN-SWA |
| 20–22 | MGE | ISEScan, IntegronFinder, ICEberg3 |
| 23–26 | mobility (4 defaults) | mob_typer (rep), MOBscan (mob), oriTfinder2 (oriT), CONJScan (MPF) |
| 27 | mobility alt | PlasmidFinder, Pfam Rep HMM etc. |
| 28 | 5-tier classify | pCONJ/pdCONJ/pMOB/pOriT/pNT |
| 30 | clustering | Track 1 (Camargo pOTU) + Track 2 (COPLA) + Track 3 (combined indirect PTU) |
| 31 | AcCNET | protein community NMI vs Track 1 (internal) + Track 3 (external) |
| 40 | quantification | minimap2 + CoverM TPM |
| 50 | host prediction | Track A (community MAG cctyper) + C (iPHoP) + D (PLSDB) |
| 60 | PLSDB lookup | Mash + NUCmer per-plasmid + Simka ecosystem NMDS |

## 02_mag_track — MAG annotation + MAG-specific tools

| Step | Folder | Tool |
|---|---|---|
| 01 | `01_raw_fasta/` | symlink: dereplicated MAG genomes from `00_shared/07` (with circ/frag if complete) |
| 02 | `02_bakta/` | Bakta `--complete` on complete-circular MAG, default on rest |
| 03 | `03_master_orf/` | combined per-MAG faa + GFF |
| 04–21 | same annotation set as plasmid | Pfam ~ ICEberg3 (18 tools) |
| 22 | MAG-specific | cctyper (CRISPR-Cas + spacer) |
| 30 | MAG metabolism | METABOLIC-G |
| 40 | abundance | CoverM (per MAG, per sample TPM) |

## 03_unbinned_track — unbinned chromosomal

| Step | Folder | Tool |
|---|---|---|
| 01 | `01_raw_fasta/` | symlink: unbinned contigs (split into circ/frag if known) |
| 02 | `02_bakta/` | Bakta |
| 03 | `03_master_orf/` | combined master |
| 04–21 | same annotation set as MAG/plasmid | Pfam ~ ICEberg3 (18 tools) |
| 22 | UB-specific | MMseqs2 LCA contig-level taxonomy |

## 04_cross_track — analyses spanning all 3 tracks

| Step | Folder | Tool |
|---|---|---|
| 01 | `01_mobile_arg/` | Step 1-7: ARG×MGE distance + CooccurrenceAffinity + module extraction + cross-source BLAST + mobility pathway + ARG network |

## Naming conventions

- **All scripts are bash** (`.sh`) or Python (`.py`) self-contained — no Snakemake/Nextflow.
- **Numbered prefixes** (01_, 02_, ...) indicate ordering. Identical numbers across tracks mean the same logical step (e.g., `05_pfam/` in all 3 tracks does the same thing on different inputs).
- **Numbers can skip** (e.g., plasmid track has 30, then 40) — leaves room for future steps without renumbering.
- **`{circ,frag}` subfolders** appear in `01_raw_fasta/`, `02_bakta/`, and `03_master_orf/` of each track.
- **`config/db_paths.sh` + `config/envs.sh`** are `source`d at the top of every script.

## Input expectations

Each project should have a working directory like:

```
my_project/
├── 00_input/
│   ├── reads/                          # *.fastq.gz
│   └── circ_frag_map.tsv               # USER provides (see docs/circ_frag_list_template.md)
├── 01_qc/                              # auto-created by pipeline
├── 02_assembly/                        # auto-created
└── ...
```

All scripts are written to use **relative paths from the project's working dir**, not from the repo. The repo is a tool set; the project is where your data lives.
