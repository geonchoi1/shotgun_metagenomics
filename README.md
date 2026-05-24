# shotgun_metagenomics

End-to-end shotgun metagenomics pipeline for **MAG + Unbinned chromosomal + Plasmidome** analysis. Works on both Illumina short reads (metaSPAdes) and PacBio HiFi long reads (metaFlye). Generic — runs on any shotgun metagenome FASTA.

## Status
**Phase 1**: Repo skeleton only. Scripts are added incrementally per pipeline section.

| Phase | Section | Status |
|---|---|---|
| 1 | Repo skeleton + config + tools install | ⏳ in-progress |
| 2 | `00_shared/01_read_qc` + `02_assembly` | pending |
| 3 | `00_shared/03-05_genomad + topology` | pending |
| 4 | `00_shared/06-08_chromosomal + MAG production` | pending |
| 5 | `01_plasmid_track/01-04` (raw + drep + Bakta + master_orf) | pending |
| 6 | `01_plasmid_track/05-22` annotation | pending |
| 7 | `01_plasmid_track/23-28` mobility typing | pending |
| 8 | `01_plasmid_track/30-31` clustering + AcCNET | pending |
| 9 | `01_plasmid_track/40-60` quantification + host + PLSDB | pending |
| 10 | `02_mag_track` + `03_unbinned_track` | pending |
| 11 | `04_cross_track/01_mobile_arg` | pending |
| 12 | `docs/` + final polish | pending |

## Pipeline overview

```
Raw reads (Illumina or HiFi)
      ↓ 00_shared/01_read_qc (fastp + dehuman bowtie2 / minimap2)
      ↓ 00_shared/02_assembly (metaSPAdes / metaFlye)
Assembly contigs (with user-provided circ/frag list)
      ↓ 00_shared/03_genomad_virus (geNomad default)
      ↓ 00_shared/04_genomad_plasmid (-s 4.8 --relaxed + 5-filter)
      ↓ 00_shared/05_topology_split (uses user circ_frag_map.tsv)
      ↓ 00_shared/06_chromosomal_extract (assembly − plasmid − virus)
      ↓ 00_shared/07_mag_production (binning + dRep + MIMAG + GTDB-Tk)
      ↓ 00_shared/08_split_binned_unbinned
              │
       ┌──────┼──────────┐
       ↓      ↓          ↓
  plasmid  MAG     Unbinned chromosomal
   track   track       track
       └────────┬─────┘
                ↓
       04_cross_track (Mobile ARG: ARG×MGE + plasmid↔MAG↔UB cross-source HGT)
```

## Repository layout

See `docs/pipeline_overview.md` for full description. Top level:

- `config/` — environment-specific paths and conda env names
- `00_shared/` — pre-track-separation steps (reads → assembly → classify → bin)
- `01_plasmid_track/` — plasmid annotation, mobility typing, clustering, host, PLSDB
- `02_mag_track/` — MAG annotation, CRISPR-Cas, metabolism, abundance
- `03_unbinned_track/` — unbinned chromosomal contig annotation
- `04_cross_track/` — analyses that span all 3 tracks (Mobile ARG flow)
- `tools/` — external standalone tools (oriTfinder2 etc., installed locally via `tools/install_*.sh`)
- `docs/` — pipeline documentation, input formats, reference papers

## Quick start (once scripts exist)

```bash
# 1. Edit config
vim config/db_paths.sh        # set DB paths for your environment
vim config/envs.sh             # confirm conda env names
vim config/threads.sh          # set THREADS default

# 2. Stage input
mkdir -p ../my_project/00_input
cp /path/to/reads/*.fastq.gz ../my_project/00_input/
echo -e "contig_id\ttopology" > ../my_project/circ_frag_map.tsv
# user fills in circ_frag_map.tsv from assembly metadata (metaFlye assembly_info.txt etc.)

# 3. Run shared upstream (read type-dependent)
export READ_TYPE=hifi   # or 'illumina'
bash 00_shared/01_read_qc/${READ_TYPE}_dehuman_*.sh
bash 00_shared/02_assembly/${READ_TYPE}_*.sh
bash 00_shared/03_genomad_virus/run_genomad_default.sh
bash 00_shared/04_genomad_plasmid/*.sh
bash 00_shared/05_topology_split/split_from_user_list.sh
bash 00_shared/06_chromosomal_extract/extract_chromosomal.sh
bash 00_shared/07_mag_production/{01..10}_*.sh
bash 00_shared/08_split_binned_unbinned/*.sh

# 4. Run track pipelines
bash 01_plasmid_track/run_all.sh
bash 02_mag_track/run_all.sh
bash 03_unbinned_track/run_all.sh

# 5. Cross-track analysis
bash 04_cross_track/01_mobile_arg/run_all.sh
```

## Required databases

All databases assumed under `/mnt/nas/DB/geon/` by default (configurable via `config/db_paths.sh`):

- `genomad_db/` — geNomad reference
- `mob_suite/` — MOB-suite
- `mobscan_db/MOBfamDB` — MOBscan HMM
- `plasmidfinder_db2/` — PlasmidFinder
- `conjscan_models/` — MacSyFinder CONJScan
- `pfam/Pfam-A.hmm`, `pfam_hallmark/` — Pfam
- `ncbifam/hmm_PGAP.LIB` — NCBIfam
- `kofam_db/` — KofamScan
- `eggnog_db/` — eggNOG-mapper
- `amrfinder_db/latest` — AMRFinderPlus
- `bacmet_db/`, `vfdb_db/`, `tadb_db/`, `dbAPIS/`, `dbcan_db/` — annotation DBs
- `defense-finder-models-v3.1/` — DefenseFinder
- `gtdbtk_db/release232/` — GTDB-Tk
- `checkm2_db/` — CheckM2
- `iphop_db/Jun_2025_pub_rw/` — iPHoP
- `plsdb/` — PLSDB
- `bakta_db/db/` — Bakta
- `rfam/` — Rfam CMs (for 5-filter F5 rRNA scan)
- `iceberg3/` — ICEberg3

See `config/db_paths.sh` for defaults.

## External standalone tools

Installed via `tools/install_*.sh` (not committed to git):

- `tools/oriTfinder2_linux/` — install via `tools/install_oritfinder2.sh`
- `tools/ares_arroyo_oriT/` — install via `tools/install_ares_arroyo.sh`

## License

TBD
