# circ_frag_map.tsv format

User-provided 2-column TSV mapping each assembly contig to its topology. This file is required by `00_shared/05_topology_split/split_from_user_list.sh` to split plasmid contigs into `circ/` (circular complete) vs `frag/` (linear fragmented).

## Format

```
contig_id<TAB>topology
```

- `contig_id`: must match the contig identifiers in the assembly FASTA (and downstream geNomad output).
- `topology`: one of `circ` or `linear`.

## Example

```
IN|contig_11059	circ
IN|contig_11447	circ
IN|contig_42177	linear
Anaerobic|contig_10318	linear
Anaerobic|contig_35686	circ
EF|contig_37903	linear
...
```

## How to generate

### For PacBio HiFi (metaFlye)

`assembly_info.txt` (one per sample) has a `circ.` column (Y/N). One-liner:

```bash
for s in IN Anaerobic Anoxic Oxic RAS EF; do
  tail -n +2 assembly_metaflye/$s/assembly_info.txt | awk -v s=$s -v OFS='\t' '
    { topo = ($4 == "Y") ? "circ" : "linear"; print s"|"$1, topo }'
done > circ_frag_map.tsv
```

### For Illumina (metaSPAdes)

metaSPAdes rarely emits circular flags. Most contigs will be `linear`. If you have specific circularization evidence (e.g., from `scaffolds.fasta` with `_loop` markers, or from external circularization tools), use it; otherwise mark all as `linear`:

```bash
grep "^>" assembly/scaffolds.fasta | sed 's/^>//;s/ .*//' | awk -v OFS='\t' '{print $1, "linear"}' > circ_frag_map.tsv
```

### General

```bash
grep "^>" assembly/contigs.fasta | sed 's/^>//;s/ .*//' | awk -v OFS='\t' '
  { topo = (FILENAME ~ /circ/) ? "circ" : "linear"; print $1, topo }' > circ_frag_map.tsv
```

## Notes

- Contig IDs must be **exact match** to what geNomad and downstream tools see. If your assembly was renamed (e.g., per-sample prefixed with `Sample|`), the IDs in this file must match too.
- This file is referenced by `05_topology_split/split_from_user_list.sh` to write `circ/all.fna` and `frag/all.fna` outputs that all per-track `01_raw_fasta/{circ,frag}/` symlink to.
- Contigs in the FASTA but not in this list are treated as `linear` by default (warning emitted).
