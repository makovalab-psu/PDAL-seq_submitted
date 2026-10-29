# Process_alignments workflow

Converts the 8-genome T2T ape Cactus alignment into BED files of one-to-one syntenic alignment
blocks, one per species, for four taxonomic groupings. These BEDs are the
`resources/alignments/{taxon}/{species}.bed` inputs of the
[PDALSeq_manuscript](../PDALSeq_manuscript/README.md) workflow.

| Grouping | Species with a block in every alignment |
|---|---|
| `hominoidae` | all 7: human, chimpanzee, bonobo, gorilla, Bornean and Sumatran orangutan, siamang |
| `hominidae` | great apes (no siamang) |
| `homininae` | human, chimpanzee, bonobo, gorilla |
| `hominini` | human, chimpanzee, bonobo |

## Inputs

| Path | Contents |
|---|---|
| `resources/alignments/8-t2t-apes-2023v2.hs1.maf` | Cactus alignment of hs1 (CHM13v2.0), hg38 and six non-human ape T2T assemblies, as MAF referenced on hs1. The file in this repository is an empty placeholder (the real file is too large for GitHub). Download: TODO. How it was made: [resources/alignments/README.md](resources/alignments/README.md) |
| `resources/alignments/Rename_species.txt`, `Rename_chromosome.txt` | Map assembly accessions and chromosome names to the names used in the PDAL-Seq workflow |
| `resources/chr_sizes/{species}.chrom.sizes` | Chromosome sizes for the 7 assemblies |

## Requirements

The rules do not use conda environments, so these must be on `PATH`: Snakemake, Python 3 with
pandas, R with tidyverse, and awk.

## Running

From this directory:

```bash
snakemake -n                 # dry run
snakemake --cores 8
```

Final outputs are `data/16_Remove_sex_chromosome_blocks/{grouping}/{species}.bed`. Copy them to
`PDALSeq_manuscript/resources/alignments/{grouping}/` to use them in the main workflow.

## Steps

| Output folder | Rule | What it does |
|---|---|---|
| `data/01_remove_hg38` | `remove_hg38` | Drop hg38 rows so human is represented only by hs1 |
| `data/02_prune_blocks` | `prune_ape_blocks` | Keep blocks ≥ 100 bp that contain every species of the grouping |
| `data/03_dup_filter` | `dup_filter` | Where a species has several sequences in a block, keep the one most similar to hs1 |
| `data/04_filter_quality` | `filter_quality` | Drop blocks in which any sequence has < 50% identity to hs1 |
| `data/05`–`07` | `rename_species`, `rename_chr` | Convert accessions and chromosome names |
| `data/08_split_by_ref_chromosomes` | `split_by_human_ref_chromosomes` | Split by hs1 chromosome |
| `data/09_sort` | `sort` | Sort blocks by hs1 position |
| `data/10`–`13` | `quantify_synteny`, `Identify_synteny`, `enforce_synteny`, `check_synteny` | For each hs1 chromosome and species, rank that species' chromosomes by number of blocks and keep the top ones that together hold 80% of blocks (`identify_synteny.R`); remove blocks that map to any other chromosome |
| `data/14_merge_maf` | `merge_maf` | Merge adjacent blocks |
| `data/15_maf_to_bg` | `maf_to_bg` | Convert blocks to per-species BED coordinates |
| `data/16_Remove_sex_chromosome_blocks` | `remove_sex_chromosome_blocks` | Remove blocks on X and Y (final output) |
| `data/check_*` | `check_alignment_stats`, `check_alignment_lengths`, `plot_alignment_lengths`, `median_alignment_lengths` | QC: block counts, length distributions, median lengths |

`data/11_identify_syteny`, `data/check_alignment_stats` and `data/check_block_lengths` are
committed so the QC plots can be viewed without rerunning the workflow.

## Minimal file checklist

Paths are relative to `Process_alignments/`. Only the MAF is missing from the repository.

```
resources/
├── alignments/
│   ├── 8-t2t-apes-2023v2.hs1.maf     # download (TODO); repository copy is an empty placeholder
│   ├── Rename_species.txt            # included
│   └── Rename_chromosome.txt         # included
└── chr_sizes/{Species}.chrom.sizes   # included, all 7 species
```
