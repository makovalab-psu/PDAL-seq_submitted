# PDAL-Seq analysis code

Analysis code for:

> Sieg J, Zeng H, Pálová H, Smeds L, Mohanty S, Lahnsteiner A, Chiaromonte F, Makova KD.
> **Single-stranded and non-canonical DNA formation in human and other ape cells with
> telomere-to-telomere genomes.** *Proceedings of the National Academy of Sciences* (2026).
> DOI: TODO

PDAL-Seq maps single-stranded DNA (ssDNA), including ssDNA in non-canonical (non-B) DNA
structures, in living cells. ssDNA is fixed in cells by permanganate oxidation, footprinted in
vitro with S1 nuclease, and the resulting breaks are sequenced. This repository contains the
code used to process PDAL-Seq data from 15 cell lines of seven ape species mapped to their
telomere-to-telomere (T2T) genomes, and to make the figures and tables in the manuscript.

The tag [`v1.0-submitted`](../../tree/v1.0-submitted) marks the code as originally
submitted. `main` is the final version, including the analyses added during peer review.

## Repository layout

```
.
├── Process_alignments/        Workflow 1: T2T ape Cactus alignment → syntenic-block BEDs
│   ├── workflow/              Snakefile, rules/, scripts/
│   └── resources/             chromosome sizes, renaming tables, alignment provenance
└── PDALSeq_manuscript/        Workflow 2: PDAL-Seq processing and all manuscript analyses
    ├── PDALseq_samples.txt    sample sheet (one row per sequencing library)
    ├── resources/Epigenetic_datasets/BG4_Rloop_CHIP_CutTag_experiments.txt
    └── workflow/              Snakefile, rules/ (one file per analysis step),
                               scripts/ (one folder per step), env/ (conda environments)
```

The two workflows are run separately. `Process_alignments` writes one BED of syntenic
alignment blocks per species and taxonomic grouping. `PDALSeq_manuscript` reads these as
`resources/alignments/{taxon}/{species}.bed`.

| Workflow | What it does | Documentation |
|---|---|---|
| `Process_alignments` | Filters the 8-genome T2T ape Cactus alignment (MAF) to one-to-one syntenic blocks for four taxonomic groupings (hominoidae, hominidae, homininae, hominini) | [Process_alignments/README.md](Process_alignments/README.md) |
| `PDALSeq_manuscript` | Trimming, mapping and deduplication of PDAL-Seq reads; comparison with non-B DNA motifs; ape and human Gaussian HMMs; functional, repeat, satellite, replication-origin, TSS, methylation and cell-cycle analyses | [PDALSeq_manuscript/README.md](PDALSeq_manuscript/README.md), [figure and dataset map](PDALSeq_manuscript/FIGURES.md) |

## Quick start

Both workflows use [Snakemake](https://snakemake.readthedocs.io) with conda and were run on a
Linux server. They need the reference genomes and annotations listed in
[PDALSeq_manuscript/README.md](PDALSeq_manuscript/README.md#external-inputs), about 290 GB in total.

```bash
git clone https://github.com/makovalab-psu/PDAL-seq_submitted.git
cd PDAL-seq_submitted/PDALSeq_manuscript
# put inputs under resources/ and raw reads under PDALseq_data/Raw_fastq/ (see README)
snakemake -n --use-conda                      # dry run
snakemake --use-conda --cores 20 --resources mem_mb=120000
```

## Data availability

- PDAL-Seq and TruSeq sequencing reads from this study: TODO (SRA BioProject accession)
- Public datasets reanalyzed here are listed with their accessions in
  [PDALSeq_manuscript/README.md](PDALSeq_manuscript/README.md#external-inputs).

## Citation

If you use this code, please cite the paper above. GitHub's "Cite this repository" button
uses [CITATION.cff](CITATION.cff).

## License

[MIT](LICENSE)
