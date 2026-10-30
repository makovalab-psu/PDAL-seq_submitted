# Step 32: siamang chr12/chr19 fix

Two naming schemes exist for the siamang T2T assembly. The sequences are identical; only the chr12 and chr19 labels differ:

| Naming | chr12_hap1 | chr19_hap1 | Used by |
|---|---|---|---|
| v2.0 (GCA_028878055.2) | 98,029,838 bp (CM054520.2) | 146,713,733 bp (CM054527.2) | UCSC hub GCA_028878055.2, RepeatMasker `mSymSyn1.pri.cur.20231205` |
| v2.1 | 146,713,733 bp | 98,029,838 bp | our FASTA, BAMs, alignment BEDs, non-B calls (`mSymSyn1.pri.cur.20240514`), CenSat, centromeres |

All analyses use v2.1 names except the two RepeatMasker-based annotations, which are now swapped to v2.1 names
(`workflow/scripts/swap_chrom_pair.awk`):

| Script / rule | Change |
|---|---|
| `17_.../Get_SST1_SATR1_SATR2_LSAU_annotations_from_rm.out.sh` | siamang SST1/SATR1/SATR2/LSAU BEDs swapped |
| rule `Get_SST1_SATR1_SATR2_LSAU_annotations_from_rm_out` | now calls the script above (it called `intersect_censat_nonB.sh`) |
| `04_.../hsat_2_3_annotation.sh` | siamang HSat2/3 BED swapped (the siamang BED is not used in any figure) |
| rule `download_ucsc_genome_files_31` | swaps chr12/chr19 in the siamang chromAlias so all hub tracks land on the right chromosome |

## What gets remade

Only Fig. S14 (siamang SST1/SATR1/SATR2/LSAU rows), from `results/17_PDALseq_reads_in_CenSat/SST1_SATR1_SATR2_LSAU_enrichment.csv`
with `new_plots/Figure_S18_Othe_satellites.R`.

## Instructions (run from the PDALseq_manuscript folder on the server)

1. Pull the updated code.

2. Move the stale siamang outputs to `backup/pre_chr12_chr19_fix/` and rebuild the 4 siamang BEDs
   (needs `bedtools` on the PATH):

   ```bash
   bash workflow/scripts/32_siamang_chr12_chr19_fix/prepare_siamang_rerun.sh
   ```

3. Dry run. `workflow/siamang_fix.smk` is the main Snakefile plus step 17 (commented out in the Snakefile).
   `--rerun-triggers input` stops Snakemake from rerunning other species' jobs because the new siamang BEDs are newer.

   ```bash
   snakemake -s workflow/siamang_fix.smk -n --use-conda --rerun-triggers input \
     results/17_PDALseq_reads_in_CenSat/SST1_SATR1_SATR2_LSAU_enrichment.csv
   ```

   Expected: `Select_sat_class_enrichment` 12 jobs (4 satellites x 3 siamang experiments) and
   `Compile_select_sat_class_enrichment` 1 job. Stop if other species or `Get_SST1_SATR1_SATR2_LSAU_annotations_from_rm_out` appear.

4. Run it: the same command without `-n`, with your usual `--cores` / profile options.

5. Remake the figure:

   ```bash
   Rscript workflow/scripts/new_plots/Figure_S18_Othe_satellites.R   # new_plots/Figure_S18_select_satellite_subclasses.svg
   ```

## Track hub

Delete the siamang chromAlias so it is downloaded again with the swap; the siamang hub files then rebuild:

```bash
rm data/31_UCSC_track_hubs/00_ucsc_genome_files/GCA_028878055.2.chromAlias.txt
```

## Notes

- The read enrichment for Fig. S14 subsamples and shuffles windows at random, so the siamang values will also move a
  little from sampling alone.
- Old outputs stay in `backup/pre_chr12_chr19_fix/` for before/after comparison; delete it when no longer needed.
