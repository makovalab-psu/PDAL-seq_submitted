#!/bin/bash
#Prepares the siamang chr12/chr19 rerun of Fig. S14 (run from the PDALseq_manuscript folder, with bedtools on the PATH)
#1. Moves the siamang outputs built from the RepeatMasker annotation into a backup folder, keeping their paths,
#   so Snakemake rebuilds them and the old versions stay available for comparison
#2. Rebuilds the siamang SST1/SATR1/SATR2/LSAU BEDs made by the multi-species script
#   Get_SST1_SATR1_SATR2_LSAU_annotations_from_rm.out.sh, so that rule (and the other species' files) is not rerun
#$1 backup folder (default backup/pre_chr12_chr19_fix)

set -euo pipefail
shopt -s nullglob

backup=${1:-backup/pre_chr12_chr19_fix}
swap="awk -v pair=chr12,chr19 -f workflow/scripts/swap_chrom_pair.awk"
rm_sats="SST1 SATR1 SATR2 LSAU"

stale=(
    data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/*/Ssyndactylus.bed
    data/17_PDALseq_reads_in_CenSat/06_Read_enrichment_in_specific_rm_catagories/*/Ssyndactylus_*.csv
    results/17_PDALseq_reads_in_CenSat/SST1_SATR1_SATR2_LSAU_enrichment.csv
)

####1. Move stale outputs to the backup folder####

for p in "${stale[@]}"; do
    if [ -e "$backup/$p" ]; then
        echo "ERROR: $backup/$p already exists; choose another backup folder" >&2
        exit 1
    fi
done

for p in "${stale[@]}"; do
    if [ ! -e "$p" ]; then
        echo "not found, skipped: $p" >&2
        continue
    fi
    mkdir -p "$backup/$(dirname "$p")"
    mv "$p" "$backup/$p"
    echo "moved $p"
done

####2. Rebuild the siamang annotation files (same commands as the multi-species script, plus the swap)####

for s in $rm_sats; do
    mkdir -p data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/$s
    { grep $s resources/repeats/linneas_annotations/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out || true; } \
        | awk '{print $5 "\t" $6  "\t" $7 "\t" $10}' | $swap | bedtools sort | bedtools merge \
        > data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/$s/Ssyndactylus.bed
done

echo "done; old outputs are in $backup"
