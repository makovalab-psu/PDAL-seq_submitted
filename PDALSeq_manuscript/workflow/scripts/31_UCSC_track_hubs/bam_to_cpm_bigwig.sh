#!/bin/bash
# bam_to_cpm_bigwig.sh
# CPM-normalized coverage bigWig from a BAM, with chromosomes renamed to UCSC native names.
#
# Arguments:
#   $1  BAM file
#   $2  UCSC chromAlias.txt for the target assembly
#   $3  UCSC chrom.sizes for the target assembly
#   $4  output bigWig
#
# Notes:
#   - CPM = coverage x 1e6 / primary mapped reads (samtools view -F 0x904).
#   - bedtools genomecov -split gives per-base coverage, so no binning.
#   - Sequences missing from chromAlias are dropped (counted in the log).

set -euo pipefail

bam=$1
alias=$2
sizes=$3
out=$4

tmp_bg="${out%.bw}.tmp.bg"
script_dir=$(dirname "$0")

total=$(samtools view -c -F 0x904 "$bam")
echo "primary mapped reads: $total"

bedtools genomecov -bg -split -ibam "$bam" \
    | awk -v scale="$(echo "$total" | awk '{printf "%.12g", 1000000 / $1}')" 'BEGIN{OFS="\t"} {$4 = $4 * scale; print}' \
    | python "$script_dir/rename_chroms.py" "$alias" \
    | LC_ALL=C sort -k1,1 -k2,2n -S 4G -T "$(dirname "$out")" > "$tmp_bg"

bedGraphToBigWig "$tmp_bg" "$sizes" "$out"
rm "$tmp_bg"
