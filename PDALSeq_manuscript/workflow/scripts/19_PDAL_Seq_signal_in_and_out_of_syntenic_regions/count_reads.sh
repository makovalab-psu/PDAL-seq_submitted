#!/bin/bash
# count_reads.sh
# Count read-aligned bases overlapping a region BED.
#
# Arguments:
#   $1  reads BED      (primary mapped reads, one row per alignment, sorted by chrom/start)
#   $2  region BED
#   $3  output text file
#
# Output (comma-separated, one line):
#   total_bases,bases_in_region,region_nt
#
# total_bases     = sum of aligned reference bases across all primary mapped reads
# bases_in_region = sum of read-aligned bases that fall inside the region BED
# region_nt       = nucleotides spanned by the region BED
#
# Each read contributes its aligned reference span across regions in proportion
# to how many of its bases fall in each region. Σ over the 16 partitioning
# regions = total_bases, so per-experiment Σ fraction = 1 by construction.

set -euo pipefail

reads=$1
bed=$2
output=$3

total_bases=$(awk '{sum += $3 - $2} END {print sum+0}' "$reads")

if [ -s "$bed" ]; then
    bases_in_region=$(bedtools intersect -a "$reads" -b "$bed" -sorted \
        | awk '{sum += $3 - $2} END {print sum+0}')
    region_nt=$(awk '{sum += $3 - $2} END {print sum+0}' "$bed")
else
    bases_in_region=0
    region_nt=0
fi

echo "${total_bases},${bases_in_region},${region_nt}" > "$output"
