#!/bin/bash
# make_strand_bigwigs.sh
# Generate forward and reverse strand bigwig files from a BAM using deeptools bamCoverage.
#
# Arguments:
#   $1  BAM file (position-sorted, indexed)
#   $2  output forward bigwig (+ strand reads)
#   $3  output reverse bigwig (- strand reads)
#   $4  number of threads (default: 4)
#
# Notes:
#   - bamCoverage --samFlagExclude 16 selects reads mapping to + genomic strand
#     (flag 16 NOT set).
#   - bamCoverage --samFlagInclude 16 selects reads mapping to - genomic strand
#     (flag 16 set).
#   - bamCoverage writes bigwig directly — no intermediate bedgraph, no sorting step,
#     no bedGraphToBigWig. This avoids the bedtools genomecov -strand - malformed output
#     and out-of-order chromosome issues seen with some species BAMs.
#   - Coverage normalized to CPM (counts per million mapped reads) for comparability
#     across experiments with different library sizes.
#   - binSize 10 bp matches the computeMatrix --binSize in downstream profiling.

set -euo pipefail

bam=$1
out_forward=$2
out_reverse=$3
threads=${4:-4}

mkdir -p "$(dirname "$out_forward")"

# Forward strand: reads with flag 16 NOT set (maps to + genomic strand)
bamCoverage \
    --bam "$bam" \
    --outFileName "$out_forward" \
    --outFileFormat bigwig \
    --samFlagExclude 16 \
    --binSize 10 \
    --normalizeUsing CPM \
    -p "$threads"

# Reverse strand: reads with flag 16 set (maps to - genomic strand)
bamCoverage \
    --bam "$bam" \
    --outFileName "$out_reverse" \
    --outFileFormat bigwig \
    --samFlagInclude 16 \
    --binSize 10 \
    --normalizeUsing CPM \
    -p "$threads"
