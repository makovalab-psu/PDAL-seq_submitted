#!/bin/bash
# prepare_nonb_bed.sh
# Prepare a non-B DNA annotation BED for strand-specific analysis.
#
# Arguments:
#   $1  input BED (raw from resources/nonb_annotations/)
#   $2  output BED
#   $3  optional pair of chromosomes to swap, e.g. chr12,chr19 (see swap_chrom_pair.awk)
#
# Behaviour:
#   - Filters out sex-chromosome features (chrX, chrY, chrW, chrZ and
#     any chromosome whose name starts with those patterns).
#   - If the input has ≥ 6 columns, outputs BED6 (preserving strand).
#   - If the input has < 6 columns, outputs BED3 (no strand info).
#   - An empty output file is valid (deeptools handles it gracefully).

set -euo pipefail

input=$1
output=$2
swap=${3:-}

mkdir -p "$(dirname "$output")"

if [ ! -s "$input" ]; then
    touch "$output"
    exit 0
fi

# Detect number of columns from first data line
ncols=$(awk 'NF>0{print NF; exit}' "$input")

if [ "$ncols" -ge 6 ]; then
    # BED6: preserve strand column
    awk 'NF>=6 {
        if ($1 ~ /^chrX/ || $1 ~ /^chrY/ || $1 ~ /^chrW/ || $1 ~ /^chrZ/) next
        print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6
    }' "$input" | awk -v pair="$swap" -f workflow/scripts/swap_chrom_pair.awk > "$output"
else
    # BED3: no strand info
    awk 'NF>=3 {
        if ($1 ~ /^chrX/ || $1 ~ /^chrY/ || $1 ~ /^chrW/ || $1 ~ /^chrZ/) next
        print $1"\t"$2"\t"$3
    }' "$input" | awk -v pair="$swap" -f workflow/scripts/swap_chrom_pair.awk > "$output"
fi
