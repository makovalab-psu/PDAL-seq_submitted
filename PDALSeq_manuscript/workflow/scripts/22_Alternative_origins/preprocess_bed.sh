#!/bin/bash
# Normalize an origin BED before liftOver:
#   - skip line 1 if first field is "Chr" or "chr" (header row)
#   - prepend "chr" to chromosome names that lack it (e.g. "1" -> "chr1")
#   - write exactly 3 columns: chr, start, end
# Usage: preprocess_bed.sh input.bed output.bed

input=$1
output=$2

awk 'BEGIN{OFS="\t"}
     NR==1 && $1~/^[Cc]hr$/ { next }
     { if ($1 !~ /^chr/) $1 = "chr"$1; print $1, $2, $3 }' "$input" > "$output"
