#!/bin/bash
# compute_strand_profile.sh
# Compute average strand-specific bigwig coverage around non-B DNA features.
#
# Arguments:
#   $1  forward bigwig
#   $2  reverse bigwig
#   $3  features BED (BED3 or BED6); assumed sorted
#   $4  output TSV
#   $5  number of threads (default: 4)
#
# BED3 (unstranded) features:
#   All features treated as + strand. Output columns: forward, reverse.
#   Represents + and − strand read coverage at each bin position.
#
# BED6 (strand-aware) features (e.g. GQ):
#   Features split into + and − strand subsets. Each subset is accumulated
#   separately, then combine_strands merges them with correct template/nontemplate
#   assignment:
#     + strand feature: template = reverse bigwig (displaced ssDNA on − strand)
#                       nontemplate = forward bigwig (G4-forming strand)
#     − strand feature: template = forward bigwig (displaced ssDNA on + strand)
#                       nontemplate = reverse bigwig (G4-forming strand)
#   deeptools has already flipped the window 5'→3' for − strand features, so bin
#   positions are comparable across strands before recoding. Output columns:
#   template, nontemplate.
#   Features with strand '.' are silently skipped (not assigned to either subset).
#
# Memory strategy: BED is split into 200k-feature chunks. Each chunk runs
# computeMatrix independently (~1-2 GB peak RAM), then profile_from_matrix.py
# streams the chunk matrix into a running-sum JSON without holding the full
# N_features × N_bins array.

set -euo pipefail

forward_bw=$1
reverse_bw=$2
bed=$3
output=$4
threads=${5:-4}

SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"
CHUNK_SIZE=200000

mkdir -p "$(dirname "$output")"

if [ ! -s "$bed" ]; then
    touch "$output"
    exit 0
fi

WORK_DIR=$(mktemp -d)
trap "rm -rf '$WORK_DIR'" EXIT INT TERM

CHUNK_MATRIX="$WORK_DIR/chunk.mat.gz"

# Detect if BED is strand-aware (≥6 columns)
ncols=$(awk 'NF>0{print NF; exit}' "$bed")

# run_accumulate <subset_bed> <state_file> <chunk_prefix>
# Splits <subset_bed> into CHUNK_SIZE chunks, runs computeMatrix on each, and
# accumulates the result into <state_file>. No-op if <subset_bed> is empty.
run_accumulate() {
    local subset_bed=$1
    local state_file=$2
    local prefix=$3
    [ ! -s "$subset_bed" ] && return 0

    split -l "$CHUNK_SIZE" -a 3 "$subset_bed" "${WORK_DIR}/${prefix}_"
    for chunk in "${WORK_DIR}/${prefix}_"*; do
        [ -f "$chunk" ] || continue
        computeMatrix reference-point \
            --referencePoint TSS \
            -S "$forward_bw" "$reverse_bw" \
            -R "$chunk" \
            --upstream 1000 \
            --downstream 1000 \
            --binSize 10 \
            --skipZeros \
            --missingDataAsZero \
            --samplesLabel forward reverse \
            -p "$threads" \
            -o "$CHUNK_MATRIX"
        python "$SCRIPT_DIR/profile_from_matrix.py" "$CHUNK_MATRIX" "$state_file"
        rm -f "$CHUNK_MATRIX"
    done
}

if [ "$ncols" -ge 6 ]; then
    STATE_PLUS="$WORK_DIR/state_plus.json"
    STATE_MINUS="$WORK_DIR/state_minus.json"

    awk '$6 == "+"' "$bed" > "$WORK_DIR/plus.bed"
    awk '$6 == "-"' "$bed" > "$WORK_DIR/minus.bed"

    run_accumulate "$WORK_DIR/plus.bed"  "$STATE_PLUS"  "plus_chunk"
    run_accumulate "$WORK_DIR/minus.bed" "$STATE_MINUS" "minus_chunk"

    python "$SCRIPT_DIR/profile_from_matrix.py" \
        --combine-strands "$STATE_PLUS" "$STATE_MINUS" "$output"
else
    STATE="$WORK_DIR/state.json"

    split -l "$CHUNK_SIZE" -a 3 "$bed" "$WORK_DIR/chunk_"
    for chunk in "$WORK_DIR"/chunk_*; do
        [ -f "$chunk" ] || continue
        computeMatrix reference-point \
            --referencePoint TSS \
            -S "$forward_bw" "$reverse_bw" \
            -R "$chunk" \
            --upstream 1000 \
            --downstream 1000 \
            --binSize 10 \
            --skipZeros \
            --missingDataAsZero \
            --samplesLabel forward reverse \
            -p "$threads" \
            -o "$CHUNK_MATRIX"
        python "$SCRIPT_DIR/profile_from_matrix.py" "$CHUNK_MATRIX" "$STATE"
        rm -f "$CHUNK_MATRIX"
    done

    python "$SCRIPT_DIR/profile_from_matrix.py" --finalize "$STATE" "$output"
fi
