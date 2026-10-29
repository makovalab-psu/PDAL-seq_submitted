#!/bin/bash
# make_regions.sh
# Generate 16 mutually-exclusive BED region files per species for syntenic analysis.
#
# Arguments:
#   $1  alignment BED     resources/alignments/hominoidae/{species}.bed
#   $2  promoter BED      resources/functional_annotations/Saswats_annotations/{species}/promoter_protCode_singleTrans.bed
#   $3  utr5 BED          ...utr5_protCode_singleTrans.bed
#   $4  exon BED          ...exon_protCode_singleTrans.bed
#   $5  intron BED        ...intron_protCode_singleTrans.bed
#   $6  utr3 BED          ...utr3_protCode_singleTrans.bed
#   $7  censat BED        resources/repeats/linneas_annotations/*CenSat*.bed
#   $8  all_repeats BED   data/19_.../01_bed_regions/{species}/all_repeats.bed
#   $9  chrom.sizes       resources/genomes/{species}.chrom.sizes
#   $10 output directory  data/19_.../01_bed_regions/{species}/
#
# Output BED files (written to output directory):
#   syntenic_{promoter,utr5,exon,intron,utr3,censat,nongenic_repeats,nongenic_nonrepetitive}.bed
#   nonsyntenic_{promoter,utr5,exon,intron,utr3,censat,nongenic_repeats,nongenic_nonrepetitive}.bed
#
# Priority order (highest → lowest):
#   censat > promoter > utr5 > utr3 > exon > intron > nongenic_repeats > nongenic_nonrepetitive
# Each base position falls into exactly one of the 8 categories per syntenic/nonsyntenic
# half — so all 16 output BEDs are mutually exclusive at the nucleotide level.

set -euo pipefail

alignment=$1
promoter=$2
utr5=$3
exon=$4
intron=$5
utr3=$6
censat=$7
repeats=$8
chrom_sizes=$9
outdir=${10}

mkdir -p "$outdir"

# ---------------------------------------------------------------------------
# Helper: sort+merge a BED, or produce an empty file if input is absent/empty
# ---------------------------------------------------------------------------
safe_sort_merge () {
    local infile=$1
    local outfile=$2
    if [ -s "$infile" ]; then
        bedtools sort -i "$infile" | bedtools merge > "$outfile"
    else
        touch "$outfile"
    fi
}

# ---------------------------------------------------------------------------
# 1. Prepare base regions
# ---------------------------------------------------------------------------

# Merged syntenic and its complement (non-syntenic).
# Must sort by genome file order (-g) so bedtools complement matches chrom.sizes order.
bedtools sort -g "$chrom_sizes" -i "$alignment" | bedtools merge > "${outdir}/tmp_syntenic.bed"
bedtools complement -i "${outdir}/tmp_syntenic.bed" -g "$chrom_sizes" > "${outdir}/tmp_nonsyntenic.bed"

# Merged genic union (for subtraction in non-genic categories)
all_genic_inputs=""
for f in "$promoter" "$utr5" "$exon" "$intron" "$utr3"; do
    [ -s "$f" ] && all_genic_inputs="$all_genic_inputs $f"
done
if [ -n "$all_genic_inputs" ]; then
    cat $all_genic_inputs | bedtools sort -i stdin | bedtools merge > "${outdir}/tmp_genic_union.bed"
else
    touch "${outdir}/tmp_genic_union.bed"
fi

# Merged CenSat and merged repeats
safe_sort_merge "$censat"  "${outdir}/tmp_censat.bed"
safe_sort_merge "$repeats" "${outdir}/tmp_repeats.bed"

# ---------------------------------------------------------------------------
# 2. Generate syntenic regions
# ---------------------------------------------------------------------------

for base in syntenic nonsyntenic; do
    base_bed="${outdir}/tmp_${base}.bed"

    # --- CenSat (highest priority): base ∩ censat ---
    if [ -s "$base_bed" ] && [ -s "${outdir}/tmp_censat.bed" ]; then
        bedtools intersect -a "$base_bed" -b "${outdir}/tmp_censat.bed" \
            | bedtools sort -i stdin \
            | bedtools merge \
            > "${outdir}/${base}_censat.bed"
    else
        touch "${outdir}/${base}_censat.bed"
    fi

    # --- Genic features in priority order: promoter > utr5 > utr3 > exon > intron ---
    # Each = (base ∩ feature) − everything already assigned (censat + higher-priority genic).
    # `assigned` accumulates the running union of all categories written so far in this base.
    assigned="${outdir}/tmp_${base}_assigned.bed"
    cp "${outdir}/${base}_censat.bed" "$assigned"

    for feat_name in promoter utr5 utr3 exon intron; do
        case $feat_name in
            promoter) feat_bed="$promoter" ;;
            utr5)     feat_bed="$utr5" ;;
            utr3)     feat_bed="$utr3" ;;
            exon)     feat_bed="$exon" ;;
            intron)   feat_bed="$intron" ;;
        esac
        out_bed="${outdir}/${base}_${feat_name}.bed"

        if [ -s "$feat_bed" ] && [ -s "$base_bed" ]; then
            bedtools intersect -a "$base_bed" -b "$feat_bed" \
                | bedtools sort -i stdin \
                | bedtools subtract -a stdin -b "$assigned" \
                | bedtools sort -i stdin \
                | bedtools merge \
                > "$out_bed"
        else
            touch "$out_bed"
        fi

        # Fold this feature's output into the cumulative mask. Write to a
        # sibling file then mv, since cat reads $assigned while > would truncate it.
        cat "$assigned" "$out_bed" \
            | bedtools sort -i stdin | bedtools merge \
            > "${outdir}/tmp_${base}_assigned_new.bed"
        mv "${outdir}/tmp_${base}_assigned_new.bed" "$assigned"
    done

    # --- Non-genic repeats: (base ∩ repeats) - genic_union - censat ---
    if [ -s "$base_bed" ] && [ -s "${outdir}/tmp_repeats.bed" ]; then
        bedtools intersect -a "$base_bed" -b "${outdir}/tmp_repeats.bed" \
            | bedtools sort -i stdin \
            | bedtools subtract -a stdin -b "${outdir}/tmp_genic_union.bed" \
            | bedtools sort -i stdin \
            | bedtools subtract -a stdin -b "${outdir}/tmp_censat.bed" \
            | bedtools sort -i stdin \
            | bedtools merge \
            > "${outdir}/${base}_nongenic_repeats.bed"
    else
        touch "${outdir}/${base}_nongenic_repeats.bed"
    fi

    # --- Non-genic non-repetitive: base - genic_union - censat - repeats ---
    if [ -s "$base_bed" ]; then
        bedtools subtract -a "$base_bed" -b "${outdir}/tmp_genic_union.bed" \
            | bedtools sort -i stdin \
            | bedtools subtract -a stdin -b "${outdir}/tmp_censat.bed" \
            | bedtools sort -i stdin \
            | bedtools subtract -a stdin -b "${outdir}/tmp_repeats.bed" \
            | bedtools sort -i stdin \
            | bedtools merge \
            > "${outdir}/${base}_nongenic_nonrepetitive.bed"
    else
        touch "${outdir}/${base}_nongenic_nonrepetitive.bed"
    fi

done

# ---------------------------------------------------------------------------
# 3. Clean up temp files
# ---------------------------------------------------------------------------
rm -f "${outdir}"/tmp_*.bed
