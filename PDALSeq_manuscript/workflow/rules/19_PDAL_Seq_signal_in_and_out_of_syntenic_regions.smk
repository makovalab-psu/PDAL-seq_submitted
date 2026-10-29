
# Analysis of PDAL-Seq signal inside and outside syntenic regions
# Addresses reviewer comment: "how much PDAL-seq signal exists in non-syntenic regions?"
# For each of 6 species, computes read counts in 16 genomic region categories:
#   (syntenic | non-syntenic) × (promoter | utr5 | exon | intron | utr3 | censat | nongenic_repeats | nongenic_nonrepetitive)

# Species included (those with alignment BED files in hominoidae/ and PDAL-Seq data)
# Note: functional annotation files (promoter, utr5, etc.) may be absent for Pabelii and Ppygmaeus.
# If so, the bash script will treat those genic categories as empty (0 reads, 0 nt).
SYNTENIC_SPECIES_19 = ["Ggorilla", "Hsapien", "Ppaniscus", "Ptroglodytes", "Pabelii", "Ppygmaeus"]

# Species → RepeatMasker .out file
SPECIES_RM_OUT_19 = {
    "Hsapien":     "resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out",
    "Ggorilla":    "resources/repeats/linneas_annotations/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out",
    "Ppaniscus":   "resources/repeats/linneas_annotations/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out",
    "Ptroglodytes":"resources/repeats/linneas_annotations/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out",
    "Pabelii":     "resources/repeats/linneas_annotations/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out",
    "Ppygmaeus":   "resources/repeats/linneas_annotations/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out",
}

# Species → CenSat BED file
SPECIES_CENSAT_19 = {
    "Hsapien":     "resources/repeats/linneas_annotations/chm13v2.0_censat_v2.0.bed",
    "Ggorilla":    "resources/repeats/linneas_annotations/mGorGor1_v2.0_CenSat_v2.0.bed",
    "Ppaniscus":   "resources/repeats/linneas_annotations/mPanPan1_v2.0_CenSat_v2.0.bed",
    "Ptroglodytes":"resources/repeats/linneas_annotations/mPanTro3_v2.0_CenSat_v2.0.bed",
    "Pabelii":     "resources/repeats/linneas_annotations/mPonAbe1_v2.0_CenSat_v2.0.bed",
    "Ppygmaeus":   "resources/repeats/linneas_annotations/mPonPyg2_v2.0_CenSat_v2.0.bed",
}

# Output region names (syntenic × category and nonsyntenic × category)
SYNTENIC_REGIONS_19 = [
    "syntenic_promoter", "syntenic_utr5", "syntenic_exon",
    "syntenic_intron", "syntenic_utr3", "syntenic_censat",
    "syntenic_nongenic_repeats", "syntenic_nongenic_nonrepetitive",
    "nonsyntenic_promoter", "nonsyntenic_utr5", "nonsyntenic_exon",
    "nonsyntenic_intron", "nonsyntenic_utr3", "nonsyntenic_censat",
    "nonsyntenic_nongenic_repeats", "nonsyntenic_nongenic_nonrepetitive",
]

# Build species → experiment list mapping (using DATA from Snakefile)
_syntenic_df = df[df['Species'].isin(SYNTENIC_SPECIES_19)]
SPECIES_TO_EXP_19 = {
    sp: list(set(_syntenic_df[_syntenic_df['Species'] == sp]['Experiment']))
    for sp in SYNTENIC_SPECIES_19
}

# Flat lists for zip-expand in compile rule
_exp_sp_pairs = [
    (exp, sp)
    for sp in SYNTENIC_SPECIES_19
    for exp in SPECIES_TO_EXP_19[sp]
]
_EXP_LIST_19 = [p[0] for p in _exp_sp_pairs]
_SP_LIST_19  = [p[1] for p in _exp_sp_pairs]


rule rm_out_to_bed_19:
    """Convert RepeatMasker .out file to a merged BED of all repeat elements."""
    input:
        lambda wildcards: SPECIES_RM_OUT_19[wildcards.species]
    output:
        "data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/01_bed_regions/{species}/all_repeats.bed"
    log:
        "logs/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/01_rm_out_to_bed/{species}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=4000
    shell:
        "awk 'NR>3 {{print $5\"\\t\"$6\"\\t\"$7}}' {input} "
        "| bedtools sort -i stdin "
        "| bedtools merge "
        "> {output} 2> {log}"


rule make_syntenic_regions_19:
    """
    Generate 16 mutually-exclusive region BED files per species:
      (syntenic | non-syntenic) ×
      (promoter | utr5 | exon | intron | utr3 | censat | nongenic_repeats | nongenic_nonrepetitive)

    Priority order (highest → lowest):
      censat > promoter > utr5 > utr3 > exon > intron > nongenic_repeats > nongenic_nonrepetitive
    Each base position is assigned to exactly one of the 8 categories within each
    syntenic/nonsyntenic half — all 16 outputs are mutually exclusive at the
    nucleotide level.
    """
    input:
        alignment  = "resources/alignments/hominoidae/{species}.bed",
        promoter   = "resources/functional_annotations/Saswats_annotations/{species}/promoter_protCode_singleTrans.bed",
        utr5       = "resources/functional_annotations/Saswats_annotations/{species}/utr5_protCode_singleTrans.bed",
        exon       = "resources/functional_annotations/Saswats_annotations/{species}/exon_protCode_singleTrans.bed",
        intron     = "resources/functional_annotations/Saswats_annotations/{species}/intron_protCode_singleTrans.bed",
        utr3       = "resources/functional_annotations/Saswats_annotations/{species}/utr3_protCode_singleTrans.bed",
        censat     = lambda wildcards: SPECIES_CENSAT_19[wildcards.species],
        repeats    = "data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/01_bed_regions/{species}/all_repeats.bed",
        chrom_sizes= "resources/genomes/{species}.chrom.sizes",
    output:
        expand(
            "data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/01_bed_regions/{{species}}/{region}.bed",
            region=SYNTENIC_REGIONS_19
        )
    log:
        "logs/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/01_make_regions/{species}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=8000
    shell:
        "bash workflow/scripts/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/make_regions.sh "
        "{input.alignment} {input.promoter} {input.utr5} {input.exon} {input.intron} {input.utr3} "
        "{input.censat} {input.repeats} {input.chrom_sizes} "
        "data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/01_bed_regions/{wildcards.species}/ "
        "2> {log}"


rule bam_to_reads_bed_19:
    """
    Convert BAM to a sorted BED of primary mapped read reference intervals.
    Used as a single-shot input to all 16 region count jobs per experiment,
    so each base of each read is counted in exactly one region (Σ fraction = 1).
    Flag filter: -F 0x904 = exclude unmapped (0x4), secondary (0x100),
    supplementary (0x800).
    Marked temp(): Snakemake deletes the per-experiment reads BED once all 16
    region count jobs for that experiment have completed.
    """
    input:
        bam = "data/01_preprocess_data/All_reads/{experiment}.bam",
        bai = "data/01_preprocess_data/All_reads/{experiment}.bam.bai",
    output:
        temp("data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/reads_bed/{experiment}.bed")
    log:
        "logs/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/reads_bed/{experiment}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=4000
    shell:
        "samtools view -b -F 0x904 {input.bam} "
        "| bedtools bamtobed -i stdin "
        "| sort -k1,1 -k2,2n "
        "> {output} 2> {log}"


rule count_reads_in_syntenic_region:
    """Count read-aligned bases overlapping a specific region BED."""
    input:
        reads = "data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/reads_bed/{experiment}.bed",
        bed   = "data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/01_bed_regions/{species}/{region}.bed",
    output:
        "data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/02_read_counts/{species}/{region}/{experiment}.txt"
    log:
        "logs/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/02_read_counts/{species}/{region}/{experiment}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=2000
    shell:
        "bash workflow/scripts/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/count_reads.sh "
        "{input.reads} {input.bed} {output} 2> {log}"


rule compile_syntenic_analysis:
    """Compile all per-region read counts into a single CSV with metadata."""
    input:
        counts = [
            f"data/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/02_read_counts/{sp}/{region}/{exp}.txt"
            for sp, exp in zip(_SP_LIST_19, _EXP_LIST_19)
            for region in SYNTENIC_REGIONS_19
        ],
        samples = "PDALseq_samples.txt",
    output:
        "results/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/PDAL_Seq_signal_in_and_out_of_syntenic_regions.csv"
    log:
        "logs/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/compile.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/19_PDAL_Seq_signal_in_and_out_of_syntenic_regions/compile_results.R "
        "{output} 2> {log}"
