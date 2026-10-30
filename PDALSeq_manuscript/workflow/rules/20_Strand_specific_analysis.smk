
# Strand-specific PDAL-Seq coverage around non-B DNA motifs
# Addresses reviewer questions about strand specificity:
#   - Reviewer 1: do forward vs. reverse reads show asymmetry at non-B structures?
#   - Reviewer 2: does the coverage pattern match the expected structural mechanism
#     (e.g., G4 on + strand → elevated - strand coverage on the complementary ssDNA)?
#
# For each experiment and non-B type, computes mean forward-strand and reverse-strand
# bigwig coverage in 10 bp bins across ±1000 bp around each feature.
# BED6 features (e.g. GQ) get strand-aware window orientation; BED3 features are
# treated as unstranded.

NONB_TYPES_20 = ["APR", "CRU", "DR", "GQ", "TRI", "Z"]

# All species present in PDALseq_samples.txt
SPECIES_20 = list(SPECIES)  # defined in Snakefile from PDALseq_samples.txt

# Species → raw assembly name prefix used in resources/nonb_annotations/
SPECIES_ASSEMBLY_20 = {
    "Hsapien":     "chm13v2.0",
    "Ggorilla":    "mGorGor1.pri.cur.20231122",
    "Ppaniscus":   "mPanPan1.pri.cur.20231122",
    "Ptroglodytes":"mPanTro3.pri.cur.20231122",
    "Pabelii":     "mPonAbe1.pri.cur.20231205",
    "Ppygmaeus":   "mPonPyg2.pri.cur.20231122",
    "Ssyndactylus":"mSymSyn1.pri.cur.20240514",
}

# Subset of experiments to run for step 20. One PDALseq + one matched TruSeq per species.
# Restricting here avoids generating bigwigs and profiles for every
# experiment in PDALseq_samples.txt, which is the main source of runtime cost.
EXPERIMENTS_20 = [
    "Hsapien_GM24385_PDALseq",      "Hsapien_GM24385_TruSeq",
    "Ggorilla_KB3781_PDALseq",      "Ggorilla_KB3781_TruSeq",
    "Ptroglodytes_AG18354_PDALseq", "Ptroglodytes_AG18354_TruSeq",
    "Ppaniscus_KB8711_PDALseq",     "Ppaniscus_KB8711_TruSeq",
    "Ppygmaeus_AG05252_PDALseq",    "Ppygmaeus_AG05252_TruSeq",
    "Pabelii_AG06213_PDALseq",      "Pabelii_AG06213_TruSeq",
    "Ssyndactylus_JAMBI_PDALseq",   "Ssyndactylus_JAMBI_TruSeq",
]

# Build species → experiment list mapping (DATA defined in Snakefile)
_df20 = df[df["Experiment"].isin(EXPERIMENTS_20)].drop_duplicates(subset=["Experiment", "Species"])
SPECIES_TO_EXP_20 = {
    sp: list(set(_df20[_df20["Species"] == sp]["Experiment"]))
    for sp in SPECIES_20
    if sp in set(_df20["Species"])
}

# Flat list of all valid (species, nonb, experiment) triples for expand
_sp_nonb_exp_triples_20 = [
    (sp, nonb, exp)
    for sp in SPECIES_20
    for exp in SPECIES_TO_EXP_20.get(sp, [])
    for nonb in NONB_TYPES_20
]
_SP_LIST_20   = [t[0] for t in _sp_nonb_exp_triples_20]
_NONB_LIST_20 = [t[1] for t in _sp_nonb_exp_triples_20]
_EXP_LIST_20  = [t[2] for t in _sp_nonb_exp_triples_20]


rule bam_to_strand_bigwigs_20:
    """
    Generate forward (+) and reverse (-) strand bigwig files from a merged BAM.
    Forward bigwig: reads mapping to the + genomic strand (flag 16 NOT set).
    Reverse bigwig: reads mapping to the - genomic strand (flag 16 set).
    Uses deeptools bamCoverage (CPM-normalized, 10 bp bins) instead of
    bedtools genomecov + bedGraphToBigWig — the bedtools -strand - flag produced
    malformed/unsorted bedgraphs for several species BAMs.
    These are kept (not temp) because all 10 non-B type profile rules depend on them.
    """
    input:
        bam = "data/01_preprocess_data/All_reads/{experiment}.bam",
        bai = "data/01_preprocess_data/All_reads/{experiment}.bam.bai",
    output:
        forward_bw = "data/20_Strand_specific_analysis/01_bigwigs/{experiment}_forward.bigwig",
        reverse_bw = "data/20_Strand_specific_analysis/01_bigwigs/{experiment}_reverse.bigwig",
    log:
        "logs/20_Strand_specific_analysis/01_bigwigs/{experiment}.txt"
    conda:
        "../env/20_strand_analysis.yml"
    threads: 4
    resources:
        mem_mb=8000
    shell:
        "bash workflow/scripts/20_Strand_specific_analysis/make_strand_bigwigs.sh "
        "{input.bam} {output.forward_bw} {output.reverse_bw} {threads} "
        "2> {log}"


rule prepare_nonb_bed_20:
    """
    Prepare raw non-B annotation BED for strand-specific profiling:
      - Filters sex-chromosome features (chrX, chrY, chrW, chrZ).
      - Outputs BED6 if input has ≥ 6 columns (preserving strand); BED3 otherwise.
    Only GQ annotations are expected to be BED6 (strand-specific G4 predictions).
    """
    input:
        lambda wildcards: (
            "resources/nonb_annotations/{}.{}.bed".format(
                SPECIES_ASSEMBLY_20[wildcards.species], wildcards.nonb
            )
        ),
    output:
        "data/20_Strand_specific_analysis/02_nonb_beds/{species}/{nonb}.bed"
    params:
        swap = lambda wildcards: NONB_CHROM_SWAP.get(wildcards.species, "")
    log:
        "logs/20_Strand_specific_analysis/02_nonb_beds/{species}/{nonb}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=2000
    shell:
        "bash workflow/scripts/20_Strand_specific_analysis/prepare_nonb_bed.sh "
        "{input} {output} '{params.swap}' 2> {log}"


rule compute_strand_profile_20:
    """
    Compute mean forward/reverse strand coverage in 10 bp bins across ±1000 bp
    around each non-B feature. compute_strand_profile.sh chunks the BED into
    200k-feature batches; each batch runs computeMatrix independently (keeping
    peak RAM at ~1-2 GB per batch), then profile_from_matrix.py streams the
    batch matrix into a running-sum JSON without ever holding the full
    N_features × N_bins array in memory. mem_mb=8000 comfortably covers one
    200k-feature computeMatrix call; the Python accumulate step is negligible.
    For BED6 features (e.g. GQ), deeptools automatically orients the window 5'→3'
    relative to each feature's strand — the coverage strand labels are NOT swapped.
    """
    input:
        forward_bw = "data/20_Strand_specific_analysis/01_bigwigs/{experiment}_forward.bigwig",
        reverse_bw = "data/20_Strand_specific_analysis/01_bigwigs/{experiment}_reverse.bigwig",
        bed        = "data/20_Strand_specific_analysis/02_nonb_beds/{species}/{nonb}.bed",
    output:
        "data/20_Strand_specific_analysis/03_profiles/{species}/{nonb}/{experiment}.tsv"
    log:
        "logs/20_Strand_specific_analysis/03_profiles/{species}/{nonb}/{experiment}.txt"
    conda:
        "../env/20_strand_analysis.yml"
    threads: 4
    resources:
        mem_mb=8000
    shell:
        "bash workflow/scripts/20_Strand_specific_analysis/compute_strand_profile.sh "
        "{input.forward_bw} {input.reverse_bw} {input.bed} {output} {threads} "
        "2> {log}"


rule compile_strand_profiles_20:
    """Compile all per-experiment strand profiles into a single CSV with metadata."""
    input:
        profiles = [
            f"data/20_Strand_specific_analysis/03_profiles/{sp}/{nonb}/{exp}.tsv"
            for sp, nonb, exp in _sp_nonb_exp_triples_20
        ],
        samples = "PDALseq_samples.txt",
    output:
        "results/20_Strand_specific_analysis/Strand_specific_coverage.csv"
    log:
        "logs/20_Strand_specific_analysis/compile.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/20_Strand_specific_analysis/compile_results.R "
        "{output} 2> {log}"
