####################################################################
# Step 22 — Alternative replication origin analysis (reviewer response)
#
# Tests PDALseq enrichment at replication origins from an independent
# method and group: Guilbaud et al. 2022 (Nucleic Acids Res. 50:7436)
# using ini-seq2 (GSM5658908, GSM6127777) and SNS-seq
# (GSM5658909, GSM5658910).
#
# The GEO BED files are in hg38 coordinates. Rule Alt_liftover_origins
# converts them to chm13v2.0 (UCSC name: hs1, chr1/chr2/... nomenclature)
# using the UCSC liftOver tool. The chain file must be placed at:
#     resources/liftover/hg38ToHs1.over.chain.gz
# (download from UCSC goldenPath hg38/liftOver/).
#
# Mirrors step 11's unscrambled / scrambled design and HMM state
# stratification. Includes the original 4 Hsapien PDALseq experiments
# (matching step 11) and K562 experiments from step 21.
#
# All outputs are under data/22_Alternative_origins/ and
# results/22_Alternative_origins/. Pre-existing bg files from steps 04
# and 21 are referenced via ancient() and never regenerated. State BEDs
# from step 10 are also ancient().
####################################################################

ALT_ORIGIN = [
    "GSM5658908_Ini-seq2.called.replication.origins",
    "GSM5658909_SNS.seq.fractions.0.5.2.kb_called.replication.origins",
    "GSM5658910_SNS.seq.fractions.2.4.kb_called.replication.origins",
    "GSM6127777_Ini-seq2.called.replication.origin.15.min",
]

HSAPIEN_4_EXPERIMENTS_22 = [
    "Hsapien_CCL86-Raji_PDALseq",
    "Hsapien_GM24385_PDALseq",
    "Hsapien_GM24143_PDALseq",
    "Hsapien_HEK293T_PDALseq",
]

K562_EXPERIMENTS_22 = ["Hsapien_K562_PDALseq"]

ALL_EXPERIMENTS_22 = HSAPIEN_4_EXPERIMENTS_22 + K562_EXPERIMENTS_22


def alt_origin_bg(wildcards):
    """Route the 1000-nt bg file to step 04 (existing) or step 21 (K562)."""
    if wildcards.experiment in K562_EXPERIMENTS_22:
        return ancient(
            f"data/21_add_K562_to_analysis/04_nonB/04_genome_window_coverage"
            f"/1000_nucleotides/{wildcards.experiment}.bg"
        )
    return ancient(
        f"data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage"
        f"/1000_nucleotides/{wildcards.experiment}.bg"
    )


# Convert hg38 → chm13v2.0 (hs1) using the dedicated liftover env so
# that no existing conda environment hash is changed.
rule Alt_liftover_origins:
    input:
        bed   = "resources/functional_annotations/alternative_annotations/{origin}.bed",
        chain = "resources/liftover/hg38ToHs1.over.chain.gz"
    output:
        lifted   = "data/22_Alternative_origins/liftover/{origin}.bed",
        unmapped = "data/22_Alternative_origins/liftover/{origin}.unmapped"
    log:
        "logs/22_Alternative_origins/liftover/{origin}.txt"
    conda:
        "../env/liftover.yml"
    shell:
        "bash workflow/scripts/22_Alternative_origins/preprocess_bed.sh "
        "{input.bed} {output.lifted}.tmp 2> {log} && "
        "liftOver {output.lifted}.tmp {input.chain} {output.lifted} {output.unmapped} 2>> {log} && "
        "rm {output.lifted}.tmp"


# Strip to 3 columns and sort so downstream bedtools operations are safe.
rule Alt_Get_origin_annotations:
    input:
        "data/22_Alternative_origins/liftover/{origin}.bed"
    output:
        "data/22_Alternative_origins/unscrambled/state_all/{origin}.bed"
    log:
        "logs/22_Alternative_origins/unscrambled/state_all/{origin}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "cut -f1-3 {input} | bedtools sort -i stdin > {output} 2> {log}"


rule Alt_Filter_by_state_annotations:
    wildcard_constraints:
        state = r"\d+"
    input:
        origin = "data/22_Alternative_origins/unscrambled/state_all/{origin}.bed",
        state  = ancient("data/10_Methylation_analysis/01_state_annotations/hominoidae/8/Hsapien/state_{state}.bed")
    output:
        "data/22_Alternative_origins/unscrambled/state_{state}/{origin}.bed"
    log:
        "logs/22_Alternative_origins/unscrambled/state_{state}/{origin}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools intersect -a {input.origin} -b {input.state} > {output} 2> {log}"


rule Alt_scramble_origin_annotations:
    input:
        "data/22_Alternative_origins/unscrambled/state_{state}/{origin}.bed"
    output:
        "data/22_Alternative_origins/scrambled/state_{state}/{origin}.bed"
    log:
        "logs/22_Alternative_origins/scrambled/state_{state}/{origin}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools shuffle -i {input} -g resources/genomes/Hsapien.chrom.sizes > {output} 2> {log}"


rule Alt_Mean_signal:
    input:
        bg         = alt_origin_bg,
        annotation = "data/22_Alternative_origins/{shuffle}/state_{state}/{origin}.bed"
    output:
        "data/22_Alternative_origins/Mean_signal/{experiment}/{shuffle}/state_{state}/{origin}.csv"
    log:
        "logs/22_Alternative_origins/Mean_signal/{experiment}/{shuffle}/state_{state}/{origin}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/11_Origin_analysis_analysis/PDAL_seq_around_origins.R "
        "{input.bg} {input.annotation} 100 {output} 2> {log}"


rule Alt_Compile_mean_signal:
    input:
        expand(
            "data/22_Alternative_origins/Mean_signal/{experiment}/{shuffle}/state_{state}/{origin}.csv",
            experiment=ALL_EXPERIMENTS_22,
            shuffle=["unscrambled", "scrambled"],
            state=["all", "4"],
            origin=ALT_ORIGIN,
        )
    output:
        "results/22_Alternative_origins/PDALseq_signal_around_alternative_origins.csv"
    log:
        "logs/22_Alternative_origins/Compile_mean_signal.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/22_Alternative_origins/Compile_alternative_origins.R {output} 2> {log}"
