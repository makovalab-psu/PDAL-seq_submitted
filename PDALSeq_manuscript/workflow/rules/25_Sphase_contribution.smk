####################################################################
# Step 25 — S-phase contribution analysis
#
# Addresses reviewer question: are PDAL-seq hits driven by S-phase
# ssDNA rather than non-B DNA? GM24143 and GM24385 cells were
# synchronized with 1.5% DMSO for 96 h to deplete S-phase, then
# PDAL-seq was repeated.
#
# Phase 1: per-sample preprocessing + technical-replicate correlation
# against each other and against matched bulk (step 03) samples.
# No merging until correlation plots are reviewed.
#
# All outputs are under data/25_Sphase_contribution/ and
# results/25_Sphase_contribution/. Step-03 bulk coverage files are
# referenced via ancient() and never regenerated.
####################################################################

DMSO_FASTQ_DIR = "Revisions_data/Raw_fastq/260512_VH00707_312_AAJ7GVWM5/fastq"

DMSO_SAMPLES = [
    "Hsa_GM24143_DMSO_PDALSeq_js4012d20mM_S4_R1_001",
    "Hsa_GM24143_DMSO_PDALSeq_js4012e20mM_S5_R1_001",
    "Hsa_GM24143_DMSO_PDALSeq_js4012f20mM_S6_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012a20mM_S1_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012b20mM_S2_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012c20mM_S3_R1_001",
]

DMSO_WINDOWS = ["1000000", "100000", "10000"]

# Bulk (pre-DMSO) samples from step 03 used for cross-condition correlation.
# Names use the Hsapien_ prefix and are stored under
# data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/.
DMSO_BULK_GM24143 = [
    "Hsapien_GM24143_LCL_PDALseq_js2007d20mM_S2_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007e20mM_S3_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007f40mM_S4_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007g40mM_S5_R1_001",
]

DMSO_BULK_GM24385 = [
    "Hsapien_GM24385_LCL_PDALseq_js1014A20_1_S2_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A20_2_S3_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A40_1_S4_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A40_2_S5_R1_001",
    "Hsapien_GM24385_TruSeq_SRR12898317_forward",
    "Hsapien_GM24385_TruSeq_SRR12898319_forward",
    "Hsapien_GM24385_TruSeq_SRR12898320_forward",
]


############################
# Preprocessing (mirrors step 21 K562 chain, R1 only, hard-coded Hsapien genome)
############################

rule DMSO_trim_fastq:
    input:
        DMSO_FASTQ_DIR + "/{sample}.fastq.gz"
    output:
        temp("data/25_Sphase_contribution/01_preprocess/Trimmed_fastq/{sample}.fastq.gz")
    log:
        "logs/25_Sphase_contribution/01_preprocess/01_trim_fastq/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=1000
    shell:
        "fastp --thread 1 -i {input} -o {output} -p -P 50 -y -x -3 2> {log}"


rule DMSO_bwa_mem:
    input:
        fastq  = "data/25_Sphase_contribution/01_preprocess/Trimmed_fastq/{sample}.fastq.gz",
        genome = ancient("resources/genomes/Hsapien.fasta"),
        amb    = ancient("resources/genomes/Hsapien.fasta.amb"),
        ann    = ancient("resources/genomes/Hsapien.fasta.ann"),
        pac    = ancient("resources/genomes/Hsapien.fasta.pac"),
        bwt    = ancient("resources/genomes/Hsapien.fasta.bwt"),
        sa     = ancient("resources/genomes/Hsapien.fasta.sa"),
    output:
        temp("data/25_Sphase_contribution/01_preprocess/Mapped_sam/{sample}.sam")
    log:
        "logs/25_Sphase_contribution/01_preprocess/02_bwa_mem/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=16000
    threads:
        5
    shell:
        "bwa mem -t {threads} -o {output} {input.genome} {input.fastq} 2> {log}"


rule DMSO_sort_sam:
    input:
        "data/25_Sphase_contribution/01_preprocess/Mapped_sam/{sample}.sam"
    output:
        temp("data/25_Sphase_contribution/01_preprocess/Mapped_bam/{sample}.bam")
    log:
        "logs/25_Sphase_contribution/01_preprocess/03_sort_sam/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=8000
    shell:
        "samtools sort {input} -T temp/DMSO_{wildcards.sample} -O bam -o {output} 2> {log}"


rule DMSO_remove_dup:
    input:
        "data/25_Sphase_contribution/01_preprocess/Mapped_bam/{sample}.bam"
    output:
        "data/25_Sphase_contribution/01_preprocess/Mapped_deduplicated_bam/{sample}.bam"
    log:
        "logs/25_Sphase_contribution/01_preprocess/04_remove_dup/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=8000
    shell:
        "samtools rmdup -s {input} {output} 2> {log}"


############################
# Technical replicate correlation
############################

rule DMSO_sort_sample_bam:
    input:
        "data/25_Sphase_contribution/01_preprocess/Mapped_deduplicated_bam/{sample}.bam"
    output:
        temp("data/25_Sphase_contribution/03_correlation/sorted_bam/{sample}.bam")
    log:
        "logs/25_Sphase_contribution/03_correlation/sorted_bam/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "samtools sort -o {output} {input} 2> {log}"


rule DMSO_make_genome_windows:
    input:
        ancient("resources/genomes/Hsapien.chrom.sizes")
    output:
        "data/25_Sphase_contribution/03_correlation/01_genome_windows/{size}_nucleotides/Hsapien.bed"
    log:
        "logs/25_Sphase_contribution/03_correlation/01_genome_windows/{size}_nucleotides.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/make_windows.sh {input} {output} {wildcards.size} 2> {log}"


rule DMSO_genome_window_coverage:
    input:
        bam     = "data/25_Sphase_contribution/03_correlation/sorted_bam/{sample}.bam",
        windows = "data/25_Sphase_contribution/03_correlation/01_genome_windows/{size}_nucleotides/Hsapien.bed"
    output:
        "data/25_Sphase_contribution/03_correlation/03_genome_window_coverage/{size}_nucleotides/{sample}.bed"
    log:
        "logs/25_Sphase_contribution/03_correlation/03_genome_window_coverage/{size}_nucleotides/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "megadepth {input.bam} --threads 1 --annotation {input.windows} --op sum > {output} 2> {log}"


rule DMSO_compile_correlation:
    input:
        dmso_cov = expand(
            "data/25_Sphase_contribution/03_correlation/03_genome_window_coverage/{size}_nucleotides/{sample}.bed",
            size=DMSO_WINDOWS, sample=DMSO_SAMPLES,
        ),
        bulk_gm24143 = [
            ancient(f"data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/{size}_nucleotides/{sample}.bed")
            for size in DMSO_WINDOWS
            for sample in DMSO_BULK_GM24143
        ],
        bulk_gm24385 = [
            ancient(f"data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/{size}_nucleotides/{sample}.bed")
            for size in DMSO_WINDOWS
            for sample in DMSO_BULK_GM24385
        ],
    output:
        "results/25_Sphase_contribution/Technical_replicate_correlation.csv"
    log:
        "logs/25_Sphase_contribution/compile_correlation.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/25_Sphase_contribution/compile_DMSO_correlation.R 2> {log}"


############################
# Phase 2 — 1 kb window coverage for S-phase contribution analysis
#
# Quantify each technical replicate (6 DMSO + 8 bulk LCL PDALseq) in
# 1 kb windows, normalize to RPM, and compile a wide-format coverage
# matrix for downstream statistical comparison.
#
# TruSeq bulk samples are excluded: different library prep makes them
# incompatible with a direct DMSO vs. PDALseq-bulk comparison.
#
# Source BAMs:
#   DMSO  — data/25_Sphase_contribution/01_preprocess/Mapped_deduplicated_bam/
#   Bulk  — PDALseq_data/Mapped_deduplicated_bam/  (ancient, no reruns)
############################

BULK_SAMPLES_1KB = [
    "Hsapien_GM24143_LCL_PDALseq_js2007d20mM_S2_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007e20mM_S3_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007f40mM_S4_R1_001",
    "Hsapien_GM24143_LCL_PDALseq_js2007g40mM_S5_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A20_1_S2_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A20_2_S3_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A40_1_S4_R1_001",
    "Hsapien_GM24385_LCL_PDALseq_js1014A40_2_S5_R1_001",
]


rule DMSO_make_1kb_windows:
    input:
        ancient("resources/genomes/Hsapien.chrom.sizes")
    output:
        "data/25_Sphase_contribution/04_1kb_coverage/windows/Hsapien.bed"
    log:
        "logs/25_Sphase_contribution/04_1kb_coverage/make_1kb_windows.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/make_windows.sh "
        "{input} {output} 1000 2> {log}"


rule DMSO_1kb_sort_bam:
    input:
        "data/25_Sphase_contribution/01_preprocess/Mapped_deduplicated_bam/{sample}.bam"
    output:
        temp("data/25_Sphase_contribution/04_1kb_coverage/dmso_sorted_bam/{sample}.bam")
    log:
        "logs/25_Sphase_contribution/04_1kb_coverage/dmso_sorted_bam/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "samtools sort -o {output} {input} 2> {log}"


rule DMSO_1kb_coverage:
    input:
        bam     = "data/25_Sphase_contribution/04_1kb_coverage/dmso_sorted_bam/{sample}.bam",
        windows = "data/25_Sphase_contribution/04_1kb_coverage/windows/Hsapien.bed"
    output:
        "data/25_Sphase_contribution/04_1kb_coverage/dmso/{sample}.bed"
    log:
        "logs/25_Sphase_contribution/04_1kb_coverage/dmso/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=4000
    shell:
        "megadepth {input.bam} --threads 1 --annotation {input.windows} --op sum > {output} 2> {log}"


rule Bulk_1kb_sort_bam:
    input:
        ancient("PDALseq_data/Mapped_deduplicated_bam/{sample}.bam")
    output:
        temp("data/25_Sphase_contribution/04_1kb_coverage/bulk_sorted_bam/{sample}.bam")
    log:
        "logs/25_Sphase_contribution/04_1kb_coverage/bulk_sorted_bam/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "samtools sort -o {output} {input} 2> {log}"


rule Bulk_1kb_coverage:
    input:
        bam     = "data/25_Sphase_contribution/04_1kb_coverage/bulk_sorted_bam/{sample}.bam",
        windows = "data/25_Sphase_contribution/04_1kb_coverage/windows/Hsapien.bed"
    output:
        "data/25_Sphase_contribution/04_1kb_coverage/bulk/{sample}.bed"
    log:
        "logs/25_Sphase_contribution/04_1kb_coverage/bulk/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=4000
    shell:
        "megadepth {input.bam} --threads 1 --annotation {input.windows} --op sum > {output} 2> {log}"


rule DMSO_compile_1kb_coverage:
    input:
        dmso_beds = expand(
            "data/25_Sphase_contribution/04_1kb_coverage/dmso/{sample}.bed",
            sample=DMSO_SAMPLES,
        ),
        bulk_beds = expand(
            "data/25_Sphase_contribution/04_1kb_coverage/bulk/{sample}.bed",
            sample=BULK_SAMPLES_1KB,
        ),
    output:
        "results/25_Sphase_contribution/Sphase_1kb_coverage.csv"
    log:
        "logs/25_Sphase_contribution/compile_1kb_coverage.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=16000
    shell:
        "Rscript workflow/scripts/25_Sphase_contribution/compile_1kb_coverage.R 2> {log}"


rule DMSO_sphase_analysis:
    input:
        "results/25_Sphase_contribution/Sphase_1kb_coverage.csv"
    output:
        sphase_143 = "results/25_Sphase_contribution/GM24143/Sphase_increases_ssDNA_windows.bed",
        g0g1_143   = "results/25_Sphase_contribution/GM24143/G0_G1_arrest_increases_ssDNA_windows.bed",
        sphase_385 = "results/25_Sphase_contribution/GM24385/Sphase_increases_ssDNA_windows.bed",
        g0g1_385   = "results/25_Sphase_contribution/GM24385/G0_G1_arrest_increases_ssDNA_windows.bed",
        csv        = "results/25_Sphase_contribution/Cell_cycle_effects_on_PDALSeq_1kbp_windows.csv",
    log:
        "logs/25_Sphase_contribution/sphase_analysis.txt"
    conda:
        "../env/25_Sphase_analysis.yml"
    resources:
        mem_mb=32000
    shell:
        "mkdir -p results/25_Sphase_contribution/GM24143 results/25_Sphase_contribution/GM24385 && "
        "Rscript workflow/scripts/25_Sphase_contribution/Sphase_analysis.R 2> {log}"


############################
# Phase 3 — Enrichment of S-phase / G0-G1 windows in ape HMM states
#
# Merges adjacent significant 1 kb windows within each of the 4 cell-cycle
# BED files, then tests each merged BED for enrichment in each of the 8
# hominoidae HMM states using the same bootstrap permutation approach as
# step 07 (Enrichment.sh: 10 000 subsamples, 100 simulations).
#
# State BEDs from step 07 are referenced via ancient() so this rule never
# triggers step-07 reruns.  Wildcards {state25} and {sphase_ann} are
# step-25-specific to avoid collision with step-07's {state}/{ann}.
############################

SPHASE_RAW_BEDS = {
    "Sphase_GM24143": "results/25_Sphase_contribution/GM24143/Sphase_increases_ssDNA_windows.bed",
    "G0G1_GM24143":   "results/25_Sphase_contribution/GM24143/G0_G1_arrest_increases_ssDNA_windows.bed",
    "Sphase_GM24385": "results/25_Sphase_contribution/GM24385/Sphase_increases_ssDNA_windows.bed",
    "G0G1_GM24385":   "results/25_Sphase_contribution/GM24385/G0_G1_arrest_increases_ssDNA_windows.bed",
}
SPHASE_ANN_NAMES = list(SPHASE_RAW_BEDS.keys())
SPHASE_STATES    = ["0", "1", "2", "3", "4", "5", "6", "7"]


rule Sphase_merge_windows:
    input:
        lambda wc: SPHASE_RAW_BEDS[wc.sphase_ann],
    output:
        "data/25_Sphase_contribution/05_merged_windows/{sphase_ann}.bed"
    log:
        "logs/25_Sphase_contribution/merge_windows/{sphase_ann}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    wildcard_constraints:
        sphase_ann = r"[A-Za-z0-9_]+"
    shell:
        "bedtools merge -i {input} > {output} 2> {log}"


rule Sphase_HMM_state_enrichment:
    input:
        state      = ancient("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/"
                             "hominoidae/8/state_{state25}.bed"),
        annotation = "data/25_Sphase_contribution/05_merged_windows/{sphase_ann}.bed",
    output:
        "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/"
        "State_functional_element_enrichment/state_{state25}/{sphase_ann}_enrichment.bed"
    log:
        "logs/25_Sphase_contribution/HMM_enrichment/state_{state25}/{sphase_ann}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=1000
    wildcard_constraints:
        state25    = r"[0-7]",
        sphase_ann = r"[A-Za-z0-9_]+"
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.state} {input.annotation} {output} 10000 100 2> {log}"


rule Sphase_compile_HMM_enrichment:
    input:
        expand(
            "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/"
            "State_functional_element_enrichment/state_{state25}/{sphase_ann}_enrichment.bed",
            state25=SPHASE_STATES, sphase_ann=SPHASE_ANN_NAMES,
        )
    output:
        "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/"
        "State_functional_element_enrichment/Sphase_enrichment.csv"
    log:
        "logs/25_Sphase_contribution/HMM_enrichment/compile.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=2000
    shell:
        "Rscript workflow/scripts/25_Sphase_contribution/compile_Sphase_HMM_enrichment.R 2> {log}"


############################
# Phase 4 — Non-B DNA and functional element enrichment/coverage
#
# For each of the 4 merged sphase window BEDs × 17 query annotations,
# compute:
#   (a) bootstrap enrichment (Enrichment.sh, 10 000 subsamples × 100 sims)
#   (b) bp coverage (bedtools intersect | awk sum)
#
# Note: genic annotations live in resources/functional_annotations/Saswats_annotations/Hsapien/,
# not resources/alignments/hominoidae/ (that directory contains multi-species alignment BEDs).
# RepeatMasker .out is converted to a merged BED before use.
# Reference is the merged sphase window BED (05_merged_windows), consistent with Phase 3.
############################

SPHASE_FUNC_ANNS = {
    "promoter": "resources/functional_annotations/Saswats_annotations/Hsapien/promoter_protCode_singleTrans.bed",
    "utr5":     "resources/functional_annotations/Saswats_annotations/Hsapien/utr5_protCode_singleTrans.bed",
    "exon":     "resources/functional_annotations/Saswats_annotations/Hsapien/exon_protCode_singleTrans.bed",
    "intron":   "resources/functional_annotations/Saswats_annotations/Hsapien/intron_protCode_singleTrans.bed",
    "utr3":     "resources/functional_annotations/Saswats_annotations/Hsapien/utr3_protCode_singleTrans.bed",
    "censat":   "resources/repeats/linneas_annotations/chm13v2.0_censat_v2.0.bed",
    "repeats":  "data/25_Sphase_contribution/06_nonb_func_enrichment/chm13v2.0_RepeatMasker.bed",
    "APR":      "resources/nonb_annotations/chm13v2.0.APR.bed",
    "CRU":      "resources/nonb_annotations/chm13v2.0.CRU.bed",
    "DR":       "resources/nonb_annotations/chm13v2.0.DR.bed",
    "GQ":       "resources/nonb_annotations/chm13v2.0.GQ.bed",
    "IR":       "resources/nonb_annotations/chm13v2.0.IR.bed",
    "MR":       "resources/nonb_annotations/chm13v2.0.MR.bed",
    "SLS":      "resources/nonb_annotations/chm13v2.0.SLS.bed",
    "STR":      "resources/nonb_annotations/chm13v2.0.STR.bed",
    "TRI":      "resources/nonb_annotations/chm13v2.0.TRI.bed",
    "Z":            "resources/nonb_annotations/chm13v2.0.Z.bed",
    "coreorigins":  "resources/functional_annotations/Saswats_annotations/Hsapien/core_origins_hglft_genome_e0596_1cc9e0.bed",
    "stochorigins": "resources/functional_annotations/Saswats_annotations/Hsapien/stochastic_origins_hglft_genome_2ff9ed_1ccfd0.bed",
}
SPHASE_FEAT_NAMES = list(SPHASE_FUNC_ANNS.keys())


rule Sphase_RM_to_bed:
    input:
        ancient("resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out")
    output:
        "data/25_Sphase_contribution/06_nonb_func_enrichment/chm13v2.0_RepeatMasker.bed"
    log:
        "logs/25_Sphase_contribution/06_nonb_func_enrichment/RM_to_bed.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=4000
    shell:
        "awk 'NR > 3 && NF >= 7 {{print $5 \"\\t\" ($6-1) \"\\t\" $7}}' {input} "
        "| bedtools sort -i stdin | bedtools merge -i stdin > {output} 2> {log}"


rule Sphase_nonb_func_enrichment:
    input:
        reference  = "data/25_Sphase_contribution/05_merged_windows/{sphase_ann}.bed",
        annotation = lambda wc: SPHASE_FUNC_ANNS[wc.sphase_feat],
    output:
        "data/25_Sphase_contribution/06_nonb_func_enrichment/{sphase_ann}/{sphase_feat}_enrichment.bed"
    log:
        "logs/25_Sphase_contribution/06_nonb_func_enrichment/{sphase_ann}/{sphase_feat}_enrichment.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=2000
    wildcard_constraints:
        sphase_ann  = r"[A-Za-z0-9_]+",
        sphase_feat = r"[A-Za-z0-9]+",
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.reference} {input.annotation} {output} 10000 100 2> {log}"


rule Sphase_compile_nonb_func_enrichment:
    input:
        expand(
            "data/25_Sphase_contribution/06_nonb_func_enrichment/{sphase_ann}/{sphase_feat}_enrichment.bed",
            sphase_ann=SPHASE_ANN_NAMES, sphase_feat=SPHASE_FEAT_NAMES,
        )
    output:
        "results/25_Sphase_contribution/Sphase_nonb_func_enrichment.csv"
    log:
        "logs/25_Sphase_contribution/06_nonb_func_enrichment/compile_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/25_Sphase_contribution/compile_Sphase_nonb_func_enrichment.R 2> {log}"


rule Sphase_nonb_func_coverage:
    input:
        reference  = "data/25_Sphase_contribution/05_merged_windows/{sphase_ann}.bed",
        annotation = lambda wc: SPHASE_FUNC_ANNS[wc.sphase_feat],
    output:
        "data/25_Sphase_contribution/06_nonb_func_enrichment/{sphase_ann}/{sphase_feat}_bp_overlap.txt"
    log:
        "logs/25_Sphase_contribution/06_nonb_func_enrichment/{sphase_ann}/{sphase_feat}_coverage.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    wildcard_constraints:
        sphase_ann  = r"[A-Za-z0-9_]+",
        sphase_feat = r"[A-Za-z0-9]+",
    shell:
        "bedtools intersect -a {input.reference} -b {input.annotation} "
        "| awk '{{sum += $3-$2}} END {{print sum+0}}' > {output} 2> {log}"


rule Sphase_compile_nonb_func_coverage:
    input:
        overlaps = expand(
            "data/25_Sphase_contribution/06_nonb_func_enrichment/{sphase_ann}/{sphase_feat}_bp_overlap.txt",
            sphase_ann=SPHASE_ANN_NAMES, sphase_feat=SPHASE_FEAT_NAMES,
        ),
        references = expand(
            "data/25_Sphase_contribution/05_merged_windows/{sphase_ann}.bed",
            sphase_ann=SPHASE_ANN_NAMES,
        ),
    output:
        "results/25_Sphase_contribution/Sphase_nonb_func_coverage.csv"
    log:
        "logs/25_Sphase_contribution/06_nonb_func_enrichment/compile_coverage.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/25_Sphase_contribution/compile_Sphase_nonb_func_coverage.R 2> {log}"


############################
# Phase 3b — Base-pair coverage of ape HMM states by S-phase / G0-G1 windows
#
# For each of the 32 state × annotation combinations, intersect the merged
# annotation BED with the state BED and sum the overlapping base pairs.
# The compile rule also reads state BED total sizes to compute the fraction
# of each state covered.
############################

rule Sphase_HMM_state_coverage:
    input:
        state      = ancient("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/"
                             "hominoidae/8/state_{state25}.bed"),
        annotation = "data/25_Sphase_contribution/05_merged_windows/{sphase_ann}.bed",
    output:
        "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/"
        "State_functional_element_enrichment/state_{state25}/{sphase_ann}_bp_overlap.txt"
    log:
        "logs/25_Sphase_contribution/HMM_coverage/state_{state25}/{sphase_ann}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    wildcard_constraints:
        state25    = r"[0-7]",
        sphase_ann = r"[A-Za-z0-9_]+"
    shell:
        "bedtools intersect -a {input.state} -b {input.annotation} "
        "| awk '{{sum += $3-$2}} END {{print sum+0}}' > {output} 2> {log}"


rule Sphase_compile_HMM_coverage:
    input:
        overlaps = expand(
            "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/"
            "State_functional_element_enrichment/state_{state25}/{sphase_ann}_bp_overlap.txt",
            state25=SPHASE_STATES, sphase_ann=SPHASE_ANN_NAMES,
        ),
        states = ancient(expand(
            "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state25}.bed",
            state25=SPHASE_STATES,
        )),
    output:
        "data/25_Sphase_contribution/Multivariate_guassian_HMM/hominoidae/8/"
        "State_functional_element_enrichment/Sphase_HMM_coverage.csv"
    log:
        "logs/25_Sphase_contribution/HMM_coverage/compile.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=2000
    shell:
        "Rscript workflow/scripts/25_Sphase_contribution/compile_Sphase_HMM_coverage.R 2> {log}"
