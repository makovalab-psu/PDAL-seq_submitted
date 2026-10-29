####################################################################
# Step 21 — Add K562 cell line (reviewer response)
#
# Self-contained chain that processes the K562 PDALseq fastqs from
# Revisions_data/Raw_fastq/260415_VH00707_304_AAJ7FM5M5/fastq/.
#
# Phase 1: technical-replicate correlation across K562 samples and
#          existing Hsapien samples (Technical_replicate_correlation.csv).
# Phase 2: per-experiment BAM merge (PDS-treated 20mM/40mM only) plus
#          K562-augmented versions of the step 04 / 11 / 12 / 16 CSVs.
#
# All outputs are written under data/21_add_K562_to_analysis/ and
# results/21_add_K562_to_analysis/. Existing data/0X_*/ trees are
# referenced via ancient() and never regenerated.
# PDALseq_samples.txt is not modified.
####################################################################

K562_FASTQ_DIR = "Revisions_data/Raw_fastq/260415_VH00707_304_AAJ7FM5M5/fastq"

K562_SAMPLES_ALL = [
    "Hsa_K562_LCL_PDALSeq_js4003_0mM_S1_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001",
]

K562_SAMPLES_FOR_CORR = [
    "Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001",
]

# Experiment names use lowercase "PDALseq" to match the existing
# convention in PDALseq_samples.txt (e.g. Hsapien_HEK293T_PDALseq).
# Sample (fastq stem) values keep capital "PDALSeq" because they are
# literal filenames.
K562_EXPERIMENTS = {
    "Hsapien_K562_PDALseq_0mM": [
        "Hsa_K562_LCL_PDALSeq_js4003_0mM_S1_R1_001",
    ],
    "Hsapien_K562_PDALseq_20mM": [
        "Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",
        "Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001",
    ],
    "Hsapien_K562_PDALseq_40mM": [
        "Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",
        "Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001",
    ],
}

# The 20mM and 40mM replicates are highly correlated and are merged into a
# single experiment for all downstream analyses. The 0mM control is
# processed through preprocessing only.
K562_EXPERIMENTS_DOWNSTREAM = ["Hsapien_K562_PDALseq"]

K562_WINDOWS = ["1000000", "100000", "10000"]
K562_NONB = ["all", "APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]
K562_NONB_NO_ALL = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]
K562_ANNOTATIONS = ANNOTATIONS  # the 17-entry global list from Snakefile
K562_ORIGIN = ["core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0"]
K562_TSS = ["promoter_protCode_singleTrans"]

# Hsapien samples already in PDALseq_samples.txt (loaded in Snakefile as `df`).
HSAPIEN_EXISTING_SAMPLES = list(df.loc[df['Species'] == 'Hsapien', 'Sequencing_sample'])

# Existing experiments referenced by the K562-augmented compile scripts.
ALL_EXISTING_EXPERIMENTS = list(set(df['Experiment']))
HSAPIEN_4_EXPERIMENTS = [
    "Hsapien_CCL86-Raji_PDALseq",
    "Hsapien_GM24385_PDALseq",
    "Hsapien_GM24143_PDALseq",
    "Hsapien_HEK293T_PDALseq",
]
# 21-experiment hardcoded list from the original step-16 compile script.
EXISTING_GLM_EXPERIMENTS_21 = [
    "Ggorilla_F7549_PDALseq", "Ppygmaeus_AG05252_TruSeq", "Pabelii_AG06213_TruSeq",
    "Hsapien_GM24143_PDALseq", "Hsapien_GM24385_TruSeq", "Hsapien_HEK293T_PDALseq",
    "Pabelii_AG06213_PDALseq", "Hsapien_CCL86-Raji_PDALseq", "Ptroglodytes_AG18358_PDALseq",
    "Ssyndactylus_JAMBI_TruSeq", "Ppygmaeus_AG05252_PDALseq", "Ppaniscus_KB14571_PDALseq",
    "Ptroglodytes_AG18354_PDALseq", "Ggorilla_KB3781_PDALseq", "Ggorilla_KB3781_TruSeq",
    "Ppaniscus_KB8711_PDALseq", "Hsapien_GM24385_PDALseq", "Ssyndactylus_JAMBI_PDALseq",
    "Ptroglodytes_AG18354_TruSeq", "Ssyndactylus_Karenina_PDALseq", "Ppaniscus_KB8711_TruSeq",
]


############################
# Preprocessing (mirrors step 01 logic, R1 only, hard-coded Hsapien genome)
############################

rule K562_trim_fastq:
    input:
        K562_FASTQ_DIR + "/{sample}.fastq.gz"
    output:
        temp("data/21_add_K562_to_analysis/01_preprocess/Trimmed_fastq/{sample}.fastq.gz")
    log:
        "logs/21_add_K562_to_analysis/01_preprocess/01_trim_fastq/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=1000
    shell:
        "fastp --thread 1 -i {input} -o {output} -p -P 50 -y -x -3 2> {log}"


rule K562_bwa_mem:
    input:
        fastq = "data/21_add_K562_to_analysis/01_preprocess/Trimmed_fastq/{sample}.fastq.gz",
        genome = ancient("resources/genomes/Hsapien.fasta"),
        amb = ancient("resources/genomes/Hsapien.fasta.amb"),
        ann = ancient("resources/genomes/Hsapien.fasta.ann"),
        pac = ancient("resources/genomes/Hsapien.fasta.pac"),
        bwt = ancient("resources/genomes/Hsapien.fasta.bwt"),
        sa  = ancient("resources/genomes/Hsapien.fasta.sa"),
    output:
        temp("data/21_add_K562_to_analysis/01_preprocess/Mapped_sam/{sample}.sam")
    log:
        "logs/21_add_K562_to_analysis/01_preprocess/02_bwa_mem/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=16000
    threads:
        5
    shell:
        "bwa mem -t {threads} -o {output} {input.genome} {input.fastq} 2> {log}"


rule K562_sort_sam:
    input:
        "data/21_add_K562_to_analysis/01_preprocess/Mapped_sam/{sample}.sam"
    output:
        temp("data/21_add_K562_to_analysis/01_preprocess/Mapped_bam/{sample}.bam")
    log:
        "logs/21_add_K562_to_analysis/01_preprocess/03_sort_sam/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=8000
    shell:
        "samtools sort {input} -T temp/K562_{wildcards.sample} -O bam -o {output} 2> {log}"


rule K562_remove_dup:
    input:
        "data/21_add_K562_to_analysis/01_preprocess/Mapped_bam/{sample}.bam"
    output:
        "data/21_add_K562_to_analysis/01_preprocess/Mapped_deduplicated_bam/{sample}.bam"
    log:
        "logs/21_add_K562_to_analysis/01_preprocess/04_remove_dup/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=8000
    shell:
        "samtools rmdup -s {input} {output} 2> {log}"


rule K562_remove_low_MapQ:
    input:
        "data/21_add_K562_to_analysis/01_preprocess/Mapped_deduplicated_bam/{sample}.bam"
    output:
        "data/21_add_K562_to_analysis/01_preprocess/Mapped_deduplicated_high_MapQ_bam/{sample}.bam"
    log:
        "logs/21_add_K562_to_analysis/01_preprocess/05_remove_low_MapQ/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=8000
    shell:
        "samtools view -b --min-MQ 20 {input} > {output} 2> {log}"


############################
# Technical replicate correlation
############################

rule K562_sort_sample_bam:
    input:
        "data/21_add_K562_to_analysis/01_preprocess/Mapped_deduplicated_bam/{sample}.bam"
    output:
        temp("data/21_add_K562_to_analysis/03_correlation/sorted_bam/{sample}.bam")
    log:
        "logs/21_add_K562_to_analysis/03_correlation/sorted_bam/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "samtools sort -o {output} {input} 2> {log}"


rule K562_make_genome_windows:
    input:
        ancient("resources/genomes/Hsapien.chrom.sizes")
    output:
        "data/21_add_K562_to_analysis/03_correlation/01_genome_windows/{size}_nucleotides/Hsapien.bed"
    log:
        "logs/21_add_K562_to_analysis/03_correlation/01_genome_windows/{size}_nucleotides.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/make_windows.sh {input} {output} {wildcards.size} 2> {log}"


rule K562_genome_window_coverage:
    input:
        bam = "data/21_add_K562_to_analysis/03_correlation/sorted_bam/{sample}.bam",
        windows = "data/21_add_K562_to_analysis/03_correlation/01_genome_windows/{size}_nucleotides/Hsapien.bed"
    output:
        "data/21_add_K562_to_analysis/03_correlation/03_genome_window_coverage/{size}_nucleotides/{sample}.bed"
    log:
        "logs/21_add_K562_to_analysis/03_correlation/03_genome_window_coverage/{size}_nucleotides/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "megadepth {input.bam} --threads 1 --annotation {input.windows} --op sum > {output} 2> {log}"


rule K562_compile_correlation:
    input:
        k562 = expand(
            "data/21_add_K562_to_analysis/03_correlation/03_genome_window_coverage/{size}_nucleotides/{sample}.bed",
            size = K562_WINDOWS, sample = K562_SAMPLES_FOR_CORR
        ),
        hsapien_existing = [
            ancient(f"data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/{size}_nucleotides/{sample}.bed")
            for size in K562_WINDOWS
            for sample in HSAPIEN_EXISTING_SAMPLES
        ],
        samples_txt = ancient("PDALseq_samples.txt"),
    output:
        "results/21_add_K562_to_analysis/Technical_replicate_correlation.csv"
    log:
        "logs/21_add_K562_to_analysis/04_compile_correlation.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/21_add_K562_to_analysis/compile_K562_correlation.R 2> {log}"


############################
# Phase 2 — Per-experiment BAM merge (PDS-treated only)
############################

def k562_per_sample_bams(wildcards):
    return [
        f"data/21_add_K562_to_analysis/01_preprocess/Mapped_deduplicated_bam/{s}.bam"
        for s in K562_EXPERIMENTS[wildcards.experiment]
    ]


# Per-rule wildcard_constraints (instead of a file-global directive) so the
# constraint cannot leak into other .smk files that also use {experiment}.
rule K562_merge_bam_all_reads:
    wildcard_constraints:
        experiment = "Hsapien_K562_PDALseq_(0mM|20mM|40mM)"
    input:
        k562_per_sample_bams
    output:
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/{experiment}.bam"
    log:
        "logs/21_add_K562_to_analysis/01_preprocess/06_merge_all_reads/{experiment}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=8000
    shell:
        "samtools merge -f -o {output} {input} 2> {log}"


rule K562_index_bam_all_reads:
    input:
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/{experiment}.bam"
    output:
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/{experiment}.bam.bai"
    log:
        "logs/21_add_K562_to_analysis/01_preprocess/07_index_all_reads/{experiment}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=4000
    shell:
        "samtools index {input} 2> {log}"


############################
# Phase 2 — Merge 20mM + 40mM into a single experiment BAM
############################

rule K562_merge_concentrations:
    input:
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq_20mM.bam",
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq_40mM.bam"
    output:
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq.bam"
    log:
        "logs/21_add_K562_to_analysis/01_preprocess/08_merge_concentrations.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=8000
    shell:
        "samtools merge -f -o {output} {input} 2> {log}"


############################
# Phase 2 — Step-04 chain (nonB_coverage_versus_read.csv)
############################

# Reuses workflow/scripts/04_*/genome_window_coverage.sh — its species
# inference (first underscore-delimited token of the sample name) yields
# "Hsapien" for Hsapien_K562_PDALseq_*, which matches the cached
# data/03_*/01_genome_windows/{size}/Hsapien.bed referenced inside the
# script. The existing step-03 windows file is required (ancient).
rule K562_window_coverage:
    input:
        bam = "data/21_add_K562_to_analysis/01_preprocess/All_reads/{experiment}.bam",
        windows = ancient("data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/{size}_nucleotides/Hsapien.bed")
    output:
        "data/21_add_K562_to_analysis/04_nonB/04_genome_window_coverage/{size}_nucleotides/{experiment}.bg"
    log:
        "logs/21_add_K562_to_analysis/04_nonB/04_genome_window_coverage/{size}_nucleotides/{experiment}.txt"
    conda:
        "../env/04_Genome_wide_comparison_to_nonB_motifs.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/04_Genome_wide_comparison_to_nonB_motifs/genome_window_coverage.sh "
        "{input.bam} {wildcards.size} {output} {wildcards.experiment} 2> {log}"


rule K562_compile_reads_vs_nonb:
    input:
        k562 = expand(
            "data/21_add_K562_to_analysis/04_nonB/04_genome_window_coverage/{size}_nucleotides/{experiment}.bg",
            size = K562_WINDOWS, experiment = K562_EXPERIMENTS_DOWNSTREAM,
        ),
        existing_cov = [
            ancient(f"data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/{size}_nucleotides/{experiment}.bg")
            for size in K562_WINDOWS
            for experiment in ALL_EXISTING_EXPERIMENTS
        ],
        existing_nonb_density = [
            ancient(f"data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/{size}_nucleotides/{species}/{nb}.bg")
            for size in K562_WINDOWS
            for species in SPECIES
            for nb in K562_NONB
        ],
        samples_txt = ancient("PDALseq_samples.txt"),
    output:
        coef = "results/21_add_K562_to_analysis/nonB_correlation_coefficient.csv",
        cov  = "results/21_add_K562_to_analysis/nonB_coverage_versus_read.csv",
    log:
        "logs/21_add_K562_to_analysis/04_nonB/05_compile_reads_vs_nonb.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=8000
    shell:
        "Rscript workflow/scripts/21_add_K562_to_analysis/compile_K562_reads_vs_nonb.R 2> {log}"


############################
# Phase 2 — Step-16 chain (glm_unique_regions_1000_nucleotides.csv)
############################

# Reuses workflow/scripts/16_*/glm_coverage.sh. Same species-inference
# behavior as step-04: yields "Hsapien" for K562 experiment names. The
# windows BED at data/16_*/03_remove_sex_chromosomes_from_windows/... is
# referenced via ancient() (it pre-exists from prior step-16 runs).
rule K562_glm_coverage:
    input:
        bam = "data/21_add_K562_to_analysis/01_preprocess/All_reads/{experiment}.bam",
        bai = "data/21_add_K562_to_analysis/01_preprocess/All_reads/{experiment}.bam.bai",
        windows = ancient("data/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/Hsapien.bed"),
    output:
        "data/21_add_K562_to_analysis/16_glm/05_glm_coverage/1000_nucleotides/{experiment}.bg"
    log:
        "logs/21_add_K562_to_analysis/16_glm/05_glm_coverage/1000_nucleotides/{experiment}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=5000
    shell:
        "bash workflow/scripts/16_lm_analysis_of_nonbmotifs/glm_coverage.sh "
        "{input.bam} {output} {wildcards.experiment} 2> {log}"


# glm_coverage_matrix.R extracts species from the experiment name (first
# underscore token) and reads the species-level nonB density / GC content
# from data/16_*/04_* and data/16_*/06_*. Those paths are species-keyed
# and pre-exist from the original step-16 run, so the reused script
# resolves them correctly without any wrapper.
rule K562_glm_matrix:
    input:
        nonb_density = [
            ancient(f"data/16_lm_analysis_of_nonbmotifs/04_nonb_densisty_in_windows/1000_nucleotides/Hsapien/{nb}.bed")
            for nb in K562_NONB_NO_ALL
        ],
        gc = ancient("data/16_lm_analysis_of_nonbmotifs/06_glm_GC_content/1000_nucleotides/Hsapien.txt"),
        windows = ancient("data/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/Hsapien.bed"),
        reads = "data/21_add_K562_to_analysis/16_glm/05_glm_coverage/1000_nucleotides/{experiment}.bg",
    output:
        "data/21_add_K562_to_analysis/16_glm/07_glm_matrix/1000_nucleotides/{experiment}.csv"
    log:
        "logs/21_add_K562_to_analysis/16_glm/07_glm_matrix/1000_nucleotides/{experiment}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/glm_coverage_matrix.R "
        "{input.reads} {output} {wildcards.experiment} 2> {log}"


rule K562_glm_unique_regions:
    input:
        "data/21_add_K562_to_analysis/16_glm/07_glm_matrix/1000_nucleotides/{experiment}.csv"
    output:
        "data/21_add_K562_to_analysis/16_glm/07_glm_unique_regions/1000_nucleotides/{experiment}.csv"
    log:
        "logs/21_add_K562_to_analysis/16_glm/07_glm_unique_regions/1000_nucleotides/{experiment}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Reads_vrs_nonb_glm_unique_mapping_regions.R "
        "{input} {output} 2> {log}"


rule K562_compile_glm_unique_regions:
    input:
        k562 = expand(
            "data/21_add_K562_to_analysis/16_glm/07_glm_unique_regions/1000_nucleotides/{experiment}.csv",
            experiment = K562_EXPERIMENTS_DOWNSTREAM,
        ),
        existing = [
            ancient(f"data/16_lm_analysis_of_nonbmotifs/07_glm_unique_regions/1000_nucleotides/{exp}.csv")
            for exp in EXISTING_GLM_EXPERIMENTS_21
        ],
    output:
        "results/21_add_K562_to_analysis/glm_unique_regions_1000_nucleotides.csv"
    log:
        "logs/21_add_K562_to_analysis/16_glm/08_compile_glm_unique_regions.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/21_add_K562_to_analysis/Compile_K562_glm_unique_mapping_regions.R "
        "{output} 2> {log}"


############################
# Phase 2 — Step-11 chain (PDALseq_signal_around_origins.csv)
############################

# The unscrambled / scrambled origin BEDs at
# data/11_*/{shuffle}/state_{state}/{origin}.bed are sample-independent
# and pre-exist from the prior step-11 run; reuse via ancient().
rule K562_Mean_signal_origins:
    input:
        bg = "data/21_add_K562_to_analysis/04_nonB/04_genome_window_coverage/1000_nucleotides/{experiment}.bg",
        annotation = ancient("data/11_Origin_analysis_analysis/{shuffle}/state_{state}/{origin}.bed"),
    output:
        "data/21_add_K562_to_analysis/11_origins/Mean_signal/{experiment}/{shuffle}/state_{state}/{origin}.csv"
    log:
        "logs/21_add_K562_to_analysis/11_origins/Mean_signal/{experiment}/{shuffle}/state_{state}/{origin}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/11_Origin_analysis_analysis/PDAL_seq_around_origins.R "
        "{input.bg} {input.annotation} 100 {output} 2> {log}"


rule K562_compile_origin_signal:
    input:
        k562 = expand(
            "data/21_add_K562_to_analysis/11_origins/Mean_signal/{experiment}/{shuffle}/state_{state}/{origin}.csv",
            experiment = K562_EXPERIMENTS_DOWNSTREAM,
            shuffle = ["unscrambled", "scrambled"],
            state = ["all", "4"],
            origin = K562_ORIGIN,
        ),
        existing = [
            ancient(f"data/11_Origin_analysis_analysis/Mean_signal/{exp}/{sh}/state_{st}/{org}.csv")
            for exp in HSAPIEN_4_EXPERIMENTS
            for sh in ["unscrambled", "scrambled"]
            for st in ["all", "4"]
            for org in K562_ORIGIN
        ],
    output:
        "results/21_add_K562_to_analysis/PDALseq_signal_around_origins.csv"
    log:
        "logs/21_add_K562_to_analysis/11_origins/06_compile_origin_signal.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/21_add_K562_to_analysis/Compile_K562_PDAL_seq_around_origins.R "
        "{output} 2> {log}"


############################
# Phase 2 — Step-12 TSS chain (PDALseq_signal_around_TSS.csv)
############################

rule K562_Mean_TSS_signal:
    input:
        bg = "data/21_add_K562_to_analysis/04_nonB/04_genome_window_coverage/1000_nucleotides/{experiment}.bg",
        annotation = ancient("data/12_TSS_analysis_analysis/{shuffle}/state_{state}/{tss}.bed"),
    output:
        "data/21_add_K562_to_analysis/12_TSS/Mean_signal/{experiment}/{shuffle}/state_{state}/{tss}.csv"
    log:
        "logs/21_add_K562_to_analysis/12_TSS/Mean_signal/{experiment}/{shuffle}/state_{state}/{tss}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/12_TSS_analysis_analysis/PDAL_seq_around_TSS.R "
        "{input.bg} {input.annotation} 100 {output} 2> {log}"


rule K562_compile_TSS_signal:
    input:
        k562 = expand(
            "data/21_add_K562_to_analysis/12_TSS/Mean_signal/{experiment}/{shuffle}/state_{state}/{tss}.csv",
            experiment = K562_EXPERIMENTS_DOWNSTREAM,
            shuffle = ["unscrambled", "scrambled"],
            state = ["all", "0", "1", "2", "3", "4", "5", "6", "7"],
            tss = K562_TSS,
        ),
        existing = [
            ancient(f"data/12_TSS_analysis_analysis/Mean_signal/{exp}/{sh}/state_{st}/{tss}.csv")
            for exp in HSAPIEN_4_EXPERIMENTS
            for sh in ["unscrambled", "scrambled"]
            for st in ["all", "0", "1", "2", "3", "4", "5", "6", "7"]
            for tss in K562_TSS
        ],
    output:
        "results/21_add_K562_to_analysis/PDALseq_signal_around_TSS.csv"
    log:
        "logs/21_add_K562_to_analysis/12_TSS/06_compile_TSS_signal.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/21_add_K562_to_analysis/Compile_K562_PDAL_seq_around_TSS.R "
        "{output} 2> {log}"


############################
# Phase 2 — Step-12 read enrichment chain (Read_enrichment.csv)
############################

rule K562_read_enrichment:
    input:
        bam = "data/21_add_K562_to_analysis/01_preprocess/All_reads/{experiment}.bam",
        bed = ancient("resources/functional_annotations/Saswats_annotations/Hsapien/{annotation}.bed"),
    output:
        "data/21_add_K562_to_analysis/12_read_enrichment/Read_enrichment/{experiment}/{annotation}.csv"
    log:
        "logs/21_add_K562_to_analysis/12_read_enrichment/Read_enrichment/{experiment}/{annotation}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/12_PDALseq_reads_human_functional_analysis/Read_enrichment.sh "
        "{input.bam} {input.bed} {output} 10000 100 2> {log}"


rule K562_compile_read_enrichment:
    input:
        k562 = expand(
            "data/21_add_K562_to_analysis/12_read_enrichment/Read_enrichment/{experiment}/{annotation}.csv",
            experiment = K562_EXPERIMENTS_DOWNSTREAM, annotation = K562_ANNOTATIONS,
        ),
        existing = [
            ancient(f"data/12_PDALseq/Read_enrichment/{exp}/{ann}.csv")
            for exp in HSAPIEN_4_EXPERIMENTS
            for ann in K562_ANNOTATIONS
        ],
    output:
        "results/21_add_K562_to_analysis/Read_enrichment.csv"
    log:
        "logs/21_add_K562_to_analysis/12_read_enrichment/06_compile_read_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/21_add_K562_to_analysis/Compile_K562_read_enrichment.R "
        "{output} 2> {log}"
