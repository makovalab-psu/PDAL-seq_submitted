####################################################################
# Step 30 — Update SI tables with K562 and DMSO revisions data
#
# Produces new rows for Tables S1, S2, and S3 for the two
# revision datasets added during peer review:
#   (1) K562 cell line  — 5 tech reps, preprocessed in step 21
#   (2) DMSO-synchronized GM24143/GM24385 — 6 tech reps, step 25
#
# All step 21/25 BAMs are referenced via ancient().
# No existing data/ or results/ trees are modified.
# Outputs are isolated under:
#   data/30_update_SI_tables_with_revisions_data/
#   results/30_update_SI_tables_with_revisions_data/
#   logs/30_update_SI_tables_with_revisions_data/
#
# Tables produced:
#   Table_S1_new_samples.csv — per-tech-rep preprocessing stats
#   Table_S2_new_samples.csv — per-experiment mapping stats (All_reads
#                               + High_MapQ for K562; All_reads only
#                               for DMSO — step 25 has no MapQ filter)
#   Table_S3_new_samples.csv — unique maps inside/outside Hsapien
#                               blacklist (per merged experiment)
####################################################################

# ---- Sample / experiment lists ---------------------------------- #

S30_K562_FASTQ_DIR = "Revisions_data/Raw_fastq/260415_VH00707_304_AAJ7FM5M5/fastq"
S30_DMSO_FASTQ_DIR = "Revisions_data/Raw_fastq/260512_VH00707_312_AAJ7GVWM5/fastq"

S30_K562_SAMPLES = [
    "Hsa_K562_LCL_PDALSeq_js4003_0mM_S1_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",
    "Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001",
]

S30_DMSO_SAMPLES = [
    "Hsa_GM24143_DMSO_PDALSeq_js4012d20mM_S4_R1_001",
    "Hsa_GM24143_DMSO_PDALSeq_js4012e20mM_S5_R1_001",
    "Hsa_GM24143_DMSO_PDALSeq_js4012f20mM_S6_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012a20mM_S1_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012b20mM_S2_R1_001",
    "Hsa_GM24385_DMSO_PDALSeq_js4012c20mM_S3_R1_001",
]

S30_ALL_SAMPLES = S30_K562_SAMPLES + S30_DMSO_SAMPLES

# Sample → experiment groupings (mirrors step 21 K562_EXPERIMENTS / step 25 DMSO design)
S30_K562_EXPERIMENTS = {
    "Hsapien_K562_PDALseq_0mM":  ["Hsa_K562_LCL_PDALSeq_js4003_0mM_S1_R1_001"],
    "Hsapien_K562_PDALseq_20mM": ["Hsa_K562_LCL_PDALSeq_js4003_20mM1_S2_R1_001",
                                   "Hsa_K562_LCL_PDALSeq_js4003_20mM2_S3_R1_001"],
    "Hsapien_K562_PDALseq_40mM": ["Hsa_K562_LCL_PDALSeq_js4003_40mM1_S4_R1_001",
                                   "Hsa_K562_LCL_PDALSeq_js4003_40mM2_S5_R1_001"],
}
# Ordered list: per-concentration first, then the final merged experiment
S30_K562_EXP_ALL = [
    "Hsapien_K562_PDALseq_0mM",
    "Hsapien_K562_PDALseq_20mM",
    "Hsapien_K562_PDALseq_40mM",
    "Hsapien_K562_PDALseq",
]

S30_DMSO_EXPERIMENTS = {
    "Hsapien_GM24143_DMSO_PDALseq": ["Hsa_GM24143_DMSO_PDALSeq_js4012d20mM_S4_R1_001",
                                      "Hsa_GM24143_DMSO_PDALSeq_js4012e20mM_S5_R1_001",
                                      "Hsa_GM24143_DMSO_PDALSeq_js4012f20mM_S6_R1_001"],
    "Hsapien_GM24385_DMSO_PDALseq": ["Hsa_GM24385_DMSO_PDALSeq_js4012a20mM_S1_R1_001",
                                      "Hsa_GM24385_DMSO_PDALSeq_js4012b20mM_S2_R1_001",
                                      "Hsa_GM24385_DMSO_PDALSeq_js4012c20mM_S3_R1_001"],
}
S30_DMSO_EXP_ALL = list(S30_DMSO_EXPERIMENTS.keys())

# All experiments for Table S2
S30_ALL_EXP = S30_K562_EXP_ALL + S30_DMSO_EXP_ALL

# Experiments for Table S3 (blacklist analysis uses final merged BAMs only)
S30_BLACKLIST_EXP = ["Hsapien_K562_PDALseq"] + S30_DMSO_EXP_ALL

# ---- Input BAM paths for All_reads merged experiments ----------- #
# K562 All_reads BAMs were created by step 21 (ancient); DMSO All_reads
# BAMs are created by S30_DMSO_merge_bam below.
S30_ALL_READS_BAM = {
    "Hsapien_K562_PDALseq_0mM":
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq_0mM.bam",
    "Hsapien_K562_PDALseq_20mM":
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq_20mM.bam",
    "Hsapien_K562_PDALseq_40mM":
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq_40mM.bam",
    "Hsapien_K562_PDALseq":
        "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq.bam",
    "Hsapien_GM24143_DMSO_PDALseq":
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/All_reads/Hsapien_GM24143_DMSO_PDALseq.bam",
    "Hsapien_GM24385_DMSO_PDALseq":
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/All_reads/Hsapien_GM24385_DMSO_PDALseq.bam",
}

def _s30_raw_fastq(wc):
    s = wc.s30sample
    if s in S30_K562_SAMPLES:
        return f"{S30_K562_FASTQ_DIR}/{s}.fastq.gz"
    return f"{S30_DMSO_FASTQ_DIR}/{s}.fastq.gz"

def _s30_dedup_bam(wc):
    s = wc.s30sample
    if s in S30_K562_SAMPLES:
        return ancient(f"data/21_add_K562_to_analysis/01_preprocess/Mapped_deduplicated_bam/{s}.bam")
    return ancient(f"data/25_Sphase_contribution/01_preprocess/Mapped_deduplicated_bam/{s}.bam")

def _s30_all_reads_bam(wc):
    path = S30_ALL_READS_BAM[wc.s30exp]
    if "K562" in wc.s30exp:
        return ancient(path)
    return path


############################
# Section A — Raw fastq read counts (Table S1)
############################

rule S30_count_raw_fastq:
    input:
        _s30_raw_fastq
    output:
        "data/30_update_SI_tables_with_revisions_data/Raw_fastq/{s30sample}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Raw_fastq/{s30sample}.txt"
    wildcard_constraints:
        s30sample = "|".join(S30_ALL_SAMPLES)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=500
    shell:
        "bash workflow/scripts/02_preprocessing_statistics/read_count_fastq.gz.sh {input} {output} 2> {log}"


############################
# Section B — Per-sample dedup/high-MapQ BAM flag stats (Table S1)
############################

rule S30_dedup_bam_flags:
    input:
        _s30_dedup_bam
    output:
        "data/30_update_SI_tables_with_revisions_data/Mapped_deduplicated_bam/{s30sample}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Mapped_deduplicated_bam/{s30sample}.txt"
    wildcard_constraints:
        s30sample = "|".join(S30_ALL_SAMPLES)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=500
    shell:
        "bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"


rule S30_high_mapq_bam_flags_per_sample:
    input:
        ancient("data/21_add_K562_to_analysis/01_preprocess/Mapped_deduplicated_high_MapQ_bam/{s30k562sample}.bam")
    output:
        "data/30_update_SI_tables_with_revisions_data/Mapped_deduplicated_high_MapQ_bam/{s30k562sample}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Mapped_deduplicated_high_MapQ_bam/{s30k562sample}.txt"
    wildcard_constraints:
        s30k562sample = "|".join(S30_K562_SAMPLES)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=500
    shell:
        "bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"


############################
# Section C — Merged BAMs for Table S2
#
# K562 All_reads merged BAMs already exist from step 21 (ancient).
# Created here:
#   • DMSO All_reads  — one merged BAM per cell line (GM24143, GM24385)
#   • K562 High_MapQ  — one merged BAM per concentration + final combined
############################

rule S30_DMSO_merge_bam:
    input:
        lambda wc: ancient(expand(
            "data/25_Sphase_contribution/01_preprocess/Mapped_deduplicated_bam/{sample}.bam",
            sample=S30_DMSO_EXPERIMENTS[wc.s30dmsoexp]
        ))
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/All_reads/{s30dmsoexp}.bam"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_bam/All_reads/{s30dmsoexp}.txt"
    wildcard_constraints:
        s30dmsoexp = "|".join(S30_DMSO_EXPERIMENTS.keys())
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=4000
    shell:
        "samtools merge -f {output} {input} 2> {log}"


rule S30_K562_merge_high_mapq_conc:
    input:
        lambda wc: ancient(expand(
            "data/21_add_K562_to_analysis/01_preprocess/Mapped_deduplicated_high_MapQ_bam/{sample}.bam",
            sample=S30_K562_EXPERIMENTS[wc.s30k562conc]
        ))
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/{s30k562conc}.bam"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/{s30k562conc}.txt"
    wildcard_constraints:
        s30k562conc = "|".join(S30_K562_EXPERIMENTS.keys())
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=4000
    shell:
        "samtools merge -f {output} {input} 2> {log}"


rule S30_K562_merge_high_mapq_combined:
    input:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/Hsapien_K562_PDALseq_20mM.bam",
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/Hsapien_K562_PDALseq_40mM.bam",
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/Hsapien_K562_PDALseq.bam"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/Hsapien_K562_PDALseq.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=4000
    shell:
        "samtools merge -f {output} {input} 2> {log}"


############################
# Section D — Stats on All_reads merged BAMs (Table S2)
############################

rule S30_all_reads_bam_flags:
    input:
        _s30_all_reads_bam
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam_flags/All_reads/{s30exp}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_bam_flags/All_reads/{s30exp}.txt"
    wildcard_constraints:
        s30exp = "|".join(S30_ALL_EXP)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=500
    shell:
        "bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"


rule S30_all_reads_mapping_quality:
    input:
        _s30_all_reads_bam
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_mapping_quality/All_reads/{s30exp}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_mapping_quality/All_reads/{s30exp}.txt"
    wildcard_constraints:
        s30exp = "|".join(S30_ALL_EXP)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=500
    shell:
        "samtools view {input} | cut -f 5 | sort -T temp/ | uniq -c > {output} 2> {log}"


rule S30_all_reads_primary_vs_supplemental:
    input:
        _s30_all_reads_bam
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_primary_vs_supplemental/All_reads/{s30exp}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_primary_vs_supplemental/All_reads/{s30exp}.txt"
    wildcard_constraints:
        s30exp = "|".join(S30_ALL_EXP)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=2000
    shell:
        "bash workflow/scripts/02_preprocessing_statistics/Estimate_unique_maps.sh {input} > {output} 2> {log}"


############################
# Section E — Stats on K562 High_MapQ merged BAMs (Table S2)
############################

rule S30_high_mapq_bam_flags_merged:
    input:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/{s30k562exp}.bam"
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam_flags/High_MapQ/{s30k562exp}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_bam_flags/High_MapQ/{s30k562exp}.txt"
    wildcard_constraints:
        s30k562exp = "|".join(S30_K562_EXP_ALL)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=500
    shell:
        "bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"


rule S30_high_mapq_mapping_quality:
    input:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/{s30k562exp}.bam"
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_mapping_quality/High_MapQ/{s30k562exp}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_mapping_quality/High_MapQ/{s30k562exp}.txt"
    wildcard_constraints:
        s30k562exp = "|".join(S30_K562_EXP_ALL)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=500
    shell:
        "samtools view {input} | cut -f 5 | sort -T temp/ | uniq -c > {output} 2> {log}"


rule S30_high_mapq_primary_vs_supplemental:
    input:
        "data/30_update_SI_tables_with_revisions_data/Merged_bam/High_MapQ/{s30k562exp}.bam"
    output:
        "data/30_update_SI_tables_with_revisions_data/Merged_primary_vs_supplemental/High_MapQ/{s30k562exp}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Merged_primary_vs_supplemental/High_MapQ/{s30k562exp}.txt"
    wildcard_constraints:
        s30k562exp = "|".join(S30_K562_EXP_ALL)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=2000
    shell:
        "bash workflow/scripts/02_preprocessing_statistics/Estimate_unique_maps.sh {input} > {output} 2> {log}"


############################
# Section F — Blacklist filtering for Table S3
#
# Uses the Hsapien blacklist from step 03 (ancient).
# filter_reads_by_blacklist.sh infers species from the sample name, which
# fails for "Hsa_" prefix samples. Instead the bedtools -L flag is called
# directly here with the explicit Hsapien blacklist path.
############################

rule S30_filter_blacklist:
    input:
        bam = lambda wc: (
            ancient(S30_ALL_READS_BAM[wc.s30blexp])
            if "K562" in wc.s30blexp
            else S30_ALL_READS_BAM[wc.s30blexp]
        ),
        blacklist = ancient("results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/Hsapien.bed")
    output:
        inside  = temp("data/30_update_SI_tables_with_revisions_data/Blacklist/inside/{s30blexp}.bam"),
        outside = temp("data/30_update_SI_tables_with_revisions_data/Blacklist/outside/{s30blexp}.bam"),
    log:
        "logs/30_update_SI_tables_with_revisions_data/Blacklist/{s30blexp}.txt"
    wildcard_constraints:
        s30blexp = "|".join(S30_BLACKLIST_EXP)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=2000
    shell:
        "samtools view -b -L {input.blacklist} -U {output.outside} {input.bam} > {output.inside} 2> {log}"


rule S30_estimate_unique_maps_blacklist:
    input:
        "data/30_update_SI_tables_with_revisions_data/Blacklist/{s30partition}/{s30blexp}.bam"
    output:
        "data/30_update_SI_tables_with_revisions_data/Blacklist_unique_maps/{s30partition}/{s30blexp}.txt"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Blacklist_unique_maps/{s30partition}/{s30blexp}.txt"
    wildcard_constraints:
        s30partition = "inside|outside",
        s30blexp     = "|".join(S30_BLACKLIST_EXP)
    conda:
        "../env/01_preprocess_data.yml"
    resources:
        mem_mb=2000
    shell:
        "bash workflow/scripts/02_preprocessing_statistics/Estimate_unique_maps.sh {input} > {output} 2> {log}"


############################
# Section G — Compile rules
############################

rule S30_compile_Table_S1:
    input:
        raw      = expand("data/30_update_SI_tables_with_revisions_data/Raw_fastq/{sample}.txt",
                          sample=S30_ALL_SAMPLES),
        dedup    = expand("data/30_update_SI_tables_with_revisions_data/Mapped_deduplicated_bam/{sample}.txt",
                          sample=S30_ALL_SAMPLES),
        high_mapq= expand("data/30_update_SI_tables_with_revisions_data/Mapped_deduplicated_high_MapQ_bam/{sample}.txt",
                          sample=S30_K562_SAMPLES),
    output:
        "results/30_update_SI_tables_with_revisions_data/Table_S1_new_samples.csv"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Compile_Table_S1.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/30_update_SI_tables_with_revisions_data/Compile_Table_S1.R {output} 2> {log}"


rule S30_compile_Table_S2:
    input:
        ar_flags = expand("data/30_update_SI_tables_with_revisions_data/Merged_bam_flags/All_reads/{exp}.txt",
                          exp=S30_ALL_EXP),
        ar_qual  = expand("data/30_update_SI_tables_with_revisions_data/Merged_mapping_quality/All_reads/{exp}.txt",
                          exp=S30_ALL_EXP),
        ar_umr   = expand("data/30_update_SI_tables_with_revisions_data/Merged_primary_vs_supplemental/All_reads/{exp}.txt",
                          exp=S30_ALL_EXP),
        hq_flags = expand("data/30_update_SI_tables_with_revisions_data/Merged_bam_flags/High_MapQ/{exp}.txt",
                          exp=S30_K562_EXP_ALL),
        hq_qual  = expand("data/30_update_SI_tables_with_revisions_data/Merged_mapping_quality/High_MapQ/{exp}.txt",
                          exp=S30_K562_EXP_ALL),
        hq_umr   = expand("data/30_update_SI_tables_with_revisions_data/Merged_primary_vs_supplemental/High_MapQ/{exp}.txt",
                          exp=S30_K562_EXP_ALL),
    output:
        "results/30_update_SI_tables_with_revisions_data/Table_S2_new_samples.csv"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Compile_Table_S2.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/30_update_SI_tables_with_revisions_data/Compile_Table_S2.R {output} 2> {log}"


rule S30_compile_Table_S3:
    input:
        inside  = expand("data/30_update_SI_tables_with_revisions_data/Blacklist_unique_maps/inside/{exp}.txt",
                         exp=S30_BLACKLIST_EXP),
        outside = expand("data/30_update_SI_tables_with_revisions_data/Blacklist_unique_maps/outside/{exp}.txt",
                         exp=S30_BLACKLIST_EXP),
    output:
        "results/30_update_SI_tables_with_revisions_data/Table_S3_new_samples.csv"
    log:
        "logs/30_update_SI_tables_with_revisions_data/Compile_Table_S3.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/30_update_SI_tables_with_revisions_data/Compile_Table_S3.R {output} 2> {log}"
