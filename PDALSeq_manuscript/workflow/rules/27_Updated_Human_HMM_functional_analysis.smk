####################################################################
# Step 27 — Updated Human HMM functional analysis
#
# Re-runs the step-09 functional annotation pipeline using the
# K562-inclusive 10-state Human GHMM fitted in step 26.
#
# Input GHMM:
#   data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/10/
#
# Coordinate reference for state BEDs (same 1 kb window set as step 09):
#   data/09_Human_HMM/01_Remove_sex_chromosomes/Hsapien_GM24385_PDALseq.bg
#
# Baseline RM intersect reused from step 12 with ancient() — same windows:
#   data/12_PDALseq/baseline_human_state_rm_intersect_bp.csv
#
# RM annotation BEDs reused from step 12 with ancient():
#   data/12_PDALseq/rm_annotation/{rm_class}.bed
#
# RM family BEDs reused from step 14 with ancient():
#   data/14_Repeat_masker_class_enrichment/01_compile_annotations/{family}.bed
#
# Requires steps 12 and 14 to be included in the Snakefile so that
# RM_CLASS and RM_FAMILY are defined.
#
# All upstream files from steps 26, 09, 12, and 14: ancient().
# All outputs are under data/27_*/ or results/27_*/.
####################################################################

_STATES_27 = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]

_ANNOTATIONS_27 = [
    "cds_protCode_singleTrans", "ngnr_protCode_singleTrans", "chrG_cpgislands",
    "promoter_protCode_singleTrans", "enhancer_protCode_singleTrans",
    "repeats_protCode_singleTrans", "exon_protCode_singleTrans",
    "rna_protCode_singleTrans", "gene_nonprotCode", "utr3_protCode_singleTrans",
    "gene_protCode_singleTrans", "utr5_protCode_singleTrans",
    "intron_protCode_singleTrans", "Recombination",
    "core_stochastic_origins_chm13v2.0_hglft.ucsc",
    "core_origins_hglft_genome_e0596_1cc9e0",
    "stochastic_origins_hglft_genome_2ff9ed_1ccfd0",
]


############################
# State annotation
############################

rule annotate_human_state_27:
    input:
        bed    = ancient("data/09_Human_HMM/01_Remove_sex_chromosomes/Hsapien_GM24385_PDALseq.bg"),
        states = ancient("data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/10/GHMM_states.csv"),
    output:
        expand("data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed", state=_STATES_27)
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/annotate_human_state.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/27_Updated_Human_HMM_functional_analysis/Annotate_state_27.R "
        "{input.bed} {input.states} 2> {log}"

rule count_bp_in_human_state_27:
    input:
        expand("data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed", state=_STATES_27)
    output:
        "data/27_Updated_Human_HMM_functional_analysis/BP_in_state.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/BP_in_state.txt"
    conda:
        "../env/hmm.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/27_Updated_Human_HMM_functional_analysis/BP_in_state_27.sh > {output} 2> {log}"

rule mean_human_PDALseq_in_states_27:
    input:
        ancient("data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/10/GHMM_states.csv")
    output:
        "data/27_Updated_Human_HMM_functional_analysis/Average_PDALseq_in_GHMM_states.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/mean_human_PDALseq_in_states.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Mean_PDALseq_signal_in_states.R {input} {output} 2> {log}"


############################
# Non-B DNA motif enrichment
############################

rule human_nonB_enrichment_in_state_27:
    input:
        state = "data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed",
        nonB  = ancient("results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Hsapien/{nonb}.bed"),
    output:
        "data/27_Updated_Human_HMM_functional_analysis/State_nonB_enrichment/state_{state}/{nonb}_enrichment.bed"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/State_nonB_enrichment/state_{state}/{nonb}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.state} {input.nonB} {output} 10000 100 2> {log}"

rule human_compile_nonB_enrichment_in_state_27:
    input:
        expand(
            "data/27_Updated_Human_HMM_functional_analysis/State_nonB_enrichment/state_{state}/{nonb}_enrichment.bed",
            state=_STATES_27, nonb=NONB,
        )
    output:
        "data/27_Updated_Human_HMM_functional_analysis/State_nonB_enrichment/hsapiens.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/State_nonB_enrichment/hsapiens.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/27_Updated_Human_HMM_functional_analysis/compile_nonB_enrichment_in_state_27.R 2> {log}"


############################
# Functional element enrichment
############################

rule human_functional_element_enrichment_in_state_27:
    input:
        state              = "data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed",
        functional_element = ancient("resources/functional_annotations/Saswats_annotations/Hsapien/{ann}.bed"),
    output:
        "data/27_Updated_Human_HMM_functional_analysis/State_functional_element_enrichment/state_{state}/{ann}_enrichment.bed"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/State_functional_element_enrichment/state_{state}/{ann}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.state} {input.functional_element} {output} 10000 100 2> {log}"

rule human_compile_functional_element_enrichment_in_state_27:
    input:
        expand(
            "data/27_Updated_Human_HMM_functional_analysis/State_functional_element_enrichment/state_{state}/{ann}_enrichment.bed",
            state=_STATES_27, ann=_ANNOTATIONS_27,
        )
    output:
        "data/27_Updated_Human_HMM_functional_analysis/State_functional_element_enrichment/hsapiens.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/State_functional_element_enrichment/hsapiens.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/27_Updated_Human_HMM_functional_analysis/compile_functional_element_enrichment_in_state_27.R 2> {log}"


############################
# Methylation enrichment
############################

rule human_methylation_in_state_27:
    input:
        state = "data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed",
        m5C   = ancient("resources/methylationData/Hsapien.bed"),
    output:
        "data/27_Updated_Human_HMM_functional_analysis/methylation/state_{state}.txt"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/methylation/state_{state}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/methylation_in_annotation.sh "
        "{input.state} {input.m5C} {output} 2> {log}"

rule human_compile_methylation_enrichment_in_state_27:
    input:
        expand("data/27_Updated_Human_HMM_functional_analysis/methylation/state_{state}.txt", state=_STATES_27)
    output:
        "data/27_Updated_Human_HMM_functional_analysis/methylation/hsapiens.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/methylation/hsapiens.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/27_Updated_Human_HMM_functional_analysis/compile_methylation_in_state_27.R 2> {log}"


############################
# RNA expression enrichment
############################

rule human_RNA_expression_in_state_27:
    input:
        state = "data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed",
        RNA   = ancient("resources/expression/RNAseq_default.bg"),
    output:
        "data/27_Updated_Human_HMM_functional_analysis/State_RNA_expression/state_{state}.txt"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/State_RNA_expression/state_{state}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/BG_enrichment.sh "
        "{input.state} {input.RNA} {output} 1000000 100 2> {log}"

rule human_compile_RNA_expression_in_state_27:
    input:
        expand("data/27_Updated_Human_HMM_functional_analysis/State_RNA_expression/state_{state}.txt", state=_STATES_27)
    output:
        "data/27_Updated_Human_HMM_functional_analysis/State_RNA_expression/hsapiens.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/State_RNA_expression/hsapiens.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/27_Updated_Human_HMM_functional_analysis/compile_RNA_expression_in_state_27.R 2> {log}"


############################
# RepeatMasker bp intersect per state
############################

rule human_state_rm_intersect_bp_27:
    input:
        state = "data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed",
        ann   = ancient("data/12_PDALseq/rm_annotation/{rm_class}.bed"),
    output:
        "data/27_Updated_Human_HMM_functional_analysis/human_state_rm_intersect_bp/{rm_class}/{state}.bed"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/human_state_rm_intersect_bp/{rm_class}/{state}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/12_PDALseq_reads_human_functional_analysis/Intersect.sh "
        "{input.state} {input.ann} {output} 2> {log}"

rule compile_human_state_rm_intersect_bp_27:
    input:
        expand(
            "data/27_Updated_Human_HMM_functional_analysis/human_state_rm_intersect_bp/{rm_class}/{state}.bed",
            state=_STATES_27, rm_class=RM_CLASS,
        )
    output:
        "data/27_Updated_Human_HMM_functional_analysis/Human_state_rm_intersect_bp.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/Human_state_rm_intersect_bp.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/27_Updated_Human_HMM_functional_analysis/Compile_human_state_rm_intersect_bp_27.R {output} 2> {log}"


############################
# RepeatMasker family enrichment per state
############################

rule rm_family_enrichment_in_state_27:
    input:
        state              = "data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed",
        functional_element = ancient("data/14_Repeat_masker_class_enrichment/01_compile_annotations/{family}.bed"),
    output:
        "data/27_Updated_Human_HMM_functional_analysis/RM_family_enrichment/state_{state}/{family}_enrichment.bed"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/RM_family_enrichment/state_{state}/{family}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.state} {input.functional_element} {output} 10000 100 2> {log}"

rule compile_rm_family_enrichment_in_state_27:
    input:
        expand(
            "data/27_Updated_Human_HMM_functional_analysis/RM_family_enrichment/state_{state}/{family}_enrichment.bed",
            state=_STATES_27, family=RM_FAMILY,
        )
    output:
        "results/27_Updated_Human_HMM_functional_analysis/RM_family_enrichment.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/RM_family_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/27_Updated_Human_HMM_functional_analysis/compile_RM_family_enrichment_27.R {output} 2> {log}"


############################
# Intersect human states with ape HMM states (step 07, hominoidae 8-state)
############################

rule intersect_human_and_ape_hmm_27:
    input:
        expand("data/27_Updated_Human_HMM_functional_analysis/state_{state}.bed", state=_STATES_27),
        ancient(expand(
            "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed",
            state=["0", "1", "2", "3", "4", "5", "6", "7"],
        )),
    output:
        "data/27_Updated_Human_HMM_functional_analysis/intersect_results.csv"
    log:
        "logs/27_Updated_Human_HMM_functional_analysis/intersect_results.txt"
    conda:
        "../env/Pyutility_plus_bedtools.yml"
    resources:
        mem_mb=1000
    shell:
        "python workflow/scripts/27_Updated_Human_HMM_functional_analysis/intersect_state_beds_27.py 2> {log}"
