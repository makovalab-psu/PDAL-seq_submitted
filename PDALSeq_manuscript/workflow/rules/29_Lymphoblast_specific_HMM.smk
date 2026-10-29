
# Cell-type-specific HMM — Lymphoblasts (Phase 1: state-number selection)
# Lymphoblasts: H. sapiens GM24385/GM24143, P. troglodytes AG18354/AG18358,
#               S. syndactylus JAMBI/Karenina

wildcard_constraints:
    states29 = r"\d+",
    state29  = r"\d+"

_STATES10_29 = [str(i) for i in range(10)]
_NONB_29     = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]
_ANNOTS_29   = [
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

rule Lymphoblast_extract_columns:
    input:
        ancient("data/07_PhyloHGMP_model/alignment_coverage/hominoidae/alignment_coverage_matrix.csv")
    output:
        "data/29_Lymphoblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/Extract_lymphoblast_columns.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/29_Lymphoblast_specific_HMM/Extract_lymphoblast_columns.R {input} {output} 2> {log}"

rule Lymphoblast_alignment_PCA:
    input:
        "data/29_Lymphoblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        csv   = "data/29_Lymphoblast_specific_HMM/PDALseq_PCA_in_alignment.csv",
        scree = "data/29_Lymphoblast_specific_HMM/PDALseq_PCA_in_scree_plot.svg"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/alignment_PCA.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/PDALseq_PCA_in_alignments.R {input} {output.csv} {output.scree} hominoidae 2> {log}"

rule Lymphoblast_alignment_Kmeans:
    input:
        "data/29_Lymphoblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        "data/29_Lymphoblast_specific_HMM/PDALseq_kmeans_scree_plot_in_alignment.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/alignment_Kmeans.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/PDALseq_kmeans_in_alignment.R {input} {output} hominoidae 2> {log}"

rule Lymphoblast_GHMM:
    input:
        "data/29_Lymphoblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        csv = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/{states29}/GHMM_states.csv",
        pkl = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/{states29}/GHMM.pkl"
    conda:
        "../env/hmm.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/GHMM/{states29}.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/GHMM.py {input} {wildcards.states29} {output.csv} {output.pkl} 2> {log}"

rule Lymphoblast_GHMM_BIC:
    input:
        expand("data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/{states}/GHMM.pkl", states=HMM_state),
        data = "data/29_Lymphoblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/BIC.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/GHMM/BIC.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/29_Lymphoblast_specific_HMM/calc_BIC.py data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/ {input.data} {output} 2> {log}"

rule Lymphoblast_GHMM_BIC_plot:
    input:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/BIC.csv"
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/BIC_plot.svg"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/GHMM/BIC_plot.txt"
    resources:
        mem_mb=500
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_BIC.R {input} {output} 2> {log}"

####Phase 2 — Functional annotation of the 10-state lymphoblast HMM####

rule Lymphoblast_annotate_state:
    """Write per-state BED files for the 10-state lymphoblast HMM."""
    input:
        bed    = ancient("resources/alignments/hominoidae/Hsapien.bed"),
        states = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/GHMM_states.csv",
    output:
        expand("data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state}.bed",
               state=_STATES10_29)
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/annotate_state.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/29_Lymphoblast_specific_HMM/Annotate_state.R "
        "{input.bed} {input.states} "
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10 "
        "2> {log}"

rule Lymphoblast_mean_PDALseq_in_states:
    input:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/GHMM_states.csv"
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/Average_PDALseq_in_GHMM_states.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/mean_PDALseq_in_states.txt"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Mean_PDALseq_signal_in_states.R "
        "{input} {output} 2> {log}"

rule Lymphoblast_conservation_in_states:
    input:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/GHMM_states.csv"
    output:
        mean       = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/mean_state_conservation.csv",
        individual = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/individual_state_conservation.csv",
        plot       = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/individual_state_conservation.svg",
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/conservation.txt"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Conservation_score.R "
        "{input} {output.mean} {output.individual} {output.plot} 2> {log}"

rule Lymphoblast_BP_in_state:
    input:
        expand("data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state}.bed",
               state=_STATES10_29)
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/BP_in_state.csv"
    log:
        "logs/29_Lymphoblast_specific_HMM/BP_in_state.txt"
    resources:
        mem_mb=1000
    shell:
        """
        for i in {{0..9}}; do
            total=$(awk '{{sum += $3 - $2 + 1}} END {{print sum+0}}' \
                data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/state_${{i}}.bed)
            echo "${{i}},${{total}}"
        done > {output} 2> {log}
        """

rule Lymphoblast_nonB_enrichment_in_state:
    input:
        state = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state29}.bed",
        nonB  = "results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Hsapien/{nonb29}.bed",
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_nonB_enrichment/state_{state29}/{nonb29}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/State_nonB_enrichment/state_{state29}/{nonb29}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.state} {input.nonB} {output} 10000 100 2> {log}"

rule Lymphoblast_compile_nonB_enrichment:
    input:
        expand("data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_nonB_enrichment/state_{state}/{nonb}_enrichment.bed",
               state=_STATES10_29, nonb=_NONB_29)
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_nonB_enrichment/compiled.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/compile_nonB_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/29_Lymphoblast_specific_HMM/compile_nonB_enrichment_in_state.R "
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10 {output} 2> {log}"

rule Lymphoblast_functional_element_enrichment_in_state:
    input:
        state = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state29}.bed",
        ann   = "resources/functional_annotations/Saswats_annotations/Hsapien/{ann29}.bed",
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_functional_element_enrichment/state_{state29}/{ann29}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/State_functional_element_enrichment/state_{state29}/{ann29}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.state} {input.ann} {output} 10000 100 2> {log}"

rule Lymphoblast_compile_functional_element_enrichment:
    input:
        expand("data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_functional_element_enrichment/state_{state}/{ann}_enrichment.bed",
               state=_STATES10_29, ann=_ANNOTS_29)
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_functional_element_enrichment/compiled.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/compile_functional_element_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/29_Lymphoblast_specific_HMM/compile_functional_element_enrichment_in_state.R "
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10 {output} 2> {log}"

rule Lymphoblast_methylation_in_state:
    input:
        state = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state29}.bed",
        m5C   = "resources/methylationData/Hsapien.bed",
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/methylation/state_{state29}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/methylation/state_{state29}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/methylation_in_annotation.sh "
        "{input.state} {input.m5C} {output} 2> {log}"

rule Lymphoblast_compile_methylation:
    input:
        expand("data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/methylation/state_{state}.txt",
               state=_STATES10_29)
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/methylation/compiled.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/compile_methylation.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/29_Lymphoblast_specific_HMM/compile_methylation_in_state.R "
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10 {output} 2> {log}"

rule Lymphoblast_RNA_expression_in_state:
    input:
        state = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state29}.bed",
        RNA   = "resources/expression/RNAseq_default.bg",
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_RNA_expression/state_{state29}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/State_RNA_expression/state_{state29}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/BG_enrichment.sh "
        "{input.state} {input.RNA} {output} 1000000 100 2> {log}"

rule Lymphoblast_compile_RNA_expression:
    input:
        expand("data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_RNA_expression/state_{state}.txt",
               state=_STATES10_29)
    output:
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_RNA_expression/compiled.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/compile_RNA_expression.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/29_Lymphoblast_specific_HMM/compile_RNA_expression_in_state.R "
        "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10 {output} 2> {log}"

rule Lymphoblast_SI_figure:
    """Assemble SI figure for the 10-state lymphoblast HMM."""
    input:
        avg_pdal     = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/Average_PDALseq_in_GHMM_states.csv",
        bp_state     = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/BP_in_state.csv",
        conservation = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/mean_state_conservation.csv",
        nonb         = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_nonB_enrichment/compiled.csv",
        ann          = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_functional_element_enrichment/compiled.csv",
        rna          = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/State_RNA_expression/compiled.csv",
        methylation  = "data/29_Lymphoblast_specific_HMM/Multivariate_guassian_HMM/10/methylation/compiled.csv",
    output:
        "results/29_Lymphoblast_specific_HMM/Figure_SI_Lymphoblast_HMM.svg"
    conda:
        "../env/25_Sphase_analysis.yml"
    log:
        "logs/29_Lymphoblast_specific_HMM/Figure_SI.txt"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/29_Lymphoblast_specific_HMM/Figure_SI_Lymphoblast_HMM.R 2> {log}"
