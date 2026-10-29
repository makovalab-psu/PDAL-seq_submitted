
# Cell-type-specific HMM — Fibroblasts (Phase 1: state-number selection)
# Fibroblasts: P. paniscus KB8711/KB14571, G. gorilla KB3781/F7549,
#              P. pygmaeus AG05252, P. abelii AG06213

wildcard_constraints:
    states28 = r"\d+",
    state28  = r"\d+"

_STATES10_28 = [str(i) for i in range(10)]
_NONB_28     = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]
_ANNOTS_28   = [
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

rule Fibroblast_extract_columns:
    input:
        ancient("data/07_PhyloHGMP_model/alignment_coverage/hominoidae/alignment_coverage_matrix.csv")
    output:
        "data/28_Fibroblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/Extract_fibroblast_columns.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/28_Fibroblast_specific_HMM/Extract_fibroblast_columns.R {input} {output} 2> {log}"

rule Fibroblast_alignment_PCA:
    input:
        "data/28_Fibroblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        csv   = "data/28_Fibroblast_specific_HMM/PDALseq_PCA_in_alignment.csv",
        scree = "data/28_Fibroblast_specific_HMM/PDALseq_PCA_in_scree_plot.svg"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/alignment_PCA.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/PDALseq_PCA_in_alignments.R {input} {output.csv} {output.scree} hominoidae 2> {log}"

rule Fibroblast_alignment_Kmeans:
    input:
        "data/28_Fibroblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        "data/28_Fibroblast_specific_HMM/PDALseq_kmeans_scree_plot_in_alignment.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/alignment_Kmeans.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/PDALseq_kmeans_in_alignment.R {input} {output} hominoidae 2> {log}"

rule Fibroblast_GHMM:
    input:
        "data/28_Fibroblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        csv = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/{states28}/GHMM_states.csv",
        pkl = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/{states28}/GHMM.pkl"
    conda:
        "../env/hmm.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/GHMM/{states28}.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/GHMM.py {input} {wildcards.states28} {output.csv} {output.pkl} 2> {log}"

rule Fibroblast_GHMM_BIC:
    input:
        expand("data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/{states}/GHMM.pkl", states=HMM_state),
        data = "data/28_Fibroblast_specific_HMM/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/BIC.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/GHMM/BIC.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/28_Fibroblast_specific_HMM/calc_BIC.py data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/ {input.data} {output} 2> {log}"

rule Fibroblast_GHMM_BIC_plot:
    input:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/BIC.csv"
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/BIC_plot.svg"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/GHMM/BIC_plot.txt"
    resources:
        mem_mb=500
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_BIC.R {input} {output} 2> {log}"

####Phase 2 — Functional annotation of the 10-state fibroblast HMM####

rule Fibroblast_annotate_state:
    """Write per-state BED files for the 10-state fibroblast HMM."""
    input:
        bed    = ancient("resources/alignments/hominoidae/Hsapien.bed"),
        states = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/GHMM_states.csv",
    output:
        expand("data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state}.bed",
               state=_STATES10_28)
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/annotate_state.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/28_Fibroblast_specific_HMM/Annotate_state.R "
        "{input.bed} {input.states} "
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10 "
        "2> {log}"

rule Fibroblast_mean_PDALseq_in_states:
    input:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/GHMM_states.csv"
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/Average_PDALseq_in_GHMM_states.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/mean_PDALseq_in_states.txt"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Mean_PDALseq_signal_in_states.R "
        "{input} {output} 2> {log}"

rule Fibroblast_conservation_in_states:
    input:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/GHMM_states.csv"
    output:
        mean       = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/mean_state_conservation.csv",
        individual = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/individual_state_conservation.csv",
        plot       = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/individual_state_conservation.svg",
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/conservation.txt"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Conservation_score.R "
        "{input} {output.mean} {output.individual} {output.plot} 2> {log}"

rule Fibroblast_BP_in_state:
    input:
        expand("data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state}.bed",
               state=_STATES10_28)
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/BP_in_state.csv"
    log:
        "logs/28_Fibroblast_specific_HMM/BP_in_state.txt"
    resources:
        mem_mb=1000
    shell:
        """
        for i in {{0..9}}; do
            total=$(awk '{{sum += $3 - $2 + 1}} END {{print sum+0}}' \
                data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/state_${{i}}.bed)
            echo "${{i}},${{total}}"
        done > {output} 2> {log}
        """

rule Fibroblast_nonB_enrichment_in_state:
    input:
        state = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state28}.bed",
        nonB  = "results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Hsapien/{nonb28}.bed",
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_nonB_enrichment/state_{state28}/{nonb28}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/State_nonB_enrichment/state_{state28}/{nonb28}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.state} {input.nonB} {output} 10000 100 2> {log}"

rule Fibroblast_compile_nonB_enrichment:
    input:
        expand("data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_nonB_enrichment/state_{state}/{nonb}_enrichment.bed",
               state=_STATES10_28, nonb=_NONB_28)
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_nonB_enrichment/compiled.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/compile_nonB_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/28_Fibroblast_specific_HMM/compile_nonB_enrichment_in_state.R "
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10 {output} 2> {log}"

rule Fibroblast_functional_element_enrichment_in_state:
    input:
        state = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state28}.bed",
        ann   = "resources/functional_annotations/Saswats_annotations/Hsapien/{ann28}.bed",
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_functional_element_enrichment/state_{state28}/{ann28}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/State_functional_element_enrichment/state_{state28}/{ann28}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh "
        "{input.state} {input.ann} {output} 10000 100 2> {log}"

rule Fibroblast_compile_functional_element_enrichment:
    input:
        expand("data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_functional_element_enrichment/state_{state}/{ann}_enrichment.bed",
               state=_STATES10_28, ann=_ANNOTS_28)
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_functional_element_enrichment/compiled.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/compile_functional_element_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/28_Fibroblast_specific_HMM/compile_functional_element_enrichment_in_state.R "
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10 {output} 2> {log}"

rule Fibroblast_methylation_in_state:
    input:
        state = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state28}.bed",
        m5C   = "resources/methylationData/Hsapien.bed",
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/methylation/state_{state28}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/methylation/state_{state28}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/methylation_in_annotation.sh "
        "{input.state} {input.m5C} {output} 2> {log}"

rule Fibroblast_compile_methylation:
    input:
        expand("data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/methylation/state_{state}.txt",
               state=_STATES10_28)
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/methylation/compiled.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/compile_methylation.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/28_Fibroblast_specific_HMM/compile_methylation_in_state.R "
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10 {output} 2> {log}"

rule Fibroblast_RNA_expression_in_state:
    input:
        state = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/state_{state28}.bed",
        RNA   = "resources/expression/RNAseq_default.bg",
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_RNA_expression/state_{state28}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/State_RNA_expression/state_{state28}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/BG_enrichment.sh "
        "{input.state} {input.RNA} {output} 1000000 100 2> {log}"

rule Fibroblast_compile_RNA_expression:
    input:
        expand("data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_RNA_expression/state_{state}.txt",
               state=_STATES10_28)
    output:
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_RNA_expression/compiled.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/compile_RNA_expression.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/28_Fibroblast_specific_HMM/compile_RNA_expression_in_state.R "
        "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10 {output} 2> {log}"

rule Fibroblast_SI_figure:
    """Assemble SI figure for the 10-state fibroblast HMM."""
    input:
        avg_pdal     = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/Average_PDALseq_in_GHMM_states.csv",
        bp_state     = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/BP_in_state.csv",
        conservation = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/mean_state_conservation.csv",
        nonb         = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_nonB_enrichment/compiled.csv",
        ann          = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_functional_element_enrichment/compiled.csv",
        rna          = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/State_RNA_expression/compiled.csv",
        methylation  = "data/28_Fibroblast_specific_HMM/Multivariate_guassian_HMM/10/methylation/compiled.csv",
    output:
        "results/28_Fibroblast_specific_HMM/Figure_SI_Fibroblast_HMM.svg"
    conda:
        "../env/25_Sphase_analysis.yml"
    log:
        "logs/28_Fibroblast_specific_HMM/Figure_SI.txt"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/28_Fibroblast_specific_HMM/Figure_SI_Fibroblast_HMM.R 2> {log}"
