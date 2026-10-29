rule Figure_2:
    input:
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonB_correlation_coeficient.csv",
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv",
        expand("resources/centromeres/{species}.bed", species = SPECIES),
        expand("results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/{species}.bed", species = SPECIES),
        "results/04_Genome_wide_comparison_to_nonB_motifs/Ssyndactylus_subTerm_monhorSF4.bed",
        "results/06_reference_free_repeat_analysis/Short_tandem_repeats_by_technical_replicate.csv"
    output:
        "plots/Figure_2_Genome_wide_plots/Figure_2.svg"
    log:
        "logs/Figures_Tables/Figure_2_Genome_wide_comparison_to_nonb_motifs.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_2_Genome_wide_comparison_to_nonb_motifs.R 2> {log}"

rule Figure_3:
    input:
        "data/07_PhyloHGMP_model/alignment_clustering/hominoidae/PDALseq_correllation_in_alignment.csv",    #Correllation
        "data/07_PhyloHGMP_model/alignment_clustering/hominoidae/PDALseq_PCA_in_alignment.csv",    #PCA
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Average_PDALseq_in_GHMM_states.csv",    #Average PDALseq
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/BP_in_state.csv",    #BP
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/mean_state_conservation.csv",    #Conservation
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_nonB_enrichment/hsapiens.csv",    #nonB annotations
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/hsapiens.csv", #Functional annotations
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/hsapiens.csv", #RNA expression
        "results/10_Methylation_analysis/Methylation_by_state_and_species.csv" #Methylation,
        "results/11_Origin_analysis_analysis/PDALseq_signal_around_origins.csv" #Origin analysis
    output:
        "plots/Figure_3_Phylogenetic_PDALseq/Figure_3.svg"
    log:
        "logs/Figure_3_Phylogenetic_PDALseq/Figure_3.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_3_Phylogenetic_analysis_of_PDAL-seq_reads.R 2> {log}"

rule Figure_4:
    input:
        "data/09_Human_HMM/04_Multivariate_guassian/10/Average_PDALseq_in_GHMM_states.csv",
        "data/09_Human_HMM/04_Multivariate_guassian/10/BP_in_state.csv",
        "data/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/hsapiens.csv",
        "data/09_Human_HMM/04_Multivariate_guassian/10/State_functional_element_enrichment/hsapiens.csv",
        "data/09_Human_HMM/04_Multivariate_guassian/10/intersect_results.csv"
    output:
        "plots/Figure_4_human_HMM_PDALseq/Figure_4.svg"
    log:
        "logs/Figure_4_human_HMM_PDALseq/Figure_4.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_4_Analysis_of_PDAL-seq_reads.R 2> {log}"

rule Table_S1:
    input:
        "results/02_preprocessing_statistics/Bam_flags_by_technical_replicate.csv"
    output:
        "tables/Table_S1_preprocessing_statistics.csv"
    log:
        "logs/Figures_Tables/Table_S1_preprocessing_statistics.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Table_S1_preprocessing_statistics.R 2> {log}"

rule Table_S2:
    input:
        "results/02_preprocessing_statistics/Mapping_statistics_by_biological_sample.csv"
    output:
        "tables/Table_S2_preprocessing_statistics_by_sample.csv"
    log:
        "logs/Figures_Tables/Table_S2_preprocessing_statistics_by_sample.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Table_S2_preprocessing_statistics_by_sample.R 2> {log}"

rule Table_S3:
    input:
        "results/03_Correlation_and_black_list_unmappable_regions/Unique_maps_inside_and_outside_blacklist.csv"
    output:
        "tables/Table_S3_Unique_maps_inside_and_outside_blacklist.csv"
    log:
        "logs/Figures_Tables/Table_S3_Unique_maps_inside_and_outside_blacklist.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Table_S3_unique_maps_inside_and_outside_blacklist.R 2> {log}"

rule Table_S4:
    input:
        "results/05_nonB_motif_cluster_analysis/linear_models.csv"
    output:
        "tables/Table_S4_Linear_model_parameters_and_test_for_multicolinearity.csv"
    log:
        "logs/Figures_Tables/Table_S4_Linear_model_parameters_and_test_for_multicolinearity.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Table_S4_Linear_modeling_results_and_test_for_multicolinearity.R 2> {log}"

rule Table_S7:
    input:
        expand("data/07_PhyloHGMP_model/05_Unique_maps_inside_and_outside_alignment/{taxon}/Unique_maps_inside_and_outside_alignment.csv", taxon = TAXON)
    output:
        "tables/Table_S7_Unique_maps_inside_alignment.csv"
    log:
        "logs/Figures_Tables/Table_S7_unique_mapping_rates_in_alingments.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Table_S7_unique_mapping_rates_in_alingments.R 2> {log}"

rule Table_SX_Mcfaddens_R2_glm_for_reads_vs_nonB_1Mbp_windows:
    input:
        "results/16_lm_analysis_of_nonbmotifs/Genome_wide_glm.csv"
    output:
        "tables/Table_SX_Mcfaddens_R2_glm_for_reads_vs_nonB_1Mbp_windows.csv"
    log:
        "logs/Figures_Tables/Table_SX_Mcfaddens_R2_glm_for_reads_vs_nonB_1Mbp_windows.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Table_SX_Mcfaddens_R2_glm_for_reads_vs_nonB_1Mbp_windows.R {input} {output} 2> {log}"


rule Figure_S1:
    input:
        "results/03_Correlation_and_black_list_unmappable_regions/Technical_replicate_correlation.csv"
    output:
        "plots/Figure_S1_technical_replicate_correlation/Figure_S1.svg"
    log:
        "logs/Figures_Tables/Figure_S1_technical_replicate_correlation.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S1_technical_replicate_correlation.R 2> {log}"

rule Figure_S2:
    input:
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv",
        expand("resources/centromeres/{species}.bed", species = SPECIES),
        expand("results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/{species}.bed", species = SPECIES),
        "results/04_Genome_wide_comparison_to_nonB_motifs/Ssyndactylus_subTerm_monhorSF4.bed"
    output:
        expand("plots/Figure_S2_Genome_wide_plots/{species}.svg", species = SPECIES)
    log:
        "logs/Figures_Tables/Figure_S2_Genome_wide_plots.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S2_Genome_wide_plots.R 2> {log}"

rule Figure_S3:
    input:
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv"
    output:
        expand("plots/Figure_S3_Genome_wide_correlation_plots/{data}.svg", data = DATA)
    log:
        "logs/Figures_Tables/Figure_S3_Genome_wide_correlation_plots.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S3_Genome_wide_correlation_plots.R 2> {log}"

rule Figure_S4:
    input:
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv",
        "resources/centromeres/Ssyndactylus.bed",
        "results/04_Genome_wide_comparison_to_nonB_motifs/Ssyndactylus_subTerm_monhorSF4.bed"
    output:
        "plots/Figure_S4_Siamang_subterm_alpha_sat/Figure_S4.svg"
    log:
        "logs/Figures_Tables/Figure_S4_Siamang_subterm_alpha_sat.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S4_Analysis_of_PDALseq_reads_in_subterminal_Siamang_alpha_sats.R {output} 2> {log}"


rule Figure_S5:
    input:
        "results/08_CHIP_CutRun_Benchmarking/Genome_wide_coverage.csv"
    output:
        "plots/Figure_S5_CHIP_CutRun_Benchmarking/Figure_S5.svg"
    log:
        "logs/Figure_S5_CHIP_CutRun_Benchmarking/Figure_S5.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S5_CHIP_CutRun_Benchmarking.R {output} 2> {log}"

rule Figure_S6:
    input:
        "results/05_nonB_motif_cluster_analysis/Cluster_count_vs_cluster_size.txt",
        "results/05_nonB_motif_cluster_analysis/Filtered_cluster_count_vs_cluster_size.txt"
    output:
        "plots/Figure_S6_Number_of_non-B_DNA_cluster_motifs_versus_clustering_distance/Figure_S6.svg"
    log:
        "logs/Figures_Tables/Figure_S6_Number_of_non-B_DNA_cluster_motifs_versus_clustering_distance.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S6_Number_of_non-B_DNA_cluster_motifs_versus_clustering_distance.R 2> {log}"

rule Figure_S7:
    input:

        "results/05_nonB_motif_cluster_analysis/linear_models.csv"
    output:
        "plots/Figure_S7_linear_model_analysis_of_nonB_motif_clusters/Figure_S7.svg"
    log:
        "logs/Figure_S7_linear_model_analysis_of_nonB_motif_clusters/Figure_S7.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S7_linear_model_analysis_of_nonB_motif_clusters.R {output} 2> {log}"

rule Figure_S8:
    input:
        expand("data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_PCA_in_alignment.csv", taxon = TAXON)
    output:
        "plots/Figure_S8_PDALseq_PCA_different_taxon/Figure_S8.svg"
    log:
        "logs/Figure_S8_PDALseq_PCA_different_taxon/Figure_S8.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S8_PDALseq_PCA_different_taxon.R {output} 2> {log}"

rule Figure_S9:
    input:
        expand("data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_kmeans_scree_plot_in_alignment.csv", taxon = TAXON),
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/BIC.csv", taxon = TAXON)
    output:
        "plots/Figure_S9_Scree_plots_different_taxon/Figure_S9.svg"
    log:
        "logs/Figure_S9_Scree_plots_different_taxon/Figure_S9.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S9_PDALseq_screeplot_different_taxon.R {output} 2> {log}"

rule Figure_S10:
    input:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/BIC.csv",
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/{states}/Average_PDALseq_in_GHMM_states.csv", states = ["2", "4", "8", "18", "20"])
    output:
        "plots/Figure_S10_Scree_plots_explanation/Figure_S10.svg"
    log:
        "logs/Figure_S10_Scree_plots_explanation/Figure_S10.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S10_BIC_versus_hidden_states_interpratation.R {output} 2> {log}"

rule Figure_S11:
    input:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Transition_probabilities.csv",
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/Average_PDALseq_in_GHMM_states.csv"
    output:
        "plots/Figure_S11_HMM_transition_probabilities/Figure_S11.svg"
    log:
        "logs/Figure_S11_HMM_transition_probabilities/Figure_S11.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S11_HMM_transition_probabilities.R {output} 2> {log}"

rule Figure_S12:
    input:
        "data/09_Human_HMM/02_Normalize_human_window_coverage/Normalized_human_window_coverage.csv"
    output:
        "plots/Figure_S12_Human_PDALseq_PCA/Figure_S12.svg"
    log:
        "logs/Figure_S12_Human_PDALseq_PCA/Figure_S12.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S12_Human_PDALseq_PCA.R {output} 2> {log}"

rule Figure_S13:
    input:
        "data/09_Human_HMM/04_Multivariate_guassian/BIC.csv"
    output:
        "plots/Figure_S13_Human_MGHMM_BIC/Figure_S13.svg"
    log:
        "logs/Figure_S13_Human_MGHMM_BIC/Figure_S13.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/Figures_Tables/Figure_S13_Human_MGHMM_BIC.R {output} 2> {log}"

