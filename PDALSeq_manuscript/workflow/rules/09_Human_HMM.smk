
HUMAN_PDALseq = ["Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq"]

rule Remove_human_black_list:
    input:
        bg = "data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000_nucleotides/{sample}.bg",
        blacklist = "results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/Hsapien.bed"
    output:
        "data/09_Human_HMM/00_Remove_human_black_list/{sample}.bg"
    log:
        "logs/09_Human_HMM/00_Remove_human_black_list/{sample}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools intersect -v -a {input.bg} -b {input.blacklist} > {output} 2> {log}"

rule Remove_human_sex_chromosomes:
    input:
        "data/09_Human_HMM/00_Remove_human_black_list/{sample}.bg"
    output:
        "data/09_Human_HMM/01_Remove_sex_chromosomes/{sample}.bg"
    log:
        "logs/09_Human_HMM/01_Remove_sex_chromosomes/{sample}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "cat {input} | grep -v chrX | grep -v chrY | grep -v chrM > {output} 2> {log}"
        
rule Normalize_human_window_coverage:
    input:
        expand("data/09_Human_HMM/01_Remove_sex_chromosomes/{sample}.bg", sample = HUMAN_PDALseq)
    output:
        "data/09_Human_HMM/02_Normalize_human_window_coverage/Normalized_human_window_coverage.csv"
    log:
        "logs/09_Human_HMM/02_Normalize_human_window_coverage.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/09_Human_HMM/Normalize_human_window_coverage.R 2> {log}"

rule Human_Kmeans_scree:
    input:
        "data/09_Human_HMM/02_Normalize_human_window_coverage/Normalized_human_window_coverage.csv"
    output:
        "data/09_Human_HMM/03_Human_Kmeans_scree/Kmeans_scree_data.csv",
        "data/09_Human_HMM/03_Human_Kmeans_scree/Kmeans_scree_plot.svg"
    log:
        "logs/09_Human_HMM/03_Human_Kmeans_scree.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/09_Human_HMM/Human_Kmeans_scree.R 2> {log}"

rule Human_multivariate_guassian_HMM:
    input:
        "data/09_Human_HMM/02_Normalize_human_window_coverage/Normalized_human_window_coverage.csv"       
    output:
        csv = "data/09_Human_HMM/04_Multivariate_guassian/{states}/GHMM_states.csv",
        pkl = "data/09_Human_HMM/04_Multivariate_guassian/{states}/GHMM.pkl"
    conda:
        "../env/hmm.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/{states}/GHMM.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/GHMM.py {input} {wildcards.states} {output.csv} {output.pkl} 2> {log}"

rule extract_human_GHMM_parameters:
    input:
        pkl = "data/09_Human_HMM/04_Multivariate_guassian/{states}/GHMM.pkl",
        data =  "data/09_Human_HMM/02_Normalize_human_window_coverage/Normalized_human_window_coverage.csv"
    output:
        Mean = "data/09_Human_HMM/04_Multivariate_guassian/{states}/Means.csv",
        Var = "data/09_Human_HMM/04_Multivariate_guassian/{states}/Varience.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/{states}/extract_GHMM_parameters.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/extract_GHMM_parameters.py {input.pkl} {input.data} {output.Mean} {output.Var}  2> {log}"
  
rule mean_human_PDALseq_in_states:
    input:
        "data/09_Human_HMM/04_Multivariate_guassian/{states}/GHMM_states.csv"
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/{states}/Average_PDALseq_in_GHMM_states.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/{states}/Average_PDALseq_in_GHMM_states.txt"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Mean_PDALseq_signal_in_states.R {input} {output} 2> {log}"

rule plot_human_HMM_vs_real_means:
    input:
        PDALseq = "data/09_Human_HMM/04_Multivariate_guassian/{states}/Average_PDALseq_in_GHMM_states.csv",
        Model = "data/09_Human_HMM/04_Multivariate_guassian/{states}/Means.csv"
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/{states}/plot_HMM_vs_real_means.png"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/{states}/plot_HMM_vs_real_means.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Compare_GHMM_model_to_means.R {input.PDALseq} {input.Model} {output} 2> {log}"

rule plot_to_test_human_Guassian:
    input:
        PDALseq = "data/09_Human_HMM/04_Multivariate_guassian/{states}/GHMM_states.csv",
        Mean = "data/09_Human_HMM/04_Multivariate_guassian/{states}/Means.csv",
        Var = "data/09_Human_HMM/04_Multivariate_guassian/{states}/Varience.csv"
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/{states}/plot_to_test_Guassian.pdf"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/{states}/plot_to_test_Guassian.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_to_test_Guassian.R {input.PDALseq} {input.Mean} {input.Var} {output} 2> {log}"

rule extract_human_GHMM_BIC:
    input:
        expand("data/09_Human_HMM/04_Multivariate_guassian/{states}/GHMM.pkl", states = ["2","4","6","8","9","10","12","14","16","18","20"]),
        data = "data/09_Human_HMM/02_Normalize_human_window_coverage/Normalized_human_window_coverage.csv"
    output:
        BIC = "data/09_Human_HMM/04_Multivariate_guassian/BIC.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/BIC.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/09_Human_HMM/calc_BIC.py data/09_Human_HMM/04_Multivariate_guassian/ {input.data} {output.BIC} 2> {log}"
        
rule plot_human_GHMM_BIC:
    input:
        "data/09_Human_HMM/04_Multivariate_guassian/BIC.csv"
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/BIC_plot.png"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/BIC_plot.txt"
    resources:
        mem_mb=500
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_BIC.R {input} {output} 2> {log}"


####Determine the functional significance of HMM states####

rule annotate_human_state:
    input:
        bed = "data/09_Human_HMM/01_Remove_sex_chromosomes/Hsapien_GM24385_PDALseq.bg",
        states = "data/09_Human_HMM/04_Multivariate_guassian/10/GHMM_states.csv"
    output:
        expand("data/09_Human_HMM/04_Multivariate_guassian/10/state_{state}.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"])
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/Annotate_state.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/09_Human_HMM/Annotate_state.R {input.bed} {input.states} 2> {log}"

#Count the bp in each state

rule Count_BP_in_human_state:
    input:
        expand("data/09_Human_HMM/04_Multivariate_guassian/10/state_{state}.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"])     
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/BP_in_state.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/BP_in_state.txt"
    resources:
        mem_mb=5000
    shell:
        "bash workflow/scripts/09_Human_HMM/BP_in_state.sh > {output} 2> {log}"

#Non-B enrichment#
     
rule human_nonB_enrichment_in_state:
    input:
        state = "data/09_Human_HMM/04_Multivariate_guassian/10/state_{state}.bed",
        nonB = "results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Hsapien/{nonb}.bed"
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/state_{state}/{nonb}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/state_{state}/{nonb}_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {input.state} {input.nonB} {output} 10000 100 2> {log}"

rule human_compile_nonB_enrichment_in_state:
    input:
        expand("data/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/state_{state}/{nonb}_enrichment.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"], nonb = NONB)
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/State_nonB_enrichment/hsapiens.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/09_Human_HMM/compile_nonB_enrichment_in_state.R 2> {log}"

#Functional element enrichment in states#

rule human_functional_element_enrichment_in_state:
    input:
        state = "data/09_Human_HMM/04_Multivariate_guassian/10/state_{state}.bed",
        functional_element = "resources/functional_annotations/Saswats_annotations/Hsapien/{ann}.bed"
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/State_functional_element_enrichment/state_{state}/{ann}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/State_functional_element_enrichmentstate_{state}/{ann}_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {input.state} {input.functional_element} {output} 10000 100 2> {log}"        
        
rule human_compile_functional_element_enrichment_in_state:
    input:
        expand("data/09_Human_HMM/04_Multivariate_guassian/10/State_functional_element_enrichment/state_{state}/{ann}_enrichment.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"], ann = ANNOTATIONS)
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/State_functional_element_enrichment/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/State_functional_element_enrichment/hsapiens.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/09_Human_HMM/compile_functional_element_enrichment_in_state.R 2> {log}"   

#Methylation enrichment in states#

rule human_methylation_in_state:
    input:
        state = "data/09_Human_HMM/04_Multivariate_guassian/10/state_{state}.bed",
        m5C = "resources/methylationData/Hsapien.bed"
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/methylation/state_{state}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/methylation/state_{state}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/methylation_in_annotation.sh {input.state} {input.m5C} {output} 2> {log}"        
        
rule human_compile_methylation_enrichment_in_state:
    input:
        expand("data/09_Human_HMM/04_Multivariate_guassian/10/methylation/state_{state}.txt", state = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"])
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/methylation/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/methylation/hsapiens.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/09_Human_HMM/compile_methylation_in_state.R 2> {log}"

#RNA expression analysis

rule  human_RNA_expression_in_state:
    input:
        state = "data/09_Human_HMM/04_Multivariate_guassian/10/state_{state}.bed",
        RNA = "resources/expression/RNAseq_default.bg"
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/State_RNA_expression/state_{state}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/State_RNA_expression/state_{state}..txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/BG_enrichment.sh {input.state} {input.RNA} {output} 1000000 100 2> {log}"
        
rule human_compile_RNA_expression_in_state:
    input:
        expand("data/09_Human_HMM/04_Multivariate_guassian/10/State_RNA_expression/state_{state}.txt", state = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"])
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/State_RNA_expression/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/State_RNA_expression/hsapiens.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/09_Human_HMM/compile_RNA_expression_in_state.R 2> {log}"   

rule intersect_human_and_ape_hmm:
    input:
        expand("data/09_Human_HMM/04_Multivariate_guassian/10/State_RNA_expression/state_{state}.txt", state = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"]),
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7"])
    output:
        "data/09_Human_HMM/04_Multivariate_guassian/10/intersect_results.csv"
    conda:
        "../env/Pyutility_plus_bedtools.yml"
    log:
        "logs/09_Human_HMM/04_Multivariate_guassian/10/intersect_results.txt"
    resources:
        mem_mb=1000
    shell:
        "python workflow/scripts/09_Human_HMM/intersect_state_beds.py 2> {log}"   

