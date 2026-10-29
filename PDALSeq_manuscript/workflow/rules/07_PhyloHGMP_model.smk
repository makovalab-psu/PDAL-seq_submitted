
#Compute and compile PDAL-seq coverage in alignments

rule alignment_coverage:
    input:
        expand("resources/alignments/{{taxon}}/{species}.bed", species = SPECIES),
        bai = "data/01_preprocess_data/All_reads/{sample}.bam.bai",
        bam = "data/01_preprocess_data/All_reads/{sample}.bam"
    output:
        "data/07_PhyloHGMP_model/01_alignment_coverage/{taxon}/{sample}.bg"
    log:
        "logs/07_PhyloHGMP_model/01_alignment_coverage/{taxon}/{sample}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=20000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/alignment_coverage.sh {input.bam} {wildcards.taxon} {output} {wildcards.sample} 2> {log}"


rule Normalize_alignment_coverage:
    input:
        bg = "data/07_PhyloHGMP_model/01_alignment_coverage/{taxon}/{sample}.bg"   
    output:
        "data/07_PhyloHGMP_model/02_Normalize_alignment_coverage/{taxon}/{sample}.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/02_Normalize_alignment_coverage/{taxon}/{sample}.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Normalize_alignment_coverage.R {input.bg} {output} 2> {log}"


rule compile_alignment_coverage:
    input:
        expand("data/07_PhyloHGMP_model/02_Normalize_alignment_coverage/{{taxon}}/{sample}.csv", sample = DATA)    
    output:
        bed = "data/07_PhyloHGMP_model/alignment_coverage/{taxon}/alignment_coverage_bed.csv",
        mtrx = "data/07_PhyloHGMP_model/alignment_coverage/{taxon}/alignment_coverage_matrix.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/02_compile_alignment_coverage/{taxon}.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/compile_alignment_coverage.R {wildcards.taxon} {output.bed} {output.mtrx} 2> {log}"

#Determine unique mapping rates inside and outside of alignments

rule filter_reads_by_alignment:
    input:
        expand("resources/alignments/{{taxon}}/{species}.bed", species = SPECIES),
        bam = "data/01_preprocess_data/All_reads/{sample}.bam"
    output:
        inside = temp("data/07_PhyloHGMP_model/03_filter_reads_by_alignment/{taxon}/inside_alignment/{sample}.bam"),
        outside = temp("data/07_PhyloHGMP_model/03_filter_reads_by_alignment/{taxon}/outside_alignment/{sample}.bam")
    log:
        "logs/07_PhyloHGMP_model/03_filter_reads_by_alignment/{taxon}/{sample}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/filter_reads_by_alignment.sh {input.bam} {output.outside} {output.inside} {wildcards.sample} {wildcards.taxon} 2> {log}"

rule estimate_unique_in_alignment:
        input:
                "data/07_PhyloHGMP_model/03_filter_reads_by_alignment/{taxon}/{alignment}/{sample}.bam"
        output:
                "data/07_PhyloHGMP_model/04_estimate_unique_maps_in_alignment/{taxon}/{alignment}/{sample}.txt"
        log:
                "logs/07_PhyloHGMP_model/04_estimate_unique_maps_in_alignment/{taxon}/{alignment}/{sample}/{sample}.txt"
        conda:
                "../env/07_PhyloHGMP_model.yml"
        resources:
                mem_mb=1000
        shell:
                "bash workflow/scripts/02_preprocessing_statistics/Estimate_unique_maps.sh {input} > {output} 2> {log}"

ALIGNMENT = ["inside_alignment", "outside_alignment"]

rule compile_unique_map_estimation_in_alignment:
        input:
                expand("data/07_PhyloHGMP_model/04_estimate_unique_maps_in_alignment/{{taxon}}/{alignment}/{sample}.txt", alignment = ALIGNMENT, sample = DATA)
        output:
                "data/07_PhyloHGMP_model/05_Unique_maps_inside_and_outside_alignment/{taxon}/Unique_maps_inside_and_outside_alignment.csv"
        log:
                "logs/07_PhyloHGMP_model/05_Unique_maps_inside_and_outside_alignment/{taxon}.txt"
        conda:
                "../env/Rutility.yml"
        resources:
                mem_mb=1000
        shell:
                "Rscript workflow/scripts/07_PhyloHGMP_model/compile_unique_map_estimation.R {wildcards.taxon} {output} 2> {log}"

####Preliminary analysis####

rule alignment_correllation:
    input:
        "data/07_PhyloHGMP_model/alignment_coverage/{taxon}/alignment_coverage_matrix.csv"    
    output:
        "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_correllation_in_alignment.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/06_alignment_correllation/{taxon}.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/PDALseq_read_correllation_in_alignments.R {input} {output} {wildcards.taxon} 2> {log}"

rule alignment_PCA:
    input:
        "data/07_PhyloHGMP_model/alignment_coverage/{taxon}/alignment_coverage_matrix.csv"    
    output:
        csv = "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_PCA_in_alignment.csv",
        scree = "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_PCA_in_scree_plot.svg"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/07_alignment_PCA/{taxon}.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/PDALseq_PCA_in_alignments.R {input} {output.csv} {output.scree} {wildcards.taxon} 2> {log}"

rule plot_alignment_PCA:
    input:
         "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_PCA_in_alignment.csv"   
    output:
        "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_PCA_plot.svg"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/07_plot_alignment_PCA/{taxon}.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_PDALseq_PCA_in_alignments.R {input} {output} {wildcards.taxon} 2> {log}"

rule alignment_Kmeans_cluster:
    input:
        "data/07_PhyloHGMP_model/alignment_coverage/{taxon}/alignment_coverage_matrix.csv"       
    output:
        "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_kmeans_scree_plot_in_alignment.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/08_alignment_kmeans/{taxon}.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/PDALseq_kmeans_in_alignment.R {input} {output} {wildcards.taxon} 2> {log}"

rule plot_alignment_Kmeans_cluster:
    input:
        "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_kmeans_scree_plot_in_alignment.csv"       
    output:
        "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_kmeans_scree_plot_in_alignment.png"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/08_plot_alignment_kmeans/{taxon}.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_PDALseq_kmeans_in_alignment.R {input} {output} {wildcards.taxon} 2> {log}"

####Fit the data to HMMs and assess HMMS####

rule Extract_PDALseq_columns:
    input:
        "data/07_PhyloHGMP_model/alignment_coverage/{taxon}/alignment_coverage_matrix.csv"       
    output:
        "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_only_alignment_coverage_matrix.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/09_Extract_PDALseq_columns/{taxon}.txt"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Extract_PDALseq_columns.R {input} {output} {wildcards.taxon} 2> {log}"

rule Multivariate_guassian_HMM:
    input:
        "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_only_alignment_coverage_matrix.csv"       
    output:
        csv = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/GHMM_states.csv",
        pkl = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/GHMM.pkl"
    conda:
        "../env/hmm.yml"
    log:
        "logs/07_PhyloHGMP_model/10_GHMM/{taxon}/{states}.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/GHMM.py {input} {wildcards.states} {output.csv} {output.pkl} 2> {log}"

rule extract_GHMM_BIC:
    input:
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{{taxon}}/{states}/GHMM.pkl", states = HMM_state),
        data = "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        BIC = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/BIC.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/07_PhyloHGMP_model/10_GHMM/{taxon}/BIC.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/calc_BIC.py data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{wildcards.taxon}/ {input.data} {output.BIC} 2> {log}"
        
rule plot_GHMM_BIC:
    input:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/BIC.csv"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/BIC_plot.svg"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/10_GHMM_plot/{taxon}/BIC.txt"
    resources:
        mem_mb=500
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_BIC.R {input} {output} 2> {log}"

rule extract_GHMM_parameters:
    input:
        pkl = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/GHMM.pkl",
        data = "data/07_PhyloHGMP_model/alignment_clustering/{taxon}/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        Mean = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Means.csv",
        Var = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Varience.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/07_PhyloHGMP_model/10_GHMM/{taxon}/{states}/extract_GHMM_parameters.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/extract_GHMM_parameters.py {input.pkl} {input.data} {output.Mean} {output.Var}  2> {log}"

rule extract_GHMM_transistion_probabilities:
    input:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/GHMM.pkl"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Transition_probabilities.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Transition_probabilities.txt"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/Extract_transition_matrix.py {input} {output} 2> {log}"
  

rule mean_PDALseq_in_states:
    input:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/GHMM_states.csv"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Average_PDALseq_in_GHMM_states.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/11_mean_PDALseq_in_states/{taxon}/{states}/Average_PDALseq_in_GHMM_states.txt"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Mean_PDALseq_signal_in_states.R {input} {output} 2> {log}"

rule plot_HMM_vs_real_means:
    input:
        PDALseq = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Average_PDALseq_in_GHMM_states.csv",
        Model = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Means.csv"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/plot_HMM_vs_real_means.png"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/12_Multivariate_guassian_HMM/{taxon}/{states}/plot_HMM_vs_real_means.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Compare_GHMM_model_to_means.R {input.PDALseq} {input.Model} {output} 2> {log}"

rule plot_to_test_Guassian:
    input:
        PDALseq = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/GHMM_states.csv",
        Mean = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Means.csv",
        Var = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/Varience.csv"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/plot_to_test_Guassian.pdf"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/12_plot_HMM_vs_real_means/{taxon}/{states}/plot_to_test_Guassian.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_to_test_Guassian.R {input.PDALseq} {input.Mean} {input.Var} {output} 2> {log}"

rule PDALseq_conservation_in_states:
    input:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/GHMM_states.csv"
    output:
        mean = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/mean_state_conservation.csv",
        individual = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/individual_state_conservation.csv",
        plot = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/individual_state_conservation.svg"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/{taxon}/{states}/State_conservation.txt"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Conservation_score.R {input} {output.mean} {output.individual} {output.plot} 2> {log}"

####Determine the functional significance of HMM states####

rule annotate_state:
    input:
        bed = "resources/alignments/hominoidae/Hsapien.bed",
        states = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/GHMM_states.csv"
    output:
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7"])
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/13_annotate_state/hominoidae/8/annotate_state.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Annotate_state.R {input.bed} {input.states} 2> {log}"

#Count the bp in each state

rule Count_BP_in_state:
    input:
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7"])     
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/BP_in_state.csv"
    conda:
        "../env/hmm.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/.txt"
    resources:
        mem_mb=5000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/BP_in_state.sh > {output} 2> {log}"

#Non-B enrichment#
     
rule nonB_enrichment_in_state:
    input:
        state = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed",
        nonB = "results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/Hsapien/{nonb}.bed"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_nonB_enrichment/state_{state}/{nonb}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_nonB_enrichment/state_{state}/{nonb}_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {input.state} {input.nonB} {output} 10000 100 2> {log}"

rule compile_nonB_enrichment_in_state:
    input:
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_nonB_enrichment/state_{state}/{nonb}_enrichment.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7"], nonb = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"])
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_nonB_enrichment/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/14_compile_nonB_enrichment_in_state.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/compile_nonB_enrichment_in_state.R 2> {log}"

#Functional element enrichment in states#

rule annotate_recombination_elements:
    input:
        "resources/functional_annotations/Refseq/genomic_ucsc_names.gff"
    output:
        "resources/functional_annotations/Saswats_annotations/Hsapien/Recombination.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/functional_annotations/Saswats_annotations/Hsapien/Recombination.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Get_recombination_features.sh {input} {output} 2> {log}"  

rule functional_element_enrichment_in_state:
    input:
        state = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed",
        functional_element = "resources/functional_annotations/Saswats_annotations/Hsapien/{ann}.bed"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/state_{state}/{ann}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/state_{state}/{ann}_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {input.state} {input.functional_element} {output} 10000 100 2> {log}"        
     
rule compile_functional_element_enrichment_in_state:
    input:
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/state_{state}/{ann}_enrichment.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7"], ann = ANNOTATIONS)
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichment/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_element_enrichmen/hsapiens.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/compile_functional_element_enrichment_in_state.R 2> {log}"   

#Refseq annotations

rule subset_ref_seq:
    input:
        "resources/functional_annotations/Refseq/genomic_ucsc_names.gff"
    output:
        "resources/functional_annotations/Refseq/{ann}.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/07_PhyloHGMP_model/functional_annotations/Refseq/{ann}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/subset_ref_seq.sh {input} {wildcards.ann} {output} 2> {log}"

rule ref_seq_functional_element_enrichment:
    input:
        state = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed",
        functional_element = "resources/functional_annotations/Refseq/{ann}.bed"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_refseq_enrichment/state_{state}/{ann}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_refseq_enrichment/state_{state}/{ann}_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {input.state} {input.functional_element} {output} 10000 1000 2> {log}"

rule compile_Refseq_element_enrichment_in_state:
    input:
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_functional_refseq_enrichment/state_{state}/{ann}_enrichment.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7"], ann = ["meiotic_recombination_region", "mitotic_recombination_region", "non_allelic_homologous_recombination_region", "recombination_feature"])
    output:
        "data/07_PhyloHGMP_model/alignment_clustering/hominoidae/8/State_functional_refseq_enrichment/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/alignment_clustering/hominoidae/8/State_functional_refseq_enrichment/hsapiens.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/compile_Refseq_enrichment_in_state.R 2> {log}"  

#Methylation enrichment in states#

rule methylation_in_state:
    input:
        state = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed",
        m5C = "resources/methylationData/Hsapien.bed"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/methylation/state_{state}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/methylation/state_{state}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/methylation_in_annotation.sh {input.state} {input.m5C} {output} 2> {log}"        
        
rule compile_methylation_enrichment_in_state:
    input:
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/methylation/state_{state}.txt", state = ["0", "1", "2", "3", "4", "5", "6", "7"])
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/methylation/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/alignment_clustering/hominoidae/8/methylation/hsapiens.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/compile_methylation_in_state.R 2> {log}" 

#RNA expression analysis

rule  RNA_expression_in_state:
    input:
        state = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed",
        RNA = "resources/expression/RNAseq_default.bg"
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/state_{state}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/state_{state}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/BG_enrichment.sh {input.state} {input.RNA} {output} 1000000 100 2> {log}"
        
rule compile_RNA_expression_in_state:
    input:
        expand("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/state_{state}.txt", state = ["0", "1", "2", "3", "4", "5", "6", "7"])
    output:
        "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/State_RNA_expression/hsapiens.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/compile_RNA_expression_in_state.R 2> {log}"        
