TSS = ["promoter_protCode_singleTrans"]

rule Get_TSS_annotations:
    input:
        "resources/functional_annotations/Saswats_annotations/Hsapien/{tss}.bed"
    output:
        "data/12_TSS_analysis_analysis/unscrambled/state_all/{tss}.bed"
    log:
        "logs/12_TSS_analysis_analysis/unscrambled/state_all/{tss}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "cp {input} {output} 2> {log}"

rule Filter_TSS_by_state_annotations:
    input:
        origin = "resources/functional_annotations/Saswats_annotations/Hsapien/promoter_protCode_singleTrans.bed",
        state = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed"
    output:
        "data/12_TSS_analysis_analysis/unscrambled/state_{state}/promoter_protCode_singleTrans.bed"
    log:
        "logs/12_TSS_analysis_analysis/unscrambled/state_{state}/promoter_protCode_singleTrans.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools intersect -a {input.origin} -b {input.state} > {output} 2> {log}"

rule scramble_TSS_annotations:
    input:
        "data/12_TSS_analysis_analysis/unscrambled/state_{state}/{tss}.bed"
    output:
        "data/12_TSS_analysis_analysis/scrambled/state_{state}/{tss}.bed"
    log:
        "data/12_TSS_analysis_analysis/scrambled/state_{state}/{tss}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools shuffle -i {input} -g resources/genomes/Hsapien.chrom.sizes > {output} 2> {log}"

rule Mean_TSS_signal:
    input:
        bg = "data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000_nucleotides/{sample}.bg",
        annotation = "data/12_TSS_analysis_analysis/{shuffle}/state_{state}/{tss}.bed"
    output:
        "data/12_TSS_analysis_analysis/Mean_signal/{sample}/{shuffle}/state_{state}/{tss}.csv"
    log:
        "logs/12_TSS_analysis_analysis/Mean_signal/{sample}/{shuffle}/state_{state}/{tss}.txt"
    conda:
        "../env/Rutility.yml"
    shell:
        "Rscript workflow/scripts/12_TSS_analysis_analysis/PDAL_seq_around_TSS.R {input.bg} {input.annotation} 100 {output} 2> {log}"
        
rule Compile_TSS_mean_signal:
    input:
        expand("data/12_TSS_analysis_analysis/Mean_signal/{sample}/{shuffle}/state_{state}/{tss}.csv", sample = ["Hsapien_CCL86-Raji_PDALseq", "Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq"], shuffle = ["unscrambled", "scrambled"], state = ["all", "0", "1", "2", "3", "4", "5", "6","7"], tss = TSS)
    output:
        "results/12_TSS_analysis_analysis/PDALseq_signal_around_TSS.csv"
    log:
        "results/12_TSS_analysis_analysis/PDALseq_signal_around_origins.txt"
    conda:
        "../env/Rutility.yml"
    shell:
        "Rscript workflow/scripts/12_TSS_analysis_analysis/Compile_PDAL_seq_around_TSS.R {output} 2> {log}"
