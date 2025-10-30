ORIGIN = ["core_origins_hglft_genome_e0596_1cc9e0", "stochastic_origins_hglft_genome_2ff9ed_1ccfd0"]

rule Get_core_annotations:
    input:
        "resources/functional_annotations/Saswats_annotations/Hsapien/{origin}.bed"
    output:
        "data/11_Origin_analysis_analysis/unscrambled/state_all/{origin}.bed"
    log:
        "logs/11_Origin_analysis_analysis/unscrambled/state_all/{origin}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "cp {input} {output} 2> {log}"

rule Filter_by_state_annotations:
    input:
        origin = "data/11_Origin_analysis_analysis/unscrambled/state_all/{origin}.bed",
        state = "data/10_Methylation_analysis/01_state_annotations/hominoidae/8/Hsapien/state_{state}.bed"
    output:
        "data/11_Origin_analysis_analysis/unscrambled/state_{state}/{origin}.bed"
    log:
        "logs/11_Origin_analysis_analysis/unscrambled/state_{state}/{origin}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools intersect -a {input.origin} -b {input.state} > {output} 2> {log}"

rule scramble_origin_annotations:
    input:
        "data/11_Origin_analysis_analysis/unscrambled/state_{state}/{origin}.bed"
    output:
        "data/11_Origin_analysis_analysis/scrambled/state_{state}/{origin}.bed"
    log:
        "data/11_Origin_analysis_analysis/scrambled/state_{state}/{origin}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools shuffle -i {input} -g resources/genomes/Hsapien.chrom.sizes > {output} 2> {log}"

rule Mean_signal:
    input:
        bg = "data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000_nucleotides/{sample}.bg",
        annotation = "data/11_Origin_analysis_analysis/{shuffle}/state_{state}/{origin}.bed"
    output:
        "data/11_Origin_analysis_analysis/Mean_signal/{sample}/{shuffle}/state_{state}/{origin}.csv"
    log:
        "logs/11_Origin_analysis_analysis/Mean_signal/{sample}/{shuffle}/state_{state}/{origin}.txt"
    conda:
        "../env/Rutility.yml"
    shell:
        "Rscript workflow/scripts/11_Origin_analysis_analysis/PDAL_seq_around_origins.R {input.bg} {input.annotation} 100 {output} 2> {log}"
        
rule Compile_mean_signal:
    input:
        expand("data/11_Origin_analysis_analysis/Mean_signal/{sample}/{shuffle}/state_{state}/{origin}.csv", sample = ["Hsapien_CCL86-Raji_PDALseq", "Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq"], shuffle = ["unscrambled", "scrambled"], state = ["all","4"], origin = ORIGIN)
    output:
        "results/11_Origin_analysis_analysis/PDALseq_signal_around_origins.csv"
    log:
        "results/11_Origin_analysis_analysis/PDALseq_signal_around_origins.txt"
    conda:
        "../env/Rutility.yml"
    shell:
        "Rscript workflow/scripts/11_Origin_analysis_analysis/Compile_PDAL_seq_around_origins.R {output} 2> {log}"
