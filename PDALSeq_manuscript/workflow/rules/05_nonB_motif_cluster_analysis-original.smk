rule Add_nonB_name:
    input:
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/{species}/{nonb}.bed"
    output:
        temp("data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{species}/{nonb}.bed")
    log:
        "logs/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{species}/{nonb}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bash workflow/scripts/05_nonB_motif_cluster_analysis/add_nonB_type_DNA_clusters.sh {input} {wildcards.nonb} {output} 2> {log}"

rule Compile_and_sort:
    input:
        expand("data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{{species}}/{nonb}.bed", nonb = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"])
    output:
        temp("data/05_nonB_motif_cluster_analysis/02_compiled_and_sorted_nonB/{species}.bed")
    log:
        "logs/05_nonB_motif_cluster_analysis/02_compiled_and_sorted_nonB/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    resources:
         mem_mb=8000
    shell:
        "cat data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{wildcards.species}/* | sort -T temp/ -V -k1,1 -k2,2 > {output} 2> {log}"

   
rule Cluster_nonB:
    input:
        "data/05_nonB_motif_cluster_analysis/02_compiled_and_sorted_nonB/{species}.bed"
    output:
        temp("data/05_nonB_motif_cluster_analysis/03_Cluster_nonB/{size}/{species}.bed")
    log:
        "logs/05_nonB_motif_cluster_analysis/03_Cluster_nonB/{size}/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    resources:
         mem_mb=8000
    shell:
        "cat {input} | bedtools cluster -d {wildcards.size} > {output} 2> {log}"

rule nonB_cluster_matrix:
    input:
        "data/05_nonB_motif_cluster_analysis/03_Cluster_nonB/{size}/{species}.bed"
    output:
        "data/05_nonB_motif_cluster_analysis/04_nonB_cluster_matrix/{size}/{species}.bed"
    log:
        "logs/05_nonB_motif_cluster_analysis/04_nonB_cluster_matrix/{size}/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    resources:
         mem_mb=1000
    shell:
        "bash workflow/scripts/05_nonB_motif_cluster_analysis/nonB_cluster_matrix.sh {input} {output} 2> {log}"

rule Count_clusters:
    input:
        expand("data/05_nonB_motif_cluster_analysis/04_nonB_cluster_matrix/{size}/{species}.bed", size = CLUSTERSIZE, species = SPECIES)
    output:
        "results/05_nonB_motif_cluster_analysis/Cluster_count_vs_cluster_size.txt"
    log:
        "logs/05_nonB_motif_cluster_analysis/05_Cluster_count_vs_cluster_size.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    resources:
         mem_mb=1000
    shell:
        "wc -l data/05_nonB_motif_cluster_analysis/04_nonB_cluster_matrix/*/*.bed > {output} 2> {log}"

rule sort_nonB_cluster_matrix:
    input:
        matrix = "data/05_nonB_motif_cluster_analysis/04_nonB_cluster_matrix/{size}/{species}.bed",
        fai = "resources/genomes/{species}.chrom.sizes"
    output:
        "data/05_nonB_motif_cluster_analysis/06_sorted_nonB_cluster_matrix/{size}/{species}.bed"
    log:
        "logs/05_nonB_motif_cluster_analysis/06_sorted_nonB_cluster_matrix/{size}/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools sort -faidx {input.fai} -i {input.matrix} > {output} 2> {log}"

rule Remove_black_list:
    input:
        clusters = "data/05_nonB_motif_cluster_analysis/06_sorted_nonB_cluster_matrix/{size}/{species}.bed",
        blacklist = "results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/{species}.bed"
    output:
        "data/05_nonB_motif_cluster_analysis/07_blacklist_removed_sorted_nonB_cluster_matrix/{size}/{species}.bed"
    log:
        "logs/05_nonB_motif_cluster_analysis/07_blacklist_removed_sorted_nonB_cluster_matrix/{size}/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools intersect -v -a {input.clusters} -b {input.blacklist} > {output} 2> {log}"

rule Remove_sex_chromosomes:
    input:
        "data/05_nonB_motif_cluster_analysis/07_blacklist_removed_sorted_nonB_cluster_matrix/{size}/{species}.bed"
    output:
        "data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/{size}/{species}.bed"
    log:
        "logs/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/{size}/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "cat {input} | grep -v chrX | grep -v chrY > {output} 2> {log}"

rule Count_filtered_clusters:
    input:
        expand("data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/{size}/{species}.bed", size = ["1000", "800", "600", "400", "200", "150", "100", "50", "10"], species = SPECIES)
    output:
        "results/05_nonB_motif_cluster_analysis/Filtered_cluster_count_vs_cluster_size.txt"
    log:
        "logs/05_nonB_motif_cluster_analysis/09_Cluster_count_vs_cluster_size.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    resources:
         mem_mb=1000
    shell:
        "wc -l data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/*/*.bed > {output} 2> {log}"

rule nonb_cluster_coverage:
    input:
        expand("data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/200/{species}.bed", species = SPECIES),
        bai = "data/01_preprocess_data/All_reads/{sample}.bam.bai",
        bam = "data/01_preprocess_data/All_reads/{sample}.bam"
    output:
        "data/05_nonB_motif_cluster_analysis/09_nonb_cluster_coverage/200/{sample}.bg"
    log:
        "logs/05_nonB_motif_cluster_analysis/09_nonb_cluster_coverage/200/{sample}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    resources:
        mem_mb=20000
    shell:
    	"bash workflow/scripts/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.sh {input.bam} 200 {output} {wildcards.sample} 2> {log}"

rule cluster_GC_content:
    input:
        bed = "data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/200/{species}.bed",
        fasta = "resources/genomes/{species}.fasta"
    output:
        "data/05_nonB_motif_cluster_analysis/10_nonB_cluster_matrix_GC_content/200/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    log:
        "logs/05_nonB_motif_cluster_analysis/09_nonb_cluster_coverage/200/{species}.txt"
    resources:
        mem_mb=1000
    shell:
    	"bash workflow/scripts/05_nonB_motif_cluster_analysis/bed_to_GC.sh {input.bed} {input.fasta} {output} 2> {log}"
		

rule compile_nonb_cluster_coverage:
    input:
        expand("data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/200/{species}.bed", species = SPECIES),
        expand("data/05_nonB_motif_cluster_analysis/10_nonB_cluster_matrix_GC_content/200/{species}.txt", species = SPECIES),
        expand("data/05_nonB_motif_cluster_analysis/09_nonb_cluster_coverage/200/{sample}.bg", sample = DATA)
    output:
        
    log:
        "logs/05_nonB_motif_cluster_analysis/11_compile_nonb_cluster_coverage.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=8000
    shell:
    	"Rscript workflow/scripts/05_nonB_motif_cluster_analysis/compile_nonb_cluster_coverage.R 2> {log}"

rule linear_regression_analysis:
    input:
        "results/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.csv"
    output:
        "results/05_nonB_motif_cluster_analysis/linear_models/{experiment}.csv"
    log:
        "logs/05_nonB_motif_cluster_analysis/12_Linear_regression_analysis/{experiment}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=50000
    shell:
    	"Rscript workflow/scripts/05_nonB_motif_cluster_analysis/Linear_regression_analysis_of_nonB_DNA_motif_clusters.R {wildcards.experiment} {output} 2> {log}"
    	

rule compile_linear_regression_analysis:
    input:
        expand("results/05_nonB_motif_cluster_analysis/linear_models/{experiment}.csv", experiment = DATA)
    output:
        "results/05_nonB_motif_cluster_analysis/linear_models.csv"
    log:
        "logs/05_nonB_motif_cluster_analysis/13_Compile_linear_regression_analysis.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=50000
    shell:
    	"Rscript workflow/scripts/05_nonB_motif_cluster_analysis/Compile_linear_models.R  2> {log}"
    	
 
rule log10_linear_regression_analysis:
    input:
        "results/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.csv"
    output:
        "results/05_nonB_motif_cluster_analysis/log10_linear_models/{experiment}.csv"
    log:
        "logs/05_nonB_motif_cluster_analysis/log10_linear_models/{experiment}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=50000
    shell:
    	"Rscript workflow/scripts/05_nonB_motif_cluster_analysis/log10_Linear_regression_analysis_of_nonB_DNA_motif_clusters.R {wildcards.experiment} {output} 2> {log}"
    	
rule compile_log10_linear_regression_analysis:
    input:
        expand("results/05_nonB_motif_cluster_analysis/log10_linear_models/{experiment}.csv", experiment = DATA)
    output:
        "results/05_nonB_motif_cluster_analysis/log10_linear_models.csv"
    log:
        "logs/05_nonB_motif_cluster_analysis/log10_linear_models.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=50000
    shell:
    	"Rscript workflow/scripts/05_nonB_motif_cluster_analysis/Compile_log10_linear_models.R  2> {log}"

rule glm_linear_regression_analysis:
    input:
        "results/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.csv"
    output:
        "results/05_nonB_motif_cluster_analysis/glm_linear_models/{experiment}.csv"
    log:
        "logs/05_nonB_motif_cluster_analysis/glm_linear_models/{experiment}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=50000
    shell:
    	"Rscript workflow/scripts/05_nonB_motif_cluster_analysis/general_inear_regression_analysis_of_nonB_DNA_motif_clusters.R {wildcards.experiment} {output} 2> {log}"
    	
rule compile_glm_linear_regression_analysis:
    input:
        expand("results/05_nonB_motif_cluster_analysis/glm_linear_models/{experiment}.csv", experiment = DATA)
    output:
        "results/05_nonB_motif_cluster_analysis/glm_linear_models.csv"
    log:
        "logs/05_nonB_motif_cluster_analysis/glm_linear_models.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=50000
    shell:
    	"Rscript workflow/scripts/05_nonB_motif_cluster_analysis/Compile_glm_linear_models.R  2> {log}"
