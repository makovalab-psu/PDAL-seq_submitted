rule genome_wide_GC_content:
    input:
        bed = "data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/1000000_nucleotides/{species}.bed",
        fasta = "resources/genomes/{species}.fasta"
    output:
        "data/16_lm_analysis_of_nonbmotifs/00_genome_wide_GC_content/1000000_nucleotides/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/00_genome_wide_GC_content/1000000_nucleotides/{species}.txt"
    resources:
        mem_mb=1000
    shell:
    	"bash workflow/scripts/05_nonB_motif_cluster_analysis/bed_to_GC.sh {input.bed} {input.fasta} {output} 2> {log}"

rule Genome_wide_GLM_analysis:
    input:
        expand("data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/1000000_nucleotides/{species}/{nonb}.bg", size = WINDOWS, species = SPECIES, nonb = NONB),
        expand("data/16_lm_analysis_of_nonbmotifs/00_genome_wide_GC_content/1000000_nucleotides/{species}.txt", species = SPECIES),
        PDALseq = "data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000000_nucleotides/{experiment}.bg"
    output:
        "data/16_lm_analysis_of_nonbmotifs/01_glm/{experiment}.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/01_glm/{experiment}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Reads_vrs_nonb_glm.R {input.PDALseq} {output} 2> {log}"
    	
rule Compile_genome_wide_GLM_analysis:
    input:
        expand("data/16_lm_analysis_of_nonbmotifs/01_glm/{experiment}.csv", experiment = DATA)
    output:
        "results/16_lm_analysis_of_nonbmotifs/Genome_wide_glm.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/Genome_wide_glm.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Compile_data.R {output} 2> {log}"
    	
rule glm_genome_wide_GLM_individual_parameters:
    input:
        expand("data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/1000000_nucleotides/{species}/{nonb}.bg", size = WINDOWS, species = SPECIES, nonb = NONB),
        expand("data/16_lm_analysis_of_nonbmotifs/00_genome_wide_GC_content/1000000_nucleotides/{species}.txt", species = SPECIES),
        PDALseq = "data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000000_nucleotides/{experiment}.bg"
    output:
        "data/16_lm_analysis_of_nonbmotifs/01_glm_individual_parameters/{experiment}.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/01_glm_individual_parameters/{experiment}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Reads_vrs_nonb_glm_individual_fit.R {input.PDALseq} {output} 2> {log}"

rule compile_genome_wide_GLM_individual_parameters:
    input:
        expand("data/16_lm_analysis_of_nonbmotifs/01_glm_individual_parameters/{experiment}.csv", experiment = DATA)
    output:
        "results/16_lm_analysis_of_nonbmotifs/Genome_wide_glm_individual_parameters.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/Genome_wide_glm_individual_parameters.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Compile_genome_wide_glm_individual_fit.R {output} 2> {log}"
    	
rule Remove_glm_black_list:
    input:
        bg = "data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/1000_nucleotides/{species}.bed",
        blacklist = "results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/{species}.bed"
    output:
        "data/16_lm_analysis_of_nonbmotifs/02_remove_blacklist_from_windows/1000_nucleotides/{species}.bed"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/02_remove_blacklist_from_windows/1000_nucleotides/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "bedtools intersect -v -a {input.bg} -b {input.blacklist} > {output} 2> {log}"

rule Remove_glm_sex_chromosomes:
    input:
        "data/16_lm_analysis_of_nonbmotifs/02_remove_blacklist_from_windows/1000_nucleotides/{species}.bed"
    output:
        "data/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/{species}.bed"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
        "cat {input} | grep -v chrX | grep -v chrY > {output} 2> {log}"

rule glm_nonB_rates:
    input:
        window = "data/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/{species}.bed",
        nonb = "data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{species}/{nonb}.bed"
    output:
        "data/16_lm_analysis_of_nonbmotifs/04_nonb_densisty_in_windows/1000_nucleotides/{species}/{nonb}.bed"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/04_nonb_densisty_in_windows/1000_nucleotides/{species}/{nonb}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    shell:
         "bedtools coverage -a {input.window} -b {input.nonb} > {output} 2> {log}"

rule glm_coverage:
    input:
        expand("data/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/{species}.bed", species = SPECIES),
        bai = "data/01_preprocess_data/All_reads/{sample}.bam.bai",
        bam = "data/01_preprocess_data/All_reads/{sample}.bam"
    output:
        "data/16_lm_analysis_of_nonbmotifs/05_glm_coverage/1000_nucleotides/{sample}.bg"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/05_glm_coverage/1000_nucleotides/{sample}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=5000
    shell:
    	"bash workflow/scripts/16_lm_analysis_of_nonbmotifs/glm_coverage.sh {input.bam} {output} {wildcards.sample} 2> {log}"

rule glm_GC_content:
    input:
        bed = "data/16_lm_analysis_of_nonbmotifs/03_remove_sex_chromosomes_from_windows/1000_nucleotides/{species}.bed",
        fasta = "resources/genomes/{species}.fasta"
    output:
        "data/16_lm_analysis_of_nonbmotifs/06_glm_GC_content/1000_nucleotides/{species}.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/06_glm_GC_content/1000_nucleotides/{species}.txt"
    resources:
        mem_mb=1000
    shell:
    	"bash workflow/scripts/05_nonB_motif_cluster_analysis/bed_to_GC.sh {input.bed} {input.fasta} {output} 2> {log}"

rule glm_matrix:
    input:
        expand("data/16_lm_analysis_of_nonbmotifs/04_nonb_densisty_in_windows/1000_nucleotides/{species}/{nonb}.bed", species = SPECIES, nonb = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]),
        expand("data/16_lm_analysis_of_nonbmotifs/06_glm_GC_content/1000_nucleotides/{species}.txt", species = SPECIES),
        reads = "data/16_lm_analysis_of_nonbmotifs/05_glm_coverage/1000_nucleotides/{sample}.bg"
    output:
        "data/16_lm_analysis_of_nonbmotifs/07_glm_matrix/1000_nucleotides/{sample}.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/07_glm_matrix/1000_nucleotides/{sample}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/glm_coverage_matrix.R {input.reads} {output} {wildcards.sample} 2> {log}"

rule glm_unique_regions:
    input:
        "data/16_lm_analysis_of_nonbmotifs/07_glm_matrix/1000_nucleotides/{sample}.csv"
    output:
        "data/16_lm_analysis_of_nonbmotifs/07_glm_unique_regions/1000_nucleotides/{sample}.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/07_glm_unique_regions/1000_nucleotides/{sample}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Reads_vrs_nonb_glm_unique_mapping_regions.R {input} {output} 2> {log}"

rule compile_glm_unique_regions:
    input:
        expand("data/16_lm_analysis_of_nonbmotifs/07_glm_unique_regions/1000_nucleotides/{sample}.csv", sample = DATA)
    output:
        "results/16_lm_analysis_of_nonbmotifs/glm_unique_regions_1000_nucleotides.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/glm_unique_regions_1000_nucleotides.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Compile_glm_unique_mapping_regions.R {output} 2> {log}"

rule glm_unique_regions_individual_parameters:
    input:
        "data/16_lm_analysis_of_nonbmotifs/07_glm_matrix/1000_nucleotides/{sample}.csv"
    output:
        "data/16_lm_analysis_of_nonbmotifs/08_glm_unique_regions_individual_parameters/1000_nucleotides/{sample}.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/08_glm_unique_regions_individual_parameters/1000_nucleotides/{sample}.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Reads_vrs_nonb_glm_unique_mapping_regions_individual_fit.R {input} {output} 2> {log}"

rule compile_glm_unique_regions_individual_parameters:
    input:
        expand("data/16_lm_analysis_of_nonbmotifs/08_glm_unique_regions_individual_parameters/1000_nucleotides/{sample}.csv", sample = DATA)
    output:
        "results/16_lm_analysis_of_nonbmotifs/glm_unique_regions_1000_nucleotides_individual_parameters.csv"
    log:
        "logs/16_lm_analysis_of_nonbmotifs/glm_unique_regions_1000_nucleotides_individual_parameters.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
    	"Rscript workflow/scripts/16_lm_analysis_of_nonbmotifs/Compile_glm_unique_mapping_regions_individual_fit.R {output} 2> {log}"

