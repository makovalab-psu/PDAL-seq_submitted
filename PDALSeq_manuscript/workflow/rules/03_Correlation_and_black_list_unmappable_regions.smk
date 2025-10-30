
####Technical rep correllation#####

rule sort_sample_bam:
    input:
       "PDALseq_data/Mapped_deduplicated_bam/{sample}.bam"
    output:
        temp("data/03_Correlation_and_black_list_unmappable_regions/sorted_bam/{sample}.bam")
    log:
        "logs/03_Correlation_and_black_list_unmappable_regions/sorted_bam/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
         mem_mb=1000
    shell:
        "samtools sort -o {output} {input} 2> {log}"

rule make_genome_windows:
    input:
        "resources/genomes/{genome}.chrom.sizes"
    output:
        "data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/{size}_nucleotides/{genome}.bed"
    log:
        "logs/03_Correlation_and_black_list_unmappable_regions/{size}_nucleotides/{genome}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/make_windows.sh {input} {output} {wildcards.size} 2> {log}"


rule genome_window_coverage:
    input:
        expand("data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/{size}_nucleotides/{genome}.bed", size = WINDOWS, genome = SPECIES),
        bam = "data/03_Correlation_and_black_list_unmappable_regions/sorted_bam/{sample}.bam"
    output:
        "data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/{size}_nucleotides/{sample}.bed"
    log:
        "logs/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/{size}_nucleotides/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/genome_window_coverage.sh {input.bam} {wildcards.size} {output} {wildcards.sample} 2> {log}"
        
rule compile_correlation:
    input:
    	expand("data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/{size}_nucleotides/{sample}.bed", size = WINDOWS, sample = SAMPLES)
    output:
    	"results/03_Correlation_and_black_list_unmappable_regions/Technical_replicate_correlation.csv"
    log:
    	"logs/03_Correlation_and_black_list_unmappable_regions/04_compile_correlation.txt"
    conda:
    	"../env/Rutility.yml"
    resources:
                mem_mb=1000
    shell:
    	"Rscript workflow/scripts/03_Correlation_and_black_list_unmappable_regions/compile_correlation.R 2> {log}"


####Blacklist####


TRUSEQ = ["Ggorilla_KB3781_TruSeq",  "Ppygmaeus_AG05252_TruSeq", "Hsapien_GM24385_TruSeq",  "Ptroglodytes_AG18354_TruSeq", "Pabelii_AG06213_TruSeq",  "Ssyndactylus_JAMBI_TruSeq", "Ppaniscus_KB8711_TruSeq"]

rule window_coverage_for_blacklist:
    input:
        expand("data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/10000_nucleotides/{genome}.bed", genome = SPECIES),
        bam = "data/01_preprocess_data/{mapq}/{sample}.bam"
    output:
        "data/03_Correlation_and_black_list_unmappable_regions/04_genome_window_coverage_for_blacklist/{mapq}/{sample}.bg"
    log:
        "logs/03_Correlation_and_black_list_unmappable_regions/04_genome_window_coverage_for_blacklist/{mapq}/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
    	"bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/genome_window_coverage_for_blacklist.sh {input.bam} 10000 {output} {wildcards.sample} 2> {log}"


rule make_raw_blacklist_files:
    input:
        expand("data/03_Correlation_and_black_list_unmappable_regions/04_genome_window_coverage_for_blacklist/{mapq}/{sample}.bg", mapq = MAPQ, sample = TRUSEQ)
    output:
        expand("data/03_Correlation_and_black_list_unmappable_regions/05_raw_blacklist_files/{species}.bed", species = SPECIES)
    log:
        "data/03_Correlation_and_black_list_unmappable_regions/05_raw_blacklist_files.txt"
    conda:
    	"../env/Rutility.yml"
    resources:
        mem_mb=8000
    shell:
    	"Rscript workflow/scripts/03_Correlation_and_black_list_unmappable_regions/black_list_windows.R 1.1 2> {log}"

rule merge_blacklist:
    input:
        "data/03_Correlation_and_black_list_unmappable_regions/05_raw_blacklist_files/{species}.bed"
    output:
        "data/03_Correlation_and_black_list_unmappable_regions/06_merged_blacklist_files/{species}.bed"
    log:
        "logs/03_Correlation_and_black_list_unmappable_regions/06_merge_blacklist_files/{species}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
    	"bedtools merge -d 1 -i {input} > {output} 2> {log}"

rule add_sex_chromosomes_to_blacklist:
    input:
        bl = "data/03_Correlation_and_black_list_unmappable_regions/06_merged_blacklist_files/{species}.bed",
        chr_sizes = "resources/genomes/{species}.chrom.sizes"
    output:
        "results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/{species}.bed"
    log:
        "logs/03_Correlation_and_black_list_unmappable_regions/08_add_sex_chromosomes_to_blacklist/{species}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
    	"bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/add_sex_chromosomes_to_blacklist.sh {input.bl} {input.chr_sizes} {output} 2> {log}"

####Black list stats####

rule filter_reads_by_blacklist:
    input:
        expand("results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/{species}.bed", species = SPECIES),
        bam = "data/01_preprocess_data/All_reads/{sample}.bam"
    output:
        inside = temp("data/03_Correlation_and_black_list_unmappable_regions/07_filter_reads_by_blacklist/inside_blacklist/{sample}.bam"),
        outside = temp("data/03_Correlation_and_black_list_unmappable_regions/07_filter_reads_by_blacklist/outside_blacklist/{sample}.bam")
    log:
        "logs/03_Correlation_and_black_list_unmappable_regions/07_blacklist_excleded_bam/{sample}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
    	"bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/filter_reads_by_blacklist.sh {input.bam} {output.outside} {output.inside} {wildcards.sample} 2> {log}"

rule estimate_unique_maps:
	input:
		"data/03_Correlation_and_black_list_unmappable_regions/07_filter_reads_by_blacklist/{blacklist}/{sample}.bam"
	output:
		"data/03_Correlation_and_black_list_unmappable_regions/08_estimate_unique_maps_by_blacklist/{blacklist}/{sample}.txt"
	log:
		"logs/03_Correlation_and_black_list_unmappable_regions/08_estimate_unique_maps_by_blacklist/{blacklist}/{sample}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/Estimate_unique_maps.sh {input} > {output} 2> {log}"

BLACKLIST = ["inside_blacklist", "outside_blacklist"]

rule compile_unique_map_estimation:
	input:
		expand("data/03_Correlation_and_black_list_unmappable_regions/08_estimate_unique_maps_by_blacklist/{blacklist}/{sample}.txt", blacklist = BLACKLIST, sample = DATA)
	output:
		"results/03_Correlation_and_black_list_unmappable_regions/Unique_maps_inside_and_outside_blacklist.csv"
	log:
		"logs/03_Correlation_and_black_list_unmappable_regions/Unique_maps_inside_and_outside_blacklist.txt"
	conda:
		"../env/Rutility.yml"
	resources:
		mem_mb=1000
	shell:
		"Rscript workflow/scripts/03_Correlation_and_black_list_unmappable_regions/compile_unique_map_estimation.R 2> {log}"

   

