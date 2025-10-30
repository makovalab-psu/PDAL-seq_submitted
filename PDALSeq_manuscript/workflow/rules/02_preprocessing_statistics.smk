
rule Count_raw_fastq_reads:
	input:
		"PDALseq_data/Raw_fastq/{sample}.fastq.gz"
	output:
		"data/02_preprocessing_statistics/Raw_fastq/{sample}.txt"
	log:
		"logs/02_preprocessing_statistic/Raw_fastq/{sample}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/read_count_fastq.gz.sh {input} {output} 2> {log}"

rule Count_trimmed_fastq_reads:
	input:
		"PDALseq_data/Trimmed_fastq/{sample}.fastq.gz"
	output:
		"data/02_preprocessing_statistics/Trimmed_fastq/{sample}.txt"
	log:
		"logs/02_preprocessing_statistic/Trimmed_fastq/{sample}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/read_count_fastq.gz.sh {input} {output} 2> {log}"

rule Mapped_bam_flags:
	input:
		"PDALseq_data/Mapped_bam/{sample}.bam"
	output:
		"data/02_preprocessing_statistics/Mapped_bam/{sample}.txt"
	log:
		"logs/02_preprocessing_statistic/Mapped_bam/{sample}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"

rule Mapped_deduplicated_bam_flags:
	input:
		"PDALseq_data/Mapped_deduplicated_bam/{sample}.bam"
	output:
		"data/02_preprocessing_statistics/Mapped_deduplicated_bam/{sample}.txt"
	log:
		"logs/02_preprocessing_statistic/Mapped_deduplicated_bam/{sample}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"

rule Mapped_deduplicated_high_MapQ_bam_flags:
	input:
		"PDALseq_data/Mapped_deduplicated_high_MapQ_bam/{sample}.bam"
	output:
		"data/02_preprocessing_statistics/Mapped_deduplicated_high_MapQ_bam/{sample}.txt"
	log:
		"logs/02_preprocessing_statistic/Mapped_deduplicated_high_MapQ_bam/{sample}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"

Preprocess_steps = ["Raw_fastq", "Trimmed_fastq", "Mapped_bam", "Mapped_deduplicated_bam", "Mapped_deduplicated_high_MapQ_bam"]

rule Compile_preprocessing_statistics:
	input:
		expand("data/02_preprocessing_statistics/{step}/{sample}.txt", step = Preprocess_steps, sample = SAMPLES) ,
		samples = "PDALseq_samples.txt"
	output:
		"results/02_preprocessing_statistics/Bam_flags_by_technical_replicate.csv"
	log:
		"logs/02_preprocessing_statistic/Compile_preprocessing_statistics.txt"
	conda:
		"../env/Rutility.yml"
	resources:
		mem_mb=1000
	shell:
		"Rscript workflow/scripts/02_preprocessing_statistics/Compile_preprocessing_statistics.R {input.samples} {output} 2> {log}"

rule Merged_bam_flags:
	input:
		"data/01_preprocess_data/{data}.bam"
	output:
		"data/02_preprocessing_statistics/Merged/bam_flags/{data}.txt"
	log:
		"logs/02_preprocessing_statistics/Merged/bam_flags/{data}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"

rule score_primary_versus_supplemental_alignment:
	input:
		"data/01_preprocess_data/{data}.bam"
	output:
		"data/02_preprocessing_statistics/Merged/primary_versus_supplemental/{data}.txt"
	log:
		"logs/02_preprocessing_statistics/Merged/primary_versus_supplemental/{data}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/Estimate_unique_maps.sh {input} > {output} 2> {log}"

rule bam_quality:
	input:
		"data/01_preprocess_data/{data}.bam"
	output:
		"data/02_preprocessing_statistics/Merged/mapping_quality/{data}.txt"
	log:
		"logs/02_preprocessing_statistics/Merged/mapping_quality/{data}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"samtools view {input} | cut -f 5 | sort -T temp/ | uniq -c > {output} 2> {log}"


rule Compile_mapping_statistics_by_biolgical_sample:
	input:
		expand("data/02_preprocessing_statistics/Merged/mapping_quality/{mapq}/{data}.txt", mapq = MAPQ, data = DATA),
		expand("data/02_preprocessing_statistics/Merged/bam_flags/{mapq}/{data}.txt", mapq = MAPQ, data = DATA),
		expand("data/02_preprocessing_statistics/Merged/primary_versus_supplemental/{mapq}/{data}.txt", mapq = MAPQ, data = DATA),
		samples = "PDALseq_samples.txt"
	output:
		"results/02_preprocessing_statistics/Mapping_statistics_by_biological_sample.csv"
	log:
		"logs/02_preprocessing_statistic/Compile_mapping_statistics_by_biolgical_sample.txt"
	conda:
		"../env/Rutility.yml"
	resources:
		mem_mb=1000
	shell:
		"Rscript workflow/scripts/02_preprocessing_statistics/Compile_mapping_statistics_by_biolgical_sample.R {input.samples} {output} 2> {log}"

