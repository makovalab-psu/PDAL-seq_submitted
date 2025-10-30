rule fastqgz_to_fasta:
	input:
		"PDALseq_data/Trimmed_fastq/{data}.fastq.gz"
	output:
		temp("data/06_reference_free_repeat_analysis/02_fasta/{data}.fasta")
	log:
		"logs/06_reference_free_repeat_analysis/01_fastqgz_to_fasta/{data}.txt"
	conda:
		"../env/06_reference_free_repeat_analysis.yml"
	resources:
		mem_mb=8000
	shell:
		"bash workflow/scripts/06_reference_free_repeat_analysis/fastqgz_to_fasta.sh {input} {output} 2> {log}"

rule trf:
	input:
		"data/06_reference_free_repeat_analysis/02_fasta/{data}.fasta"
	output:
		temp("data/06_reference_free_repeat_analysis/02_trf/{data}.dat")
	log:
		"logs/06_reference_free_repeat_analysis/02_trf/{data}.txt"
	conda:
		"../env/06_reference_free_repeat_analysis.yml"
	resources:
		mem_mb=8000
	shell:
		"workflow/scripts/06_reference_free_repeat_analysis/trf409.legacylinux64 {input} 2 7 7 80 10 50 2000 -l 6 -f -d -h -ngs > {output} 2> {log}"

rule parse_trf:
	input:
		"data/06_reference_free_repeat_analysis/02_trf/{data}.dat"
	output:
		temp("data/06_reference_free_repeat_analysis/03_trf_parsed/{data}.dat")
	log:
		"logs/06_reference_free_repeat_analysis/03_trf_parsed/{data}.txt"
	conda:
		"../env/06_reference_free_repeat_analysis.yml"
	resources:
		mem_mb=8000
	shell:
		"python workflow/scripts/06_reference_free_repeat_analysis/parseTRFngsKeepHeader.py {input} {output} 2> {log}"

rule raw_count:
	input:
		"data/06_reference_free_repeat_analysis/03_trf_parsed/{data}.dat"
	output:
		temp("data/06_reference_free_repeat_analysis/04_raw_count/{data}.dat")
	log:
		"logs/06_reference_free_repeat_analysis/04_raw_count/{data}.txt"
	conda:
		"../env/06_reference_free_repeat_analysis.yml"
	resources:
		mem_mb=8000
	shell:
		"bash workflow/scripts/06_reference_free_repeat_analysis/rawcounts.sh {input} {output}"

rule raw_length:
	input:
		"data/06_reference_free_repeat_analysis/03_trf_parsed/{data}.dat"
	output:
		temp("data/06_reference_free_repeat_analysis/05_raw_length/{data}.dat")
	log:
		"logs/06_reference_free_repeat_analysis/05_raw_length/{data}.txt"
	conda:
		"../env/06_reference_free_repeat_analysis.yml"
	resources:
		mem_mb=8000
	shell:
		"bash workflow/scripts/06_reference_free_repeat_analysis/rawlength.sh {input} {output}"

rule filter_raw:
	input:
		counts = "data/06_reference_free_repeat_analysis/04_raw_count/{data}.dat",
		lengths = "data/06_reference_free_repeat_analysis/05_raw_length/{data}.dat"
	output:
		temp("data/06_reference_free_repeat_analysis/06_filtered/{data}.dat")
	log:
		"logs/06_reference_free_repeat_analysis/06_filter_raw/{data}.txt"
	conda:
		"../env/06_reference_free_repeat_analysis.yml"
	resources:
		mem_mb=8000
	shell:
		"bash workflow/scripts/06_reference_free_repeat_analysis/filter_raw_files.sh {input.counts} {input.lengths} {output} 2> {log}"	

rule kbp_per_Mbp:
	input:
		fasta = "data/06_reference_free_repeat_analysis/02_fasta/{data}.fasta",
		filtered = "data/06_reference_free_repeat_analysis/06_filtered/{data}.dat"
	output:
		"data/06_reference_free_repeat_analysis/07_kbp_per_Mbp/{data}.dat"
	log:
		"logs/06_reference_free_repeat_analysis/07_kbp_per_Mbp/{data}.dat"
	conda:
		"../env/06_reference_free_repeat_analysis.yml"
	resources:
		mem_mb=8000
	shell:
		"bash workflow/scripts/06_reference_free_repeat_analysis/kbp_per_Mbp.sh {input.fasta} {input.filtered} {output} 2> {log}"	

rule Compile_kbp_per_Mbp:
	input:
		expand("data/06_reference_free_repeat_analysis/07_kbp_per_Mbp/{data}.dat", data = SAMPLES),
		samples = "PDALseq_samples.txt"
	output:
		"results/06_reference_free_repeat_analysis/Short_tandem_repeats_by_technical_replicate.csv"
	log:
		"logs/06_reference_free_repeat_analysis/08_Compile_kbp_per_Mbp.txt"
	conda:
		"../env/Rutility.yml"
	resources:
	 	 mem_mb=1000
	shell:
		"Rscript workflow/scripts/06_reference_free_repeat_analysis/Tidy_csv.R {input.samples} {output} 2> {log}"	


