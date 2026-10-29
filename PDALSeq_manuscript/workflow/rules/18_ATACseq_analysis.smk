
rule ATACseq_trim_fastq:
    input:
        R1 = "resources/ATAC_Seq/{sra}_1.fastq.gz",
        R2 = "resources/ATAC_Seq/{sra}_2.fastq.gz"
    output:
        R1 = temp("data/18_ATACseq_analysis/01_trim_fastq/{sra}_1.fastq.gz"),
        R2 = temp("data/18_ATACseq_analysis/01_trim_fastq/{sra}_2.fastq.gz")
    log:
        "logs/18_ATACseq_analysis/01_trim_fastq/{sra}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=1000
    threads:
        5
    shell:
        "fastp --thread {threads} --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} -p -P 50 -y -x -3 2> {log}"

rule ATACseq_bwa_mem:
    input:
        expand("resources/genomes/Hsapien.fasta.{ifile}", ifile = ["amb", "ann", "pac", "bwt", "sa"]),
	fasta = "resources/genomes/Hsapien.fasta",
	R1 = "data/18_ATACseq_analysis/01_trim_fastq/{sra}_1.fastq.gz",
        R2 = "data/18_ATACseq_analysis/01_trim_fastq/{sra}_2.fastq.gz"
    output:
        temp("data/18_ATACseq_analysis/02_bwa_mem/{sra}.sam")
    log:
        "logs/18_ATACseq_analysis/02_bwa_mem/{sra}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=16000
    threads:
        5
    shell:
        "bwa mem -t {threads} -o {output} {input.fasta} {input.R1} {input.R2} 2> {log}"

rule ATACseq_sort_sam:
    input:
         "data/18_ATACseq_analysis/02_bwa_mem/{sra}.sam"
    output:
         temp("data/18_ATACseq_analysis/03_sort_sam/{sra}.bam")
    log:
        "logs/18_ATACseq_analysis/03_sort_sam/{sra}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools sort {input} -T temp/{wildcards.sra} -O bam -o {output} 2> {log}"

rule ATACseq_remove_dup:
    input:
         "data/18_ATACseq_analysis/03_sort_sam/{sra}.bam"
    output:
         "data/18_ATACseq_analysis/04_remove_dup/{sra}.bam"
    log:
        "logs/18_ATACseq_analysis/04_remove_dup/{sra}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools rmdup -s {input} {output} 2> {log}"

rule ATACseq_bam_flags:
	input:
		"data/18_ATACseq_analysis/04_remove_dup/{sra}.bam"
	output:
		"data/18_ATACseq_analysis/05_bam_flags/{sra}.txt"
	log:
		"logs/18_ATACseq_analysis/05_bam_flags/{sra}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"

rule ATACseq_bg:
	input:
		"data/18_ATACseq_analysis/04_remove_dup/{sra}.bam"
	output:
		"data/18_ATACseq_analysis/06_bedgraph/{sra}.bg"
	log:
		"logs/18_ATACseq_analysis/06_bedgraph/{sra}.txt"
	conda:
		"../env/01_preprocess_data.yml"
	resources:
		mem_mb=1000
	shell:
		"bedtools genomecov -ibam {input} -bg > {output} 2> {log}"

rule  ATACseq_in_state:
    input:
        state = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/state_{state}.bed",
    	bg = "data/18_ATACseq_analysis/06_bedgraph/{sra}.bg"
    output:
    	"data/18_ATACseq_analysis/07_ape_state_ATACseq/state_{state}/{sra}.csv"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/18_ATACseq_analysis/07_ape_state_ATACseq/state_{state}/{sra}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/BG_enrichment.sh {input.state} {input.bg} {output} 1000000 100 2> {log}"
        
rule compile_ATACseq_in_state:
    input:
    	expand("data/18_ATACseq_analysis/07_ape_state_ATACseq/state_{state}/{sra}.csv", state = ["0", "1", "2", "3", "4", "5", "6", "7"], sra = ATAC_SRA)
    output:
    	"results/18_ATACseq_analysis/ATACseq_read_enrichment.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/18_ATACseq_analysis/08_compile_ATAC_seq_read_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/18_ATACseq_analysis/compile_ATACseq_signal_in_state.R 2> {log}"       


