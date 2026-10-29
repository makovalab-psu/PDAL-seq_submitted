
rule benchmark_trim_fastq:
    input:
        R1 = "resources/Epigenetic_datasets/{sra}_1.fastq.gz",
        R2 = "resources/Epigenetic_datasets/{sra}_2.fastq.gz"
    output:
        R1 = temp("data/08_CHIP_CutRun_Benchmarking/01_trim_fastq/{sra}_1.fastq.gz"),
        R2 = temp("data/08_CHIP_CutRun_Benchmarking/01_trim_fastq/{sra}_2.fastq.gz")
    log:
        "logs/08_CHIP_CutRun_Benchmarking/01_trim_fastq/{sra}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=1000
    threads:
        5
    shell:
        "fastp --thread {threads} --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} -p -P 50 -y -x -3 2> {log}"

rule benchmark_bwa_mem:
    input:
        expand("resources/genomes/Hsapien.fasta.{ifile}", ifile = ["amb", "ann", "pac", "bwt", "sa"]),
        fasta = "resources/genomes/Hsapien.fasta",
        R1 = "data/08_CHIP_CutRun_Benchmarking/01_trim_fastq/{sra}_1.fastq.gz",
        R2 = "data/08_CHIP_CutRun_Benchmarking/01_trim_fastq/{sra}_2.fastq.gz"
    output:
        temp("data/08_CHIP_CutRun_Benchmarking/02_bwa_mem/{sra}.sam")
    log:
        "logs/08_CHIP_CutRun_Benchmarking/02_bwa_mem/{sra}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=16000
    threads:
        5
    shell:
        "bwa mem -t {threads} -o {output} {input.fasta} {input.R1} {input.R2} 2> {log}"

rule benchmark_sort_sam:
    input:
         "data/08_CHIP_CutRun_Benchmarking/02_bwa_mem/{sra}.sam"
    output:
         temp("data/08_CHIP_CutRun_Benchmarking/03_sort_sam/{sra}.bam")
    log:
        "logs/08_CHIP_CutRun_Benchmarking/03_sort_sam/{sra}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools sort {input} -T temp/{wildcards.sra} -O bam -o {output} 2> {log}"

rule benchmark_remove_dup:
    input:
         "data/08_CHIP_CutRun_Benchmarking/03_sort_sam/{sra}.bam"
    output:
         "data/08_CHIP_CutRun_Benchmarking/04_remove_dup/{sra}.bam"
    log:
        "logs/08_CHIP_CutRun_Benchmarking/04_remove_dup/{sra}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools rmdup -s {input} {output} 2> {log}"


rule benchmark_bam_flags:
        input:
                "data/08_CHIP_CutRun_Benchmarking/04_remove_dup/{sra}.bam"
        output:
                "data/08_CHIP_CutRun_Benchmarking/05_bam_flags/{sra}.txt"
        log:
                "logs/08_CHIP_CutRun_Benchmarking/05_bam_flags/{sra}..txt"
        conda:
                "../env/01_preprocess_data.yml"
        resources:
                mem_mb=1000
        shell:
                "bash workflow/scripts/02_preprocessing_statistics/bam_flags.sh {input} {output} 2> {log}"


#Reads in windows

rule benchmark_window_coverage:
    input:
        bed = "data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/1000000_nucleotides/Hsapien.bed",
        bam = "data/08_CHIP_CutRun_Benchmarking/04_remove_dup/{sra}.bam"
    output:
        "data/08_CHIP_CutRun_Benchmarking/06_genome_window_coverage/1000000_nucleotides/{sra}.bg"
    log:
        "logs/08_CHIP_CutRun_Benchmarking/06_genome_window_coverage/1000000_nucleotides/{sra}.txt"
    conda:
        "../env/04_Genome_wide_comparison_to_nonB_motifs.yml"
    resources:
        mem_mb=1000
    shell:
        "megadepth {input.bam} --threads 1 --annotation ./{input.bed} --op sum > {output} 2> {log}"

rule Compile_benchmark_window_coverage:
    input:
        expand("data/08_CHIP_CutRun_Benchmarking/06_genome_window_coverage/1000000_nucleotides/{sra}.bg", sra = BENCHMARKS), 
        "data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/1000000_nucleotides/Hsapien_HEK293T_PDALseq.bg",
        SRA = "resources/Epigenetic_datasets/BG4_Rloop_CHIP_CutTag_experiments.txt"
    output:
        "results/08_CHIP_CutRun_Benchmarking/Genome_wide_coverage.csv"
    log:
        "logs/08_CHIP_CutRun_Benchmarking/Genome_wide_coverage.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/08_CHIP_CutRun_Benchmarking/Compile_benchmark_window_coverage.R {input.SRA} {output} 2> {log}"


