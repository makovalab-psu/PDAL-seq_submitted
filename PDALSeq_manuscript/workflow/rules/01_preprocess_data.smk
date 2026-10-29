
#rule get_fastq:
#    input:
#        "PDALseq_samples.txt"
#    output:
#        "PDALseq_data/Raw_fastq/{SRA}.fastq.gz"
#    log:
#        "logs/01_preprocess_data/01_get_fasta/{SRA}.txt"
#    conda:
#        "../env/01_preprocess_data.yml"
#    resources:
#         mem_mb=1000
#    shell:
#        "bash workflow/scripts/01_preprocess_data/01_get_fastq.sh {wildcards.SRA} {output} 2> {log}"


rule trim_fastq:
    input:
        "PDALseq_data/Raw_fastq/{sample}.fastq.gz"
    output:
        fastq = temp("PDALseq_data/Trimmed_fastq/{sample}.fastq.gz")
        #report = "PDALseq_data/Trimmed_fastq/{sample}_fastp_report.html"
    log:
        "logs/01_preprocess_data/02_Trim_fastq/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=1000
    shell:
        "fastp --thread 1 -i {input} -o {output.fastq} -p -P 50 -y -x -3 2> {log}"


rule bwa_index:
    input:
        "resources/genomes/{species}.fasta"
    output:
        "resources/genomes/{species}.fasta.amb",
        "resources/genomes/{species}.fasta.ann",
        "resources/genomes/{species}.fasta.pac",
        "resources/genomes/{species}.fasta.bwt",
        "resources/genomes/{species}.fasta.sa"
    log:
        "logs/01_preprocess_data/03_bwa_index/{species}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=1000 
    shell:
        "bwa index {input} 2> {log}"

rule bwa_mem:
    input:
        expand("resources/genomes/{species}.fasta.{ifile}", species = SPECIES, ifile = ["amb", "ann", "pac", "bwt", "sa"]),
        fastq = "PDALseq_data/Trimmed_fastq/{sample}.fastq.gz"
    output:
        temp("PDALseq_data/Mapped_sam/{sample}.sam")
    log:
        "logs/01_preprocess_data/04_bwa_mem/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=16000
    threads:
        5
    shell:
        "bash workflow/scripts/01_preprocess_data/map_reads.sh {input.fastq} {wildcards.sample} {output} {threads} 2> {log}"

rule sort_sam:
    input:
         "PDALseq_data/Mapped_sam/{sample}.sam"
    output:
         temp("PDALseq_data/Mapped_bam/{sample}.bam")
    log:
        "logs/01_preprocess_data/05_sort_sam/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools sort {input} -T temp/{wildcards.sample} -O bam -o {output} 2> {log}"

rule remove_dup:
    input:
         "PDALseq_data/Mapped_bam/{sample}.bam"
    output:
         "PDALseq_data/Mapped_deduplicated_bam/{sample}.bam"
    log:
        "logs/01_preprocess_data/06_remove_dup/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools rmdup -s {input} {output} 2> {log}"

rule Remove_low_MapQ:
    input:
         "PDALseq_data/Mapped_deduplicated_bam/{sample}.bam"
    output:
         "PDALseq_data/Mapped_deduplicated_high_MapQ_bam/{sample}.bam"
    log:
        "logs/01_preprocess_data/07_remove_low_MapQ/{sample}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools view -b --min-MQ 20 {input} > {output} 2> {log}"
      

rule Make_merge_list_all_reads:
    input:
        expand("PDALseq_data/Mapped_deduplicated_bam/{sample}.bam", sample=SAMPLES),
        Samples = "PDALseq_samples.txt"
    output:
        "data/01_preprocess_data/All_reads/{experiment}.txt"
    log:
        "logs/01_preprocess_data/08_Make_merge_list/All_reads/{experiment}.txt"
    conda:
        "../env/Pyutility.yml"
    resources:
         mem_mb=1000
    shell:
        "python workflow/scripts/01_preprocess_data/Make_merge_list_all_reads.py {input.Samples} {wildcards.experiment} {output} 2> {log}"

rule Make_merge_list_high_MapQ_reads:
    input:
        expand("PDALseq_data/Mapped_deduplicated_high_MapQ_bam/{sample}.bam", sample=SAMPLES),
        Samples = "PDALseq_samples.txt"
    output:
        "data/01_preprocess_data/High_MapQ/{experiment}.txt"
    log:
        "logs/01_preprocess_data/08_Make_merge_list/High_MapQ/{experiment}.txt"
    conda:
        "../env/Pyutility.yml"
    resources:
         mem_mb=1000
    shell:
        "python workflow/scripts/01_preprocess_data/Make_merge_list_High_MAPQ_reads.py {input.Samples} {wildcards.experiment} {output} 2> {log}"

 
rule Merge_mapped_deduplicated_bam:
    input:
        expand("PDALseq_data/Mapped_deduplicated_bam/{sample}.bam", sample=SAMPLES),
        bamlist = "data/01_preprocess_data/All_reads/{experiment}.txt"
    output:
        "data/01_preprocess_data/All_reads/{experiment}.bam"
    log:
        "logs/01_preprocess_data/09_Merge_mapped_deduplicated_bam/{experiment}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools merge -b {input.bamlist} -o {output} 2> {log}"

rule Merge_mapped_deduplicated_high_MapQ_bam:
    input:
        expand("PDALseq_data/Mapped_deduplicated_high_MapQ_bam/{sample}.bam", sample=SAMPLES),
        bamlist = "data/01_preprocess_data/High_MapQ/{experiment}.txt"
    output:
        "data/01_preprocess_data/High_MapQ/{experiment}.bam"
    log:
        "logs/01_preprocess_data/10_Merge_mapped_deduplicated_high_MapQ_bam/{experiment}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools merge -b {input.bamlist} -o {output} 2> {log}"
        
rule index_bam:
    input:
        "data/01_preprocess_data/{experiment}.bam"
    output:
        "data/01_preprocess_data/{experiment}.bam.bai"
    log:
        "logs/01_preprocess_data/Index_bam/{experiment}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "samtools index {input} 2> {log}"
        
rule bam_to_bg:
    input:
       "data/01_preprocess_data/{data}.bam",
    output:
        temp("data/01_preprocess_data/{data}.bg")
    log:
        "logs/01_preprocess_data/11_bam_to_bg/{data}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "bedtools genomecov -bg -split -ibam {input} > {output} 2> {log}"

rule chr_size:
    input:
        "resources/genomes/{genome}.fasta"
    output:
        "resources/genomes/{genome}.chrom.sizes"
    log:
        "logs/01_preprocess_data/12_chr_size/{genome}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "faSize -detailed {input} > {output} 2> {log}"

rule bg_to_bigwig:
    input:
        expand("resources/genomes/{species}.chrom.sizes", species = SPECIES),
        bg = "data/01_preprocess_data/{data}.bg"
    output:
        "data/01_preprocess_data/{data}.bigwig"
    log:
        "logs/01_preprocess_data/13_bg_to_bigwig/{data}.txt"
    conda:
        "../env/01_preprocess_data.yml"
    resources:
         mem_mb=8000
    shell:
        "bash workflow/scripts/01_preprocess_data/bg_to_bigwig.sh {input.bg} {wildcards.data} {output} 2> {log}"


