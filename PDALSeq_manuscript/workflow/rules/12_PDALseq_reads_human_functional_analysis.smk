
rule read_enrichment:
    input:
        bam = "data/01_preprocess_data/All_reads/{experiment}.bam",
        bed = "resources/functional_annotations/Saswats_annotations/Hsapien/{annotation}.bed"
    output:
        "data/12_PDALseq/Read_enrichment/{experiment}/{annotation}.csv"
    log:
        "logs/12_PDALseq/Read_enrichment/{experiment}/{annotation}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/12_PDALseq_reads_human_functional_analysis/Read_enrichment.sh {input.bam} {input.bed} {output} 10000 100 2> {log}"
        
rule compile_read_enrichment:
    input:
        expand("data/12_PDALseq/Read_enrichment/{experiment}/{annotation}.csv", experiment = ["Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq", "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq"], annotation = ANNOTATIONS)
    output:
        "data/12_PDALseq/Read_enrichment.csv"
    log:
        "logs/12_PDALseq/Read_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/12_PDALseq_reads_human_functional_analysis/Compile_read_enrichment.R {output} 2> {log}"

RM_CLASS = ["DNA","DNA?","DNA/Crypton","DNA/Crypton-A","DNA/hAT","DNA/hAT?","DNA/hAT-Ac","DNA/hAT-Blackjack","DNA/hAT-Charlie","DNA/hAT-hAT19","DNA/hAT-Tag1","DNA/hAT-Tip100","DNA/hAT-Tip100?","DNA?/hAT-Tip100?","DNA/Kolobok","DNA/Merlin","DNA/MULE-MuDR","DNA/PIF-Harbinger","DNA/PiggyBac","DNA/TcMar","DNA/TcMar-Mariner","DNA/TcMar-Pogo","DNA/TcMar-Tc1","DNA/TcMar-Tc2","DNA/TcMar-Tigger","LINE/CR1","LINE/Dong-R4","LINE/I-Jockey","LINE/L1","LINE/L1-Tx1","LINE/L2","LINE/Penelope","LINE/RTE-BovB","LINE/RTE-X","Low_complexity","LTR","LTR/ERV1","LTR/ERVK","LTR/ERVL","LTR/ERVL-MaLR","LTR/Gypsy","RC/Helitron","Retroposon/SVA","rRNA","Satellite","Satellite/acro","Satellite/centr","Satellite/subtelo","scRNA","Simple_repeat","SINE/5S-Deu-L2","SINE/Alu","SINE/MIR","SINE/tRNA","SINE/tRNA-Deu","SINE/tRNA-RTE","snRNA","srpRNA","tRNA","Unknown"]
        
rule get_rm_annotation:
    input:
        "resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out"
    output:
        "data/12_PDALseq/rm_annotation/{rm_class}.bed"
    log:
        "logs/12_PDALseq/rm_annotation/{rm_class}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/bed_generate.sh {input} {wildcards.rm_class} > {output} 2> {log}"
	    
rule human_state_rm_enrichment:
    input:
        state = "data/09_Human_HMM/04_Multivariate_guassian/10/{state}.bed",
        ann = "data/12_PDALseq/rm_annotation/{rm_class}.bed"
    output:
        "data/12_PDALseq/human_state_rm_enrichment/{rm_class}/{state}.bed"
    log:
        "logs/12_PDALseq/human_state_rm_enrichment/{rm_class}/{state}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {input.state} {input.ann} {output} 10000 100 2> {log}"

rule compile_human_state_rm_enrichment:
    input:
        expand("data/12_PDALseq/human_state_rm_enrichment/{rm_class}/{state}.bed", state = ["state_0", "state_1", "state_2", "state_3", "state_4", "state_5", "state_6", "state_7", "state_8", "state_9"], rm_class = RM_CLASS)
    output:
        "data/12_PDALseq/Human_state_rm_enrichment.csv"
    log:
        "logs/12_PDALseq/Human_state_rm_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/12_PDALseq_reads_human_functional_analysis/Compile_human_state_rm_enrichment.R {output} 2> {log}"
	    
rule human_state_rm_intersect_bp:
    input:
        state = "data/09_Human_HMM/04_Multivariate_guassian/10/{state}.bed",
        ann = "data/12_PDALseq/rm_annotation/{rm_class}.bed"
    output:
        "data/12_PDALseq/human_state_rm_intersect_bp/{rm_class}/{state}.bed"
    log:
        "logs/12_PDALseq/human_state_rm_intersect_bp/{rm_class}/{state}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    shell:
        "bash workflow/scripts//12_PDALseq_reads_human_functional_analysis/Intersect.sh {input.state} {input.ann} {output} 2> {log}"	    
	    
rule compile_human_state_rm_intersect_bp:
    input:
        expand("data/12_PDALseq/human_state_rm_intersect_bp/{rm_class}/{state}.bed", state = ["state_0", "state_1", "state_2", "state_3", "state_4", "state_5", "state_6", "state_7", "state_8", "state_9"], rm_class = RM_CLASS)
    output:
        "data/12_PDALseq/Human_state_rm_intersect_bp.csv"
    log:
        "logs/12_PDALseq/Human_state_rm_intersect_bp.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/12_PDALseq_reads_human_functional_analysis/Compile_human_state_rm_intersect_bp.R {output} 2> {log}"  
           
rule baseline_intersect:
    input:
        state = "data/09_Human_HMM/01_Remove_sex_chromosomes/Hsapien_GM24385_PDALseq.bg",
        ann = "data/12_PDALseq/rm_annotation/{rm_class}.bed"
    output:
        "data/12_PDALseq/baseline_human_state_rm_intersect_bp/{rm_class}.csv"
    log:
        "logs/12_PDALseq/baseline_human_state_rm_intersect_bp/{rm_class}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    shell:
        "bash workflow/scripts/12_PDALseq_reads_human_functional_analysis/Intersect.sh {input.state} {input.ann} {output} 2> {log}"

rule compile_baseline_intersect:
    input:
        expand("data/12_PDALseq/baseline_human_state_rm_intersect_bp/{rm_class}.csv", rm_class = RM_CLASS)
    output:
        "data/12_PDALseq/baseline_human_state_rm_intersect_bp.csv"
    log:
        "logs/12_PDALseq/baseline_human_state_rm_intersect_bp.txt"
    conda:
        "../env/Rutility.yml"
    shell:
        "Rscript workflow/scripts/12_PDALseq_reads_human_functional_analysis/Compile_baseline_human_state_rm_intersect_bp.R {output} 2> {log}"



  
