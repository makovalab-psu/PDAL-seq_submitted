
RM_FAMILY = ["Psuedogene", "LINE", "SINE", "LTR", "DNA_transposons", "Satelitte", "Simple_repeat", "Low_complexity"]

rule get_rm_class:
    input:
        expand("data/12_PDALseq/rm_annotation/{rm_class}.bed", rm_class = RM_CLASS)
    output:
        expand("data/14_Repeat_masker_class_enrichment/01_compile_annotations/{family}.bed", family = RM_FAMILY)
    log:
        "log/14_Repeat_masker_class_enrichment/01_compile_annotations.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    shell:
        "bash workflow/scripts/14_Repeat_masker_class_enrichment_in_states/get_rm_class.sh 2> {log}"
        
rule rm_family_enrichment_in_state:
    input:
    	state = "data/09_Human_HMM/04_Multivariate_guassian/10/state_{state}.bed",
    	functional_element = "data/14_Repeat_masker_class_enrichment/01_compile_annotations/{family}.bed"
    output:
    	"data/14_Repeat_masker_class_enrichment/02_annotations_enrichment/state_{state}/{family}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/14_Repeat_masker_class_enrichment/02_annotations_enrichment/state_{state}/{family}_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {input.state} {input.functional_element} {output} 10000 100 2> {log}"


rule compile_rm_family_enrichment_in_state:
    input:
    	expand("data/14_Repeat_masker_class_enrichment/02_annotations_enrichment/state_{state}/{family}_enrichment.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"], family = RM_FAMILY)
    output:
    	"results/14_Repeat_masker_class_enrichment/03_RM_damily_enrichment.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/14_Repeat_masker_class_enrichment/03_RM_damily_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/14_Repeat_masker_class_enrichment_in_states/compile_RM_family_enrichment.R {output} 2> {log}"

