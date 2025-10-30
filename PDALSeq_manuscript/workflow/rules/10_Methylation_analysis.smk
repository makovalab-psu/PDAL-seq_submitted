
rule annotate_state_for_methylation:
    input:
    	bed = "resources/alignments/hominoidae/{species}.bed",
    	states = "data/07_PhyloHGMP_model/Multivariate_guassian_HMM/hominoidae/8/GHMM_states.csv"
    output:
    	expand("data/10_Methylation_analysis/01_state_annotations/hominoidae/8/{{species}}/state_{state}.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7"])
    conda:
        "../env/Rutility.yml"
    log:
        "logs/10_Methylation_analysis/01_state_annotations/hominoidae/8/{species}.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/10_Methylation_analysis/Annotate_state.R {input.bed} {input.states} data/10_Methylation_analysis/01_state_annotations/hominoidae/8/{wildcards.species} 2> {log}"

rule methylation_in_ape_states:
    input:
    	state = "data/10_Methylation_analysis/01_state_annotations/hominoidae/8/{species}/state_{state}.bed",
    	m5C = "resources/methylationData/{species}.bed"
    output:
    	"data/10_Methylation_analysis/02_state_annotations/hominoidae/8/{species}/state_{state}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/10_Methylation_analysis/02_state_annotations/hominoidae/8/{species}/state_{state}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/methylation_in_annotation.sh {input.state} {input.m5C} {output} 2> {log}"        
        
rule compile_methylation_enrichment_in_ape_state:
    input:
    	expand("data/10_Methylation_analysis/02_state_annotations/hominoidae/8/{species}/state_{state}.txt", species = ["Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii"], state = ["0", "1", "2", "3", "4", "5", "6", "7"])
    output:
    	"results/10_Methylation_analysis/Methylation_by_state_and_species.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/10_Methylation_analysis/Methylation_by_state_and_species.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/10_Methylation_analysis/compile_methylation_in_state.R {output} 2> {log}" 
        

rule nonb_methylation_in_ape_states:
    input:
    	state = "data/10_Methylation_analysis/01_state_annotations/hominoidae/8/{species}/state_{state}.bed",
    	m5C = "resources/methylationData/{species}.bed",
    	nonB = "data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{species}/{nonb}.bed"
    output:
    	"data/10_Methylation_analysis/03_nonB_methylation_in_state/hominoidae/8/{species}/{nonb}/state_{state}.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/10_Methylation_analysis/03_nonB_methylation_in_state/hominoidae/8/{species}/{nonb}/state_{state}.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/methylation_in_nonB_annotation.sh {input.state} {input.m5C} {input.nonB} {output} 2> {log}"        

rule compile_nonB_methylation_in_ape_state:
    input:
    	expand("data/10_Methylation_analysis/03_nonB_methylation_in_state/hominoidae/8/{species}/{nonb}/state_{state}.txt", species = ["Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii"], state = ["0", "1", "2", "3", "4", "5", "6", "7"], nonb = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"])
    output:
    	"results/10_Methylation_analysis/Methylation_by_state_nonB_motif_and_species.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/10_Methylation_analysis/Methylation_by_state_nonB_motif_and_species.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/10_Methylation_analysis/compile_methylation_in_nonB_motifs_and_state.R {output} 2> {log}" 

