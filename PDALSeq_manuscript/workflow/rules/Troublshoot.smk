

####Determine the functional significance of HMM states####

rule annotate_state:
    input:
    	bed = "resources/alignments/hominoidea/Hsapien.bed",
    	states = "results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/GHMM_states.csv"
    output:
    	expand("results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_annotations/hsapiens/state_{state}.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"])
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/11_annotate_state/hominoidea/hsapiens/annotate_state.txt"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Annotate_state.R {input.bed} {input.states} 2> {log}"

rule functional_element_enrichment_in_state:
    input:
    	state = "results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_annotations/hsapiens/state_{state}.bed",
    	functional_element = "resources/functional_annotations/Saswats_annotations/Hsapien/{ann}.bed"
    output:
    	"results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_functional_element_enrichment/hsapiens/state_{state}/{ann}_enrichment.bed"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    log:
        "logs/07_PhyloHGMP_model/13_functional_element_enrichment_in_state/hominoidea/10/State_nonB_enrichment/hsapiens/state_{state}/{ann}_enrichment.txt"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/Enrichment.sh {input.state} {input.functional_element} {output} 10000 100 2> {log}"        
        
       
rule compile_nonB_enrichment_in_state:
    input:
    	expand("results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_nonB_enrichment/hsapiens/state_{state}/{nonb}_enrichment.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"], nonb = NONB)
    output:
    	"results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_nonB_enrichment/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/14_compile_nonB_enrichment_in_state.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/compile_nonB_enrichment_in_state.R 2> {log}"
        
rule compile_functional_element_enrichment_in_state:
    input:
    	expand("results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_functional_element_enrichment/hsapiens/state_{state}/{ann}_enrichment.bed", state = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"], ann = ANNOTATIONS)
    output:
    	"results/07_PhyloHGMP_model/alignment_clustering/hominoidea/10/State_functional_element_enrichment/hsapiens.csv"
    conda:
        "../env/Rutility.yml"
    log:
        "logs/07_PhyloHGMP_model/15_compile_functional_element_enrichment_in_state.txt"
    resources:
        mem_mb=1000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/compile_functional_element_enrichment_in_state.R 2> {log}"   
