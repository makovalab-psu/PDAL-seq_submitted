
#Remove alignments that come from hg38
rule remove_hg38:
    input:
        "resources/alignments/{alignments}.maf"
    output:
        "data/01_remove_hg38/{alignments}.maf"
    shell:
        "grep -v hg38 {input} > {output}"

#Prune blocks that do not contain all primates and are shorter than 10 bp
rule prune_ape_blocks:
    input:
        "data/01_remove_hg38/{alignments}.maf"
    output:
        temp("data/02_prune_blocks/{class}/{alignments}.maf")
    shell:
        "python workflow/scripts/Prune_{wildcards.class}_alignment.py {input} {output} 100"

#Keep only the alignment per species in a block. If duplicates are present keeps the sequence with the highest bit score match with the reference 
rule dup_filter:
    input:
        "data/02_prune_blocks/{class}/{alignments}.maf"
    output:
        temp("data/03_dup_filter/{class}/{alignments}.maf")
    shell:
        "python workflow/scripts/dup_filter.py {input} {output}"

#Filter out alignments that match the reference less than 50% based on bit scores
rule filter_quality:
    input:
        "data/03_dup_filter/{class}/{alignments}.maf"
    output:
        temp("data/04_filter_quality/{class}/{alignments}.maf")
    shell:
        "python workflow/scripts/Filter_quality.py {input} {output} 0.5 {wildcards.class}"

#Determine the original species and chromosomes names so I can rename them to match the PDALseq_manuscript pipeline
#rule original_species_chr_names:
#    input:
#        "data/04_filter_quality/{alignments}.maf"
#    output:
#        species = "data/05_original_species_chr_names/species_{alignments}.txt",
#        chr = "data/05_original_species_chr_names/chr_{alignments}.txt"
#    shell:
#        "python workflow/scripts/original_species_chr_names.py {input} {output.species} {output.chr}"

#Rename the species to match the PDALseq_manuscript pipeline
rule rename_species:
    input:
        #"data/05_original_species_chr_names/species_{alignments}.txt",
        maf = "data/04_filter_quality/{class}/{alignments}.maf",
        names = "resources/alignments/Rename_species.txt"
    output:
        temp("data/06_rename_species/{class}/{alignments}.maf")
    shell:
        "bash workflow/scripts/find_and_replace.sh {input.maf} {input.names} {output}"

#Rename the chromosome to match the PDALseq_manuscript pipeline
rule rename_chr:
    input:
        #"data/05_original_species_chr_names/chr_{alignments}.txt",
        maf = "data/06_rename_species/{class}/{alignments}.maf",
        names = "resources/alignments/Rename_chromosome.txt"
    output:
        temp("data/07_rename_chr/{class}/{alignments}.maf")
    shell:
        "bash workflow/scripts/find_and_replace.sh {input.maf} {input.names} {output}"

#Split the alignment into one file per human reference chromosome
rule split_by_human_ref_chromosomes:
    input:
        maf = "data/07_rename_chr/{class}/8-t2t-apes-2023v2.hs1.maf",
        chr_sizes = "resources/chr_sizes/Hsapien.chrom.sizes"
    output:
        temp(expand("data/08_split_by_ref_chromosomes/Hsapien/{{class}}/{chr}.maf", chr = Hsapien_chr))
    shell:
        "python workflow/scripts/split_reference_chr.py {input.maf} {input.chr_sizes} data/08_split_by_ref_chromosomes/Hsapien/{wildcards.class}/"

#Sort the alignment by start position in the human reference chromosome
rule sort:
    input:
        "data/08_split_by_ref_chromosomes/Hsapien/{class}/{chr}.maf"
    output:
        "data/09_sort/Hsapien/{class}/{chr}.maf"
    shell:
        "python workflow/scripts/sort_alignments_by_reference_position.py {input} {output}"

#Determine what chromosomes from other species are aligning to the reference
rule quantify_synteny:
    input:
        "data/09_sort/Hsapien/{class}/{chr}.maf"
    output:
        "data/10_quantify_syteny/Hsapien/{class}/{chr}.txt"
    shell:
        "python workflow/scripts/quantify_synteny.py {input} {output}"

#Identify the syntetic relationships based on abundance. Identify which chromosomes in gorilla make up 97.5% of the alinments to human chr1 
rule Identify_synteny:
    input:
        "data/10_quantify_syteny/Hsapien/{class}/{chr}.txt"
    output:
        "data/11_identify_syteny/Hsapien/{class}/{chr}.txt"
    shell:
        "Rscript workflow/scripts/identify_synteny.R {input} 0.8 {output}"

#Enforces the sytenic relationship observed in the previous step. Gets rid of blocks containing off alignments.
rule enforce_synteny:
    input:
        syn = "data/11_identify_syteny/Hsapien/{class}/{chr}.txt",
        maf = "data/09_sort/Hsapien/{class}/{chr}.maf"
    output:
        temp("data/12_enforce_syteny/Hsapien/{class}/{chr}.maf")
    shell:
        "python workflow/scripts/enforce_synteny.py {input.syn} {input.maf} {output}"

#Check how well the previous step worked
rule quantify_synteny_after_enforcement:
    input:
        "data/12_enforce_syteny/Hsapien/{class}/{chr}.maf"
    output:
        "data/13_quantify_syteny_after_enforcement/Hsapien/{class}/{chr}.txt"
    shell:
        "python workflow/scripts/quantify_synteny.py {input} {output}"

#Plots how well the previous step work. Plots should match the known sytenic relationship
rule check_synteny:
    input:
        "data/13_quantify_syteny_after_enforcement/Hsapien/{class}/{chr}.txt"
    output:
        "data/13_quantifycheck_synteny/Hsapien/{class}/{chr}.png"
    shell:
        "Rscript workflow/scripts/check_synteny.R {input} {output}"

#Merge into one maf file
rule merge_maf:
    input:
        expand("data/12_enforce_syteny/Hsapien/{{class}}/{chr}.maf", chr = Hsapien_chr)
    output:
        "data/14_merge_maf/{class}/Hsapien.maf"
    shell:
        "python workflow/scripts/merge_maf.py data/12_enforce_syteny/Hsapien/{wildcards.class}/ {output}"

#Convert to a 1 indexed bed file for calcualting read coverage
rule maf_to_bg:
    input:
        "data/14_merge_maf/{class}/Hsapien.maf"
    output:
        "data/15_maf_to_bg/{class}/{species}.bed"
    shell:
        "grep {wildcards.species} {input} | awk -f workflow/scripts/maf_to_bed.awk > {output}"

#Remove sex chromosome blocks
rule remove_sex_chromosome_blocks:
    input:
        ref = "data/15_maf_to_bg/{class}/Hsapien.bed",
        bed = "data/15_maf_to_bg/{class}/{species}.bed"
    output:
        "data/16_Remove_sex_chromosome_blocks/{class}/{species}.bed"
    shell:
        "bash workflow/scripts/Remove_sex_chromosome_blocks.sh {input.ref} {input.bed} {output}"

#Check alignment statisitics (how many alignments are there for each species, how much do they overlap, how much of the genome is relected?)
rule check_alignment_stats:
    input:
        bed = "data/16_Remove_sex_chromosome_blocks/{class}/{species}.bed",
        chr_sizes = "resources/chr_sizes/{species}.chrom.sizes"
    output:
        "data/check_alignment_stats/{class}/{species}.txt"
    shell:
        "bash workflow/scripts/check_alignment_stats.sh {input.bed} {input.chr_sizes} > {output}"


#Determine the length in basepairs of each alignment and count
rule check_alignment_lengths:
    input:
        "data/16_Remove_sex_chromosome_blocks/{class}/{species}.bed"
    output:
        "data/check_block_lengths/{class}/{species}.hist"
    shell:
        "cat {input} | awk -f workflow/scripts/bed_annotation_length.awk | sort | uniq -c > {output}"

#Make a histogram plot of alignment length
rule plot_alignment_lengths:
    input:
        expand("data/check_block_lengths/{{class}}/{species}.hist", species = SPECIES)
    output:
        "data/check_block_lengths/{class}/all.png"
    shell:
        "Rscript workflow/scripts/Plot_length.R data/check_block_lengths/{wildcards.class}/ {output}"

#Make a histogram plot of alignment length
rule median_alignment_lengths:
    input:
        expand("data/check_block_lengths/{{class}}/{species}.hist", species = SPECIES)
    output:
        "data/check_block_lengths/{class}/Meadian_alignment_lengths_by_species.txt"
    shell:
        "Rscript workflow/scripts/Median_allinment_length.R data/check_block_lengths/{wildcards.class}/ {output}"