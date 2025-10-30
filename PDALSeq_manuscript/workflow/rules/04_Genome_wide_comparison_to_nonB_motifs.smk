#Clean nonB annotations. I want to switch the old names with new names. I also want to remove scores for the G4 predictions.

OLD = ["chm13v2.0", "mPanTro3.pri.cur.20231122", "mPanPan1.pri.cur.20231122", "mGorGor1.pri.cur.20231122", "mPonAbe1.pri.cur.20231205", "mPonPyg2.pri.cur.20231122", "mSymSyn1.pri.cur.20240514"]

rule Clean_nonb:
    input:
        expand("resources/nonb_annotations/{old}.{nonb}.bed", old = OLD, nonb = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"])
    output:
        expand("results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/{species}/{nonb}.bed", species = SPECIES, nonb = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"])
    log:
        "logs/04_Genome_wide_comparison_to_nonB_motifs/01_Clean_nonb.txt"
    conda:
        "../env/04_Genome_wide_comparison_to_nonB_motifs.yml"
    shell:
        "bash workflow/scripts/04_Genome_wide_comparison_to_nonB_motifs/clean_nonb.sh 2> {log}"

rule Combine_nonB:
    input:
        expand("results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/{{species}}/{nonb}.bed", nonb = ["APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"])
    output:
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/{species}/all.bed"
    log:
        "logs/04_Genome_wide_comparison_to_nonB_motifs/02_Combine_nonb/{species}/all.txt"
    conda:
        "../env/04_Genome_wide_comparison_to_nonB_motifs.yml"
    shell:
        "bash workflow/scripts/04_Genome_wide_comparison_to_nonB_motifs/all.sh results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/{wildcards.species}/ {output} 2> {log}"

#Non-B DNA motif density

rule nonb_density:
    input:
        window = "data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/{size}_nucleotides/{species}.bed",
        nonb = "results/04_Genome_wide_comparison_to_nonB_motifs/nonb_annotations/{species}/{nonb}.bed"
    output:
        "data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/{size}_nucleotides/{species}/{nonb}.bg"
    log:
        "logs/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/{size}_nucleotides/{species}/{nonb}.txt"
    conda:
        "../env/04_Genome_wide_comparison_to_nonB_motifs.yml"
    shell:
        "bedtools coverage -a {input.window} -b {input.nonb} > {output} 2> {log}"

#Reads in windows

rule window_coverage:
    input:
        expand("data/03_Correlation_and_black_list_unmappable_regions/01_genome_windows/{{size}}_nucleotides/{genome}.bed", genome = SPECIES),
        bam = "data/01_preprocess_data/All_reads/{sample}.bam"
    output:
        "data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/{size}_nucleotides/{sample}.bg"
    log:
        "logs/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/{size}_nucleotides/{sample}.txt"
    conda:
        "../env/04_Genome_wide_comparison_to_nonB_motifs.yml"
    resources:
        mem_mb=1000
    shell:
    	"bash workflow/scripts/04_Genome_wide_comparison_to_nonB_motifs/genome_window_coverage.sh {input.bam} {wildcards.size} {output} {wildcards.sample} 2> {log}"

#Compile data

rule compile_reads_vs_nonb:
    input:
        expand("data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/{size}_nucleotides/{species}/{nonb}.bg", size = WINDOWS, species = SPECIES, nonb = NONB),
        expand("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/{size}_nucleotides/{sample}.bg", size = WINDOWS, sample = DATA)
    output:
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonB_correlation_coeficient.csv",
        "results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv"
    log:
        "logs/04_Genome_wide_comparison_to_nonB_motifs/05_compile_reads_vs_nonb.txt"
    conda:
    	"../env/Rutility.yml"
    resources:
        mem_mb=1000
    shell:
    	"Rscript workflow/scripts/04_Genome_wide_comparison_to_nonB_motifs/compile_reads_vs_nonb.R 2> {log}"

#Generate satallite annotations

rule hsat_2_3_annotation:
    input:
        "resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out",
        "resources/repeats/linneas_annotations/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out",
        "resources/repeats/linneas_annotations/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out",
        "resources/repeats/linneas_annotations/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out",
        "resources/repeats/linneas_annotations/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out",
        "resources/repeats/linneas_annotations/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out",
        "resources/repeats/linneas_annotations/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out"
    output:
    	expand("results/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotations/{species}.bed", species = SPECIES)
    log:
        "logs/04_Genome_wide_comparison_to_nonB_motifs/06_hsat_2_3_annotations.txt"
    conda:
        "../env/04_Genome_wide_comparison_to_nonB_motifs.yml"
    resources:
        mem_mb=1000
    shell:
    	"bash workflow/scripts/04_Genome_wide_comparison_to_nonB_motifs/hsat_2_3_annotation.sh 2> {log}"

rule siamang_subterminal_alphasat_annotation:
    input:
        "resources/repeats/linneas_annotations/mSymSyn1_v2.1_CenSat_v2.0.bed"
    output:
    	"results/04_Genome_wide_comparison_to_nonB_motifs/Ssyndactylus_subTerm_monhorSF4.bed"
    log:
        "logs/04_Genome_wide_comparison_to_nonB_motifs/07_Ssyndactylus_subTerm_mon.txt"
    conda:
        "../env/04_Genome_wide_comparison_to_nonB_motifs.yml"
    resources:
        mem_mb=1000
    shell:
       "bash workflow/scripts/04_Genome_wide_comparison_to_nonB_motifs/Siamang_subterminal_alphasat_annotation.sh 2> {log}"


