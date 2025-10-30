
####Analysis of reads in parent family###

CENSAT = ["active_hor", "bSat", "cenSat", "ct", "dhor", "GAP", "GAP,HSat3", "gSat", "hor", "HSat1A", "HSat1B", "HSat2", "HSat2_3", "HSat3", "mixedSuperFamily", "mon", "mon-hor", "rDNA", "subTerm", "subTerm_mon", "subTerm_mon-hor"]

rule Catagorize_censat_by_class:
    input:
        "resources/repeats/linneas_annotations/chm13v2.0_censat_v2.0.bed",
        "resources/repeats/linneas_annotations/mGorGor1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPanPan1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPanTro3_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPonAbe1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPonPyg2_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mSymSyn1_v2.1_CenSat_v2.0.bed"
    output:
        expand("data/17_PDALseq_reads_in_CenSat/01_censat_classes/{species}/{censat}.bed", species = SPECIES, censat = CENSAT)
    log:
        "logs/16_lm_analysis_of_nonbmotifs/Catagorize_censat_by_class.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
    	"Rscript workflow/scripts/17_PDALseq_reads_in_CenSat/Catagorize_CenSat.R 2> {log}"

rule CenSat_read_enrichment:
    input:
        expand("data/17_PDALseq_reads_in_CenSat/01_censat_classes/{species}/{censat}.bed", species = SPECIES, censat = CENSAT),
        bam = "data/01_preprocess_data/All_reads/{experiment}.bam"
    output:
        "data/17_PDALseq_reads_in_CenSat/02_censat_enrichment_classes/{censat}/{experiment}.csv"
    log:
        "logs/17_PDALseq_reads_in_CenSat/02_censat_enrichment_classes/{censat}/{experiment}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/17_PDALseq_reads_in_CenSat/Read_enrichment.sh {input.bam} {output} 1000 100 {wildcards.experiment} {wildcards.censat} 2> {log}"

rule Compile_CenSat_enrichment:
    input:
        expand("data/17_PDALseq_reads_in_CenSat/02_censat_enrichment_classes/{censat}/{experiment}.csv", experiment = DATA, censat = CENSAT)
    output:
        "results/17_PDALseq_reads_in_CenSat/Censat_class_enrichment.csv"
    log:
        "logs/17_PDALseq_reads_in_CenSat/Censat_class_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/17_PDALseq_reads_in_CenSat/Compile_CenSat_Class_enrichment.R {output} 2> {log}"

rule calculate_nonB_motif_censat_intersect:
    input:
        expand("data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{species}/{nonb}.bed", species = SPECIES, nonb = ["all", "APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]),
        expand("data/17_PDALseq_reads_in_CenSat/01_censat_classes/{species}/{censat}.bed", species = SPECIES, censat = CENSAT)
    output:
        "results/17_PDALseq_reads_in_CenSat/Censat_nonB_intersect.txt"
    log:
        "logs/17_PDALseq_reads_in_CenSat/Censat_nonB_intersect.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "python3 workflow/scripts/17_PDALseq_reads_in_CenSat/Censat_nonB_intersect.py {output} 2> {log}"

####Analysis of reads in subfamillies####

df_subclass = pd.read_csv("data/17_PDALseq_reads_in_CenSat/00_scratch/Filtered_subclasses_for_analysis.csv", sep=',') 
SUBCLASS = df_subclass['file_path']

rule Catagorize_censat_by_subclass:
    input:
        "resources/repeats/linneas_annotations/chm13v2.0_censat_v2.0.bed",
        "resources/repeats/linneas_annotations/mGorGor1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPanPan1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPanTro3_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPonAbe1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPonPyg2_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mSymSyn1_v2.1_CenSat_v2.0.bed"
    output:
        expand("data/17_PDALseq_reads_in_CenSat/03_censat_subclasses/{species}/{subclass}.bed", species = SPECIES, subclass = SUBCLASS)
    log:
        "logs/16_lm_analysis_of_nonbmotifs/Catagorize_censat_by_subclass.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
    	"Rscript workflow/scripts/17_PDALseq_reads_in_CenSat/Catagorize_CenSat_by_subclass.R 2> {log}"

rule CenSat_subclass_read_enrichment:
    input:
        expand("data/17_PDALseq_reads_in_CenSat/03_censat_subclasses/{species}/{subclass}.bed", species = SPECIES, subclass = SUBCLASS),
        bam = "data/01_preprocess_data/All_reads/{experiment}.bam"
    output:
        "data/17_PDALseq_reads_in_CenSat/04_censat_enrichment_subclasses/{censat}/{experiment}.csv"
    log:
        "logs/17_PDALseq_reads_in_CenSat/02_censat_enrichment_subclasses/{censat}/{experiment}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/17_PDALseq_reads_in_CenSat/Read_enrichment.sh {input.bam} {output} 1000 100 {wildcards.experiment} {wildcards.censat} 2> {log}"

rule Compile_CenSat_subclass_enrichment:
    input:
        expand("data/17_PDALseq_reads_in_CenSat/04_censat_enrichment_subclasses/{censat}/{experiment}.csv", experiment = DATA, censat = SUBCLASS)
    output:
        "results/17_PDALseq_reads_in_CenSat/Censat_subclass_enrichment.txt"
    log:
        "logs/17_PDALseq_reads_in_CenSat/Censat_subclass_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/17_PDALseq_reads_in_CenSat/Compile_CenSat_subclass_enrichment.R {output} 2> {log}"


rule calculate_nonB_motif_censat_subclass_intersect:
    input:
        expand("data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{species}/{nonb}.bed", species = SPECIES, nonb = ["all", "APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]),
         expand("data/17_PDALseq_reads_in_CenSat/03_censat_subclasses/{species}/{subclass}.bed", species = SPECIES, subclass = SUBCLASS)
    output:
        "results/17_PDALseq_reads_in_CenSat/Censat_subclass_nonB_intersect.txt"
    log:
        "logs/17_PDALseq_reads_in_CenSat/Censat_subclass_nonB_intersect.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "python3 workflow/scripts/17_PDALseq_reads_in_CenSat/Censat_subclass_nonB_intersect.py {output} 2> {log}"

rule Get_SST1_SATR1_SATR2_LSAU_annotations_from_rm_out:
    input:
        "resources/repeats/linneas_annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out",
        "resources/repeats/linneas_annotations/mGorGor1.pri.cur.20231122.RepeatMasker-combo3.out",
        "resources/repeats/linneas_annotations/mPanPan1.pri.cur.20231122.RepeatMasker-combo3.out",
        "resources/repeats/linneas_annotations/mPanTro3.pri.cur.20231122.RepeatMasker-combo.out",
        "resources/repeats/linneas_annotations/mPonAbe1.pri.cur.20231205.RepeatMasker-combo.out",
        "resources/repeats/linneas_annotations/mPonPyg2.pri.cur.20231122.RepeatMasker-combo.out",
        "resources/repeats/linneas_annotations/mSymSyn1.pri.cur.20231205.RepeatMasker-combo.out"
    output:
        expand("data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/{sat}/{species}.bed", sat = ["SST1", "SATR1","SATR2", "LSAU"], species = SPECIES)
    log:
        "logs/17_PDALseq_reads_in_CenSat/05_Get_SST1_SATR1_SATR2_LSAU_annotations_from_rm_out.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/17_PDALseq_reads_in_CenSat/intersect_censat_nonB.sh 2> {log}"

rule Select_sat_class_enrichment:
    input:
        expand("data/17_PDALseq_reads_in_CenSat/05_specific_rm_catagories/{sat}/{species}.bed", sat = ["SST1", "SATR1","SATR2", "LSAU"], species = SPECIES),
        bam = "data/01_preprocess_data/All_reads/{experiment}.bam"
    output:
        "data/17_PDALseq_reads_in_CenSat/06_Read_enrichment_in_specific_rm_catagories/{sat}/{experiment}.csv"
    log:
        "logs/17_PDALseq_reads_in_CenSat/06_Read_enrichment_in_specific_rm_catagories/{sat}/{experiment}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/17_PDALseq_reads_in_CenSat/Select_sat_read_enrichment.sh {input.bam} {output} 1000 100 {wildcards.experiment} {wildcards.sat} 2> {log}"

rule Compile_select_sat_class_enrichment:
    input:
        expand("data/17_PDALseq_reads_in_CenSat/06_Read_enrichment_in_specific_rm_catagories/{sat}/{experiment}.csv", experiment = DATA, sat = ["SST1", "SATR1","SATR2", "LSAU"])
    output:
        "results/17_PDALseq_reads_in_CenSat/SST1_SATR1_SATR2_LSAU_enrichment.csv"
    log:
        "logs/17_PDALseq_reads_in_CenSat/SST1_SATR1_SATR2_LSAU_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/17_PDALseq_reads_in_CenSat/Compile_select_sat_class_enrichment.R {output} 2> {log}"

Altemose = ["active_HOR", "inactive_HOR", "diverged_monomeric_HOR", "HSat1", "HSat2", "HSat3", "bSat", "Other_satellite", "rDNA", "Transition_region"]

rule Catagorize_censat_by_Altemose_paper_class:
    input:
        "resources/repeats/linneas_annotations/chm13v2.0_censat_v2.0.bed",
        "resources/repeats/linneas_annotations/mGorGor1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPanPan1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPanTro3_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPonAbe1_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mPonPyg2_v2.0_CenSat_v2.0.bed",
        "resources/repeats/linneas_annotations/mSymSyn1_v2.1_CenSat_v2.0.bed"
    output:
        expand("data/17_PDALseq_reads_in_CenSat/07_Altemose_paper_censat_classes/{species}/{censat}.bed", species = SPECIES, censat = Altemose)
    log:
        "logs/16_lm_analysis_of_nonbmotifs/Catagorize_censat_by_Altemose_paper_class.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
    	"Rscript workflow/scripts/17_PDALseq_reads_in_CenSat/Catagorize_CenSat_by_Altemose_labels.R 2> {log}"

rule Altemose_paper_names_CenSat_read_enrichment:
    input:
        expand("data/17_PDALseq_reads_in_CenSat/07_Altemose_paper_censat_classes/{species}/{censat}.bed", species = SPECIES, censat = Altemose),
        bam = "data/01_preprocess_data/All_reads/{experiment}.bam"
    output:
        "data/17_PDALseq_reads_in_CenSat/08_Altemose_paper_censat_enrichment_classes/{censat}/{experiment}.csv"
    log:
        "logs/17_PDALseq_reads_in_CenSat/08_Altemose_paper_censat_enrichment_classes/{censat}/{experiment}.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/17_PDALseq_reads_in_CenSat/Altemose_read_enrichment.sh {input.bam} {output} 1000 100 {wildcards.experiment} {wildcards.censat} 2> {log}"

rule Compile_Altemose_labels_CenSat_enrichment:
    input:
        expand("data/17_PDALseq_reads_in_CenSat/08_Altemose_paper_censat_enrichment_classes/{censat}/{experiment}.csv", experiment = DATA, censat = Altemose)
    output:
        "results/17_PDALseq_reads_in_CenSat/08_Altemose_paper_censat_enrichment_classes.csv"
    log:
        "logs/17_PDALseq_reads_in_CenSat/Censat_class_enrichment.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/17_PDALseq_reads_in_CenSat/Compile_Altemose_paper_CenSat_Class_enrichment.R {output} 2> {log}"
    

rule calculate_nonB_motif_Altemose_censat_intersect:
    input:
        expand("data/05_nonB_motif_cluster_analysis/01_Add_nonB_name/{species}/{nonb}.bed", species = SPECIES, nonb = ["all", "APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z"]),
        expand("data/17_PDALseq_reads_in_CenSat/07_Altemose_paper_censat_classes/{species}/{censat}.bed", species = SPECIES, censat = Altemose)
    output:
        "results/17_PDALseq_reads_in_CenSat/Altemose_paper_Censat_nonB_intersect.txt"
    log:
        "logs/17_PDALseq_reads_in_CenSat/Altemose_paper_Censat_nonB_intersect.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "python3 workflow/scripts/17_PDALseq_reads_in_CenSat/Altemose_paper_Censat_nonB_intersect.py {output} 2> {log}"        

