####################################################################
# Step 26 — Ape PCA and Human HMM K-means scree with K562
#
# Adds K562 (Hsapien_K562_PDALseq) to two analyses without touching
# any prior code:
#   A) Hominoidae alignment-based PCA (mirrors step 07 chain)
#   B) Human 1 kb genome-window K-means scree (mirrors step 09 chain)
#
# All existing step 07 / step 09 / step 21 files are ancient().
# All outputs are under data/26_*/ or results/26_*/.
####################################################################

_K562_BAM_26 = ancient(
    "data/21_add_K562_to_analysis/01_preprocess/All_reads/Hsapien_K562_PDALseq.bam"
)

HOMINOIDAE_19_EXPERIMENTS = [
    "Hsapien_GM24385_PDALseq", "Hsapien_GM24143_PDALseq",
    "Hsapien_HEK293T_PDALseq", "Hsapien_CCL86-Raji_PDALseq",
    "Ptroglodytes_AG18354_PDALseq", "Ptroglodytes_AG18358_PDALseq",
    "Ppaniscus_KB8711_PDALseq", "Ppaniscus_KB14571_PDALseq",
    "Ggorilla_KB3781_PDALseq", "Ggorilla_F7549_PDALseq",
    "Ppygmaeus_AG05252_PDALseq", "Pabelii_AG06213_PDALseq",
    "Ssyndactylus_JAMBI_PDALseq", "Ssyndactylus_Karenina_PDALseq",
    "Hsapien_GM24385_TruSeq", "Ptroglodytes_AG18354_TruSeq",
    "Ppaniscus_KB8711_TruSeq", "Ggorilla_KB3781_TruSeq",
    "Ppygmaeus_AG05252_TruSeq", "Pabelii_AG06213_TruSeq",
    "Ssyndactylus_JAMBI_TruSeq",
]

_HOMINOIDAE_SPECIES = [
    "Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla",
    "Ppygmaeus", "Pabelii", "Ssyndactylus",
]

_HUMAN_4 = [
    "Hsapien_GM24385_PDALseq",
    "Hsapien_GM24143_PDALseq",
    "Hsapien_HEK293T_PDALseq",
    "Hsapien_CCL86-Raji_PDALseq",
]


############################
# Part A — Alignment PCA with K562
############################

rule K562PCA_alignment_coverage:
    input:
        bam = _K562_BAM_26,
        beds = ancient(expand(
            "resources/alignments/hominoidae/{species}.bed",
            species=_HOMINOIDAE_SPECIES,
        )),
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/01_alignment_coverage/hominoidae/Hsapien_K562_PDALseq.bg"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/01_alignment_coverage.txt"
    conda:
        "../env/07_PhyloHGMP_model.yml"
    resources:
        mem_mb=20000
    shell:
        "bash workflow/scripts/07_PhyloHGMP_model/alignment_coverage.sh "
        "{input.bam} hominoidae {output} Hsapien_K562_PDALseq 2> {log}"


rule K562PCA_normalize_alignment_coverage:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/01_alignment_coverage/hominoidae/Hsapien_K562_PDALseq.bg"
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/02_normalize_alignment_coverage/hominoidae/Hsapien_K562_PDALseq.csv"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/02_normalize_alignment_coverage.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=10000
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/Normalize_alignment_coverage.R {input} {output} 2> {log}"


rule K562PCA_compile_alignment_coverage:
    input:
        k562_csv = "data/26_Ape_PCA_and_Human_HMM_with_K562/02_normalize_alignment_coverage/hominoidae/Hsapien_K562_PDALseq.csv",
        existing = ancient(expand(
            "data/07_PhyloHGMP_model/02_Normalize_alignment_coverage/hominoidae/{experiment}.csv",
            experiment=HOMINOIDAE_19_EXPERIMENTS,
        )),
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage/hominoidae/alignment_coverage_matrix.csv"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/03_compile_alignment_coverage.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage_with_K562.R 2> {log}"


rule K562PCA_alignment_PCA:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage/hominoidae/alignment_coverage_matrix.csv"
    output:
        csv   = "results/26_Ape_PCA_and_Human_HMM_with_K562/hominoidae/PDALseq_PCA_in_alignment.csv",
        scree = "results/26_Ape_PCA_and_Human_HMM_with_K562/hominoidae/PDALseq_PCA_in_alignment_scree.svg",
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/04_alignment_PCA.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=20000
    shell:
        "Rscript workflow/scripts/26_Ape_PCA_and_Human_HMM_with_K562/PDALseq_PCA_with_K562.R "
        "{input} {output.csv} {output.scree} 2> {log}"


############################
# Part B — Human K-means scree with K562
############################

rule K562PCA_make_1kb_windows:
    input:
        ancient("resources/genomes/Hsapien.chrom.sizes")
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/01_genome_windows/1000_nucleotides/Hsapien.bed"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/05_make_1kb_windows.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=1000
    shell:
        "bash workflow/scripts/03_Correlation_and_black_list_unmappable_regions/make_windows.sh "
        "{input} {output} 1000 2> {log}"


rule K562PCA_1kb_window_coverage:
    input:
        bam     = _K562_BAM_26,
        windows = "data/26_Ape_PCA_and_Human_HMM_with_K562/01_genome_windows/1000_nucleotides/Hsapien.bed",
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/02_1kb_window_coverage/Hsapien_K562_PDALseq.bg"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/06_1kb_window_coverage.txt"
    conda:
        "../env/03_Correlation_and_black_list_unmappable_regions.yml"
    resources:
        mem_mb=4000
    shell:
        "megadepth {input.bam} --threads 1 --annotation {input.windows} --op sum > {output} 2> {log}"


rule K562PCA_remove_blacklist:
    input:
        bg        = "data/26_Ape_PCA_and_Human_HMM_with_K562/02_1kb_window_coverage/Hsapien_K562_PDALseq.bg",
        blacklist = ancient("results/03_Correlation_and_black_list_unmappable_regions/blacklist_files/Hsapien.bed"),
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/03_remove_blacklist/Hsapien_K562_PDALseq.bg"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/07_remove_blacklist.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    resources:
        mem_mb=1000
    shell:
        "bedtools intersect -v -a {input.bg} -b {input.blacklist} > {output} 2> {log}"


rule K562PCA_remove_sex_chr:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/03_remove_blacklist/Hsapien_K562_PDALseq.bg"
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/04_remove_sex_chr/Hsapien_K562_PDALseq.bg"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/08_remove_sex_chr.txt"
    conda:
        "../env/05_nonB_motif_cluster_analysis.yml"
    resources:
        mem_mb=1000
    shell:
        "grep -v chrX {input} | grep -v chrY | grep -v chrM > {output} 2> {log}"


rule K562PCA_normalize_human_window_coverage:
    input:
        k562_bg  = "data/26_Ape_PCA_and_Human_HMM_with_K562/04_remove_sex_chr/Hsapien_K562_PDALseq.bg",
        existing = ancient(expand(
            "data/09_Human_HMM/01_Remove_sex_chromosomes/{sample}.bg",
            sample=_HUMAN_4,
        )),
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/05_normalized_coverage/Normalized_human_window_coverage_with_K562.csv"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/09_normalize_human_window_coverage.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/26_Ape_PCA_and_Human_HMM_with_K562/Normalize_human_window_coverage_with_K562.R 2> {log}"


rule K562PCA_Human_Kmeans_scree:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/05_normalized_coverage/Normalized_human_window_coverage_with_K562.csv"
    output:
        scree = "results/26_Ape_PCA_and_Human_HMM_with_K562/Kmeans_scree_data.csv",
        svg   = "results/26_Ape_PCA_and_Human_HMM_with_K562/Kmeans_scree_plot.svg",
        cor   = "results/26_Ape_PCA_and_Human_HMM_with_K562/Correlation_data.csv",
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/10_Human_Kmeans_scree.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/26_Ape_PCA_and_Human_HMM_with_K562/Human_Kmeans_scree_with_K562.R 2> {log}"


############################
# Part C — Ape alignment GHMM with K562
############################

_HMM_STATES_26 = ["2", "4", "6", "8", "10", "12", "14", "16", "18", "20"]


rule K562PCA_extract_PDALseq_columns:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage/hominoidae/alignment_coverage_matrix.csv"
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage/hominoidae/PDALseq_only_alignment_coverage_matrix.csv"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/11_extract_PDALseq_columns.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=4000
    shell:
        "Rscript workflow/scripts/26_Ape_PCA_and_Human_HMM_with_K562/Extract_PDALseq_columns_with_K562.R {input} {output} 2> {log}"


rule K562PCA_alignment_GHMM:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage/hominoidae/PDALseq_only_alignment_coverage_matrix.csv"
    output:
        csv = "data/26_Ape_PCA_and_Human_HMM_with_K562/Multivariate_guassian_HMM/hominoidae/{states26}/GHMM_states.csv",
        pkl = "data/26_Ape_PCA_and_Human_HMM_with_K562/Multivariate_guassian_HMM/hominoidae/{states26}/GHMM.pkl",
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/12_alignment_GHMM/{states26}.txt"
    conda:
        "../env/hmm.yml"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/GHMM.py "
        "{input} {wildcards.states26} {output.csv} {output.pkl} 2> {log}"


rule K562PCA_extract_alignment_GHMM_BIC:
    input:
        pkls = expand(
            "data/26_Ape_PCA_and_Human_HMM_with_K562/Multivariate_guassian_HMM/hominoidae/{states26}/GHMM.pkl",
            states26=_HMM_STATES_26,
        ),
        data = "data/26_Ape_PCA_and_Human_HMM_with_K562/compile_alignment_coverage/hominoidae/PDALseq_only_alignment_coverage_matrix.csv",
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/Multivariate_guassian_HMM/hominoidae/BIC.csv"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/13_alignment_GHMM_BIC.txt"
    conda:
        "../env/hmm.yml"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/calc_BIC.py "
        "data/26_Ape_PCA_and_Human_HMM_with_K562/Multivariate_guassian_HMM/hominoidae/ "
        "{input.data} {output} 2> {log}"


rule K562PCA_plot_alignment_GHMM_BIC:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/Multivariate_guassian_HMM/hominoidae/BIC.csv"
    output:
        "results/26_Ape_PCA_and_Human_HMM_with_K562/hominoidae/Alignment_GHMM_BIC_plot.svg"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/14_plot_alignment_GHMM_BIC.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=500
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_BIC.R {input} {output} 2> {log}"


############################
# Part D — Human GHMM with K562
############################

rule K562PCA_Human_GHMM:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/05_normalized_coverage/Normalized_human_window_coverage_with_K562.csv"
    output:
        csv = "data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/{states26}/GHMM_states.csv",
        pkl = "data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/{states26}/GHMM.pkl",
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/15_Human_GHMM/{states26}.txt"
    conda:
        "../env/hmm.yml"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/07_PhyloHGMP_model/GHMM.py "
        "{input} {wildcards.states26} {output.csv} {output.pkl} 2> {log}"


rule K562PCA_extract_human_GHMM_BIC:
    input:
        pkls = expand(
            "data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/{states26}/GHMM.pkl",
            states26=_HMM_STATES_26,
        ),
        data = "data/26_Ape_PCA_and_Human_HMM_with_K562/05_normalized_coverage/Normalized_human_window_coverage_with_K562.csv",
    output:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/BIC.csv"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/16_human_GHMM_BIC.txt"
    conda:
        "../env/hmm.yml"
    resources:
        mem_mb=5000
    shell:
        "python workflow/scripts/26_Ape_PCA_and_Human_HMM_with_K562/calc_BIC_human_K562.py "
        "data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/ "
        "{input.data} {output} 2> {log}"


rule K562PCA_plot_human_GHMM_BIC:
    input:
        "data/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM/BIC.csv"
    output:
        "results/26_Ape_PCA_and_Human_HMM_with_K562/Human_GHMM_BIC_plot.svg"
    log:
        "logs/26_Ape_PCA_and_Human_HMM_with_K562/17_plot_human_GHMM_BIC.txt"
    conda:
        "../env/Rutility.yml"
    resources:
        mem_mb=500
    shell:
        "Rscript workflow/scripts/07_PhyloHGMP_model/plot_BIC.R {input} {output} 2> {log}"
